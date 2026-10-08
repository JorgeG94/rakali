!! Time-stepping driver for the shallow water solver
module rki_solver
   !! Calls compute kernels and manages time advancement.
   !! Time integration: SSP-RK2 (Strong Stability Preserving Runge-Kutta)
   !!   Stage 1: U* = U^n + dt * L(U^n)
   !!   Stage 2: U^{n+1} = 0.5 * U^n + 0.5 * (U* + dt * L(U*))
   use rki_constants, only: wp, CP_WATER, VCOORD_ZSTAR_FULL
   use rki_state, only: state_t
   use rki_grid, only: hgrid_t
   use rki_boundary, only: apply_boundary_conditions, apply_nh_boundary_conditions, &
                           bc_params_t
   use rki_kernel_sponge, only: apply_sponge
   use rki_kernel_flux, only: compute_flux_hll, compute_flux_hllc
   use rki_kernel_update, only: compute_update, compute_rk2_average, &
                                compute_ssp_combine
   use rki_kernel_timestep, only: compute_timestep, compute_timestep_baroclinic
   use rki_kernel_friction, only: apply_friction, apply_coriolis
   use rki_kernel_pressure, only: compute_depth_averaged_momentum, &
                                  save_barotropic_state, &
                                  compute_layer_from_barotropic, &
                                  compute_vertical_velocity, &
                                  solve_nh_pressure, &
                                  apply_nh_correction, &
                                  compute_rk2_average_layers, &
                                  compute_nh_baro_tendency, &
                                  apply_baro_nh_tendency, &
                                  reconcile_barotropic
   use rki_ml_common, only: ML_MIN_DEPTH
   use rki_ml_dynamics, only: ml_distribute_barotropic, &
                              ml_distribute_with_remap, &
                              ml_remap_tracer, &
                              ml_distribute_with_remap_full, &
                              ml_apply_layer_bcs_full, &
                              ml_compute_total_depth, &
                              ml_compute_barotropic_momentum, &
                              ml_compute_vertical_velocity, &
                              ml_compute_vertical_exchange, &
                              ml_rk2_average_layers, &
                              ml_init_layers_from_barotropic, &
                              ml_save_layer_state, &
                              ml_apply_layer_bcs, &
                              ml_dynamics_workspaces_cleanup, &
                              ml_apply_wind_to_surface_layer
   use rki_ml_tracers, only: ml_advect_tracer, &
                             ml_advect_tracer_ppm, &
                             ml_distribute_tracer, &
                             ml_vertical_exchange_tracer, &
                             ml_save_tracer, &
                             ml_rk2_average_tracer, &
                             ml_apply_tracer_bcs, &
                             ml_clamp_tracer, &
                             ml_horizontal_diffusion_tracer, &
                             ml_apply_surface_tracer_fluxes, &
                             ml_vertical_diffusion_tracer, &
                             ml_settling_tracer, &
                             ml_erode_bed_sediment, &
                             tracer_scratch_t, &
                             tracer_scratch_ensure, &
                             tracer_scratch_destroy
   use rki_recon_weno, only: ml_advect_tracer_weno, &
                             TRACER_RECON_PPM
   use rki_ml_eos, only: ml_equation_of_state, ml_equation_of_state_shared
   use rki_ml_pressure, only: ml_baroclinic_pressure_gradient, &
                              ml_baroclinic_pressure_gradient_jacobian
   use rki_ml_vertical_viscosity, only: ml_vertical_viscosity
   use rki_ml_kpp, only: ml_kpp_boundary_layer, ml_kpp_viscosity, &
                         ml_kpp_diffusion_tracers
   use rki_ml_keps, only: ml_keps_update, ml_keps_apply_momentum, &
                          ml_keps_apply_diffusion_tracer, &
                          ml_keps_friction_velocities
   use rki_ml_friction, only: ml_apply_bottom_friction, ml_apply_coriolis
   use rki_ml_horizontal_viscosity, only: ml_horizontal_viscosity_smag, &
                                          ml_horizontal_viscosity_biharmonic, &
                                          ml_horizontal_viscosity_cleanup
   use rki_ml_truncate, only: ml_truncate_velocity
   use rki_kernel_fused, only: fused_update_stage1, fused_update_stage2, &
                               fused_update_stage_generic
   use rki_timestepper, only: timestepper_t, make_timestepper, TIME_RK2
   use rki_kernel_fused_ml, only: ml_fused_post_update, ml_fused_rk2_finalise, &
                                  ml_fused_tracer_post
   use rki_profiler, only: profiler_start, profiler_stop

   implicit none
   private

   public :: solver_step
   public :: accumulate_boundary_mass_out
      !! Exposed for the analytical mass-budget test (test_console_stats).
   public :: accumulate_mass_src
      !! Exposed for the analytical mass-budget test (test_console_stats).
   public :: solver_enter_data
   public :: solver_exit_data
   public :: solver_update_host
   public :: solver_update_forcing
   public :: solver_update_state

   ! Persistent caller-owned scratch for the `pure` structured tracer
   ! kernels (upwind / PPM advection, horizontal diffusion).  Allocated
   ! once in solver_enter_data, released in solver_exit_data.  Module-level
   ! so the step routines can reach it without threading it through every
   ! call.
   type(tracer_scratch_t), save :: tracer_scratch

contains

   subroutine solver_compute_barotropic_flux(state)
      !! Dispatch to the HLL or HLLC barotropic flux kernel based on
      !! `state%barotropic%use_hllc`.  Runs on device — both kernels
      !! consume / produce the same arrays.
      type(state_t), intent(inout) :: state

      if (state%barotropic%use_hllc) then
         call compute_flux_hllc(state%barotropic%h, state%barotropic%hu, &
                                state%barotropic%hv, state%barotropic%b, &
                                state%barotropic%flux_h, state%barotropic%flux_hu, &
                                state%barotropic%flux_hv, &
                                state%barotropic%mass_flux_x, &
                                state%barotropic%mass_flux_y, &
                                state%grid)
      else
         call compute_flux_hll(state%barotropic%h, state%barotropic%hu, &
                               state%barotropic%hv, state%barotropic%b, &
                               state%barotropic%flux_h, state%barotropic%flux_hu, &
                               state%barotropic%flux_hv, &
                               state%barotropic%mass_flux_x, &
                               state%barotropic%mass_flux_y, &
                               state%grid)
      end if
   end subroutine solver_compute_barotropic_flux

   subroutine solver_enter_data(state)
      !! Copy state arrays to GPU device. Called once at simulation start;
      !! all subsequent computation runs on device until solver_exit_data.
      type(state_t), intent(inout) :: state
         !! Solution state to place on device

      integer :: it
      logical :: want_ppm, want_upwind, want_hdiff, want_weno

      ! Map the parent state struct first so the device knows the shape of
      ! state_t; otherwise attaching a sub-struct creates a "partial
      ! presence" condition that errors out when kernels reference state
      ! from an `!$acc parallel loop`.
      !$omp target enter data map(to: state)

      !$omp target enter data map(to: state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
      !$omp&                  state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
      !$omp&                  state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
      !$omp&                  state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0)

      ! Surface forcing: value-type member, so only its allocatable arrays
      ! need explicit attaches (scalars/descriptors came with copyin(state)).
      !$omp target enter data map(to: state%forcing%wind_x, state%forcing%wind_y, &
      !$omp&                  state%forcing%atm_pressure, &
      !$omp&                  state%forcing%precipitation, &
      !$omp&                  state%forcing%heat_flux, &
      !$omp&                  state%forcing%evap_minus_precip)

      ! Chapman BC eta_old (one slot per cell along each edge). Always
      ! mapped because bc_init_chapman unconditionally allocates all four
      ! arrays; the kernel only reads them when the matching edge is
      ! BC_CHAPMAN.
      !$omp target enter data map(to: state%bc%eta_old_west, state%bc%eta_old_east, &
      !$omp&                  state%bc%eta_old_south, state%bc%eta_old_north)

      ! Coupled hydrostatic multi-layer arrays
      if (state%use_multilayer) then
         !$omp target enter data map(to: state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref, &
         !$omp&                  state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
         !$omp&                  state%multilayer%h_layer0, state%multilayer%hu_layer0, state%multilayer%hv_layer0, &
         !$omp&                  state%multilayer%w_ml, state%multilayer%hu_ml_pre, state%multilayer%hv_ml_pre, &
         !$omp&                  state%multilayer%rho_layer, state%multilayer%rho_ref)
         if (allocated(state%multilayer%vcoord_z_ref_col)) then
            ! associate hoists the leaf so the attach path is one-deep
            ! (ifx OpenMP-target chokes on deep `state%multilayer%X` chains
            ! in map clauses; enforced by openmp_portability_lint).
            associate (vref_col => state%multilayer%vcoord_z_ref_col)
               !$omp target enter data map(to: vref_col)
            end associate
         end if

         ! KPP boundary-layer workspace (allocated only when the overlay is
         ! enabled). associate-leaf for the one-deep attach path.
         if (allocated(state%multilayer%kpp_bl_depth)) then
            associate (kpp_bl => state%multilayer%kpp_bl_depth, &
                       kpp_us => state%multilayer%kpp_u_star, &
                       kpp_bf => state%multilayer%kpp_buoy_flux)
               !$omp target enter data map(to: kpp_bl, kpp_us, kpp_bf)
            end associate
         end if

         ! k-epsilon interface-shaped prognostics (allocated only when the
         ! closure is on). State, not scratch — implicit diffusion reads the
         ! previous-step nu_t — so copyin and keep device-resident.
         if (allocated(state%multilayer%keps_tke)) then
            associate (ke_tke => state%multilayer%keps_tke, &
                       ke_eps => state%multilayer%keps_eps, &
                       ke_nut => state%multilayer%keps_nut, &
                       ke_kt => state%multilayer%keps_kt, &
                       ke_usb => state%multilayer%keps_ustar_b, &
                       ke_uss => state%multilayer%keps_ustar_s)
               !$omp target enter data map(to: ke_tke, ke_eps, ke_nut, ke_kt, ke_usb, ke_uss)
            end associate
         end if

         if (allocated(state%multilayer%bed_sediment)) then
            associate (bed_sed => state%multilayer%bed_sediment)
               !$omp target enter data map(to: bed_sed)
            end associate
         end if

         ! Generic tracer registry (salinity + temperature). Two-step
         ! attach: array descriptor first (so the device knows the shape of
         ! tracers(:)), then each element's hTr/hTr0. Without the parent
         ! copyin, NVHPC stdpar can't dereference tracers(it)%hTr from a
         ! do-concurrent body.
         if (allocated(state%multilayer%tracers)) then
            call solver_tracers_enter_data(state)
         end if

         ! Persistent, caller-owned scratch for the `pure` tracer kernels,
         ! allocated once up front (never lazily in the per-step kernel).
         ! Guard each group by what the config actually runs so a PPM-only
         ! run does not allocate the upwind snapshot (and vice-versa), and
         ! the hdiff scratch only exists when a tracer diffuses.
         if (allocated(state%multilayer%tracers)) then
            ! use_ppm_tracer selects the high-order branch; tracer_recon_code
            ! then picks PPM vs WENO within it, so exactly one of the PPM /
            ! WENO scratch groups is allocated (never both).
            want_ppm = state%multilayer%use_ppm_tracer .and. &
                       state%multilayer%tracer_recon_code == TRACER_RECON_PPM
            want_weno = state%multilayer%use_ppm_tracer .and. &
                        state%multilayer%tracer_recon_code /= TRACER_RECON_PPM
            want_upwind = .not. state%multilayer%use_ppm_tracer
            ! Horizontal-diffusion scratch is needed iff the hdiff kernel can
            ! fire for any tracer.  The kernel's per-step call guard reads the
            ! PER-TRACER hdiff_kappa, but that is SYNCED from the multilayer
            ! scalar every step (see solver_step) — at enter_data time the
            ! per-tracer copy may still be zero while the scalar is already set
            ! (config sets both; a direct scalar-set, as tests do, sets only the
            ! scalar).  So gate on the SCALAR (the source of truth, fixed before
            ! enter_data) plus the static do_horizontal_diffusion flag, never the
            ! not-yet-synced per-tracer kappa.
            want_hdiff = .false.
            if (state%multilayer%hdiff_kappa > 0.0_wp) then
               do it = 1, size(state%multilayer%tracers)
                  if (state%multilayer%tracers(it)%do_horizontal_diffusion) then
                     want_hdiff = .true.
                  end if
               end do
            end if
            call tracer_scratch_ensure(tracer_scratch, &
                                       state%grid%nx_total, state%grid%ny_total, &
                                       state%multilayer%nz_ml, &
                                       want_ppm, want_upwind, want_hdiff, &
                                       want_weno=want_weno)
         end if
      end if

      ! Non-hydrostatic arrays
      if (state%use_nonhydrostatic) then
         !$omp target enter data map(to: state%nh%dsig, state%nh%u_layer, state%nh%v_layer, &
         !$omp&                  state%nh%w, state%nh%q_nh, &
         !$omp&                  state%nh%u_layer0, state%nh%v_layer0, state%nh%w0, &
         !$omp&                  state%nh%cg_r, state%nh%cg_p, state%nh%cg_Ap, &
         !$omp&                  state%nh%hu_pre, state%nh%hv_pre)
         if (state%nh%nh_preconditioner) then
            !$omp target enter data map(to: state%nh%cg_z)
         end if
      end if

      ! Mode splitting arrays
      if (allocated(state%nh%h_avg)) then
         !$omp target enter data map(alloc: state%nh%h_avg, state%nh%hu_avg, state%nh%hv_avg, &
         !$omp&                  state%nh%h_bt0, state%nh%hu_bt0, state%nh%hv_bt0, &
         !$omp&                  state%nh%baro_tend_x, state%nh%baro_tend_y)
      end if

   end subroutine solver_enter_data

   subroutine solver_exit_data(state)
      !! Copy state arrays back from GPU device. Called once at simulation
      !! end or for final output.
      type(state_t), intent(inout) :: state
         !! Solution state to retrieve from device

      !$omp target exit data map(from: state%barotropic%h, state%barotropic%hu, state%barotropic%hv)
      ! The Chapman BC's eta_old_* is PERSISTENT state (the previous step's
      ! boundary surface); it must round-trip, or every snapshot / restart
      ! write (which exits and re-enters the device data) would re-enter with
      ! a stale host copy.
      !$omp target exit data map(from: state%bc%eta_old_west, state%bc%eta_old_east, &
      !$omp&                  state%bc%eta_old_south, state%bc%eta_old_north)
      !$omp target exit data map(delete: state%barotropic%b, state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
      !$omp&                 state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
      !$omp&                 state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0)

      ! Surface forcing: only the allocatable arrays need detaching; the
      ! value-type forcing struct drops with the parent below.
      !$omp target exit data map(delete: state%forcing%wind_x, state%forcing%wind_y, &
      !$omp&                 state%forcing%atm_pressure, &
      !$omp&                 state%forcing%precipitation, &
      !$omp&                 state%forcing%heat_flux, &
      !$omp&                 state%forcing%evap_minus_precip)

      ! Coupled hydrostatic multi-layer arrays
      if (state%use_multilayer) then
         !$omp target exit data map(from: state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
         !$omp&                 state%multilayer%w_ml, state%multilayer%rho_layer)
         !$omp target exit data map(delete: state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref, &
         !$omp&                 state%multilayer%h_layer0, state%multilayer%hu_layer0, state%multilayer%hv_layer0, &
         !$omp&                 state%multilayer%hu_ml_pre, state%multilayer%hv_ml_pre, state%multilayer%rho_ref)
         if (allocated(state%multilayer%vcoord_z_ref_col)) then
            associate (vref_col => state%multilayer%vcoord_z_ref_col)
               !$omp target exit data map(delete: vref_col)
            end associate
         end if

         if (allocated(state%multilayer%kpp_bl_depth)) then
            associate (kpp_bl => state%multilayer%kpp_bl_depth, &
                       kpp_us => state%multilayer%kpp_u_star, &
                       kpp_bf => state%multilayer%kpp_buoy_flux)
               !$omp target exit data map(from: kpp_bl, kpp_us, kpp_bf)
            end associate
         end if

         if (allocated(state%multilayer%keps_tke)) then
            associate (ke_tke => state%multilayer%keps_tke, &
                       ke_eps => state%multilayer%keps_eps, &
                       ke_nut => state%multilayer%keps_nut, &
                       ke_kt => state%multilayer%keps_kt, &
                       ke_usb => state%multilayer%keps_ustar_b, &
                       ke_uss => state%multilayer%keps_ustar_s)
               ! copyout, not delete: k and epsilon are prognostic and nu_t /
               ! kappa_t carry the previous step's value, so they must survive
               ! the exit/enter round trip around a snapshot or restart write.
               !$omp target exit data map(from: ke_tke, ke_eps, ke_nut, ke_kt, ke_usb, ke_uss)
            end associate
         end if

         if (allocated(state%multilayer%bed_sediment)) then
            associate (bed_sed => state%multilayer%bed_sediment)
               !$omp target exit data map(from: bed_sed)
            end associate
         end if

         ! Generic tracer registry: per-element components first, then
         ! the array descriptor (mirror of the enter_data order).
         if (allocated(state%multilayer%tracers)) then
            call solver_tracers_exit_data(state)
         end if
      end if

      ! Persistent multilayer kernel workspaces (PPM, hdiff, Smagorinsky,
      ! vcoord remap), lazy-allocated on first call; released here so
      ! re-running at a different grid size in the same process resizes.
      call tracer_scratch_destroy(tracer_scratch)
      call ml_horizontal_viscosity_cleanup()
      call ml_dynamics_workspaces_cleanup()

      ! Non-hydrostatic arrays
      if (state%use_nonhydrostatic) then
         !$omp target exit data map(from: state%nh%u_layer, state%nh%v_layer, state%nh%w, state%nh%q_nh)
         if (state%nh%nh_preconditioner) then
            !$omp target exit data map(delete: state%nh%cg_z)
         end if
         !$omp target exit data map(delete: state%nh%dsig, &
         !$omp&                 state%nh%u_layer0, state%nh%v_layer0, state%nh%w0, &
         !$omp&                 state%nh%cg_r, state%nh%cg_p, state%nh%cg_Ap, &
         !$omp&                 state%nh%hu_pre, state%nh%hv_pre)
      end if

      ! Mode splitting arrays
      if (allocated(state%nh%h_avg)) then
         !$omp target exit data map(delete: state%nh%h_avg, state%nh%hu_avg, state%nh%hv_avg, &
         !$omp&                 state%nh%h_bt0, state%nh%hu_bt0, state%nh%hv_bt0, &
         !$omp&                 state%nh%baro_tend_x, state%nh%baro_tend_y)
      end if

      ! Parent state struct last (mirror of enter_data order).
      !$omp target exit data map(delete: state)

   end subroutine solver_exit_data

   subroutine solver_update_host(state)
      !! Copy solution fields GPU->host without removing GPU data. Uses
      !! `!$acc update self(...)` so the device copy stays mapped — the
      !! lightweight sync for callbacks (no full exit/enter round-trip).
      type(state_t), intent(inout) :: state
         !! Solution state to sync

      integer :: it

      !$omp target update from(state%barotropic%h, state%barotropic%hu, state%barotropic%hv)

      ! Coupled hydrostatic multi-layer arrays
      if (state%use_multilayer) then
         !$omp target update from(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
         !$omp&                   state%multilayer%w_ml, state%multilayer%rho_layer)

         ! Generic tracer registry: pull each element's hTr back
         if (allocated(state%multilayer%tracers)) then
            do it = 1, size(state%multilayer%tracers)
               !$omp target update from(state%multilayer%tracers(it)%hTr)
            end do
         end if

         ! Bed sediment reservoir (2D), present only with sediment on.
         if (allocated(state%multilayer%bed_sediment)) then
            !$omp target update from(state%multilayer%bed_sediment)
         end if
      end if

      ! Non-hydrostatic arrays
      if (state%use_nonhydrostatic) then
         !$omp target update from(state%nh%u_layer, state%nh%v_layer, state%nh%w, state%nh%q_nh)
      end if

   end subroutine solver_update_host

   subroutine solver_update_forcing(state)
      !! Copy forcing fields host->GPU (in-place `!$acc update device`).
      !! Called after Python/C sets wind/pressure/precip on the host
      !! mid-simulation.
      type(state_t), intent(inout) :: state
         !! Solution state whose forcing arrays to sync

      !$omp target update to(state%forcing%wind_x, state%forcing%wind_y, &
      !$omp&                 state%forcing%atm_pressure, &
      !$omp&                 state%forcing%precipitation, &
      !$omp&                 state%forcing%heat_flux, &
      !$omp&                 state%forcing%evap_minus_precip)

   end subroutine solver_update_forcing

   subroutine solver_update_state(state)
      !! Copy solution fields host->GPU (in-place `!$acc update device`).
      !! Mirror of solver_update_host (pushes the same arrays it pulls);
      !! called after Python/C injects new state mid-simulation.
      type(state_t), intent(inout) :: state
         !! Solution state whose solution arrays to push to device

      integer :: it

      !$omp target update to(state%barotropic%h, state%barotropic%hu, state%barotropic%hv)

      ! Coupled hydrostatic multi-layer arrays
      if (state%use_multilayer) then
         !$omp target update to(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
         !$omp&                 state%multilayer%w_ml, state%multilayer%rho_layer)

         ! Generic tracer registry: push each element's hTr to device
         if (allocated(state%multilayer%tracers)) then
            do it = 1, size(state%multilayer%tracers)
               !$omp target update to(state%multilayer%tracers(it)%hTr)
            end do
         end if
      end if

      ! Non-hydrostatic arrays
      if (state%use_nonhydrostatic) then
         !$omp target update to(state%nh%u_layer, state%nh%v_layer, state%nh%w, state%nh%q_nh)
      end if

   end subroutine solver_update_state

   subroutine solver_tracers_enter_data(state)
      !! Map state%multilayer%tracers(:) into device memory. Two-step:
      !! copyin the array descriptor first (so the device sees the shape of
      !! `tracers(:)`), then per-element `hTr`/`hTr0`. Without the parent
      !! map, NVHPC stdpar can't dereference `tracers(it)%hTr` from a
      !! do-concurrent body (CUDA_ERROR_ILLEGAL_ADDRESS).
      type(state_t), intent(inout) :: state

      integer :: it

      !$omp target enter data map(to: state%multilayer%tracers)
      do it = 1, size(state%multilayer%tracers)
         !$omp target enter data map(to: state%multilayer%tracers(it)%hTr, &
         !$omp&                  state%multilayer%tracers(it)%hTr0)
      end do

      ! The device map above is an element-by-element snapshot: growing
      ! the registry after this point would strand every mapping.
      ! register_passive_tracer refuses while locked.
      state%multilayer%registry_locked = .true.
   end subroutine solver_tracers_enter_data

   subroutine solver_tracers_exit_data(state)
      !! Release the tracer registry from device memory (reverse of
      !! solver_tracers_enter_data: per-element first, then the array
      !! descriptor). Copies hTr back for diagnostics/tests; hTr0 (RK2
      !! workspace) is dropped.
      type(state_t), intent(inout) :: state

      integer :: it

      do it = 1, size(state%multilayer%tracers)
         !$omp target exit data map(from: state%multilayer%tracers(it)%hTr)
         !$omp target exit data map(delete: state%multilayer%tracers(it)%hTr0)
      end do
      !$omp target exit data map(delete: state%multilayer%tracers)

      state%multilayer%registry_locked = .false.
   end subroutine solver_tracers_exit_data

   subroutine accumulate_boundary_mass_out(state, stage_weight)
      !! Accumulate the net volume that left the domain this RK stage into
      !! `state%barotropic%mass_out`.  By the divergence theorem the boundary
      !! outflux equals `-Σ_interior flux_h·area` (interior face fluxes
      !! telescope), and `flux_h` is the SAME continuity tendency the update
      !! consumes (`h += dt·flux_h`), so with the matching RK stage weight
      !! this closes the mass budget to round-off — the console `Error` then
      !! measures only genuine non-conservation (e.g. dry-cell clamps).
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: stage_weight
         !! Stage contribution to the step's h update (0.5 for each SSP-RK2
         !! stage: h^{n+1} = h^n + 0.5·dt·(flux_h^1 + flux_h^2)).
      integer :: i, j, ng, i_lo, i_hi, j_lo, j_hi
      real(wp) :: acc, area

      ng = state%grid%nghost
      i_lo = ng + 1
      i_hi = state%grid%nx_total - ng
      j_lo = ng + 1
      j_hi = state%grid%ny_total - ng
      area = state%grid%dx*state%grid%dy
      acc = 0.0_wp
      !$omp target teams distribute parallel do collapse(2)
      do j = j_lo, j_hi
      do i = i_lo, i_hi
         acc = acc + state%barotropic%flux_h(i, j)
      end do
      end do
      !$omp end target teams distribute parallel do
      ! flux_h is d h/dt; Σ·area is the advective d(mass)/dt.  Outflux
      ! (positive = leaving) is its negative, times the stage's dt weight.
      state%barotropic%mass_out = state%barotropic%mass_out &
                                  - stage_weight*state%dt*acc*area
      state%barotropic%mass_out_tracked = .true.
   end subroutine accumulate_boundary_mass_out

   subroutine accumulate_mass_src(state, stage_weight, evap_mass_sink)
      !! Accumulate the net volume ADDED to the domain this RK stage by
      !! surface sources into `state%barotropic%mass_src`.
      !!
      !! Sources tracked:
      !!   - Precipitation: `h += dt·P` → contributes `+stage_weight·dt·Σ P·area`
      !!     (positive when P > 0, i.e. rain adds mass).
      !!   - Evaporative mass sink (when `evap_mass_sink=.true.`): the kernel
      !!     removes `dt·max(E,0)` metres from the column, so this contributes
      !!     `−stage_weight·dt·Σ max(E,0)·area` (negative; reduces mass_src).
      !!     When `evap_mass_sink=.false.` (default virtual-salt-flux mode) the
      !!     evaporation does NOT change `h`, so there is nothing to track here.
      !!
      !! Sign convention (mirrors `accumulate_boundary_mass_out`):
      !!   mass_src > 0  means sources have added net volume to the domain.
      !!
      !! Forcing off (`state%forcing%active = .false.`) → no precipitation,
      !! no evaporation → accumulator is a no-op and `mass_src` stays 0,
      !! which is correct (the budget still closes to round-off).
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: stage_weight
         !! Stage contribution (0.5 per SSP-RK2 stage).
      logical, intent(in) :: evap_mass_sink
         !! Whether the evaporative mass-sink mode is active (real water
         !! removal from `h`); when `.false.` the virtual-salt-flux mode
         !! leaves `h` unchanged, so evap contributes nothing here.

      integer :: i, j, ng, i_lo, i_hi, j_lo, j_hi
      real(wp) :: acc, area

      if (.not. state%forcing%active) then
         state%barotropic%mass_src_tracked = .true.
         return
      end if

      ng = state%grid%nghost
      i_lo = ng + 1
      i_hi = state%grid%nx_total - ng
      j_lo = ng + 1
      j_hi = state%grid%ny_total - ng
      area = state%grid%dx*state%grid%dy
      acc = 0.0_wp

      if (evap_mass_sink) then
         ! Net freshwater flux = precip − evap (evap_minus_precip stores E−P).
         ! Precip  (E−P < 0): adds mass → src += dt·|P|·area (P = -(E-P))
         ! Evap    (E−P > 0): removes mass → src −= dt·E·area
         ! Combined: net h increment from freshwater = dt·(−evap_minus_precip)
         ! BUT: the kernel only acts on E (removes min(dt·E, avail)); for the
         ! budget we track the INTENDED flux = dt·(P − E) = dt·(−(E−P)) per cell.
         ! Precipitation is always added in full (no availability limit) while
         ! evaporation is clamped; we track the intended field here (same as the
         ! conservative-rounding in the actual kernel for the precip leg).
         ! → acc = Σ_interior (−evap_minus_precip(i,j)) = Σ (P − E) per cell.
         !$omp target teams distribute parallel do collapse(2)
         do j = j_lo, j_hi
         do i = i_lo, i_hi
            acc = acc - state%forcing%evap_minus_precip(i, j)
         end do
         end do
         !$omp end target teams distribute parallel do
      else
         ! Virtual-salt-flux mode: only precipitation changes h.
         !$omp target teams distribute parallel do collapse(2)
         do j = j_lo, j_hi
         do i = i_lo, i_hi
            acc = acc + state%forcing%precipitation(i, j)
         end do
         end do
         !$omp end target teams distribute parallel do
      end if

      state%barotropic%mass_src = state%barotropic%mass_src &
                                  + stage_weight*state%dt*acc*area
      state%barotropic%mass_src_tracked = .true.
   end subroutine accumulate_mass_src

   subroutine accumulate_sponge_mass_src(state, stage_weight)
      !! Accumulate the net volume change caused by the sponge h-relaxation
      !! this RK stage into `state%barotropic%mass_src`.  The sponge applies
      !! `Δh = −dt·α·(h − h_ref)`, which can add or remove mass; this routine
      !! measures the actual h-change by summing h over interior cells before
      !! and after `apply_sponge`.  No-op when no sponge edge is configured.
      !!
      !! Must be called IN PLACE OF a bare `apply_sponge` call — it calls
      !! `apply_sponge` internally.
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: stage_weight
         !! Stage contribution (0.5 per SSP-RK2 stage).

      integer :: i, j, ng, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h_before, h_after, area

      ! Fast-exit when no sponge is configured (bc%sponge_width <= 0).
      if (state%bc%sponge_width <= 0) then
         call apply_sponge(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                           state%grid, state%dt, state%bc)
         return
      end if

      ng = state%grid%nghost
      i_lo = ng + 1
      i_hi = state%grid%nx_total - ng
      j_lo = ng + 1
      j_hi = state%grid%ny_total - ng
      area = state%grid%dx*state%grid%dy

      ! Sum h over interior before the sponge update.
      h_before = 0.0_wp
      !$omp target teams distribute parallel do collapse(2)
      do j = j_lo, j_hi
      do i = i_lo, i_hi
         h_before = h_before + state%barotropic%h(i, j)
      end do
      end do
      !$omp end target teams distribute parallel do

      call apply_sponge(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                        state%grid, state%dt, state%bc)

      ! Sum h after the sponge update.
      h_after = 0.0_wp
      !$omp target teams distribute parallel do collapse(2)
      do j = j_lo, j_hi
      do i = i_lo, i_hi
         h_after = h_after + state%barotropic%h(i, j)
      end do
      end do
      !$omp end target teams distribute parallel do

      ! Δh > 0: sponge added mass; Δh < 0: sponge removed mass.
      ! stage_weight accounts for the RK2 contribution (0.5 per stage).
      state%barotropic%mass_src = state%barotropic%mass_src &
                                  + stage_weight*(h_after - h_before)*area
   end subroutine accumulate_sponge_mass_src

   subroutine solver_step(state, dt_max, time)
      !! Advance the solution by one timestep using SSP-RK2
      !!
      !! 1. Apply BCs, compute timestep (CFL)
      !! 2. Save U^n
      !! 3. Stage 1: Euler step U* = U^n + dt * L(U^n)
      !! 4. Apply BCs to U*
      !! 5. Stage 2: Euler step U** = U* + dt * L(U*)
      !! 6. Average: U^{n+1} = 0.5 * (U^n + U**)
      type(state_t), intent(inout) :: state
         !! Solution state to advance
      real(wp), intent(in), optional :: dt_max
         !! Maximum timestep cap (e.g. t_end - t_current)
      real(wp), intent(in), optional :: time
         !! Current simulation time (s), needed for tidal BCs

      ! Set time in BC params for tidal/Chapman boundaries
      if (present(time)) then
         state%bc%time = time
      end if
      state%bc%dx = state%grid%dx
      state%bc%dy = state%grid%dy

      if (state%use_multilayer) then
         call solver_step_multilayer_impl(state, dt_max)
      else if (allocated(state%nh%h_avg) .and. state%use_nonhydrostatic) then
         call solver_step_split_impl(state, dt_max)
      else
         call solver_step_impl(state, dt_max)
      end if

   end subroutine solver_step

   subroutine solver_step_impl(state, dt_max)
      !! SSP-RK2 step for unsplit (2D or NH without mode splitting)
      type(state_t), intent(inout) :: state
      real(wp), intent(in), optional :: dt_max

      logical :: use_nh
      integer :: nz
      type(timestepper_t) :: ts

      use_nh = state%use_nonhydrostatic
      nz = state%nh%nz

      ! Build the explicit-integrator coefficient table for this step.
      ! Host-side scalars only (never device-mapped); rebuilding each step
      ! is a few branch tests.  NH stays on RK2 (configure guard forbids a
      ! non-RK2 scheme unless the run is pure 2D barotropic).
      ts = make_timestepper(state%time_scheme)

      ! ---- NH: compute depth-averaged momentum from layers ----
      if (use_nh) then
         call compute_depth_averaged_momentum(state%barotropic%h, state%nh%u_layer, &
                                              state%nh%v_layer, state%nh%dsig, &
                                              state%barotropic%hu, state%barotropic%hv, &
                                              state%grid, nz)
      end if

      ! ---- Compute timestep from U^n ----
      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)
      if (use_nh) then
         call apply_nh_boundary_conditions(state%nh%q_nh, state%nh%u_layer, &
                                           state%nh%v_layer, state%nh%w, &
                                           state%barotropic%h, state%barotropic%b, &
                                           state%grid, nz, state%bc)
      end if

      if (state%dt_fixed > 0.0_wp) then
         state%dt = state%dt_fixed
         if (present(dt_max)) state%dt = min(state%dt, dt_max)
      else if (mod(state%step_count, state%cfl_interval) == 0 .or. state%dt <= 0.0_wp) then
         call compute_timestep(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                               state%grid, state%cfl, state%dt)

         ! SSP coefficient: rk43 (C=2) is stable at 2x the forward-Euler
         ! CFL dt.  Exact x1.0 for euler/rk2/rk3, so default RK2 is
         ! bit-identical here.  Applied before dt_max so the cap still wins.
         state%dt = state%dt*ts%ssp_coeff

         if (present(dt_max)) then
            state%dt = min(state%dt, dt_max)
         end if

      end if
      state%step_count = state%step_count + 1

      ! Propagate dt to BC params for Chapman radiation
      state%bc%dt = state%dt

      ! 2D barotropic fused path: save + update + friction + coriolis in one
      ! kernel per RK stage. NH/multilayer use the unfused path below (layer
      ! operations sit between update and friction).
      if (.not. use_nh .and. .not. state%use_multilayer &
          .and. ts%scheme /= TIME_RK2) then
         ! Selectable explicit integrator (euler/rk3/rk43): generic
         ! N-stage advance.  RK2 falls through to the untouched fused
         ! fast path below (its bit-identity gate).
         call solver_advance_2d_generic(state, ts)

      else if (.not. use_nh .and. .not. state%use_multilayer) then

         ! ==== Fused RK2 Stage 1 ====
         call profiler_start("boundary")
         call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                        state%grid, state%bc)
         call profiler_stop("boundary")
         call profiler_start("flux")
         call solver_compute_barotropic_flux(state)
         call accumulate_boundary_mass_out(state, 0.5_wp)
         call profiler_stop("flux")
         call profiler_start("update")
         call fused_update_stage1(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                  state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                                  state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                                  state%grid, state%dt, &
                                  state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                                  state%forcing%atm_pressure, state%forcing%precipitation, &
                                  state%barotropic%coriolis_f, &
                                  state%forcing%active)
         ! Accumulate the precipitation mass source for this stage.  The 2D
         ! barotropic path has no evap mass sink (no multilayer tracers), so
         ! evap_mass_sink=.false.
         call accumulate_mass_src(state, 0.5_wp, evap_mass_sink=.false.)
         call profiler_stop("update")

         ! Sponge: measure the h-change so it enters the mass-source budget.
         call accumulate_sponge_mass_src(state, 0.5_wp)

         ! ==== Fused RK2 Stage 2 ====
         state%bc%time = state%bc%time + state%dt

         call profiler_start("boundary")
         call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                        state%grid, state%bc)
         call profiler_stop("boundary")
         call profiler_start("flux")
         call solver_compute_barotropic_flux(state)
         call accumulate_boundary_mass_out(state, 0.5_wp)
         call profiler_stop("flux")
         call profiler_start("rk2_average")
         call fused_update_stage2(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                  state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                                  state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                                  state%grid, state%dt, &
                                  state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                                  state%forcing%atm_pressure, state%forcing%precipitation, &
                                  state%barotropic%coriolis_f, &
                                  state%forcing%active)
         ! Stage 2 precipitation source (same sign, same weight).
         call accumulate_mass_src(state, 0.5_wp, evap_mass_sink=.false.)
         call profiler_stop("rk2_average")

         ! Sponge stage 2.
         call accumulate_sponge_mass_src(state, 0.5_wp)

      else
         ! Unfused path: NH and multilayer.

         ! ---- Save U^n for the averaging step ----
         call solver_save_state(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, state%grid)
         if (use_nh) then
            call solver_save_state_layers(state%nh%u_layer, state%nh%v_layer, state%nh%w, &
                                          state%nh%u_layer0, state%nh%v_layer0, state%nh%w0, &
                                          state%grid, nz)
         end if

         ! ==== RK2 Stage 1 ====
         call solver_nh_rk_stage(state, nz, use_nh)

         ! ---- Stage 2: U** = U* + dt * L(U*) ----
         ! Advance BC time by dt for tidal boundaries at stage 2
         state%bc%time = state%bc%time + state%dt

         ! NH: recompute depth-averaged momentum from updated layers
         if (use_nh) then
            call compute_depth_averaged_momentum(state%barotropic%h, state%nh%u_layer, &
                                                 state%nh%v_layer, state%nh%dsig, &
                                                 state%barotropic%hu, state%barotropic%hv, &
                                                 state%grid, nz)
         end if

         call solver_nh_rk_stage(state, nz, use_nh)

         ! ---- Average: U^{n+1} = 0.5 * (U^n + U**) ----
         call profiler_start("rk2_average")
         call compute_rk2_average(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                  state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                                  state%grid)
         if (use_nh) then
            call compute_rk2_average_layers(state%nh%u_layer, state%nh%v_layer, state%nh%w, &
                                            state%nh%u_layer0, state%nh%v_layer0, state%nh%w0, &
                                            state%grid, nz)
         end if
         call profiler_stop("rk2_average")

      end if  ! fused vs unfused path

   end subroutine solver_step_impl

   subroutine solver_advance_2d_generic(state, ts)
      !! Generic N-stage explicit SSP advance for the 2D barotropic path
      !! (euler / rk3 / rk43).  Loops the same forward-Euler stage
      !! primitive used by the fused RK2 fast path, then applies the
      !! Shu-Osher convex blend with the start-of-step state U^n:
      !!
      !!     E      = U^(i-1) + beta(i)*dt * L(U^(i-1))   (fused stage)
      !!     U^(i)  = alpha(i)*U^n + (1 - alpha(i))*E      (combine)
      !!
      !! `state%dt` already carries the ssp_coeff scaling.  No save/average
      !! kernels are reused from RK2 — U^n is backed up once, the blend is
      !! a per-stage compute_ssp_combine, both wet/dry-positive (SSP).
      type(state_t), intent(inout) :: state
      type(timestepper_t), intent(in) :: ts

      integer :: stage
      real(wp) :: t_base, dt_stage

      ! Back up U^n once; every stage blends against this start-of-step state.
      call solver_save_state(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                             state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                             state%grid)

      t_base = state%bc%time

      do stage = 1, ts%n_stages
         dt_stage = ts%beta(stage)*state%dt

         ! Time-dependent BCs see the stage abscissa t_base + c(i)*dt.
         state%bc%time = t_base + ts%c(stage)*state%dt

         call profiler_start("boundary")
         call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                        state%grid, state%bc)
         call profiler_stop("boundary")

         call profiler_start("flux")
         call solver_compute_barotropic_flux(state)
         call profiler_stop("flux")

         call profiler_start("update")
         call fused_update_stage_generic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                         state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                                         state%grid, dt_stage, &
                                         state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                                         state%forcing%atm_pressure, state%forcing%precipitation, &
                                         state%barotropic%coriolis_f, &
                                         state%forcing%active)

         ! Convex blend with U^n (skip when alpha == 0: the stage iterate
         ! IS the new running state, e.g. euler and the first rk43 stages).
         if (ts%alpha(stage) /= 0.0_wp) then
            call compute_ssp_combine(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                     state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                                     ts%alpha(stage), state%grid)
         end if
         call profiler_stop("update")

         call apply_sponge(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                           state%grid, dt_stage, state%bc)
      end do

      ! Restore the canonical end-of-step BC time (t_base + dt).
      state%bc%time = t_base + state%dt

   end subroutine solver_advance_2d_generic

   subroutine solver_nh_rk_stage(state, nz, use_nh)
      !! One RK stage: BCs, flux, update, friction, then NH correction
      type(state_t), intent(inout) :: state
      integer, intent(in) :: nz
      logical, intent(in) :: use_nh

      integer :: nh_iter

      ! Save barotropic state before HLL update
      if (use_nh) then
         call save_barotropic_state(state%barotropic%hu, state%barotropic%hv, &
                                    state%nh%hu_pre, state%nh%hv_pre, state%grid)
      end if

      call profiler_start("boundary")
      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)
      call profiler_stop("boundary")

      call profiler_start("flux")
      call solver_compute_barotropic_flux(state)
      call profiler_stop("flux")

      call profiler_start("update")
      call compute_update(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                          state%grid, state%dt)
      call profiler_stop("update")

      ! Distribute barotropic acceleration to layers
      if (use_nh) then
         call compute_layer_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                            state%nh%hu_pre, state%nh%hv_pre, &
                                            state%nh%u_layer, state%nh%v_layer, &
                                            state%nh%dsig, state%grid, nz, state%dt)

         ! Compute vertical velocity from continuity
         call compute_vertical_velocity(state%barotropic%h, state%nh%u_layer, state%nh%v_layer, &
                                        state%barotropic%b, state%nh%dsig, state%nh%w, &
                                        state%grid, nz)

         ! Solve NH pressure (CG) and apply correction
         call profiler_start("nh_pressure")
         if (mod(state%step_count, state%nh%nh_solve_interval) == 0) then
            call solve_nh_pressure(state%barotropic%h, state%nh%u_layer, state%nh%v_layer, &
                                   state%nh%w, state%nh%q_nh, state%nh%dsig, &
                                   state%nh%cg_r, state%nh%cg_p, state%nh%cg_Ap, state%nh%cg_z, &
                                   state%grid, nz, state%dt, &
                                   state%nh%nh_max_iterations, state%nh%nh_tolerance, nh_iter, &
                                   state%nh%nh_preconditioner)
         end if

         ! NH BCs for q_nh before applying correction
         call apply_nh_boundary_conditions(state%nh%q_nh, state%nh%u_layer, &
                                           state%nh%v_layer, state%nh%w, &
                                           state%barotropic%h, state%barotropic%b, &
                                           state%grid, nz, state%bc)

         call apply_nh_correction(state%barotropic%h, state%nh%u_layer, state%nh%v_layer, &
                                  state%nh%w, state%nh%q_nh, state%nh%dsig, &
                                  state%grid, nz, state%dt)
         call profiler_stop("nh_pressure")
      end if

      call profiler_start("friction")
      call apply_friction(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%grid, state%dt, &
                          state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                          state%forcing%atm_pressure, state%forcing%precipitation, forcing_active=state%forcing%active)

      call apply_coriolis(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%grid, state%dt, state%barotropic%coriolis_f)

      call profiler_stop("friction")

      ! Sponge layer relaxation (no-op when sponge_width <= 0)
      call apply_sponge(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                        state%grid, state%dt, state%bc)

   end subroutine solver_nh_rk_stage

   ! ================================================================
   ! Coupled hydrostatic multi-layer driver
   ! ================================================================

   subroutine solver_step_multilayer_impl(state, dt_max)
      !! SSP-RK2 step for coupled hydrostatic multi-layer equations
      type(state_t), intent(inout) :: state
      real(wp), intent(in), optional :: dt_max

      integer :: nz, it

      nz = state%multilayer%nz_ml

      ! ---- Compute barotropic quantities from layers ----
      ! Skipped at nz=1: the fused per-stage kernel already leaves h, hu, hv
      ! consistent with the layers, so these recomputes would be pure copies.
      if (nz > 1) then
         call ml_compute_total_depth(state%multilayer%h_layer, state%barotropic%h, &
                                     state%grid, nz)
         call ml_compute_barotropic_momentum(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                             state%multilayer%hv_layer, state%barotropic%hu, state%barotropic%hv, &
                                             state%grid, nz)
      end if

      ! ---- Apply barotropic BCs and distribute to layer ghost cells ----
      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)
      if (state%multilayer%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
         call ml_apply_layer_bcs_full(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                      state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                      state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref_col, &
                                      state%multilayer%vcoord%zstar_h_min, &
                                      state%grid, nz)
      else
         call ml_apply_layer_bcs(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                 state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%multilayer%dsig_ml, &
                                 state%grid, nz, &
                                 bc_west=state%bc%west, bc_east=state%bc%east, &
                                 bc_south=state%bc%south, bc_north=state%bc%north)
      end if

      ! ---- Compute timestep from barotropic state ----
      if (state%dt_fixed > 0.0_wp) then
         state%dt = state%dt_fixed
         if (present(dt_max)) state%dt = min(state%dt, dt_max)
      else if (mod(state%step_count, state%cfl_interval) == 0 .or. state%dt <= 0.0_wp) then
         call compute_timestep(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                               state%grid, state%cfl, state%dt)
         if (present(dt_max)) state%dt = min(state%dt, dt_max)
      end if
      state%step_count = state%step_count + 1
      state%bc%dt = state%dt

      ! ---- Save state for RK2 ----
      ! At nz=1 the fused stage-1 kernel absorbs `solver_save_state` +
      ! `ml_save_layer_state`; nz>1 uses the standalone kernels.
      if (nz > 1) then
         call solver_save_state(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, state%grid)
         call ml_save_layer_state(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                  state%multilayer%h_layer0, state%multilayer%hu_layer0, state%multilayer%hv_layer0, &
                                  state%grid, nz)
      end if
      do it = 1, size(state%multilayer%tracers)
         call ml_save_tracer(state%multilayer%tracers(it)%hTr, &
                             state%multilayer%tracers(it)%hTr0, state%grid, nz)
      end do

      ! ==== RK2 Stage 1 ====
      call solver_multilayer_rk_stage(state, nz, is_stage1=.true.)

      ! ---- Update barotropic for BCs at stage 2 ----
      ! Redundant at nz=1 (fused stage-1 already left h, hu, hv consistent).
      state%bc%time = state%bc%time + state%dt
      if (nz > 1) then
         call ml_compute_total_depth(state%multilayer%h_layer, state%barotropic%h, &
                                     state%grid, nz)
         call ml_compute_barotropic_momentum(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                             state%multilayer%hv_layer, state%barotropic%hu, state%barotropic%hv, &
                                             state%grid, nz)
      end if

      ! ==== RK2 Stage 2 ====
      call solver_multilayer_rk_stage(state, nz, is_stage1=.false.)

      ! ---- RK2 Average + recompute barotropic in one fused kernel ----
      ! One per-cell loop averages all layer arrays and accumulates h, hu, hv;
      ! replaces ml_rk2_average_layers + ml_rk2_average_tracer ×2 +
      ! ml_compute_total_depth + ml_compute_barotropic_momentum. Same path for
      ! nz=1 and nz>1 (no absorbed kernel has neighbour reads).
      call profiler_start("ml_work")
      call ml_fused_rk2_finalise(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                 state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                 state%multilayer%h_layer0, state%multilayer%hu_layer0, state%multilayer%hv_layer0, &
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
    state%multilayer%tracers(state%multilayer%idx_salinity)%hTr0, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr0, &
                                 state%grid, nz)
      ! Fused finalise RK2-averages S/T; extra passive tracers averaged here.
      do it = 1, size(state%multilayer%tracers)
         if (it == state%multilayer%idx_salinity .or. it == state%multilayer%idx_temperature) cycle
         call ml_rk2_average_tracer(state%multilayer%tracers(it)%hTr, &
                                    state%multilayer%tracers(it)%hTr0, state%grid, nz)
      end do
      call profiler_stop("ml_work")

      ! Sediment settling — applied ONCE per step on the RK2-averaged state
      ! (not per RK stage) so suspended loss = bed gain and total sediment is
      ! conserved. Constant w_s downward sink; bottom-layer flux deposits into
      ! the bed. Positivity-preserving (flux capped at layer mass), no clamp.
      if (state%multilayer%use_sediment .and. state%multilayer%idx_sediment > 0 .and. &
          state%multilayer%sediment_w_s > 0.0_wp) then
         call ml_settling_tracer(state%multilayer%h_layer, &
                                 state%multilayer%tracers(state%multilayer%idx_sediment)%hTr, &
                                 state%multilayer%bed_sediment, &
                                 state%multilayer%sediment_w_s, state%grid, nz, state%dt)
      end if

      ! Bed erosion — bed sediment lifted into the bottom layer where the
      ! (Manning-consistent) bed stress exceeds the critical erosion stress.
      ! After settling, on the same RK2-averaged state; conservative against
      ! the bed reservoir. Net bed change = deposition minus erosion.
      if (state%multilayer%use_sediment .and. state%multilayer%idx_sediment > 0 .and. &
          state%multilayer%sediment_erosion_rate > 0.0_wp .and. &
          state%multilayer%sediment_tau_crit_erosion > 0.0_wp) then
         call ml_erode_bed_sediment(state%multilayer%h_layer, &
                                    state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                    state%multilayer%tracers(state%multilayer%idx_sediment)%hTr, &
                                    state%multilayer%bed_sediment, &
                                    state%barotropic%manning_n, state%multilayer%rho_0, &
                                    state%multilayer%sediment_tau_crit_erosion, &
                                    state%multilayer%sediment_erosion_rate, &
                                    state%multilayer%sediment_porosity, &
                                    state%grid, nz, state%dt)
      end if

      ! Refresh EOS from the rk2-averaged tracers. Skipped at nz=1 (rho_layer
      ! has no consumer) and when alpha_T=beta_S=0 (density stays at rho_0)
      ! — unless the shared nonlinear EOS is selected (&tracer_nml eos_form),
      ! which doesn't read alpha_T/beta_S and must always refresh.
      if (nz > 1 .and. (state%multilayer%eos_use_shared .or. &
                        state%multilayer%alpha_T /= 0.0_wp .or. state%multilayer%beta_S /= 0.0_wp)) then
         call profiler_start("ml_eos")
         if (state%multilayer%eos_use_shared) then
            call ml_equation_of_state_shared(state%multilayer%h_layer, &
                                             state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                             state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                             state%multilayer%rho_layer, state%multilayer%eos, &
                                             state%grid, nz)
         else
            call ml_equation_of_state(state%multilayer%h_layer, &
                                      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                      state%multilayer%rho_layer, state%multilayer%rho_0, &
                                      state%multilayer%beta_S, state%multilayer%S_ref, &
                                      state%multilayer%alpha_T, state%multilayer%T_ref, &
                                      state%grid, nz)
         end if
         call profiler_stop("ml_eos")
      end if

   end subroutine solver_step_multilayer_impl

   subroutine solver_multilayer_rk_stage(state, nz, is_stage1)
      !! One RK stage: barotropic HLL (well-balanced, handles all pressure
      !! gradients), distribute acceleration to layers preserving shear,
      !! vertical exchange. `is_stage1` (nz=1 fused fast path only): also save
      !! pre-update state into the `*_0` arrays; ignored on the nz>1 path.
      type(state_t), intent(inout) :: state
      integer, intent(in) :: nz
      logical, intent(in) :: is_stage1

      integer :: it

      ! Mirror the S/T scalar config (held on the multilayer state + BC) onto
      ! their tracer-registry entries so the generic per-tracer kernels (both
      ! fast and nz>1 paths) read config uniformly. Extra passive tracers keep
      ! their own config. Cheap host-side scalar copies, idempotent; before
      ! the nz=1 dispatch so the fast path sees the synced config.
      it = state%multilayer%idx_salinity
      if (it > 0) then
         state%multilayer%tracers(it)%kappa_bg = state%multilayer%kappa_S_bg
         state%multilayer%tracers(it)%hdiff_kappa = state%multilayer%hdiff_kappa
         state%multilayer%tracers(it)%tr_min = state%multilayer%S_min
         state%multilayer%tracers(it)%tr_max = state%multilayer%S_max
         state%multilayer%tracers(it)%tr_inflow = state%bc%inflow_salinity
      end if
      it = state%multilayer%idx_temperature
      if (it > 0) then
         state%multilayer%tracers(it)%kappa_bg = state%multilayer%kappa_T_bg
         state%multilayer%tracers(it)%hdiff_kappa = state%multilayer%hdiff_kappa
         state%multilayer%tracers(it)%tr_min = state%multilayer%T_min
         state%multilayer%tracers(it)%tr_max = state%multilayer%T_max
         state%multilayer%tracers(it)%tr_inflow = state%bc%inflow_temperature
      end if

      ! At nz=1 the post-update sequence collapses into one fused per-cell
      ! kernel (+ fused tracer surface-flux/clamp): the depth-averaged-with-
      ! tracers fast path. nz>1 falls through to the un-fused logic below.
      if (nz == 1) then
         call solver_multilayer_rk_stage_nz1(state, is_stage1)
         return
      end if

      ! Apply barotropic BCs
      call profiler_start("boundary")
      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)
      if (state%multilayer%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
         call ml_apply_layer_bcs_full(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                      state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                      state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref_col, &
                                      state%multilayer%vcoord%zstar_h_min, &
                                      state%grid, nz)
      else
         call ml_apply_layer_bcs(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                 state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%multilayer%dsig_ml, &
                                 state%grid, nz, &
                                 bc_west=state%bc%west, bc_east=state%bc%east, &
                                 bc_south=state%bc%south, bc_north=state%bc%north)
      end if
      call profiler_stop("boundary")

      ! Save barotropic momentum before HLL update
      call solver_copy_2d(state%multilayer%hu_ml_pre, state%barotropic%hu, state%grid)
      call solver_copy_2d(state%multilayer%hv_ml_pre, state%barotropic%hv, state%grid)

      ! Barotropic HLL flux and update (well-balanced)
      call profiler_start("flux")
      call solver_compute_barotropic_flux(state)
      ! Mass budget: the layers are kept consistent with the barotropic h, and
      ! the RK2 finalise recomputes h from the averaged layers, so the total
      ! mass change is 0.5·dt·(flux_h^1 + flux_h^2) — the same telescoping as
      ! the 2D path.  Fires for both stages (this routine is called twice).
      call accumulate_boundary_mass_out(state, 0.5_wp)
      call profiler_stop("flux")

      call profiler_start("update")
      call compute_update(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                          state%grid, state%dt)
      call profiler_stop("update")

      ! Tracer step ordering (mirrors unstructured solver):
      !   1. Tracer BCs on OLD h_layer (before redistribution).
      !   2. CWC tracer advection: f_layer built from the HLL mass flux so
      !      `sum_k f_layer = mass_flux_h` (the barotropic continuity
      !      divergence) — uniform tracers exactly invariant, mass conserved.
      !   3. Distribute barotropic update to layers (does NOT touch tracer
      !      mass; step 2 already produced the right hTr for the new h_layer).
      ! For zsigma (`vcoord%needs_remap()`) the remap path reorders internally
      ! and conservatively remaps tracers between old and new thicknesses.

      ! Tracer BCs (operate on OLD h_layer) — loop the registry.
      do it = 1, size(state%multilayer%tracers)
         call ml_apply_tracer_bcs(state%multilayer%tracers(it)%hTr, state%multilayer%h_layer, &
                                  state%multilayer%tracers(it)%tr_inflow, &
                                  state%grid, nz, state%bc)
      end do

      ! CWC tracer advection on OLD h_layer with NEW mass_flux (same
      ! ordering as unstructured solver_unstr_step).
      do it = 1, size(state%multilayer%tracers)
         if (.not. state%multilayer%tracers(it)%do_horizontal_advection) cycle
         if (state%multilayer%use_ppm_tracer) then
            if (state%multilayer%tracer_recon_code /= TRACER_RECON_PPM) then
               call ml_advect_tracer_weno( &
                  state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                  state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                  state%multilayer%tracers(it)%hTr, tracer_scratch%weno_tr_snap, state%grid, nz, &
                  state%dt, state%multilayer%tracer_recon_code)
            else
               call ml_advect_tracer_ppm(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                         state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                                         state%multilayer%tracers(it)%hTr, &
                                         tracer_scratch%ppm_tr_snap, tracer_scratch%ppm_ql_x, tracer_scratch%ppm_qr_x, &
                                         tracer_scratch%ppm_ql_y, tracer_scratch%ppm_qr_y, &
                                         state%grid, nz, &
                                         state%dt)
            end if
         else
            call ml_advect_tracer(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                  state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                                  state%multilayer%tracers(it)%hTr, tracer_scratch%adv_hTr_snap, state%grid, nz, &
                                  state%dt)
         end if
      end do

      ! Distribute barotropic to layers (sigma) or remap to the zsigma target
      ! layout (the remap path conservatively projects the post-advection
      ! tracer profile onto the new vcoord-prescribed thicknesses per cell).
      if (state%multilayer%vcoord%needs_remap()) then
         if (state%multilayer%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
            call ml_distribute_with_remap_full(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                               state%multilayer%hu_ml_pre, state%multilayer%hv_ml_pre, &
                                               state%multilayer%vcoord_z_ref_col, &
                                               state%multilayer%vcoord%zstar_h_min, &
                                               state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                               state%grid, nz)
         else
            call ml_distribute_with_remap(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%hu_ml_pre, state%multilayer%hv_ml_pre, &
                                          state%multilayer%vcoord%coord_type, &
                                          state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref, &
                                          state%multilayer%vcoord%depth_transition, &
                                          state%multilayer%vcoord%blend_width, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                          state%grid, nz)
         end if
         ! Conservative per-tracer remap (saved remap_dz_old -> new h_layer),
         ! looped over the whole registry.
         do it = 1, size(state%multilayer%tracers)
            call ml_remap_tracer(state%multilayer%tracers(it)%hTr, state%barotropic%h, &
                                 state%multilayer%vcoord%remap_method, &
                                 state%multilayer%h_layer, state%grid, nz)
         end do
      else
         call ml_distribute_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                       state%multilayer%hu_ml_pre, state%multilayer%hv_ml_pre, &
                                       state%multilayer%dsig_ml, &
                                       state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                       state%grid, nz)
      end if

      ! Horizontal diffusion of tracers (smooths sharp density fronts) —
      ! per-tracer hdiff_kappa (synced from the state scalar for S/T).
      do it = 1, size(state%multilayer%tracers)
         if (state%multilayer%tracers(it)%do_horizontal_diffusion .and. &
             state%multilayer%tracers(it)%hdiff_kappa > 0.0_wp) then
            call ml_horizontal_diffusion_tracer(state%multilayer%h_layer, &
                                                state%multilayer%tracers(it)%hTr, &
                                                state%multilayer%tracers(it)%hdiff_kappa, &
                                                tracer_scratch%hdiff_hTr_pre, &
                                                state%grid, nz, &
                                                state%dt)
         end if
      end do

      ! Diagnostic vertical velocity and vertical exchange (momentum + tracers).
      ! Skipped at nz=1 (no interior interface; all w terms vanish).
      if (nz > 1) then
         call ml_compute_vertical_velocity(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                           state%multilayer%hv_layer, state%multilayer%dsig_ml, &
                                           state%multilayer%w_ml, &
                                           state%grid, nz)
         call ml_compute_vertical_exchange(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                           state%multilayer%hv_layer, state%multilayer%w_ml, &
                                           state%grid, nz, state%dt)
         do it = 1, size(state%multilayer%tracers)
            if (state%multilayer%tracers(it)%do_vertical_exchange) then
               call ml_vertical_exchange_tracer(state%multilayer%h_layer, &
                                                state%multilayer%tracers(it)%hTr, state%multilayer%w_ml, &
                                                state%grid, nz, state%dt)
            end if
         end do
      end if

      ! Equation of state + Ri-dependent vertical mixing. EOS feeds only BPG
      ! and the vertical-mixing kernels (both nz>1-guarded). At alpha_T=beta_S
      ! =0 EOS is a no-op; vmix still reads the init-populated rho_layer.
      if (nz > 1) then
         if (state%multilayer%eos_use_shared) then
            ! Shared nonlinear EOS (&tracer_nml eos_form) — refresh
            ! unconditionally (alpha_T/beta_S don't gate this path).
            call ml_equation_of_state_shared(state%multilayer%h_layer, &
                                             state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                             state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                             state%multilayer%rho_layer, state%multilayer%eos, &
                                             state%grid, nz)
         else if (state%multilayer%alpha_T /= 0.0_wp .or. state%multilayer%beta_S /= 0.0_wp) then
            call ml_equation_of_state(state%multilayer%h_layer, &
                                      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                      state%multilayer%rho_layer, state%multilayer%rho_0, &
                                      state%multilayer%beta_S, state%multilayer%S_ref, &
                                      state%multilayer%alpha_T, state%multilayer%T_ref, &
                                      state%grid, nz)
         end if
         if (state%multilayer%vmix_use_keps) then
            ! k-epsilon closure: advance k/eps, recover nu_t, apply to
            ! momentum via implicit vertical-diffusion sweep. Takes precedence
            ! over KPP / PP81.
            call ml_keps_friction_velocities(state%forcing%wind_x, state%forcing%wind_y, &
                                             state%multilayer%h_layer, state%multilayer%hu_layer, &
                                             state%multilayer%hv_layer, state%barotropic%manning_n, &
                                             state%barotropic%use_linear_bottom_drag, &
                                             state%barotropic%bottom_drag_r, state%multilayer%rho_0, &
                                             state%multilayer%keps_ustar_s, state%multilayer%keps_ustar_b, &
                                             state%grid, nz)
            call ml_keps_update(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                state%multilayer%hv_layer, state%multilayer%rho_layer, &
                                state%multilayer%keps_tke, state%multilayer%keps_eps, &
                                state%multilayer%keps_nut, state%multilayer%keps_kt, &
                                state%multilayer%rho_0, state%dt, state%grid, nz, &
                                state%multilayer%keps_ustar_b, state%multilayer%keps_ustar_s, &
                                scheme=state%multilayer%keps_stability)
            call ml_keps_apply_momentum(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                        state%multilayer%hv_layer, state%multilayer%keps_nut, &
                                        state%grid, nz, state%dt)
            ! Tracer vertical mixing with k-eps scalar diffusivity kt
            ! (+ per-tracer background floor), implicit + conservative.
            do it = 1, size(state%multilayer%tracers)
               if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                  call ml_keps_apply_diffusion_tracer(state%multilayer%h_layer, &
                                                      state%multilayer%tracers(it)%hTr, &
                                                      state%multilayer%keps_kt, &
                                                      state%multilayer%tracers(it)%kappa_bg, &
                                                      state%grid, nz, state%dt)
               end if
            end do
         else if (state%multilayer%vmix_use_kpp) then
            ! KPP overlay: one BL-depth sweep feeds the momentum + tracer apply
            ! kernels (max shape-function diffusivity onto PP81 + non-local
            ! transport). Replaces the three PP81 calls below.
            call ml_kpp_boundary_layer(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                       state%multilayer%hv_layer, state%multilayer%rho_layer, &
                                       state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                       state%multilayer%rho_0, state%multilayer%alpha_T, &
                                       state%multilayer%beta_S, CP_WATER, &
                                       state%forcing%wind_x, state%forcing%wind_y, &
                                       state%forcing%heat_flux, state%forcing%evap_minus_precip, &
                                       state%multilayer%kpp_ri_crit, state%multilayer%kpp_c_vt2, &
                                       state%multilayer%kpp_u_star, state%multilayer%kpp_buoy_flux, &
                                       state%multilayer%kpp_bl_depth, state%grid, nz)
            call ml_kpp_viscosity(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                  state%multilayer%hv_layer, state%multilayer%rho_layer, &
                                  state%multilayer%rho_0, state%multilayer%vmix_nu_max, &
                                  state%multilayer%vmix_nu_bg, state%multilayer%vmix_ri_alpha, &
                                  state%multilayer%kpp_u_star, state%multilayer%kpp_buoy_flux, &
                                  state%multilayer%kpp_bl_depth, state%grid, nz, state%dt)
            call ml_kpp_diffusion_tracers(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                          state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                          state%multilayer%rho_layer, state%multilayer%rho_0, &
                                          state%multilayer%vmix_nu_max, state%multilayer%kappa_S_bg, &
                                          state%multilayer%kappa_T_bg, state%multilayer%vmix_ri_alpha, &
                                          state%multilayer%kpp_cs_nonlocal, state%forcing%heat_flux, &
                                          state%forcing%evap_minus_precip, CP_WATER, &
                                          state%multilayer%kpp_u_star, state%multilayer%kpp_buoy_flux, &
                                          state%multilayer%kpp_bl_depth, state%grid, nz, state%dt)
            ! KPP fuses S/T only; extra passive tracers fall back to the
            ! PP81 interior closure for vertical diffusion.
            do it = 1, size(state%multilayer%tracers)
               if (it == state%multilayer%idx_salinity .or. it == state%multilayer%idx_temperature) cycle
               if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                  call ml_vertical_diffusion_tracer(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                                    state%multilayer%hv_layer, state%multilayer%tracers(it)%hTr, &
                                                    state%multilayer%rho_layer, state%multilayer%rho_0, &
                                                    state%multilayer%vmix_nu_max, state%multilayer%tracers(it)%kappa_bg, &
                                                    state%multilayer%vmix_ri_alpha, &
                                                    state%grid, nz, state%dt)
               end if
            end do
         else
            call ml_vertical_viscosity(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                       state%multilayer%rho_layer, state%multilayer%rho_0, &
                                       state%multilayer%vmix_nu_max, state%multilayer%vmix_nu_bg, &
                                       state%multilayer%vmix_ri_alpha, &
                                       state%grid, nz, state%dt)
            do it = 1, size(state%multilayer%tracers)
               if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                  call ml_vertical_diffusion_tracer(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                                    state%multilayer%hv_layer, state%multilayer%tracers(it)%hTr, &
                                                    state%multilayer%rho_layer, state%multilayer%rho_0, &
                                                    state%multilayer%vmix_nu_max, state%multilayer%tracers(it)%kappa_bg, &
                                                    state%multilayer%vmix_ri_alpha, &
                                                    state%grid, nz, state%dt)
               end if
            end do
         end if
      end if

      ! Surface tracer fluxes (net heat + virtual freshwater): top-layer BC
      ! for vertical diffusion; before the clamp so overshoot is caught.
      call ml_apply_surface_tracer_fluxes(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                          state%forcing%heat_flux, state%forcing%evap_minus_precip, &
                                          state%multilayer%evap_mass_sink, &
                                          state%multilayer%evap_latent_cooling, &
                                          state%multilayer%rho_0, CP_WATER, &
                                          state%grid, nz, state%dt)

      ! Enforce physical tracer bounds (catches mode-split overshoot) —
      ! per-tracer bounds (synced from the state scalars for S/T).
      do it = 1, size(state%multilayer%tracers)
         if (state%multilayer%tracers(it)%do_clamp) then
            call ml_clamp_tracer(state%multilayer%tracers(it)%hTr, state%multilayer%h_layer, &
                                 state%multilayer%tracers(it)%tr_min, state%multilayer%tracers(it)%tr_max, &
                                 state%grid, nz)
         end if
      end do

      ! BPG: skipped at nz=1 (single-layer pressure column is constant, so
      ! the horizontal gradient is zero).
      if (nz > 1) then
         if (state%multilayer%use_bpg_jacobian) then
            call ml_baroclinic_pressure_gradient_jacobian(state%multilayer%h_layer, state%multilayer%rho_layer, &
                                                          state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                                          state%multilayer%rho_0, state%multilayer%rho_ref, &
                                                          state%multilayer%rho_ref_dz, state%multilayer%nz_ref, &
                                                          state%grid, nz, state%dt)
         else
            call ml_baroclinic_pressure_gradient(state%multilayer%h_layer, state%multilayer%rho_layer, &
                                                 state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                                 state%multilayer%rho_0, state%multilayer%rho_ref, &
                                                 state%multilayer%rho_ref_dz, state%multilayer%nz_ref, &
                                                 state%grid, nz, state%dt)
         end if
      end if

      ! Per-layer bottom friction (damps baroclinic shear that depth-averaged friction misses)
      call ml_apply_bottom_friction(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                    state%barotropic%manning_n, state%barotropic%use_linear_bottom_drag, &
                                    state%barotropic%bottom_drag_r, state%grid, nz, state%dt)

      ! Smagorinsky horizontal viscosity (damps sharp velocity gradients at fronts)
      if (state%multilayer%smagorinsky_cs > 0.0_wp) then
         call ml_horizontal_viscosity_smag(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                           state%multilayer%smagorinsky_cs, &
                                           state%grid, nz, state%dt)
      end if

      ! Biharmonic (del^4) lateral viscosity (scale-selective grid-noise damping)
      if (state%multilayer%biharmonic_nu4 > 0.0_wp) then
         call ml_horizontal_viscosity_biharmonic(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                                 state%multilayer%hv_layer, state%multilayer%biharmonic_nu4, &
                                                 state%grid, nz, state%dt)
      end if

      ! CFL-based velocity truncation (MOM6-style)
      call ml_truncate_velocity(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                state%grid, nz, state%dt, &
                                0.8_wp)

      ! Per-layer Coriolis: rotate each layer's (hu_k, hv_k) through f*dt
      ! before recomputing the depth-average — correctly rotates baroclinic
      ! shear (reduces to depth-averaged for uniform columns).
      call profiler_start("friction")
      call ml_apply_coriolis(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                             state%grid, nz, state%dt, state%barotropic%coriolis_f)
      call profiler_stop("friction")

      ! Surface-concentrated wind stress: apply to the top layer so it
      ! generates vertical shear (proper Ekman / mixed-layer dynamics),
      ! instead of the depth-uniform barotropic forcing below.  The recompute
      ! folds it into the barotropic; the wind is then suppressed in
      ! apply_friction so it is not applied twice.
      if (state%barotropic%wind_surface_layer .and. state%forcing%active) then
         call ml_apply_wind_to_surface_layer(state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                             state%multilayer%h_layer, &
                                             state%forcing%wind_x, state%forcing%wind_y, &
                                             state%grid, nz, state%dt)
      end if

      ! Recompute barotropic from layers (BPG + vertical exchange + Coriolis shift momentum)
      call ml_compute_barotropic_momentum(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, state%barotropic%hu, state%barotropic%hv, &
                                          state%grid, nz)

      ! Wind/pressure/precip on barotropic (friction + Coriolis already applied per-layer above)
      call profiler_start("friction")
      call apply_friction(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%grid, state%dt, &
                          state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                          state%forcing%atm_pressure, state%forcing%precipitation, &
                          forcing_active=state%forcing%active, &
                          include_wind=.not. state%barotropic%wind_surface_layer)
      call profiler_stop("friction")

      ! Accumulate precipitation + evap mass-source for this RK stage.
      ! apply_friction adds `h += dt·precip` (barotropic); ml_apply_surface_tracer_fluxes
      ! adds the evap mass sink on the layer thickness.  Both contribute to mass_src.
      call accumulate_mass_src(state, 0.5_wp, &
                               evap_mass_sink=state%multilayer%evap_mass_sink)

      ! Sponge: accumulate the h-change from sponge relaxation in the budget.
      call accumulate_sponge_mass_src(state, 0.5_wp)

   end subroutine solver_multilayer_rk_stage

   subroutine solver_multilayer_rk_stage_nz1(state, is_stage1)
      !! Fused RK-stage path specialised for nz_ml=1 (depth-averaged
      !! 2D with optional tracers).
      !!
      !! Reorders the tracer pipeline to run before the per-cell post-
      !! update sequence so that (compute_update + distribute_barotropic
      !! + bottom_friction + truncate_velocity + recompute_barotropic +
      !! wind / atm_pressure / precip + coriolis) collapses into a
      !! single fused kernel `ml_fused_post_update`. Mirrors what the
      !! 2D HLL path does with `fused_update_stage1/2`. At nz_ml=1 the
      !! vertical-mixing block (vert_exchange, EOS, vmix, BPG) is
      !! definitionally empty, so this path is complete; for nz_ml>1
      !! the un-fused stage above runs instead.
      !!
      !! When `is_stage1` is true the fused kernel additionally saves
      !! the pre-update barotropic and layer state into the `*_0`
      !! arrays, absorbing `solver_save_state` and `ml_save_layer_state`.
      type(state_t), intent(inout) :: state
      logical, intent(in) :: is_stage1
      integer, parameter :: nz_local = 1
      integer :: it

      ! ---- Barotropic + layer BCs ----
      call profiler_start("boundary")
      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)
      call ml_apply_layer_bcs(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                              state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%multilayer%dsig_ml, &
                              state%grid, nz_local, &
                              bc_west=state%bc%west, bc_east=state%bc%east, &
                              bc_south=state%bc%south, bc_north=state%bc%north)
      call profiler_stop("boundary")

      ! ---- Barotropic HLL flux (mass + momentum) ----
      ! The pre-update momentum needed by the fused-kernel's distribute
      ! step is captured in registers inside the kernel itself, so the
      ! un-fused path's hu_ml_pre / hv_ml_pre copies are skipped here.
      call profiler_start("flux")
      call solver_compute_barotropic_flux(state)
      call profiler_stop("flux")

      ! ---- Tracer BCs + CWC advection on OLD h_layer + mass_flux ----
      ! Reordered to run before the fused post-update kernel, since the
      ! fused kernel writes the NEW h_layer in place. The CWC property
      ! holds because the mass_flux is the same continuity flux that
      ! the barotropic update inside the fused kernel will consume.
      ! Tracer BCs on OLD h_layer, looped over the registry (per-tracer
      ! inflow synced for S/T in solver_multilayer_rk_stage above).
      do it = 1, size(state%multilayer%tracers)
         call ml_apply_tracer_bcs(state%multilayer%tracers(it)%hTr, state%multilayer%h_layer, &
                                  state%multilayer%tracers(it)%tr_inflow, &
                                  state%grid, nz_local, state%bc)
      end do

      ! CWC tracer advection on OLD h_layer with NEW mass_flux, looped
      ! over the whole registry (mirrors the nz>1 path).
      do it = 1, size(state%multilayer%tracers)
         if (.not. state%multilayer%tracers(it)%do_horizontal_advection) cycle
         if (state%multilayer%use_ppm_tracer) then
            if (state%multilayer%tracer_recon_code /= TRACER_RECON_PPM) then
               call ml_advect_tracer_weno( &
                  state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                  state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                  state%multilayer%tracers(it)%hTr, tracer_scratch%weno_tr_snap, state%grid, nz_local, &
                  state%dt, state%multilayer%tracer_recon_code)
            else
               call ml_advect_tracer_ppm(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                         state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                                         state%multilayer%tracers(it)%hTr, &
                                         tracer_scratch%ppm_tr_snap, tracer_scratch%ppm_ql_x, tracer_scratch%ppm_qr_x, &
                                         tracer_scratch%ppm_ql_y, tracer_scratch%ppm_qr_y, &
                                         state%grid, nz_local, &
                                         state%dt)
            end if
         else
            call ml_advect_tracer(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                  state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                                  state%multilayer%tracers(it)%hTr, tracer_scratch%adv_hTr_snap, state%grid, nz_local, &
                                  state%dt)
         end if
      end do

      ! ---- The big fusion: 7 per-cell kernels collapsed into one ----
      call profiler_start("update")
      call ml_fused_post_update(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                                state%multilayer%h_layer0, state%multilayer%hu_layer0, state%multilayer%hv_layer0, &
                                state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                                state%multilayer%dsig_ml, &
                                state%forcing%wind_x, state%forcing%wind_y, &
                                state%forcing%atm_pressure, state%forcing%precipitation, &
                                state%grid, nz_local, state%dt, &
                                state%barotropic%manning_n, &
                                state%barotropic%use_linear_bottom_drag, state%barotropic%bottom_drag_r, &
                                state%barotropic%coriolis_f, state%forcing%active, &
                                do_save=is_stage1)
      call profiler_stop("update")

      ! ---- Tracer pipeline tail: horizontal diffusion, surface flux,
      !      physical clamp. ----
      ! Horizontal tracer diffusion, looped over the registry with
      ! per-tracer hdiff_kappa (synced from the state scalar for S/T).
      do it = 1, size(state%multilayer%tracers)
         if (state%multilayer%tracers(it)%do_horizontal_diffusion .and. &
             state%multilayer%tracers(it)%hdiff_kappa > 0.0_wp) then
            call ml_horizontal_diffusion_tracer(state%multilayer%h_layer, state%multilayer%tracers(it)%hTr, &
                                                state%multilayer%tracers(it)%hdiff_kappa, &
                                                tracer_scratch%hdiff_hTr_pre, &
                                                state%grid, nz_local, &
                                                state%dt)
         end if
      end do

      ! Surface tracer flux + S/T clamp fused into one per-cell pass.
      ! Replaces 3 launches: ml_apply_surface_tracer_fluxes +
      ! ml_clamp_tracer(S) + ml_clamp_tracer(T).
      call ml_fused_tracer_post(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%h_layer, &
                                state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                state%forcing%heat_flux, state%forcing%evap_minus_precip, &
                                state%multilayer%evap_mass_sink, &
                                state%multilayer%evap_latent_cooling, &
                                state%multilayer%S_min, state%multilayer%S_max, &
                                state%multilayer%T_min, state%multilayer%T_max, &
                                state%multilayer%rho_0, CP_WATER, ML_MIN_DEPTH, &
                                state%grid, nz_local, state%dt)

      ! Clamp extra registered tracers — the fused tracer-post above clamps
      ! S/T only, so extras (e.g. sediment) need their own positivity clamp.
      do it = 1, size(state%multilayer%tracers)
         if (it == state%multilayer%idx_salinity .or. it == state%multilayer%idx_temperature) cycle
         if (state%multilayer%tracers(it)%do_clamp) then
            call ml_clamp_tracer(state%multilayer%tracers(it)%hTr, state%multilayer%h_layer, &
                                 state%multilayer%tracers(it)%tr_min, state%multilayer%tracers(it)%tr_max, &
                                 state%grid, nz_local)
         end if
      end do

      ! ---- Smagorinsky horizontal viscosity (optional) ----
      ! Smag modifies hu_layer/hv_layer; refresh barotropic hu/hv so
      ! the next stage's BCs see consistent state.
      if (state%multilayer%smagorinsky_cs > 0.0_wp) then
         call ml_horizontal_viscosity_smag(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer, &
                                           state%multilayer%smagorinsky_cs, &
                                           state%grid, nz_local, state%dt)
      end if
      if (state%multilayer%biharmonic_nu4 > 0.0_wp) then
         call ml_horizontal_viscosity_biharmonic(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                                 state%multilayer%hv_layer, state%multilayer%biharmonic_nu4, &
                                                 state%grid, nz_local, state%dt)
      end if
      if (state%multilayer%smagorinsky_cs > 0.0_wp .or. state%multilayer%biharmonic_nu4 > 0.0_wp) then
         call ml_compute_barotropic_momentum(state%multilayer%h_layer, state%multilayer%hu_layer, &
                                             state%multilayer%hv_layer, state%barotropic%hu, state%barotropic%hv, &
                                             state%grid, nz_local)
      end if

      ! ---- Sponge zone relaxation ----
      ! NOTE: the nz=1 ML fast path does NOT call accumulate_mass_src /
      ! accumulate_sponge_mass_src, so mass_src_tracked stays .false. and
      ! the console falls back to raw drift for this path.  This is a safe,
      ! conservative fallback: the nz=1 fused path is typically used for
      ! barotropic ML spin-up or validation, not the flagship forced runs
      ! (Sydney Harbour / coastal configs with precip/evap/sponge use nz>1).
      ! Instrumenting it requires threading the evap_mass_sink flag and the
      ! budget accumulators through the fused kernel; deferred until needed.
      call apply_sponge(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                        state%grid, state%dt, state%bc)
   end subroutine solver_multilayer_rk_stage_nz1

   ! ================================================================
   ! Barotropic-baroclinic mode splitting driver
   ! ================================================================

   subroutine solver_step_split_impl(state, dt_max)
      !! Split-step driver: N barotropic substeps + 1 baroclinic step
      type(state_t), intent(inout) :: state
      real(wp), intent(in), optional :: dt_max

      integer :: nz, n_split, n_sub
      real(wp) :: dt_fast, dt_slow, inv_n

      nz = state%nh%nz

      ! ---- Compute both timesteps ----
      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)
      call apply_nh_boundary_conditions(state%nh%q_nh, state%nh%u_layer, &
                                        state%nh%v_layer, state%nh%w, &
                                        state%barotropic%h, state%barotropic%b, &
                                        state%grid, nz, state%bc)

      call compute_depth_averaged_momentum(state%barotropic%h, state%nh%u_layer, &
                                           state%nh%v_layer, state%nh%dsig, &
                                           state%barotropic%hu, state%barotropic%hv, &
                                           state%grid, nz)

      if (state%dt_fixed > 0.0_wp) then
         dt_slow = state%dt_fixed
         dt_fast = state%dt_fixed
         if (present(dt_max)) then
            dt_slow = min(dt_slow, dt_max)
            dt_fast = min(dt_fast, dt_max)
         end if
      else if (mod(state%step_count, state%cfl_interval) == 0 .or. state%dt <= 0.0_wp) then
         call compute_timestep(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                               state%grid, state%cfl, dt_fast)
         call compute_timestep_baroclinic(state%barotropic%h, state%nh%u_layer, &
                                          state%nh%v_layer, state%nh%dsig, &
                                          state%grid, nz, state%cfl, dt_slow)

         if (present(dt_max)) then
            dt_slow = min(dt_slow, dt_max)
            dt_fast = min(dt_fast, dt_max)
         end if

      else
         dt_fast = state%dt
         dt_slow = state%dt
         if (present(dt_max)) then
            dt_slow = min(dt_slow, dt_max)
            dt_fast = min(dt_fast, dt_max)
         end if
      end if
      state%step_count = state%step_count + 1

      ! Compute split ratio
      if (state%nh%split_ratio > 1) then
         n_split = state%nh%split_ratio
      else
         n_split = max(1, ceiling(dt_slow/dt_fast))
      end if

      if (n_split > 1) then
         dt_fast = dt_slow/real(n_split, wp)
      else
         dt_slow = dt_fast
      end if

      state%dt = dt_slow
      state%bc%dt = dt_fast

      ! ---- Compute NH barotropic tendency ----
      call compute_nh_baro_tendency(state%barotropic%h, state%nh%q_nh, state%nh%dsig, &
                                    state%nh%baro_tend_x, state%nh%baro_tend_y, &
                                    state%grid, nz)

      ! ---- Save full state ----
      call solver_save_state(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                             state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, state%grid)
      call solver_save_state_layers(state%nh%u_layer, state%nh%v_layer, state%nh%w, &
                                    state%nh%u_layer0, state%nh%v_layer0, state%nh%w0, &
                                    state%grid, nz)

      ! ---- Zero accumulators ----
      call solver_zero_2d(state%nh%h_avg, state%grid)
      call solver_zero_2d(state%nh%hu_avg, state%grid)
      call solver_zero_2d(state%nh%hv_avg, state%grid)

      ! ============================================================
      ! BAROTROPIC SUBCYCLING
      ! ============================================================
      do n_sub = 1, n_split
         call solver_save_state(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                state%nh%h_bt0, state%nh%hu_bt0, state%nh%hv_bt0, state%grid)

         ! RK2 Stage 1
         call solver_barotropic_rk_stage(state, dt_fast)

         ! RK2 Stage 2
         state%bc%time = state%bc%time + dt_fast
         call solver_barotropic_rk_stage(state, dt_fast)

         ! RK2 average
         call compute_rk2_average(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                  state%nh%h_bt0, state%nh%hu_bt0, state%nh%hv_bt0, &
                                  state%grid)

         ! Accumulate
         call solver_accumulate_2d(state%nh%h_avg, state%barotropic%h, state%grid)
         call solver_accumulate_2d(state%nh%hu_avg, state%barotropic%hu, state%grid)
         call solver_accumulate_2d(state%nh%hv_avg, state%barotropic%hv, state%grid)
      end do

      ! Finalize time average
      if (state%nh%split_averaging .and. n_split > 1) then
         inv_n = 1.0_wp/real(n_split, wp)
         call solver_scale_2d(state%nh%h_avg, inv_n, state%grid)
         call solver_scale_2d(state%nh%hu_avg, inv_n, state%grid)
         call solver_scale_2d(state%nh%hv_avg, inv_n, state%grid)
      else
         call solver_copy_2d(state%nh%h_avg, state%barotropic%h, state%grid)
         call solver_copy_2d(state%nh%hu_avg, state%barotropic%hu, state%grid)
         call solver_copy_2d(state%nh%hv_avg, state%barotropic%hv, state%grid)
      end if

      ! ============================================================
      ! BAROCLINIC STEP
      ! ============================================================
      ! Restore layer state from beginning-of-step save (swap src/dst)
      call solver_save_state_layers(state%nh%u_layer0, state%nh%v_layer0, state%nh%w0, &
                                    state%nh%u_layer, state%nh%v_layer, state%nh%w, &
                                    state%grid, nz)

      ! Single baroclinic fractional step at dt_slow
      call solver_baroclinic_rk_stage(state, nz, dt_slow)

      ! ============================================================
      ! RECONCILIATION
      ! ============================================================
      call solver_copy_2d(state%barotropic%h, state%nh%h_avg, state%grid)
      call reconcile_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                state%nh%u_layer, state%nh%v_layer, state%nh%dsig, &
                                state%grid, nz)

   end subroutine solver_step_split_impl

   subroutine solver_barotropic_rk_stage(state, dt_fast)
      !! One RK stage of barotropic subcycling: BCs, flux, update, NH tendency, friction, sponge
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: dt_fast

      call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                     state%grid, state%bc)

      call solver_compute_barotropic_flux(state)

      call compute_update(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                          state%grid, dt_fast)

      call apply_baro_nh_tendency(state%barotropic%hu, state%barotropic%hv, state%barotropic%h, &
                                  state%nh%baro_tend_x, state%nh%baro_tend_y, &
                                  state%grid, dt_fast)

      call apply_friction(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%grid, dt_fast, &
                          state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                          state%forcing%atm_pressure, state%forcing%precipitation, forcing_active=state%forcing%active)

      call apply_coriolis(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                          state%grid, dt_fast, state%barotropic%coriolis_f)

      call apply_sponge(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                        state%grid, dt_fast, state%bc)

   end subroutine solver_barotropic_rk_stage

   subroutine solver_baroclinic_rk_stage(state, nz, dt_slow)
      !! One RK stage of baroclinic step: layer distribution, w, CG Poisson, NH correction
      type(state_t), intent(inout) :: state
      integer, intent(in) :: nz
      real(wp), intent(in) :: dt_slow

      integer :: nh_iter

      ! hu0/hv0 hold the pre-step barotropic momentum
      call compute_layer_from_barotropic(state%nh%h_avg, state%nh%hu_avg, state%nh%hv_avg, &
                                         state%barotropic%hu0, state%barotropic%hv0, &
                                         state%nh%u_layer, state%nh%v_layer, &
                                         state%nh%dsig, state%grid, nz, dt_slow)

      call compute_vertical_velocity(state%nh%h_avg, state%nh%u_layer, state%nh%v_layer, &
                                     state%barotropic%b, state%nh%dsig, state%nh%w, &
                                     state%grid, nz)

      if (mod(state%step_count, state%nh%nh_solve_interval) == 0) then
         call solve_nh_pressure(state%nh%h_avg, state%nh%u_layer, state%nh%v_layer, &
                                state%nh%w, state%nh%q_nh, state%nh%dsig, &
                                state%nh%cg_r, state%nh%cg_p, state%nh%cg_Ap, state%nh%cg_z, &
                                state%grid, nz, dt_slow, &
                                state%nh%nh_max_iterations, state%nh%nh_tolerance, nh_iter, &
                                state%nh%nh_preconditioner)
      end if

      call apply_nh_boundary_conditions(state%nh%q_nh, state%nh%u_layer, &
                                        state%nh%v_layer, state%nh%w, &
                                        state%nh%h_avg, state%barotropic%b, &
                                        state%grid, nz, state%bc)

      call apply_nh_correction(state%nh%h_avg, state%nh%u_layer, state%nh%v_layer, &
                               state%nh%w, state%nh%q_nh, state%nh%dsig, &
                               state%grid, nz, dt_slow)

   end subroutine solver_baroclinic_rk_stage

   ! ================================================================
   ! Helper routines
   ! ================================================================

   subroutine solver_save_state(h, hu, hv, h0, hu0, hv0, grid)
      !! Copy current 2D state into saved arrays
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: h0(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hu0(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hv0(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp target teams distribute parallel do collapse(2)
      do j = 1, ny
      do i = 1, nx
         h0(i, j) = h(i, j)
         hu0(i, j) = hu(i, j)
         hv0(i, j) = hv(i, j)
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine solver_save_state

   subroutine solver_save_state_layers(u_layer, v_layer, w, &
                                            u_layer0, v_layer0, w0, &
                                            grid, nz)
      !! Copy layer state: first 3 args (source) -> last 3 args (destination)
      !! Used for both save (layers -> layers0) and restore (layers0 -> layers)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(out) :: u_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: v_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: w0(grid%nx_total, grid%ny_total, nz + 1)

      integer :: i, j, k, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp target teams distribute parallel do collapse(3)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         u_layer0(i, j, k) = u_layer(i, j, k)
         v_layer0(i, j, k) = v_layer(i, j, k)
      end do
      end do
      end do
      !$omp end target teams distribute parallel do

      !$omp target teams distribute parallel do collapse(3)
      do k = 1, nz + 1
      do j = 1, ny
      do i = 1, nx
         w0(i, j, k) = w(i, j, k)
      end do
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine solver_save_state_layers

   subroutine solver_zero_2d(f, grid)
      !! Zero a 2D field
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: f(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp target teams distribute parallel do collapse(2)
      do j = 1, ny
      do i = 1, nx
         f(i, j) = 0.0_wp
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine solver_zero_2d

   subroutine solver_accumulate_2d(acc, f, grid)
      !! Accumulate: acc = acc + f
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: acc(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: f(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp target teams distribute parallel do collapse(2)
      do j = 1, ny
      do i = 1, nx
         acc(i, j) = acc(i, j) + f(i, j)
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine solver_accumulate_2d

   subroutine solver_scale_2d(f, sc, grid)
      !! Scale: f = f * sc
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: f(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: sc

      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp target teams distribute parallel do collapse(2)
      do j = 1, ny
      do i = 1, nx
         f(i, j) = f(i, j)*sc
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine solver_scale_2d

   subroutine solver_copy_2d(dst, src, grid)
      !! Copy: dst = src
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: dst(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: src(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp target teams distribute parallel do collapse(2)
      do j = 1, ny
      do i = 1, nx
         dst(i, j) = src(i, j)
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine solver_copy_2d

end module rki_solver
