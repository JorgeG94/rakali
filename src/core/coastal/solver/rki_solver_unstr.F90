!! SSP-RK2 driver for the unstructured central-upwind solver
module rki_solver_unstr
   !! Two-stage strong-stability-preserving Runge-Kutta time loop for
   !! the unstructured central-upwind solver. Based on the
   !! evolve step (`shallow_water_domain.py:2000-2100`):
   !!
   !!   backup Q^n
   !!   stage 1: extrapolate, flux, update(dt), protect          -> Q*
   !!   stage 2: extrapolate, flux, update(dt), protect          -> Q**
   !!   saxpy(0.5, 0.5)                                          -> Q^{n+1}
   !!   protect
   !!
   !! Two public entry points:
   !!
   !!   * `solver_unstr_step(state, mesh, dt)` -- one SSP-RK2 step at
   !!     a caller-supplied `dt`. The driver does NOT call the CFL
   !!     kernel here; this is the right granularity for the future
   !!     hybrid driver, which will own the global timestep across
   !!     structured + unstructured halves and pass it in.
   !!
   !!   * `solver_unstr_advance(state, mesh, t_end)` -- step until
   !!     `state%time >= t_end`, recomputing `dt` from CFL each step
   !!     (or every `cfl_interval` steps once that knob is plumbed
   !!     through). Convenience wrapper for tests and stand-alone runs.
   !!
   !! Both entry points assume the data has already been entered onto
   !! the device. The caller is responsible for the
   !! `mesh_unstr_enter_data` / `state_unstr_enter_data` /
   !! `..._exit_data` lifecycle. This matches the structured solver's
   !! convention (`solver_step` does not touch enter/exit_data; the
   !! API layer handles it once at simulation start and end).
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   use rki_kernel_extrapolate_unstr, only: extrapolate_second_order_edge_unstr
   use rki_kernel_flux_unstr, only: compute_fluxes_central_unstr
   use rki_kernel_update_unstr, only: backup_conserved_unstr, &
                                      apply_explicit_update_unstr, &
                                      saxpy_conserved_unstr
   use rki_kernel_protect_unstr, only: protect_unstr
   use rki_kernel_friction_unstr, only: apply_manning_friction_unstr
   use rki_kernel_sponge_unstr, only: apply_sponge_unstr
   use rki_kernel_timestep_unstr, only: compute_timestep_unstr
   use rki_kernel_fused_unstr, only: fused_update_stage1_unstr, &
                                     fused_update_stage2_unstr, &
                                     fused_save_state_unstr, &
                                     fused_update_stage_generic_unstr, &
                                     ssp_combine_unstr
   use rki_timestepper, only: timestepper_t, make_timestepper, TIME_RK2, &
                              TIME_SEMI_IMPLICIT
   use rki_si_unstr_step, only: si_unstr_step
   use rki_si_unstr_step_ml, only: si_unstr_step_ml
   use rki_ml_dynamics_unstr, only: ml_init_layers_unstr, &
                                    ml_save_barotropic_pre_unstr, &
                                    ml_distribute_barotropic_unstr, &
                                    ml_distribute_with_remap_unstr, &
                                    ml_distribute_with_remap_full_unstr, &
                                    ml_remap_tracer_unstr, &
                                    ml_compute_barotropic_momentum_unstr, &
                                    ml_save_layer_state_unstr, &
                                    ml_compute_vertical_velocity_unstr, &
                                    ml_vertical_exchange_unstr, &
                                    ml_build_vcoord_zstar_full_unstr
   use rki_ml_tracers_unstr, only: ml_init_tracer_unstr, &
                                   ml_seed_standard_tracers_unstr, &
                                   ml_save_tracer_unstr, &
                                   ml_rk2_average_tracer_unstr, &
                                   ml_distribute_tracer_unstr, &
                                   ml_clamp_tracer_unstr, &
                                   ml_advect_tracer_unstr, &
                                   ml_advect_tracer_muscl_unstr, &
                                   ml_apply_surface_tracer_fluxes_unstr, &
                                   ml_vertical_exchange_tracer_unstr, &
                                   ml_horizontal_diffusion_tracer_unstr, &
                                   ml_vertical_diffusion_tracer_unstr
   use rki_ml_eos_unstr, only: ml_equation_of_state_unstr, &
                               ml_compute_reference_density_unstr
   use rki_ml_pressure_unstr, only: ml_baroclinic_pressure_gradient_unstr, &
                                    ml_baroclinic_pressure_gradient_jacobian_unstr
   use rki_ml_vertical_viscosity_unstr, only: ml_vertical_viscosity_unstr
   use rki_ml_kpp_unstr, only: ml_kpp_boundary_layer_unstr, &
                               ml_kpp_viscosity_unstr, &
                               ml_kpp_diffusion_tracers_unstr
   use rki_ml_keps_unstr, only: ml_keps_update_unstr, &
                                ml_keps_apply_momentum_unstr, &
                                ml_keps_apply_diffusion_tracer_unstr, &
                                ml_keps_friction_velocities_unstr
   use rki_ml_horizontal_viscosity_unstr, only: ml_horizontal_viscosity_smag_unstr
   use rki_ml_friction_unstr, only: ml_apply_coriolis_unstr, &
                                    ml_apply_bottom_friction_unstr
   use rki_ml_truncate_unstr, only: ml_truncate_velocity_unstr
   implicit none
   private

   public :: solver_unstr_step
   public :: solver_unstr_advance
   public :: solver_unstr_update_host
   public :: solver_unstr_update_forcing
   public :: solver_unstr_init_multilayer

   real(wp), parameter :: SAXPY_HALF = 0.5_wp

contains

   subroutine solver_unstr_init_multilayer(state, mesh)
      !! Populate the layer fields from the current barotropic state
      !! and seed the passive tracers (S, T) from the scalar initial
      !! values carried on `state`.
      !!
      !! Must be called once after `state_unstr_enter_data` (so the
      !! barotropic IC is on the device) and before the first step.
      !!
      !! After EOS, builds the horizontally-uniform reference-density
      !! profile `rho_ref(z)` so the BPG kernels see the small density
      !! anomaly `rho - rho_ref(z)` rather than `rho - rho_0` (Phase B2).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: it

      ! For VCOORD_ZSTAR_FULL build the per-cell reference z-level
      ! pattern from the local bathymetry before laying out the
      ! initial layer profile.
      if (state%multilayer%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
         call ml_build_vcoord_zstar_full_unstr(state)
      end if
      call ml_init_layers_unstr(state)
      ! Seed salinity + temperature from their IC scalars (each dispatching
      ! uniform-vs-stratified; S_init_bottom > S_init_surface = salty bed).
      call ml_seed_standard_tracers_unstr(state)
      ! Seed any extra registered passive tracers from their tr_init so
      ! the registry "just works" — register_passive_tracer + a uniform IC
      ! without the caller touching hTr directly.
      do it = 1, size(state%multilayer%tracers)
         if (it == state%multilayer%idx_salinity .or. it == state%multilayer%idx_temperature) cycle
         call ml_init_tracer_unstr(state, state%multilayer%tracers(it)%hTr, state%multilayer%tracers(it)%tr_init)
      end do
      call ml_equation_of_state_unstr(state)
      call ml_compute_reference_density_unstr(state, mesh)
   end subroutine solver_unstr_init_multilayer

   subroutine solver_unstr_step(state, mesh, dt)
      !! Take one SSP-RK2 step at the caller-supplied timestep `dt`.
      !! Updates `(stage_cv, xmom_cv, ymom_cv)` from `t` to `t + dt`.
      !! Increments `state%time` by `dt` and `state%step_count` by 1.
      !!
      !! Multilayer mode (`state%use_multilayer`) shadows the barotropic
      !! stages: before each HLL update the barotropic momentum is
      !! snapshotted, and after each update the barotropic acceleration
      !! is distributed uniformly across layers. At the very end of the
      !! step the layers are RK2-averaged in their own right and the
      !! barotropic fields are resynced via `sum_k` so the invariant
      !! `sum_k hu_layer == xmom_cv` holds after every call.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      integer :: it

      ! Semi-implicit (Casulli) free-surface: 2D barotropic path only.  The
      ! step owns its own eta-solve + edge back-substitution + time/step
      ! bookkeeping; no SSP stage table (make_timestepper would fail loud).
      if (state%use_semi_implicit) then
         call solver_unstr_step_semi_implicit(state, mesh, dt)
         return
      end if

      ! Selectable explicit integrator (euler/rk3/rk43): 2D barotropic
      ! path only.  RK2 (default) and every multilayer run fall through to
      ! the historical two-stage fused pipeline below.
      if (.not. state%use_multilayer .and. state%time_scheme /= TIME_RK2) then
         call solver_unstr_step_generic(state, mesh, dt)
         state%time = state%time + dt
         state%step_count = state%step_count + 1
         return
      end if

      if (state%use_multilayer) then
         ! Sync the salinity/temperature scalar config (held on the
         ! multilayer state for backward compat) onto their tracer-registry
         ! entries so the generic per-tracer kernels below read per-tracer
         ! config uniformly.  Extra passive tracers carry their own config
         ! (set by register_passive_tracer).  Cheap host-side copies.
         it = state%multilayer%idx_salinity
         if (it > 0) then
            state%multilayer%tracers(it)%kappa_bg = state%multilayer%kappa_S_bg
            state%multilayer%tracers(it)%hdiff_kappa = state%multilayer%hdiff_kappa
            state%multilayer%tracers(it)%tr_min = state%multilayer%S_min
            state%multilayer%tracers(it)%tr_max = state%multilayer%S_max
            state%multilayer%tracers(it)%bc_ext = state%multilayer%bc_S_ext
         end if
         it = state%multilayer%idx_temperature
         if (it > 0) then
            state%multilayer%tracers(it)%kappa_bg = state%multilayer%kappa_T_bg
            state%multilayer%tracers(it)%hdiff_kappa = state%multilayer%hdiff_kappa
            state%multilayer%tracers(it)%tr_min = state%multilayer%T_min
            state%multilayer%tracers(it)%tr_max = state%multilayer%T_max
            state%multilayer%tracers(it)%bc_ext = state%multilayer%bc_T_ext
         end if

         call ml_save_layer_state_unstr(state)
         do it = 1, size(state%multilayer%tracers)
            call ml_save_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                      state%multilayer%tracers(it)%hTr0)
         end do
      end if

      ! ----- RK substage 1: save + Euler update + friction + protect -----
      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)
      if (state%use_multilayer) then
         call ml_save_barotropic_pre_unstr(state)
      end if
      call fused_update_stage1_unstr(state, mesh, dt)
      call apply_sponge_unstr(state, mesh, dt)
      if (state%use_multilayer) then
         ! Tracer advection first -- it reads the OLD `h_layer` to get
         ! the upwind concentration and the just-stashed `mass_flux_e`
         ! that is consistent with the barotropic continuity update we
         ! just applied. That CWC recipe keeps uniform tracers exactly
         ! uniform and tracer mass globally conserved.  Loop the whole
         ! tracer registry; each tracer carries its own inflow bc_ext.
         do it = 1, size(state%multilayer%tracers)
            if (.not. state%multilayer%tracers(it)%do_horizontal_advection) cycle
            if (state%multilayer%use_ppm_tracer) then
               call ml_advect_tracer_muscl_unstr(state, mesh, state%multilayer%tracers(it)%hTr, &
                                                 dt, state%multilayer%tracers(it)%bc_ext)
            else
               call ml_advect_tracer_unstr(state, mesh, state%multilayer%tracers(it)%hTr, &
                                           dt, state%multilayer%tracers(it)%bc_ext)
            end if
         end do
         ! Horizontal tracer diffusion (Laplacian, after advection) —
         ! per-tracer hdiff_kappa (synced from the state scalar for S/T).
         do it = 1, size(state%multilayer%tracers)
            if (state%multilayer%tracers(it)%do_horizontal_diffusion .and. &
                state%multilayer%tracers(it)%hdiff_kappa > 0.0_wp) then
               call ml_horizontal_diffusion_tracer_unstr(state, mesh, &
                                                         state%multilayer%tracers(it)%hTr, &
                                                         state%multilayer%tracers(it)%hdiff_kappa, dt)
            end if
         end do
         ! Distribute the barotropic update to layers (momentum + thickness).
         ! For non-sigma coordinates this reshapes to the vcoord target;
         ! the per-tracer conservative remap then runs over the registry.
         if (state%multilayer%vcoord%needs_remap()) then
            if (state%multilayer%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
               call ml_distribute_with_remap_full_unstr(state)
            else
               call ml_distribute_with_remap_unstr(state)
            end if
            do it = 1, size(state%multilayer%tracers)
               call ml_remap_tracer_unstr(state, state%multilayer%tracers(it)%hTr)
            end do
         else
            call ml_distribute_barotropic_unstr(state)
         end if
         ! Smagorinsky horizontal viscosity (after distribute, per-layer)
         if (state%multilayer%smagorinsky_cs > 0.0_wp) then
            call ml_horizontal_viscosity_smag_unstr(state, mesh, dt)
         end if
         ! Vertical velocity and exchange: only when baroclinic physics active.
         if (state%multilayer%beta_S > 0.0_wp .or. state%multilayer%alpha_T > 0.0_wp) then
            call ml_compute_vertical_velocity_unstr(state, mesh)
            call ml_vertical_exchange_unstr(state, dt)
            do it = 1, size(state%multilayer%tracers)
               if (state%multilayer%tracers(it)%do_vertical_exchange) then
                  call ml_vertical_exchange_tracer_unstr(state, state%multilayer%tracers(it)%hTr, dt)
               end if
            end do
            ! Vertical mixing: k-epsilon when enabled (precedence), else KPP
            ! overlay (BL-depth sweep + shape-function diffusivity + non-local
            ! transport), else PP81.  KPP fuses S/T; extra passive tracers fall
            ! back to the PP81 interior closure with their own kappa_bg.
            if (state%multilayer%vmix_use_keps) then
               ! k-epsilon: friction velocities -> advance k/eps + recover
               ! nu_t/kt -> implicit momentum + per-tracer vertical diffusion.
               call ml_keps_friction_velocities_unstr(state)
               call ml_keps_update_unstr(state, dt, .true.)
               call ml_keps_apply_momentum_unstr(state, dt)
               do it = 1, size(state%multilayer%tracers)
                  if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                     call ml_keps_apply_diffusion_tracer_unstr(state, &
                                                               state%multilayer%tracers(it)%hTr, &
                                                               state%multilayer%tracers(it)%kappa_bg, dt)
                  end if
               end do
            else if (state%multilayer%vmix_use_kpp) then
               call ml_kpp_boundary_layer_unstr(state, &
                                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr)
               call ml_kpp_viscosity_unstr(state, dt)
               call ml_kpp_diffusion_tracers_unstr(state, &
                                                   state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                                   state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                                   state%multilayer%kappa_S_bg, state%multilayer%kappa_T_bg, dt)
               do it = 1, size(state%multilayer%tracers)
                  if (it == state%multilayer%idx_salinity .or. it == state%multilayer%idx_temperature) cycle
                  if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                     call ml_vertical_diffusion_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                                             state%multilayer%tracers(it)%kappa_bg, dt)
                  end if
               end do
            else
               call ml_vertical_viscosity_unstr(state, dt)
               do it = 1, size(state%multilayer%tracers)
                  if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                     call ml_vertical_diffusion_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                                             state%multilayer%tracers(it)%kappa_bg, dt)
                  end if
               end do
            end if
         end if
         ! Surface tracer fluxes (net heat + virtual freshwater).  Runs every
         ! RK2 stage unconditionally to mirror the structured pipeline -- a
         ! passive temperature tracer (alpha_T=0) still gets warmed by Q_net.
         call ml_apply_surface_tracer_fluxes_unstr(state, &
                                                   state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                                   state%multilayer%tracers(state%multilayer%idx_temperature)%hTr)
         ! Per-layer bottom friction (Manning, on the lowest active layer):
         ! damps baroclinic shear the depth-averaged friction misses.
         call ml_apply_bottom_friction_unstr(state, dt)
         ! CFL-based velocity truncation (MOM6-style, cfl_trunc=0.8): clip
         ! runaway layer velocities so a few bad cells don't tank the step.
         call ml_truncate_velocity_unstr(state, mesh, dt, 0.8_wp)
         ! Per-layer Coriolis: rotate every layer's (hu_k, hv_k) before
         ! recomputing the depth-averaged so the barotropic mean reflects
         ! the rotated layer momenta.  Mirrors the structured slot
         ! (between truncate / final per-layer kernels and the recompute).
         call ml_apply_coriolis_unstr(state, dt)
         call ml_compute_barotropic_momentum_unstr(state)
      end if

      ! ----- RK substage 2: update + friction + RK2 average + protect -----
      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)
      if (state%use_multilayer) then
         call ml_save_barotropic_pre_unstr(state)
      end if
      call fused_update_stage2_unstr(state, mesh, dt)
      call apply_sponge_unstr(state, mesh, dt)
      if (state%use_multilayer) then
         do it = 1, size(state%multilayer%tracers)
            if (.not. state%multilayer%tracers(it)%do_horizontal_advection) cycle
            if (state%multilayer%use_ppm_tracer) then
               call ml_advect_tracer_muscl_unstr(state, mesh, state%multilayer%tracers(it)%hTr, &
                                                 dt, state%multilayer%tracers(it)%bc_ext)
            else
               call ml_advect_tracer_unstr(state, mesh, state%multilayer%tracers(it)%hTr, &
                                           dt, state%multilayer%tracers(it)%bc_ext)
            end if
         end do
         ! Horizontal tracer diffusion (stage 2) — per-tracer hdiff_kappa.
         do it = 1, size(state%multilayer%tracers)
            if (state%multilayer%tracers(it)%do_horizontal_diffusion .and. &
                state%multilayer%tracers(it)%hdiff_kappa > 0.0_wp) then
               call ml_horizontal_diffusion_tracer_unstr(state, mesh, &
                                                         state%multilayer%tracers(it)%hTr, &
                                                         state%multilayer%tracers(it)%hdiff_kappa, dt)
            end if
         end do
         ! `fused_update_stage2_unstr` already RK2-blended the
         ! barotropic state (stage_cv, xmom_cv, ymom_cv) against the
         ! `_backup` snapshots captured at step start, and
         ! the distribute propagates that blend straight into
         ! `h_layer`. Applying `ml_rk2_average_layers_unstr` here
         ! would blend a second time (effectively 0.25/0.75 weights
         ! instead of 0.5/0.5), breaking tracer-concentration
         ! consistency. So we only RK2 the tracers, which carry
         ! their OWN `_layer0` backups.
         if (state%multilayer%vcoord%needs_remap()) then
            if (state%multilayer%vcoord%coord_type == VCOORD_ZSTAR_FULL) then
               call ml_distribute_with_remap_full_unstr(state)
            else
               call ml_distribute_with_remap_unstr(state)
            end if
            do it = 1, size(state%multilayer%tracers)
               call ml_remap_tracer_unstr(state, state%multilayer%tracers(it)%hTr)
            end do
         else
            call ml_distribute_barotropic_unstr(state)
         end if
         ! Smagorinsky horizontal viscosity (stage 2)
         if (state%multilayer%smagorinsky_cs > 0.0_wp) then
            call ml_horizontal_viscosity_smag_unstr(state, mesh, dt)
         end if
         ! Vertical velocity and exchange (stage 2): only when baroclinic physics active.
         if (state%multilayer%beta_S > 0.0_wp .or. state%multilayer%alpha_T > 0.0_wp) then
            call ml_compute_vertical_velocity_unstr(state, mesh)
            call ml_vertical_exchange_unstr(state, dt)
            do it = 1, size(state%multilayer%tracers)
               if (state%multilayer%tracers(it)%do_vertical_exchange) then
                  call ml_vertical_exchange_tracer_unstr(state, state%multilayer%tracers(it)%hTr, dt)
               end if
            end do
            ! Vertical mixing (stage 2): k-epsilon when enabled (precedence),
            ! else KPP overlay, else PP81.  KPP fuses S/T; extra passive
            ! tracers fall back to PP81 vdiff.
            if (state%multilayer%vmix_use_keps) then
               ! k-epsilon: friction velocities -> advance k/eps + recover
               ! nu_t/kt -> implicit momentum + per-tracer vertical diffusion.
               call ml_keps_friction_velocities_unstr(state)
               call ml_keps_update_unstr(state, dt, .true.)
               call ml_keps_apply_momentum_unstr(state, dt)
               do it = 1, size(state%multilayer%tracers)
                  if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                     call ml_keps_apply_diffusion_tracer_unstr(state, &
                                                               state%multilayer%tracers(it)%hTr, &
                                                               state%multilayer%tracers(it)%kappa_bg, dt)
                  end if
               end do
            else if (state%multilayer%vmix_use_kpp) then
               call ml_kpp_boundary_layer_unstr(state, &
                                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr)
               call ml_kpp_viscosity_unstr(state, dt)
               call ml_kpp_diffusion_tracers_unstr(state, &
                                                   state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                                   state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                                   state%multilayer%kappa_S_bg, state%multilayer%kappa_T_bg, dt)
               do it = 1, size(state%multilayer%tracers)
                  if (it == state%multilayer%idx_salinity .or. it == state%multilayer%idx_temperature) cycle
                  if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                     call ml_vertical_diffusion_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                                             state%multilayer%tracers(it)%kappa_bg, dt)
                  end if
               end do
            else
               call ml_vertical_viscosity_unstr(state, dt)
               do it = 1, size(state%multilayer%tracers)
                  if (state%multilayer%tracers(it)%do_vertical_diffusion) then
                     call ml_vertical_diffusion_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                                             state%multilayer%tracers(it)%kappa_bg, dt)
                  end if
               end do
            end if
         end if
         ! Surface tracer fluxes (stage 2).  Must run BEFORE the RK2
         ! average so the flux contribution is time-integrated by the
         ! 0.5*(hTr0 + hTr) blend, matching the structured pipeline.
         call ml_apply_surface_tracer_fluxes_unstr(state, &
                                                   state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                                   state%multilayer%tracers(state%multilayer%idx_temperature)%hTr)
         ! RK2 average + conservative clamp over the whole registry —
         ! per-tracer bounds (synced from the state scalars for S/T).
         do it = 1, size(state%multilayer%tracers)
            call ml_rk2_average_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                             state%multilayer%tracers(it)%hTr0)
         end do
         do it = 1, size(state%multilayer%tracers)
            if (state%multilayer%tracers(it)%do_clamp) then
               call ml_clamp_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                          state%multilayer%tracers(it)%tr_min, &
                                          state%multilayer%tracers(it)%tr_max)
            end if
         end do
         call ml_equation_of_state_unstr(state)
         ! Baroclinic pressure gradient: density-driven horizontal acceleration.
         ! Uses rho_layer from EOS, modifies hu/hv_layer per layer.
         if (state%multilayer%use_bpg_jacobian) then
            call ml_baroclinic_pressure_gradient_jacobian_unstr(state, mesh, dt)
         else
            call ml_baroclinic_pressure_gradient_unstr(state, mesh, dt)
         end if
         ! Per-layer bottom friction + CFL velocity truncation (after BPG,
         ! before Coriolis) — mirrors the structured stage order.
         call ml_apply_bottom_friction_unstr(state, dt)
         call ml_truncate_velocity_unstr(state, mesh, dt, 0.8_wp)
         ! Per-layer Coriolis after BPG, before recompute (matches stage 1).
         call ml_apply_coriolis_unstr(state, dt)
         call ml_compute_barotropic_momentum_unstr(state)
      end if

      state%time = state%time + dt
      state%step_count = state%step_count + 1
   end subroutine solver_unstr_step

   subroutine solver_unstr_step_generic(state, mesh, dt)
      !! Generic N-stage explicit SSP advance for the 2D barotropic
      !! unstructured path (euler/rk3/rk43).  Loops the same KNP
      !! forward-Euler stage primitive used by the fused RK2 fast path:
      !!
      !!     E      = U^(i-1) + beta(i)*dt * L(U^(i-1))   (fused stage)
      !!     U^(i)  = alpha(i)*U^n + (1 - alpha(i))*E      (convex blend)
      !!
      !! U^n is backed up once; the blend is a per-stage ssp_combine_unstr.
      !! Both keep the protect (wet/dry positivity) so the SSP guarantee
      !! carries through.  `dt` already carries any ssp_coeff scaling.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      type(timestepper_t) :: ts
      integer :: stage
      real(wp) :: dt_stage

      ts = make_timestepper(state%time_scheme)

      ! Back up U^n once; every stage blends against this start-of-step state.
      call fused_save_state_unstr(state, mesh)

      do stage = 1, ts%n_stages
         dt_stage = ts%beta(stage)*dt

         call extrapolate_second_order_edge_unstr(state, mesh)
         call compute_fluxes_central_unstr(state, mesh)
         call fused_update_stage_generic_unstr(state, mesh, dt_stage)

         ! Convex blend with U^n (skip when alpha == 0: the stage iterate
         ! IS the new running state, e.g. euler and the first rk43 stages).
         if (ts%alpha(stage) /= 0.0_wp) then
            call ssp_combine_unstr(state, mesh, ts%alpha(stage))
         end if
      end do
   end subroutine solver_unstr_step_generic

   subroutine solver_unstr_step_semi_implicit(state, mesh, dt)
      !! One semi-implicit (Casulli) free-surface step on the triangulation.
      !! Hands the conserved variables to `si_unstr_step`, which owns the frozen-Hf predictor ->
      !! TPFA Helmholtz eta-solve -> edge back-substitution -> A-grid cell
      !! momentum reconstruction, plus the time/step bookkeeping.  Open-edge
      !! Dirichlet coupling (tidal / clamped) is read inside the step from
      !! `state%bc_type` / `state%bc_eta_ext` (the coefficient-space route);
      !! wall edges get Hf=0 (no-flux Neumann).  The caller sets `bc_eta_ext`
      !! to the prescribed surface (tests evaluate the tide at the
      !! theta-weighted time).  The step ends with barotropic stage_cv /
      !! xmom_cv / ymom_cv / height_cv valid — downstream unchanged.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      integer :: iters

      ! The ML/2D split MUST match state_unstr_enable_si's slot choice: an
      ! nz_ml = 1 "multilayer" state gets the 2D si slot (the column solve
      ! needs >= 2 layers), so dispatching on use_multilayer alone would
      ! run the ML step against unallocated column slots.
      if (state%use_multilayer .and. state%multilayer%nz_ml > 1) then
         ! P3.2: the 3D multilayer column-coupled SI step (per-column
         ! implicit vertical viscosity + Hhat/Gbar coupling into the SAME
         ! TPFA eta-solve + dsig-proportional registry tracer advection).
         call si_unstr_step_ml(state, mesh, dt, iters)
      else
         call si_unstr_step(state, mesh, dt, iters)
      end if
   end subroutine solver_unstr_step_semi_implicit

   subroutine solver_unstr_advance(state, mesh, t_end)
      !! Step until `state%time >= t_end`. Recomputes the CFL-limited
      !! timestep from the per-cell wave speeds before each step using
      !! `compute_timestep_unstr`.
      !!
      !! NOTE: the very first step needs `state%barotropic%max_speed` to be
      !! populated, which only happens after a flux pass. We bootstrap
      !! by running extrapolate + flux once before the loop starts so
      !! the CFL kernel has something to read.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: t_end

      real(wp) :: dt, ssp_coeff
      type(timestepper_t) :: ts

      ! Semi-implicit mode: no SSP stage table (make_timestepper is
      ! fail-loud on SEMI_IMPLICIT) and no gravity-wave CFL — the eta-solve
      ! is unconditionally stable in the surface-wave mode and runs at the
      ! configured dt_fixed (stored by state_unstr_enable_si).  The
      ! explicit bootstrap flux pass is skipped too: max_speed feeds only
      ! the CFL kernel.
      if (state%use_semi_implicit) then
         state%dt = state%dt_fixed
         do while (state%time < t_end)
            dt = state%dt
            ! Cap the final step so we land exactly on t_end.
            if (state%time + dt > t_end) dt = t_end - state%time
            if (dt <= 0.0_wp) exit
            call solver_unstr_step(state, mesh, dt)
         end do
         return
      end if

      ! ssp_coeff scales the CFL-limited dt for schemes stable beyond the
      ! forward-Euler limit (rk43, C=2 -> 2x dt).  Exact x1.0 for
      ! euler/rk2/rk3, so the default RK2 dt is bit-identical.
      ts = make_timestepper(state%time_scheme)
      ssp_coeff = ts%ssp_coeff

      ! Bootstrap: populate state%barotropic%max_speed for the first CFL call.
      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)
      call compute_timestep_unstr(state, mesh)
      state%dt = state%dt*ssp_coeff

      do while (state%time < t_end)
         dt = state%dt
         ! Cap the final step so we land exactly on t_end.
         if (state%time + dt > t_end) dt = t_end - state%time
         if (dt <= 0.0_wp) exit

         call solver_unstr_step(state, mesh, dt)

         ! Recompute the CFL timestep using the wave speeds the flux
         ! kernel just wrote.
         call compute_timestep_unstr(state, mesh)
         state%dt = state%dt*ssp_coeff
      end do
   end subroutine solver_unstr_advance

   subroutine solver_unstr_update_host(state)
      !! Copy solution fields from GPU to host without tearing down the
      !! device mapping.
      !!
      !! Uses `!$acc update self(...)` so the device copy remains
      !! mapped and valid. This is the lightweight sync for output
      !! callbacks and Python getters — avoids the full exit/enter
      !! round-trip.
      type(unstr_state_t), intent(inout) :: state

      integer :: it

      !$acc update self(state%barotropic%stage_cv, state%barotropic%xmom_cv, state%barotropic%ymom_cv, &
      !$acc&            state%barotropic%height_cv)

      if (state%use_multilayer) then
         !$acc update self(state%multilayer%h_layer, state%multilayer%hu_layer, state%multilayer%hv_layer)

         if (allocated(state%multilayer%tracers)) then
            do it = 1, size(state%multilayer%tracers)
               !$acc update self(state%multilayer%tracers(it)%hTr)
            end do
         end if
      end if

   end subroutine solver_unstr_update_host

   subroutine solver_unstr_update_forcing(state)
      !! Refresh the device copy of the surface-forcing arrays after
      !! the host has been updated by the C API setters.
      !!
      !! Mirror of `solver_update_forcing` on the structured side: uses
      !! `!$acc update device(...)` so the device mapping stays alive
      !! and the next step picks up the new values without an
      !! exit/enter round-trip.
      type(unstr_state_t), intent(inout) :: state

      !$acc update device(state%forcing%heat_flux, &
      !$acc&              state%forcing%evap_minus_precip, &
      !$acc&              state%forcing%wind_x, state%forcing%wind_y)
   end subroutine solver_unstr_update_forcing

end module rki_solver_unstr
