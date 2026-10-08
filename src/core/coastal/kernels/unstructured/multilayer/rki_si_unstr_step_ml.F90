!! Semi-implicit free-surface time step — unstructured 3D hydrostatic multilayer.
module rki_si_unstr_step_ml
   !! One semi-implicit (Casulli/SCHISM-style) step for the UNSTRUCTURED 3D
   !! multilayer coastal path (`time_scheme="semi_implicit"` +
   !! `use_multilayer`, `nz>1`, sigma) — built on the P3.1 TPFA eta-solve (the operator
   !! never changes; only the Hf/Gbar fills do — the "strict superset"
   !! contract) and the `(nz, ncells)` column kernels of
   !! `rki_si_unstr_column`.  The step (one implicit theta-step; the theta
   !! scheme IS the integrator):
   !!
   !!   1. save U^n edges + layer state; eta^n = stage at circumcentres;
   !!   2. explicit predictor tendencies IN PLACE on hu/hv_layer:
   !!      baroclinic PGF, Coriolis, Smagorinsky (existing unstructured
   !!      kernels), and surface wind stress on the TOP layer (k=nz) —
   !!      dt*tau/RHO_WATER from the uniform scalar wind (rides the
   !!      pure-u^n predictor into G_k; R_k / Hhat / the operator are
   !!      untouched by linearity; no momentum self-advection — the SI
   !!      target regime is advection-weak, same envelope as structured);
   !!   3. PP81 interface viscosity + per-column implicit solve
   !!      (G_k, R_k, Hhat = -sum h_k R_k, Gbar = sum h_k G_k);
   !!   4. Hf from Hhat via the P3.1 edge face-depth kernel reused verbatim
   !!      (Hhat passed as the cell "depth"; a Dirichlet edge gets
   !!      0.5*(Hhat_c + h_ghost) with h_ghost = eta_bc - b — the
   !!      structured ghost-Hhat seed, edge-wise), edge predictor from the
   !!      face-averaged Gbar, theta-blend RHS + area-weighted CG — ALL
   !!      P3.1 machinery reused verbatim;
   !!   5. per-layer edge-transport recovery (column sums telescope EXACTLY
   !!      to the 2D back-substituted edge transport => mass to round-off;
   !!      Dirichlet ghost R = -1 seed inline — the estuary-mouth lesson);
   !!   6. tracers: every registry tracer advected depth-proportionally
   !!      with the EXACT blended edge mass flux the eta update used
   !!      (uniform tracers invariant per layer; OUTER SHIM: the host
   !!      registry loop passes tracers(it)%hTr/hTr0 into the flat-impl
   !!      kernel — the array-of-derived-types device rule);
   !!   7. h^{n+1} = eta^{n+1} - b; h_layer = dsig*h; per-layer cell
   !!      momenta = least-squares edge reconstruction; barotropic momenta
   !!      = layer sums; U^{n+1} edge prognostic = column sums;
   !!   8. tracer vertical diffusion, surface fluxes, clamp, EOS — the
   !!      existing unstructured kernels reused as-is.  NO diasurface
   !!      exchange: the dsig-proportional advection makes sigma surfaces
   !!      material for the transported fields (the structured design
   !!      decision, carried over).
   !!
   !! The step ends with valid barotropic stage/xmom/ymom/height +
   !! multilayer *_layer state — downstream (output, gauges, restart)
   !! unchanged.  Single-rank, owned-cell loops.
   use rki_constants, only: wp, GRAVITY, CP_WATER
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   use rki_si_unstr_helmholtz, only: si_unstr_face_depths, si_unstr_compute_rhs, &
                                     si_unstr_cg_solve, si_unstr_apply_helmholtz
   use rki_si_unstr_step, only: SI_UNSTR_THETA, &
                                SI_UNSTR_NEWTON_MAX, SI_UNSTR_NEWTON_TOL
   use rki_si_unstr_wetdry, only: si_unstr_wetdry_fill_bsub, &
                                  si_unstr_wetdry_volume, si_unstr_wetdry_wetfrac, &
                                  si_unstr_wetdry_mass_diag, &
                                  si_unstr_wetdry_newton_residual, &
                                  si_unstr_wetdry_wetcount
   use pic_logger, only: global_logger
   use rki_si_unstr_column, only: si_unstr_build_nu_interface, &
                                  si_unstr_column_solve, &
                                  si_unstr_fill_predictor_edges, &
                                  si_unstr_recover_edge_layers, &
                                  si_unstr_fill_mass_flux_edges, &
                                  si_unstr_advect_tracer_dsig, &
                                  si_unstr_edge_layers_to_cells, &
                                  si_unstr_update_thickness_ml, &
                                  si_unstr_apply_wind_top_layer
   use rki_si_unstr_closures, only: si_unstr_fill_nu_iface_from_keps, &
                                    si_unstr_keps_copy_iface, &
                                    si_unstr_build_nu_interface_kpp, &
                                    si_unstr_kpp_diffusion_tracers
   use rki_ml_dynamics_unstr, only: ml_save_layer_state_unstr, &
                                    ml_compute_barotropic_momentum_unstr
   use rki_ml_pressure_unstr, only: ml_baroclinic_pressure_gradient_unstr, &
                                    ml_baroclinic_pressure_gradient_jacobian_unstr
   use rki_ml_friction_unstr, only: ml_apply_coriolis_unstr
   use rki_ml_horizontal_viscosity_unstr, only: ml_horizontal_viscosity_smag_unstr
   use rki_ml_tracers_unstr, only: ml_save_tracer_unstr, ml_clamp_tracer_unstr, &
                                   ml_vertical_diffusion_tracer_unstr, &
                                   ml_apply_surface_tracer_fluxes_unstr
   use rki_ml_eos_unstr, only: ml_equation_of_state_unstr
   use rki_ml_keps_unstr, only: ml_keps_friction_velocities_unstr, &
                                ml_keps_update_unstr, &
                                ml_keps_apply_diffusion_tracer_unstr
   use rki_ml_kpp_unstr, only: ml_kpp_boundary_layer_unstr
   use rki_boundary_unstr, only: MAX_BC_GROUPS, BC_WALL, BC_TIDAL, BC_CLAMPED
   implicit none
   private

   public :: si_unstr_step_ml

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.

contains

   subroutine si_unstr_step_ml(state, mesh, dt, iter_out)
      !! Advance the unstructured 3D multilayer state one semi-implicit
      !! step.  The caller (`solver_unstr_step_semi_implicit`) has already
      !! halo-exchanged; Dirichlet boundary coupling reads
      !! `state%bc_type` / `state%bc_eta_ext` (tidal/clamped) inside the
      !! edge kernels — the P3.1 coefficient-space route.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt
      integer, intent(out) :: iter_out

      integer :: nct, nco, n_bc, nz, it, ip
      real(wp) :: theta, g_theta2_dt2, g_omtheta_dt
      integer :: cg_it
      real(wp) :: mass_res
      real(wp) :: tr_bc(MAX_BC_GROUPS)

      nct = state%ncells_total
      nco = state%ncells_owned
      n_bc = size(state%bc_type)
      nz = state%multilayer%nz_ml
      theta = SI_UNSTR_THETA
      g_theta2_dt2 = GRAVITY*theta*theta*dt*dt
      g_omtheta_dt = GRAVITY*(1.0_wp - theta)*dt

      ! Mirror the explicit stage's S/T config sync onto the tracer registry
      ! (host scalars, idempotent) so the generic per-tracer kernels read
      ! per-tracer config uniformly.
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

      ! ---- 1. Save U^n edges + layer state; eta^n = stage; height ----
      call si_ml_save_and_eta(nct, nco, state%barotropic%stage_cv, &
                              state%barotropic%bed_cv, state%barotropic%height_cv, &
                              state%si%eta, state%si%u_edge, state%si%u_n_edge)
      call ml_save_layer_state_unstr(state)

      ! ---- 2. Explicit predictor tendencies IN PLACE on hu/hv_layer ----
      ! EOS FIRST (the unstructured house order — the explicit stage also
      ! runs EOS immediately before the BPG): rho_layer is zero-initialised,
      ! so a BPG before the first EOS would see a -rho_0 anomaly kick.
      call ml_equation_of_state_unstr(state)
      if (state%multilayer%use_bpg_jacobian) then
         call ml_baroclinic_pressure_gradient_jacobian_unstr(state, mesh, dt)
      else
         call ml_baroclinic_pressure_gradient_unstr(state, mesh, dt)
      end if
      call ml_apply_coriolis_unstr(state, dt)
      if (state%multilayer%smagorinsky_cs > 0.0_wp) then
         call ml_horizontal_viscosity_smag_unstr(state, mesh, dt)
      end if
      call si_unstr_apply_wind_top_layer(nz, nct, nco, state%multilayer%h_layer, &
                                         state%multilayer%hu_layer, &
                                         state%multilayer%hv_layer, &
                                         state%barotropic%wind_tau_x, &
                                         state%barotropic%wind_tau_y, dt)

      ! ---- 3. Closure nu (PP81 / k-eps / KPP) + per-column implicit solve ----
      ! The column solve consumes ANY interface viscosity through nu_iface
      ! (R = solve(-1) is nu-independent), so closures are pure FILLS:
      ! * k-eps: the UNCHANGED ml_keps_update_unstr at SI cadence — already
      !   backward-Euler/Patankar (see local_archive/prototypes/si_closures); then
      !   keps_nut copied into nu_iface.  Optional keps_si_picard > 1
      !   Picard: save pre-step k/eps, re-run from the same initial state
      !   with the refreshed nu_t production coefficient.  The IMPLICIT
      !   momentum apply is the column tridiagonal — do NOT call
      !   ml_keps_apply_momentum_unstr here.
      ! * KPP: ml_kpp_boundary_layer_unstr reused verbatim (diagnostic),
      !   then the LMD94 shape-function overlay MAX'd onto PP81 in nu_iface
      !   — the implicit tridiagonal replaces the explicit CFL-clamped flux.
      ! * default (PP81): bit-identical to the unextended step.
      if (state%multilayer%vmix_use_keps) then
         call ml_keps_friction_velocities_unstr(state)
         if (state%multilayer%keps_si_picard > 1) then
            call si_unstr_keps_copy_iface(nz + 1, nct, nco, &
                                          state%multilayer%keps_tke, &
                                          state%multilayer%keps_tke0)
            call si_unstr_keps_copy_iface(nz + 1, nct, nco, &
                                          state%multilayer%keps_eps, &
                                          state%multilayer%keps_eps0)
         end if
         do ip = 1, max(state%multilayer%keps_si_picard, 1)
            if (ip > 1) then
               ! Picard: restore pre-step k/eps; production now reads
               ! the pass-(ip-1) nu_t (the refreshed coefficient).
               call si_unstr_keps_copy_iface(nz + 1, nct, nco, &
                                             state%multilayer%keps_tke0, &
                                             state%multilayer%keps_tke)
               call si_unstr_keps_copy_iface(nz + 1, nct, nco, &
                                             state%multilayer%keps_eps0, &
                                             state%multilayer%keps_eps)
            end if
            call ml_keps_update_unstr(state, dt, .true.)
         end do
         call si_unstr_fill_nu_iface_from_keps(nz, nz + 1, nct, nco, &
                                               state%multilayer%h_layer, &
                                               state%multilayer%keps_nut, &
                                               state%si%nu_iface)
      else if (state%multilayer%vmix_use_kpp) then
         call ml_kpp_boundary_layer_unstr(state, &
                                          state%multilayer%tracers(state%multilayer%idx_salinity)%hTr)
         call si_unstr_build_nu_interface_kpp(nz, nz + 1, nct, nco, &
                                              state%multilayer%h_layer, &
                                              state%multilayer%hu_layer, &
                                              state%multilayer%hv_layer, &
                                              state%multilayer%rho_layer, &
                                              state%multilayer%rho_0, &
                                              state%multilayer%vmix_nu_max, &
                                              state%multilayer%vmix_nu_bg, &
                                              state%multilayer%vmix_ri_alpha, &
                                              state%multilayer%kpp_u_star, &
                                              state%multilayer%kpp_buoy_flux, &
                                              state%multilayer%kpp_bl_depth, &
                                              state%si%nu_iface)
      else
         call si_unstr_build_nu_interface(nz, nz + 1, nct, nco, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%multilayer%rho_layer, &
                                          state%multilayer%rho_0, &
                                          state%multilayer%vmix_nu_max, &
                                          state%multilayer%vmix_nu_bg, &
                                          state%multilayer%vmix_ri_alpha, &
                                          state%si%nu_iface)
      end if
      call si_unstr_column_solve(nz, nz + 1, nct, nco, state%multilayer%h_layer, &
                                 state%multilayer%hu_layer, &
                                 state%multilayer%hv_layer, state%si%nu_iface, &
                                 state%barotropic%manning_n, dt, &
                                 state%si%g_col_u, state%si%g_col_v, &
                                 state%si%r_col, state%si%Hhat_cell, &
                                 state%si%Gbar_cell_x, state%si%Gbar_cell_y)

      ! ---- 4. Hf from Hhat; zero dead-edge U^n (R3); edge predictor; RHS; CG ----
      ! WET/DRY (state%si%wetdry): the Hhat conveyance is GATED by the
      ! geometric max-eta/sill wet/dry rule (a face dead by geometry — the
      ! higher surface below the higher bed — carries no conveyance even where
      ! Hhat is nonzero), then a Newton outer loop wraps the eta-solve with the
      ! surface wet-fraction mass diagonal a(eta)+floor.  Everything else —
      ! the R3 dead-edge zeroing (STRUCTURAL here: the RHS consumes per-edge
      ! transports directly), the frozen-column responses (g/r/Hhat/Gbar
      ! depend on h_layer at time n, not eta), the dsig scale-together dry
      ! columns, the hold-last-wet tracer drain — is the all-wet ML machinery
      ! reused verbatim.  All-wet Hhat basins reduce the gated conveyance
      ! exactly to the P3.1 face depth => the wet/dry branch is a strict
      ! superset (the structured "strict superset" contract, edge-wise).
      if (state%si%wetdry) then
         call si_unstr_wetdry_fill_bsub(mesh, state%si%n_sub, nct, nco, &
                                        state%barotropic%bed_cv, &
                                        state%si%bsub_fill, state%si%b_sub, &
                                        state%si%b_min)
         call si_ml_wetdry_face_depths(mesh, nct, nco, n_bc, state%si%eta, &
                                       state%si%b_min, state%si%Hhat_cell, &
                                       state%bc_type, state%bc_eta_ext, &
                                       state%si%Hf_edge, state%si%eta_ghost_edge)
      else
         call si_unstr_face_depths(mesh, nct, nco, n_bc, state%si%Hhat_cell, &
                                   state%barotropic%bed_cv, state%bc_type, &
                                   state%bc_eta_ext, state%si%Hf_edge, &
                                   state%si%eta_ghost_edge)
      end if
      call si_ml_zero_dead_edges(nct, nco, state%si%Hf_edge, state%si%u_n_edge)
      call si_unstr_fill_predictor_edges(mesh, nz, nct, nco, state%si%Gbar_cell_x, &
                                         state%si%Gbar_cell_y, &
                                         state%multilayer%hu_layer0, &
                                         state%multilayer%hv_layer0, &
                                         state%si%u_n_edge, state%si%eta, &
                                         state%si%Hf_edge, state%si%eta_ghost_edge, &
                                         state%si%u_star, g_omtheta_dt)
      call si_unstr_compute_rhs(mesh, nct, nco, state%si%eta, state%si%u_star, &
                                state%si%u_n_edge, state%si%eta_r, dt, theta, &
                                state%si%Hf_edge, state%si%eta_ghost_edge, &
                                g_theta2_dt2)
      if (state%si%wetdry) then
         ! The Newton wrapper writes state%si%mass_res / newton_iters directly.
         call si_ml_wetdry_newton(state, mesh, nct, nco, g_theta2_dt2, cg_it)
      else
         call si_unstr_cg_solve(mesh, nct, nco, state%si%eta, state%si%eta_r, &
                                state%si%eta_p, state%si%eta_Ap, state%si%eta_z, &
                                state%si%eta_diag, state%si%Hf_edge, g_theta2_dt2, &
                                state%si%precond, state%si%rtol, cg_it, mass_res)
         state%si%mass_res = mass_res
      end if
      state%si%cg_iters = cg_it
      iter_out = cg_it

      ! ---- 5. Per-layer edge-transport recovery ----
      ! eta^n is still in stage_cv (untouched until step 7); eta^{n+1} in
      ! si%eta.  h_layer is still the time-n field that fed the tridiagonal
      ! (one-h-per-stage rule).
      call si_unstr_recover_edge_layers(mesh, nz, nct, nco, &
                                        state%multilayer%h_layer, &
                                        state%multilayer%hu_layer0, &
                                        state%multilayer%hv_layer0, &
                                        state%multilayer%dsig_ml, &
                                        state%barotropic%bed_cv, &
                                        state%si%g_col_u, state%si%g_col_v, &
                                        state%si%r_col, state%barotropic%stage_cv, &
                                        state%si%eta, state%si%Hf_edge, &
                                        state%si%eta_ghost_edge, dt, theta, &
                                        state%si%u_edge_layer)

      ! ---- 6. Blended mass flux + dsig tracer advection (registry loop,
      !         OUTER SHIM: plain arrays into the flat-impl kernel) ----
      call si_unstr_fill_mass_flux_edges(nz, nct, nco, state%si%u_edge_layer, &
                                         state%si%u_n_edge, state%si%Hf_edge, &
                                         theta, state%si%mass_flux_edge, &
                                         state%si%u_edge)
      do it = 1, size(state%multilayer%tracers)
         if (.not. state%multilayer%tracers(it)%do_horizontal_advection) cycle
         call ml_save_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                   state%multilayer%tracers(it)%hTr0)
         ! Pass the Dirichlet inflow table through a FRESH host copy: the
         ! tracers array-of-DTs is device-mapped at enter_data, so the
         ! per-step host-side config sync above never reaches the device
         ! copy of the inline `bc_ext` component — a kernel reading
         ! `tracers(it)%bc_ext` directly sees a STALE (init-time) device
         ! value and the flood tide injects tracer 0 (measured 29-PSU
         ! freshening on GPU, invisible on the CPU builds where host ==
         ! device).  A local array is never in the present table => the
         ! launch copies the current host values in.
         tr_bc = state%multilayer%tracers(it)%bc_ext
         call si_unstr_advect_tracer_dsig(mesh, nz, nct, nco, n_bc, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%si%mass_flux_edge, &
                                          tr_bc, &
                                          state%multilayer%tracers(it)%hTr0, &
                                          state%multilayer%tracers(it)%hTr, dt)
      end do

      ! ---- 7. Thickness + momentum closers ----
      call si_unstr_update_thickness_ml(nz, nct, nco, state%si%eta, &
                                        state%barotropic%bed_cv, &
                                        state%multilayer%dsig_ml, &
                                        state%barotropic%stage_cv, &
                                        state%barotropic%height_cv, &
                                        state%multilayer%h_layer)
      call si_unstr_edge_layers_to_cells(mesh, nz, nct, nco, &
                                         state%si%u_edge_layer, &
                                         state%multilayer%hu_layer, &
                                         state%multilayer%hv_layer)
      call ml_compute_barotropic_momentum_unstr(state)

      ! ---- 8. Tracer vertical diffusion + surface fluxes + clamp + EOS ----
      ! Closure-aware tracer vertical diffusion (mirrors the structured
      ! SI dispatch, but implicit-capable kernels only — SI step sizes):
      ! * k-eps: implicit Thomas diffusion with the k-eps scalar diffusivity
      !   kt for EVERY registry tracer (conservative).
      ! * KPP: implicit Thomas diffusion with the PP81-cubic + overlay
      !   diffusivity + non-local transport for S/T; extra passive tracers
      !   fall back to the PP81 kernel (as on the explicit path).
      ! * default: the PP81 kernel, unchanged (bit-identical).
      if (state%multilayer%vmix_use_keps) then
         do it = 1, size(state%multilayer%tracers)
            if (state%multilayer%tracers(it)%do_vertical_diffusion) then
               call ml_keps_apply_diffusion_tracer_unstr(state, &
                                                         state%multilayer%tracers(it)%hTr, &
                                                         state%multilayer%tracers(it)%kappa_bg, &
                                                         dt)
            end if
         end do
      else if (state%multilayer%vmix_use_kpp) then
         call si_unstr_kpp_diffusion_tracers(nz, nz + 1, nct, nco, &
                                             state%multilayer%h_layer, state%multilayer%hu_layer, &
                                             state%multilayer%hv_layer, &
                                             state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                             state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                             state%multilayer%rho_layer, state%multilayer%rho_0, &
                                             state%multilayer%vmix_nu_max, &
                                             state%multilayer%kappa_S_bg, state%multilayer%kappa_T_bg, &
                                             state%multilayer%vmix_ri_alpha, &
                                             state%multilayer%kpp_cs_nonlocal, &
                                             state%forcing%heat_flux, state%forcing%evap_minus_precip, &
                                             CP_WATER, state%multilayer%kpp_u_star, &
                                             state%multilayer%kpp_buoy_flux, state%multilayer%kpp_bl_depth, dt)
         ! Extra passive tracers (not S/T) fall back to PP81.
         do it = 1, size(state%multilayer%tracers)
            if (it == state%multilayer%idx_salinity .or. &
                it == state%multilayer%idx_temperature) cycle
            if (state%multilayer%tracers(it)%do_vertical_diffusion) then
               call ml_vertical_diffusion_tracer_unstr(state, &
                                                       state%multilayer%tracers(it)%hTr, &
                                                       state%multilayer%tracers(it)%kappa_bg, &
                                                       dt)
            end if
         end do
      else
         do it = 1, size(state%multilayer%tracers)
            if (state%multilayer%tracers(it)%do_vertical_diffusion) then
               call ml_vertical_diffusion_tracer_unstr(state, &
                                                       state%multilayer%tracers(it)%hTr, &
                                                       state%multilayer%tracers(it)%kappa_bg, &
                                                       dt)
            end if
         end do
      end if
      if (state%multilayer%idx_salinity > 0 .and. &
          state%multilayer%idx_temperature > 0) then
         call ml_apply_surface_tracer_fluxes_unstr(state, &
                                                   state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                                   state%multilayer%tracers(state%multilayer%idx_temperature)%hTr)
      end if
      do it = 1, size(state%multilayer%tracers)
         if (state%multilayer%tracers(it)%do_clamp) then
            call ml_clamp_tracer_unstr(state, state%multilayer%tracers(it)%hTr, &
                                       state%multilayer%tracers(it)%tr_min, &
                                       state%multilayer%tracers(it)%tr_max)
         end if
      end do
      ! (EOS runs at the TOP of the next step, before its BPG — the
      ! unstructured house order; no end-of-step call needed.)

      state%time = state%time + dt
      state%step_count = state%step_count + 1
      state%dt = dt
   end subroutine si_unstr_step_ml

   ! ------------------------------------------------------------------
   ! Flat-impl helpers (explicit-shape, do concurrent, device-resident).
   ! ------------------------------------------------------------------

   subroutine si_ml_save_and_eta(nct, nco, stage, bed, height, eta, &
                                      u_edge, u_n_edge)
      !! Save the edge transports U^n, refresh height = max(stage - b, 0),
      !! and form eta^n = stage at circumcentres (the CG warm start).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: stage(nct), bed(nct)
      real(wp), intent(inout) :: height(nct)
      real(wp), intent(out) :: eta(nct)
      real(wp), intent(in) :: u_edge(EDGES_PER_TRI, nct)
      real(wp), intent(out) :: u_n_edge(EDGES_PER_TRI, nct)

      integer :: c, i

      !$omp target teams distribute parallel do private(i)
      do c = 1, nco
         do i = 1, EDGES_PER_TRI
            u_n_edge(i, c) = u_edge(i, c)
         end do
         height(c) = max(stage(c) - bed(c), 0.0_wp)
         eta(c) = stage(c)
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_ml_save_and_eta

   subroutine si_ml_zero_dead_edges(nct, nco, Hf_edge, u_n_edge)
      !! Zero the saved U^n on dead (wall / Hf=0) edges — Casulli rule R3:
      !! a dead edge carries no (1-theta) flux into the theta-blend RHS.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(inout) :: u_n_edge(EDGES_PER_TRI, nct)

      integer :: c, i

      !$omp target teams distribute parallel do private(i)
      do c = 1, nco
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, c) <= 0.0_wp) u_n_edge(i, c) = 0.0_wp
         end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_ml_zero_dead_edges

   subroutine si_ml_wetdry_face_depths(mesh, nct, nco, n_bc, eta, b_min, &
                                            Hhat, bc_type, bc_eta_ext, &
                                            Hf_edge, eta_ghost_edge)
      !! Frozen wet/dry edge depths for the 3D column-integrated eta-solve —
      !! the triangular twin of the structured `si_ml_wetdry_face_depths`:
      !! the effective CONVEYANCE is the Hhat-based two-cell average (as the
      !! all-wet ML path uses), GATED by the geometric max-eta/sill wet/dry
      !! rule so an edge dead by geometry (higher surface below the higher
      !! bed) carries Hf=0 — a wall — regardless of Hhat.  The surface
      !! geometry decides wet/dead; where wet, the implicit column response
      !! Hhat sets the conveyance magnitude.  b_min = the geometric sill bed
      !! (deepest sub-bottom).  Dirichlet edges gate the ghost-coupled face
      !! the same way (ghost surface = the prescribed eta_bc; ghost sill
      !! mirrors the cell sill).  All-wet Hhat basins reduce this EXACTLY to
      !! the P3.1 `si_unstr_face_depths` conveyance (0.5*(Hhat_c+Hhat_nbr)),
      !! so the wet/dry branch is a strict superset.  Cell-centric
      !! double-visit `(3,ncell)`, cell-index-last.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      integer, intent(in) :: n_bc
         !! BC group table size.
      real(wp), intent(in) :: eta(nct)
         !! eta^n at circumcentres (the FROZEN surface).
      real(wp), intent(in) :: b_min(nct)
         !! Geometric sill bed (deepest sub-bottom) at cells (m).
      real(wp), intent(in) :: Hhat(nct)
         !! Column effective depth (conveyance magnitude, m).
      integer, intent(in) :: bc_type(n_bc)
      real(wp), intent(in) :: bc_eta_ext(n_bc)
      real(wp), intent(out) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(out) :: eta_ghost_edge(EDGES_PER_TRI, nct)

      integer :: c, i, kn, g, bt
      real(wp) :: geo, eta_bc

      !$omp target teams distribute parallel do private(i, kn, g, bt, geo, eta_bc)
      do c = 1, nco
         do i = 1, EDGES_PER_TRI
            eta_ghost_edge(i, c) = 0.0_wp
            Hf_edge(i, c) = 0.0_wp
            kn = mesh%neighbours(i, c)
            if (kn > 0) then
               ! Interior edge: geometric sill gate on eta^n, Hhat conveyance.
               geo = max(max(eta(c), eta(kn)) - max(b_min(c), b_min(kn)), 0.0_wp)
               if (geo > 0.0_wp) Hf_edge(i, c) = 0.5_wp*(Hhat(c) + Hhat(kn))
            else
               g = mesh%bc_group(i, c)
               bt = BC_WALL
               if (g >= 1 .and. g <= n_bc) bt = bc_type(g)
               if ((bt == BC_TIDAL .or. bt == BC_CLAMPED) .and. g >= 1) then
                  eta_bc = bc_eta_ext(g)
                  ! Ghost sill mirrors the cell sill.  Where geometrically wet,
                  ! the conveyance is the SAME 0.5*(Hhat_c + h_ghost) mean the
                  ! all-wet P3.1 path uses (Hhat interior, geometric ghost
                  ! depth eta_bc - b) — matching it exactly so the recovery
                  ! reads the same Dirichlet Hf the eta-solve used (the
                  ! estuary-mouth mass<->tracer consistency lesson).
                  geo = max(max(eta(c), eta_bc) - b_min(c), 0.0_wp)
                  if (geo > 0.0_wp) Hf_edge(i, c) = &
                     0.5_wp*(Hhat(c) + max(eta_bc - b_min(c), 0.0_wp))
                  eta_ghost_edge(i, c) = eta_bc
               end if
            end if
         end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_ml_wetdry_face_depths

   subroutine si_ml_wetdry_newton(state, mesh, nct, nco, g_theta2_dt2, iter_out)
      !! Newton outer loop for the 3D column-integrated wet/dry eta-solve on
      !! the triangulation — the triangular twin of the structured
      !! `si_ml_wetdry_newton` (and the 2D-unstr `si_unstr_advance_wetdry`'s
      !! Newton loop).  The column responses (g_col/r_col/Hhat/Gbar) and the
      !! frozen Hf are held FIXED (they depend on h_layer at time n, not eta),
      !! so only the surface mass system F = v(eta) + T*eta - rhs_d iterates,
      !! with T the Hhat-weighted TPFA stiffness and the surface wet-fraction
      !! Newton diagonal a(eta)+floor.  On entry `state%si%eta_r` holds the
      !! LINEAR RHS (leading eta^n from si_unstr_compute_rhs); this swaps the
      !! leading term for v(eta^n) once, then Newton-solves.  Host-side
      !! orchestration; flat device kernels.  Writes state%si%newton_iters
      !! and state%si%mass_res = the NONLINEAR mass residual ||F||^2 at
      !! Newton exit (the conservation diagnostic — NOT the inner-CG residual
      !! of the last delta-solve, which is an internal solver number).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nct, nco
      real(wp), intent(in) :: g_theta2_dt2
      integer, intent(out) :: iter_out
         !! Inner CG iterations summed over the Newton outers.

      integer :: k, nit, cg_it, cg_it_total, wetcount, wetcount_prev
      real(wp) :: rhs2, f2, tol2, mass_res
      logical :: converged

      associate (si => state%si)
         ! rhs_d = eta_r + (v(eta^n) - eta^n): swap the linear RHS leading
         ! term for the wetted volume (same as the 2D/structured wet/dry).
         call si_unstr_wetdry_volume(si%n_sub, nct, nco, si%eta, si%b_sub, si%v_cell)
         !$omp target teams distribute parallel do
         do k = 1, nco
            si%rhs_d(k) = si%eta_r(k) + (si%v_cell(k) - si%eta(k))
         end do
         !$omp end target teams distribute parallel do

         rhs2 = 0.0_wp
         !$omp target teams distribute parallel do
         do k = 1, nco
            rhs2 = rhs2 + si%rhs_d(k)*si%rhs_d(k)
         end do
         !$omp end target teams distribute parallel do
         tol2 = max(SI_UNSTR_NEWTON_TOL*SI_UNSTR_NEWTON_TOL*rhs2, 1.0e-60_wp)

         cg_it_total = 0
         cg_it = 0
         mass_res = 0.0_wp
         wetcount_prev = -1
         converged = .false.
         do nit = 1, SI_UNSTR_NEWTON_MAX
            call si_unstr_wetdry_wetfrac(si%n_sub, nct, nco, si%eta, si%b_sub, si%a_wet)
            call si_unstr_wetdry_volume(si%n_sub, nct, nco, si%eta, si%b_sub, si%v_cell)
            call si_unstr_wetdry_mass_diag(nct, nco, si%a_wet, si%Hf_edge, si%mass_diag)
            call si_unstr_apply_helmholtz(mesh, nct, nco, si%eta, si%newton_Jeta, &
                                          si%Hf_edge, g_theta2_dt2, &
                                          mass_diag=si%mass_diag)
            call si_unstr_wetdry_newton_residual(nct, nco, si%eta, si%newton_Jeta, &
                                                 si%mass_diag, si%v_cell, &
                                                 si%rhs_d, si%newton_F, f2)
            call si_unstr_wetdry_wetcount(si%n_sub, nct, nco, si%eta, si%b_sub, wetcount)
            if (nit > 1 .and. wetcount == wetcount_prev .and. f2 < tol2) then
               converged = .true.
               exit
            end if
            wetcount_prev = wetcount

            ! Solve J*dz = -F from zero (full-array zero so the ghost band of
            ! the search space stays 0); eta += dz.
            !$omp target teams distribute parallel do
            do k = 1, nct
               si%newton_dz(k) = 0.0_wp
            end do
            !$omp end target teams distribute parallel do
            !$omp target teams distribute parallel do
            do k = 1, nco
               si%eta_r(k) = -si%newton_F(k)
            end do
            !$omp end target teams distribute parallel do
            call si_unstr_cg_solve(mesh, nct, nco, si%newton_dz, si%eta_r, &
                                   si%eta_p, si%eta_Ap, si%eta_z, si%eta_diag, &
                                   si%Hf_edge, g_theta2_dt2, si%precond, &
                                   si%rtol, cg_it, mass_res, &
                                   mass_diag=si%mass_diag)
            cg_it_total = cg_it_total + cg_it
            !$omp target teams distribute parallel do
            do k = 1, nco
               si%eta(k) = si%eta(k) + si%newton_dz(k)
            end do
            !$omp end target teams distribute parallel do
         end do

         si%newton_iters = min(nit, SI_UNSTR_NEWTON_MAX)
         ! Report the NONLINEAR mass residual ||F||^2 (the conservation
         ! diagnostic), not `mass_res` (the inner-CG residual of the final
         ! delta-solve — an internal solver number that misleads anyone
         ! reading state%si%mass_res as a mass defect).
         si%mass_res = f2
         if (.not. converged) then
            call global_logger%warning("si_ml_wetdry_newton: Newton did not "// &
                                       "converge in the iteration cap; last ||F||^2 above tol.")
         end if
      end associate

      iter_out = cg_it_total
   end subroutine si_ml_wetdry_newton

end module rki_si_unstr_step_ml
