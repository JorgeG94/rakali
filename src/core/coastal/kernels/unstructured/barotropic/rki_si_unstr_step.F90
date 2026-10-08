!! Semi-implicit free-surface time step (unstructured 2D barotropic).
module rki_si_unstr_step
   !! One semi-implicit (Casulli) step for the 2D barotropic free surface on
   !! the UNSTRUCTURED (triangular) backend — the TPFA / edge-transport
   !! scheme.  INTERNALLY STAGGERED: the predictor / back-sub
   !! transports live on EDGES (per-cell double-visit `u_edge(3,ncell)`,
   !! transport along each edge's outward normal), matching Casulli's
   !! staggered scheme, so the eta update is an exact edge-flux
   !! finite-volume continuity step (mass conserved to CG tolerance +
   !! round-off).  eta lives at the triangle CIRCUMCENTRE.  The step:
   !!
   !!   1. save U^n (u_edge -> u_n_edge) and form eta^n = stage_cv at
   !!      circumcentres (stage = h + b already; eta == stage);
   !!   2. build per-edge frozen effective depths Hf (wall edges 0, Dirichlet
   !!      edges real) + the per-edge Dirichlet ghost eta;
   !!   3. EDGE predictor U*_e = U^n_e - g(1-theta)*dt*edge_flux(Hf,eta^n)
   !!      + wind, gated on Hf>0 (dead edges => exactly 0);
   !!   4. theta-blend RHS from the edge transports (+ Dirichlet ghost fold);
   !!   5. area-weighted CG solve of the SPD TPFA Helmholtz for eta^{n+1}
   !!      (warm-started from eta^n; the load-bearing area weight);
   !!   6. EDGE back-substitution U^{n+1}_e = U*_e - g*theta*dt*edge_flux(Hf,
   !!      eta^{n+1}) — the identical frozen Hf the operator used, so the
   !!      edge transports satisfy discrete continuity EXACTLY;
   !!   7. reconstruct the A-grid cell momenta (Green-Gauss from the edge
   !!      normals) for the KNP seam + recover h^{n+1} = eta^{n+1} - b.
   !!
   !! Single-rank, owned-cell loops (no per-CG-iteration halo).  The step
   !! ends with barotropic stage_cv/xmom_cv/ymom_cv/height_cv valid, so BCs,
   !! diagnostics, output, gauges work unchanged.
   use rki_constants, only: wp, GRAVITY, RHO_WATER
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   use rki_si_unstr_helmholtz, only: si_unstr_edge_flux, si_unstr_face_depths, &
                                     si_unstr_compute_rhs, si_unstr_cg_solve, &
                                     si_unstr_apply_helmholtz
   use rki_si_unstr_mg, only: agglom_hierarchy_unstr_t, &
                              agglom_hierarchy_unstr_refresh
   use rki_si_unstr_state, only: SI_PRECOND_MG
   use rki_si_unstr_wetdry, only: si_unstr_wetdry_fill_bsub, &
                                  si_unstr_wetdry_volume, &
                                  si_unstr_wetdry_wetfrac, &
                                  si_unstr_wetdry_mass_diag, &
                                  si_unstr_wetdry_face_depths, &
                                  si_unstr_wetdry_newton_residual, &
                                  si_unstr_wetdry_wetcount
   use pic_logger, only: global_logger
   implicit none
   private

   public :: si_unstr_step
   public :: si_unstr_reconstruct_cell_momenta
   public :: SI_UNSTR_THETA
   public :: SI_UNSTR_NEWTON_MAX, SI_UNSTR_NEWTON_TOL

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.
   real(wp), parameter :: SI_UNSTR_THETA = 0.6_wp
      !! Implicitness parameter.  theta >= 0.5 => unconditionally stable for
      !! the gravity wave; theta > 0.5 adds damping ~(theta-0.5).
   real(wp), parameter :: DELTA_EPS = 1.0e-30_wp
      !! 1/0 armour on the circumcentre-line distance.
   integer, parameter :: SI_UNSTR_NEWTON_MAX = 10
      !! Cap on the wet/dry Newton outer iterations per step (prototype
      !! measured 1-3 outer/step on triangles, max 5; fail-loud ceiling).
   real(wp), parameter :: SI_UNSTR_NEWTON_TOL = 1.0e-11_wp
      !! Newton convergence tolerance on the relative nonlinear residual,
      !! combined with a frozen wet active set (matches structured).

contains

   subroutine si_unstr_step(state, mesh, dt, iter_out)
      !! Advance the 2D barotropic free surface one semi-implicit step on the
      !! triangulation.  Reads the frozen Hf / ghost from the current
      !! surface, solves the Helmholtz eta-system, and back-substitutes the
      !! edge transports.  `iter_out` returns the CG iteration count (summed
      !! over the Newton outers in wet/dry mode).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt
      integer, intent(out) :: iter_out

      integer :: nct, nco, n_bc
      real(wp) :: theta, g_theta2_dt2, gtheta_dt, g_omtheta_dt
      integer :: cg_it
      real(wp) :: mass_res

      nct = state%ncells_total
      nco = state%ncells_owned
      n_bc = size(state%bc_type)
      theta = SI_UNSTR_THETA
      g_theta2_dt2 = GRAVITY*theta*theta*dt*dt
      gtheta_dt = GRAVITY*theta*dt
      g_omtheta_dt = GRAVITY*(1.0_wp - theta)*dt

      if (state%si%wetdry) then
         call si_unstr_advance_wetdry(state, mesh, nct, nco, n_bc, dt, theta, &
                                      g_theta2_dt2, gtheta_dt, g_omtheta_dt, &
                                      cg_it)
         mass_res = state%si%mass_res
      else
         call si_unstr_advance(mesh, nct, nco, n_bc, &
                               state%barotropic%stage_cv, state%barotropic%xmom_cv, &
                               state%barotropic%ymom_cv, state%barotropic%bed_cv, &
                               state%barotropic%height_cv, &
                               state%si%eta, state%si%eta_r, state%si%eta_p, &
                               state%si%eta_Ap, state%si%eta_z, state%si%eta_diag, &
                               state%si%u_edge, state%si%u_star, state%si%Hf_edge, &
                               state%si%eta_ghost_edge, state%si%u_n_edge, &
                               state%bc_type, state%bc_eta_ext, &
                               state%barotropic%wind_tau_x, state%barotropic%wind_tau_y, &
                               dt, theta, g_theta2_dt2, gtheta_dt, g_omtheta_dt, &
                               state%si%precond, state%si%rtol, cg_it, mass_res, &
                               mg=state%si%mg)
         state%si%mass_res = mass_res
      end if

      state%si%cg_iters = cg_it
      state%time = state%time + dt
      state%step_count = state%step_count + 1
      state%dt = dt
      iter_out = cg_it
   end subroutine si_unstr_step

   subroutine si_unstr_advance(mesh, nct, nco, n_bc, stage, xmom, ymom, bed, &
                               height, eta, eta_r, eta_p, eta_Ap, eta_z, &
                               eta_diag, u_edge, u_star, Hf_edge, &
                               eta_ghost_edge, u_n_edge, bc_type, &
                               bc_eta_ext, wind_tau_x, wind_tau_y, &
                               dt, theta, g_theta2_dt2, gtheta_dt, &
                               g_omtheta_dt, precond, rtol, iter_out, &
                               mass_res_out, mg)
      !! Flat-impl core: explicit-shape state arrays only (NVHPC device
      !! codegen — assumed-shape dummies in do-concurrent kernels emit
      !! per-launch descriptor-walk memcpys).  The edge workspaces
      !! `eta_ghost_edge` / `u_n_edge` are persistent slots on the SI state
      !! (mapped in its enter_data — no per-step allocation).
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total) — the explicit-shape dim.
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      integer, intent(in) :: n_bc
         !! BC group table size (`MAX_BC_GROUPS`).
      real(wp), intent(inout) :: stage(nct), xmom(nct), ymom(nct), height(nct)
      real(wp), intent(in) :: bed(nct)
      real(wp), intent(inout) :: eta(nct), eta_r(nct), eta_p(nct), eta_Ap(nct)
      real(wp), intent(inout) :: eta_z(nct), eta_diag(nct)
      real(wp), intent(inout) :: u_edge(EDGES_PER_TRI, nct)
      real(wp), intent(inout) :: u_star(EDGES_PER_TRI, nct)
      real(wp), intent(inout) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(inout) :: eta_ghost_edge(EDGES_PER_TRI, nct)
      real(wp), intent(inout) :: u_n_edge(EDGES_PER_TRI, nct)
      integer, intent(in) :: bc_type(n_bc)
      real(wp), intent(in) :: bc_eta_ext(n_bc)
      real(wp), intent(in) :: wind_tau_x, wind_tau_y
      integer, intent(in) :: precond
      real(wp), intent(in) :: dt, theta, g_theta2_dt2, gtheta_dt, g_omtheta_dt, rtol
      integer, intent(out) :: iter_out
      real(wp), intent(out) :: mass_res_out
      type(agglom_hierarchy_unstr_t), intent(inout), optional :: mg
         !! Aggregation-MG hierarchy (present iff precond == MG); its operator
         !! is refreshed here for this step's frozen Hf before the solve.

      integer :: k, i, kn
      real(wp) :: inv_delta, dgdn, un_e, tau_n

      ! ---- 1. Save U^n; refresh height h = max(stage - b, 0); eta^n = stage
      !         at circumcentres (stage = h + b already, so eta == stage) ----
      do concurrent(k=1:nco)
         do i = 1, EDGES_PER_TRI
            u_n_edge(i, k) = u_edge(i, k)
         end do
         height(k) = max(stage(k) - bed(k), 0.0_wp)
         eta(k) = stage(k)
      end do

      ! ---- 2. Frozen Hf + per-edge Dirichlet ghost eta ----
      call si_unstr_face_depths(mesh, nct, nco, n_bc, height, bed, bc_type, &
                                bc_eta_ext, Hf_edge, eta_ghost_edge)

      ! ---- 2b. Refresh the aggregation-MG operator for this step's frozen Hf
      !          (host build from the live Hf; structures static).  Once per
      !          step, NOT per V-cycle. ----
      if (precond == SI_PRECOND_MG .and. present(mg)) then
         !$omp target update from(Hf_edge)
         call agglom_hierarchy_unstr_refresh(mg, mesh, nco, g_theta2_dt2, &
                                             Hf_edge(:, 1:nco))
      end if

      ! ---- 3. EDGE predictor U*_e (dead edges => 0) + wind ----
      ! U^n_e = the saved edge transport; the wind stress projects onto the
      ! outward edge-normal.  Wall edges (Hf=0) stay exactly 0.
      do concurrent(k=1:nco) local(i, kn, inv_delta, dgdn, un_e, tau_n)
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, k) > 0.0_wp) then
               inv_delta = 1.0_wp/max(mesh%edge_delta(i, k), DELTA_EPS)
               kn = mesh%neighbours(i, k)
               if (kn > 0) then
                  dgdn = si_unstr_edge_flux(Hf_edge(i, k), eta(k), eta(kn), inv_delta)
               else
                  ! Dirichlet edge: gradient against the seeded ghost surface.
                  dgdn = si_unstr_edge_flux(Hf_edge(i, k), eta(k), &
                                            eta_ghost_edge(i, k), inv_delta)
               end if
               un_e = u_n_edge(i, k)
               ! Wind stress projected onto the outward edge-normal (Pa/rho).
               tau_n = (wind_tau_x*mesh%edge_n1(i, k) &
                        + wind_tau_y*mesh%edge_n2(i, k))/RHO_WATER
               u_star(i, k) = un_e - g_omtheta_dt*dgdn + dt*tau_n
            else
               u_star(i, k) = 0.0_wp
               u_n_edge(i, k) = 0.0_wp   ! R3: dead edge carries no U^n flux
            end if
         end do
      end do

      ! ---- 4. theta-blend RHS (+ Dirichlet ghost fold) ----
      call si_unstr_compute_rhs(mesh, nct, nco, eta, u_star, u_n_edge, eta_r, &
                                dt, theta, Hf_edge, eta_ghost_edge, g_theta2_dt2)

      ! ---- 5. area-weighted CG solve for eta^{n+1} (warm start eta^n) ----
      if (present(mg)) then
         call si_unstr_cg_solve(mesh, nct, nco, eta, eta_r, eta_p, eta_Ap, eta_z, &
                                eta_diag, Hf_edge, g_theta2_dt2, precond, rtol, &
                                iter_out, mass_res_out, mg=mg)
      else
         call si_unstr_cg_solve(mesh, nct, nco, eta, eta_r, eta_p, eta_Ap, eta_z, &
                                eta_diag, Hf_edge, g_theta2_dt2, precond, rtol, &
                                iter_out, mass_res_out)
      end if

      ! ---- 6. EDGE back-substitution with the SAME frozen Hf ----
      do concurrent(k=1:nco) local(i, kn, inv_delta, dgdn)
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, k) > 0.0_wp) then
               inv_delta = 1.0_wp/max(mesh%edge_delta(i, k), DELTA_EPS)
               kn = mesh%neighbours(i, k)
               if (kn > 0) then
                  dgdn = si_unstr_edge_flux(Hf_edge(i, k), eta(k), eta(kn), inv_delta)
               else
                  dgdn = si_unstr_edge_flux(Hf_edge(i, k), eta(k), &
                                            eta_ghost_edge(i, k), inv_delta)
               end if
               u_edge(i, k) = u_star(i, k) - gtheta_dt*dgdn
            else
               u_edge(i, k) = 0.0_wp
            end if
         end do
      end do

      ! ---- 7. Reconstruct A-grid cell momenta + recover h^{n+1} ----
      call si_unstr_reconstruct_cell_momenta(mesh, nct, nco, u_edge, xmom, ymom)

      do concurrent(k=1:nco)
         stage(k) = eta(k)
         height(k) = max(eta(k) - bed(k), 0.0_wp)
      end do
   end subroutine si_unstr_advance

   subroutine si_unstr_advance_wetdry(state, mesh, nct, nco, n_bc, dt, theta, &
                                      g_theta2_dt2, gtheta_dt, g_omtheta_dt, &
                                      iter_out)
      !! One Casulli wetting/drying SI step on the triangulation — the
      !! Newton wrapper around the area-weighted TPFA CG (the structured
      !! `si_advance_wetdry`, edge-wise; V/a/Newton per-cell, prototype-
      !! confirmed mesh-agnostic).  The nonlinear system is
      !!   F(eta) = v(eta) + T*eta - rhs_d = 0,  J = diag(a+floor) + T
      !! with T the SAME frozen-Hf stiffness the linear code has (gated-mean
      !! edge depths — all-wet reduces exactly to the linear edge depth) and
      !! the constant mass `1` replaced per cell by the wet fraction a.
      !! V convex => plain Newton converges monotonically (no line search).
      !!
      !! The FROZEN-Hf step (rules R1-R5):
      !!   0. fill b_sub/b_min from bed; save U^n; eta^n = stage;
      !!   1. freeze Hf by the gated-mean/sill rule (+ Dirichlet ghost);
      !!   2. edge predictor U* (dead edges => exactly 0; U^n zeroed on dead
      !!      edges — R3 is STRUCTURAL here: the RHS consumes per-edge
      !!      transports directly, so no reconstruction masking is needed);
      !!   3. rhs_d = v(eta^n) - dt*div[theta*U* + (1-theta)*U^n] (+ ghost);
      !!   4. Newton: a,v,F at eta^k; converged (frozen active set + small
      !!      residual)? else CG-solve J*dz = -F from zero, eta += dz
      !!      (dry cells settle below bed — never clamped);
      !!   5. back-sub with the SAME frozen Hf (R2);
      !!   6. stage = eta; h = max(eta - b, 0); cell momenta.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nct, nco, n_bc
      real(wp), intent(in) :: dt, theta, g_theta2_dt2, gtheta_dt, g_omtheta_dt
      integer, intent(out) :: iter_out
         !! Inner CG iterations summed over the Newton outers.

      integer :: k, i, kn, nit, cg_it, cg_it_total, wetcount, wetcount_prev
      real(wp) :: inv_delta, dgdn, tau_n
      real(wp) :: rhs2, f2, tol2, mass_res
      logical :: converged

      associate (si => state%si, baro => state%barotropic)

         ! ---- 0. Sub-bottom table; save U^n; eta^n = stage; height ----
         call si_unstr_wetdry_fill_bsub(mesh, si%n_sub, nct, nco, baro%bed_cv, &
                                        si%bsub_fill, si%b_sub, si%b_min)
         do concurrent(k=1:nco) local(i)
            do i = 1, EDGES_PER_TRI
               si%u_n_edge(i, k) = si%u_edge(i, k)
            end do
            baro%height_cv(k) = max(baro%stage_cv(k) - baro%bed_cv(k), 0.0_wp)
            si%eta(k) = baro%stage_cv(k)
         end do

         ! ---- 1. Freeze Hf (gated-mean/sill rule) + Dirichlet ghost ----
         call si_unstr_wetdry_face_depths(mesh, nct, nco, n_bc, si%eta, &
                                          si%b_min, state%bc_type, &
                                          state%bc_eta_ext, si%Hf_edge, &
                                          si%eta_ghost_edge)

         ! ---- 2. EDGE predictor (dead edges => exactly 0) + wind (R3) ----
         do concurrent(k=1:nco) local(i, kn, inv_delta, dgdn, tau_n)
            do i = 1, EDGES_PER_TRI
               if (si%Hf_edge(i, k) > 0.0_wp) then
                  inv_delta = 1.0_wp/max(mesh%edge_delta(i, k), DELTA_EPS)
                  kn = mesh%neighbours(i, k)
                  if (kn > 0) then
                     dgdn = si_unstr_edge_flux(si%Hf_edge(i, k), si%eta(k), &
                                               si%eta(kn), inv_delta)
                  else
                     dgdn = si_unstr_edge_flux(si%Hf_edge(i, k), si%eta(k), &
                                               si%eta_ghost_edge(i, k), inv_delta)
                  end if
                  tau_n = (baro%wind_tau_x*mesh%edge_n1(i, k) &
                           + baro%wind_tau_y*mesh%edge_n2(i, k))/RHO_WATER
                  si%u_star(i, k) = si%u_n_edge(i, k) - g_omtheta_dt*dgdn &
                                    + dt*tau_n
               else
                  si%u_star(i, k) = 0.0_wp
                  si%u_n_edge(i, k) = 0.0_wp
               end if
            end do
         end do

         ! ---- 3. rhs_d = v(eta^n) - dt*div[blend] (+ Dirichlet fold) ----
         ! Reuse the linear RHS (leading eta^n + ghost fold via the frozen
         ! Hf) then swap the leading term for v(eta^n).
         call si_unstr_compute_rhs(mesh, nct, nco, si%eta, si%u_star, &
                                   si%u_n_edge, si%eta_r, dt, theta, &
                                   si%Hf_edge, si%eta_ghost_edge, g_theta2_dt2)
         call si_unstr_wetdry_volume(si%n_sub, nct, nco, si%eta, si%b_sub, &
                                     si%v_cell)
         do concurrent(k=1:nco)
            si%rhs_d(k) = si%eta_r(k) + (si%v_cell(k) - si%eta(k))
         end do

         ! ---- 4. Newton outer loop ----
         rhs2 = 0.0_wp
         do concurrent(k=1:nco) reduce(+:rhs2)
            rhs2 = rhs2 + si%rhs_d(k)*si%rhs_d(k)
         end do
         tol2 = max(SI_UNSTR_NEWTON_TOL*SI_UNSTR_NEWTON_TOL*rhs2, 1.0e-60_wp)

         cg_it_total = 0
         cg_it = 0
         mass_res = 0.0_wp
         wetcount_prev = -1
         converged = .false.
         do nit = 1, SI_UNSTR_NEWTON_MAX
            call si_unstr_wetdry_wetfrac(si%n_sub, nct, nco, si%eta, si%b_sub, &
                                         si%a_wet)
            call si_unstr_wetdry_volume(si%n_sub, nct, nco, si%eta, si%b_sub, &
                                        si%v_cell)
            call si_unstr_wetdry_mass_diag(nct, nco, si%a_wet, si%Hf_edge, &
                                           si%mass_diag)
            call si_unstr_apply_helmholtz(mesh, nct, nco, si%eta, si%newton_Jeta, &
                                          si%Hf_edge, g_theta2_dt2, &
                                          mass_diag=si%mass_diag)
            call si_unstr_wetdry_newton_residual(nct, nco, si%eta, si%newton_Jeta, &
                                                 si%mass_diag, si%v_cell, &
                                                 si%rhs_d, si%newton_F, f2)
            call si_unstr_wetdry_wetcount(si%n_sub, nct, nco, si%eta, si%b_sub, &
                                          wetcount)
            if (nit > 1 .and. wetcount == wetcount_prev .and. f2 < tol2) then
               converged = .true.
               exit
            end if
            wetcount_prev = wetcount

            ! Solve J*dz = -F from zero (full-array zero so the ghost band
            ! of the search space stays 0); eta += dz.
            do concurrent(k=1:nct)
               si%newton_dz(k) = 0.0_wp
            end do
            do concurrent(k=1:nco)
               si%eta_r(k) = -si%newton_F(k)
            end do
            call si_unstr_cg_solve(mesh, nct, nco, si%newton_dz, si%eta_r, &
                                   si%eta_p, si%eta_Ap, si%eta_z, si%eta_diag, &
                                   si%Hf_edge, g_theta2_dt2, si%precond, &
                                   si%rtol, cg_it, mass_res, &
                                   mass_diag=si%mass_diag)
            cg_it_total = cg_it_total + cg_it
            do concurrent(k=1:nco)
               si%eta(k) = si%eta(k) + si%newton_dz(k)
            end do
         end do

         si%newton_iters = min(nit, SI_UNSTR_NEWTON_MAX)
         si%mass_res = mass_res
         if (.not. converged) then
            call global_logger%warning("si_unstr_advance_wetdry: Newton did "// &
                                       "not converge in the iteration cap; last ||F||^2 above tol.")
         end if

         ! ---- 5. EDGE back-substitution with the SAME frozen Hf (R2) ----
         do concurrent(k=1:nco) local(i, kn, inv_delta, dgdn)
            do i = 1, EDGES_PER_TRI
               if (si%Hf_edge(i, k) > 0.0_wp) then
                  inv_delta = 1.0_wp/max(mesh%edge_delta(i, k), DELTA_EPS)
                  kn = mesh%neighbours(i, k)
                  if (kn > 0) then
                     dgdn = si_unstr_edge_flux(si%Hf_edge(i, k), si%eta(k), &
                                               si%eta(kn), inv_delta)
                  else
                     dgdn = si_unstr_edge_flux(si%Hf_edge(i, k), si%eta(k), &
                                               si%eta_ghost_edge(i, k), inv_delta)
                  end if
                  si%u_edge(i, k) = si%u_star(i, k) - gtheta_dt*dgdn
               else
                  si%u_edge(i, k) = 0.0_wp
               end if
            end do
         end do

         ! ---- 6. Cell momenta + thickness (dry cells clamp h to 0) ----
         call si_unstr_reconstruct_cell_momenta(mesh, nct, nco, si%u_edge, &
                                                baro%xmom_cv, baro%ymom_cv)
         do concurrent(k=1:nco)
            baro%stage_cv(k) = si%eta(k)
            baro%height_cv(k) = max(si%eta(k) - baro%bed_cv(k), 0.0_wp)
         end do
      end associate

      iter_out = cg_it_total
   end subroutine si_unstr_advance_wetdry

   pure subroutine si_unstr_reconstruct_cell_momenta(mesh, nct, nco, u_edge, &
                                                     xmom, ymom)
      !! Least-squares reconstruction of the cell-centre transport vector
      !! (hu, hv) from the three edge-normal transports U_e = (hu,hv)·n_e:
      !!   M = sum_e len_e * n_e n_e^T,   rhs = sum_e len_e * U_e * n_e,
      !!   (hu, hv) = M^{-1} rhs.
      !! Exact for a constant transport field (the three normals span R^2 on
      !! a non-degenerate triangle => M SPD, invertible).  The KNP seam only
      !! needs a diagnostic cell momentum; eta/mass already closed on the
      !! edges.  Shared by the 2D step and (per layer) the ML step.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total) — the explicit-shape dim.
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: u_edge(EDGES_PER_TRI, nct)
         !! Edge-normal transports (m^2/s), + along the outward normal.
      real(wp), intent(out) :: xmom(nct), ymom(nct)
         !! Reconstructed cell transport vector (m^2/s).

      integer :: k, i
      real(wp) :: m11, m12, m22, rx, ry, det, wlen, enx, eny, ue

      do concurrent(k=1:nco) &
         local(i, m11, m12, m22, rx, ry, det, wlen, enx, eny, ue)
         m11 = 0.0_wp
         m12 = 0.0_wp
         m22 = 0.0_wp
         rx = 0.0_wp
         ry = 0.0_wp
         do i = 1, EDGES_PER_TRI
            wlen = mesh%edge_len(i, k)
            enx = mesh%edge_n1(i, k)
            eny = mesh%edge_n2(i, k)
            ue = u_edge(i, k)
            m11 = m11 + wlen*enx*enx
            m12 = m12 + wlen*enx*eny
            m22 = m22 + wlen*eny*eny
            rx = rx + wlen*ue*enx
            ry = ry + wlen*ue*eny
         end do
         det = m11*m22 - m12*m12
         ! det ~ 0 only for a degenerate cell (collinear edge normals) —
         ! excluded by the fail-loud mesh quality gate at enable_si; the
         ! zero fallback is 1/0 armour for a diagnostic-only quantity
         ! (mass/eta close on the edges regardless).
         if (abs(det) > tiny(1.0_wp)) then
            xmom(k) = (m22*rx - m12*ry)/det
            ymom(k) = (m11*ry - m12*rx)/det
         else
            xmom(k) = 0.0_wp
            ymom(k) = 0.0_wp
         end if
      end do
   end subroutine si_unstr_reconstruct_cell_momenta

end module rki_si_unstr_step
