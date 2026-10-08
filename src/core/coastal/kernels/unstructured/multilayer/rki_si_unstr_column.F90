!! Per-column + per-edge-layer kernels for the unstructured 3-D SI step.
module rki_si_unstr_column
   !! The SCHISM-inspired (Casulli 1990; Casulli & Cattani 1994; Zhang &
   !! Baptista 2008) per-column algebra of the 3D semi-implicit multilayer
   !! step, ported per-cell to the `(nz, ncells)` axis layout (k FIRST,
   !! cell LAST — the unstructured convention, OPPOSITE of the structured
   !! `(nx,ny,nz)`), plus the edge-wise recovery / tracer kernels that
   !! replace the structured face kernels:
   !!
   !! 1. `si_unstr_build_nu_interface` — PP81 Richardson interface viscosity.
   !! 2. `si_unstr_column_solve` — plain-u^n predictor responses G_k/R_k
   !!    (face/edge-based PGF split: the predictor carries NO surface PGF),
   !!    Thomas factorisation, three-RHS back-solve, Hhat/Gbar integrals.
   !!    Bottom drag: Manning quadratic on the bed row (the unstructured
   !!    barotropic state carries no linear-drag knob — Manning-only,
   !!    documented deviation from the structured `use_linear_drag` option).
   !! 3. `si_unstr_fill_predictor_edges` — the theta-blend RHS edge
   !!    predictor U*_e from the face-averaged Gbar (vector projected on
   !!    the outward edge-normal) minus the explicit (1-theta) surface PGF.
   !! 4. `si_unstr_recover_edge_layers` — post-CG per-layer EDGE transport
   !!    recovery; column sums telescope EXACTLY to the 2D back-substituted
   !!    edge transport (the mass-consistency invariant).  Dirichlet
   !!    boundary edges use the structured mouth-consistency seeds inline:
   !!    ghost G = 0, ghost R = -1, ghost h_k = dsig_k*max(eta_bc - b_c, 0).
   !! 5. `si_unstr_fill_mass_flux_edges` — the blended column edge flux the
   !!    eta update used (drives the dsig tracer split + next step's U^n).
   !! 6. `si_unstr_advect_tracer_dsig` — depth-proportional, race-free,
   !!    first-order-upwind tracer advection: per-layer edge flux EXACTLY
   !!    `dsig_k * mass_flux_e` (the only split consistent with the
   !!    dsig-slaved thickness update => uniform tracers invariant per
   !!    layer); upwind values read from the PRE-STEP snapshot only.
   !!    FLAT-IMPL: called per registry tracer through the outer shim in
   !!    the step module (the array-of-derived-types device rule).
   !! 7. `si_unstr_edge_layers_to_cells` — per-layer least-squares A-grid
   !!    cell-momentum reconstruction (the 2D reconstruction, per layer).
   !! 8. `si_unstr_update_thickness_ml` — h^{n+1} = eta^{n+1} - b and the
   !!    sigma slaving h_k = dsig_k*h.
   !! 9. `si_unstr_apply_wind_top_layer` — dt*tau/rho on the surface layer
   !!    (k = nz) transport; the uniform-scalar unstructured wind.
   !!
   !! Layer convention (load-bearing): BOTTOM-UP — k=1 bed, k=nz surface.
   !! Cell-centric double-visit throughout (each cell owns its 3 edge
   !! values; interior edges computed twice, bit-exactly antisymmetric —
   !! no scatter, no atomics).  Explicit-shape dummies in every DC kernel
   !! (NVHPC descriptor-walk rule); fixed-size local columns `NZ_STACK_MAX`.
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, GRAVITY, H_VANISHED, DRY_TOLERANCE, &
                            RHO_WATER
#else
   use rki_constants, only: NZ_STACK_MAX, wp, GRAVITY, H_VANISHED, DRY_TOLERANCE, &
                            RHO_WATER
#endif
   use rki_mesh_unstr, only: mesh_unstr_t
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: si_unstr_build_nu_interface
   public :: si_unstr_column_solve
   public :: si_unstr_fill_predictor_edges
   public :: si_unstr_recover_edge_layers
   public :: si_unstr_fill_mass_flux_edges
   public :: si_unstr_advect_tracer_dsig
   public :: si_unstr_edge_layers_to_cells
   public :: si_unstr_update_thickness_ml
   public :: si_unstr_apply_wind_top_layer

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.
   real(wp), parameter :: DELTA_EPS = 1.0e-30_wp
      !! 1/0 armour on the circumcentre-line distance (gated mesh => never
      !! fires).

contains

   subroutine si_unstr_build_nu_interface(nz, nzp1, nct, nco, h_layer, &
                                               hu_layer, hv_layer, rho_layer, &
                                               rho_0, nu_max, nu_bg, ri_alpha, &
                                               nu_iface)
      !! PP81 Richardson-number interface viscosity, per cell (mirror of the
      !! structured `si_build_nu_interface`, `(nz, ncells)` order):
      !!   nu_iface(1, c)      = 0   (bed interface)
      !!   nu_iface(k+1, c)    = nu between layers k and k+1
      !!   nu_iface(nz+1, c)   = 0   (surface interface)
      !! Ri = g*(drho/rho_0)*h_avg / |du|^2 (floored at 0);
      !! nu = nu_max/(1 + ri_alpha*Ri)^2 + nu_bg; vanished layers decouple.
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nzp1
         !! nz + 1 (the interface-array leading dim, explicit-shape).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! Layer thicknesses [m], bottom-up (k=1 bed).
      real(wp), intent(in) :: hu_layer(nz, nct)
         !! Layer x-transports [m^2/s].
      real(wp), intent(in) :: hv_layer(nz, nct)
         !! Layer y-transports [m^2/s].
      real(wp), intent(in) :: rho_layer(nz, nct)
         !! Layer densities [kg/m^3].
      real(wp), intent(in) :: rho_0
         !! Reference density [kg/m^3].
      real(wp), intent(in) :: nu_max
         !! Maximum (convective) eddy viscosity [m^2/s].
      real(wp), intent(in) :: nu_bg
         !! Background viscosity floor [m^2/s].
      real(wp), intent(in) :: ri_alpha
         !! Richardson tuning coefficient (PP81: 5.0).
      real(wp), intent(out) :: nu_iface(nzp1, nct)
         !! Interface viscosity [m^2/s], vdiff convention.

      integer :: c, k
      real(wp) :: hk, hkp1, h_avg, uk, ukp1, vk, vkp1
      real(wp) :: du, dv, drho, shear2, ri_num, denom, nu_t

      !$omp parallel do private(k, hk, hkp1, h_avg, uk, ukp1, vk, vkp1, du, dv, drho, shear2, ri_num, denom, nu_t)
      do c = 1, nco
         nu_iface(1, c) = 0.0_wp
         nu_iface(nz + 1, c) = 0.0_wp
         do k = 1, nz - 1
            hk = h_layer(k, c)
            hkp1 = h_layer(k + 1, c)
            if (hk <= H_VANISHED .or. hkp1 <= H_VANISHED) then
               nu_iface(k + 1, c) = 0.0_wp
            else
               h_avg = 0.5_wp*(hk + hkp1)
               uk = hu_layer(k, c)/hk
               ukp1 = hu_layer(k + 1, c)/hkp1
               vk = hv_layer(k, c)/hk
               vkp1 = hv_layer(k + 1, c)/hkp1
               du = ukp1 - uk
               dv = vkp1 - vk
               drho = rho_layer(k, c) - rho_layer(k + 1, c)
               shear2 = du*du + dv*dv + 1.0e-20_wp
               ri_num = GRAVITY*(drho/rho_0)*h_avg/shear2
               if (ri_num < 0.0_wp) ri_num = 0.0_wp
               denom = 1.0_wp + ri_alpha*ri_num
               nu_t = nu_max/(denom*denom) + nu_bg
               nu_iface(k + 1, c) = nu_t
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_build_nu_interface

   subroutine si_unstr_column_solve(nz, nzp1, nct, nco, h_layer, hu_layer, &
                                         hv_layer, nu_iface, manning_n, dt, &
                                         g_col_u, g_col_v, r_col, Hhat_cell, &
                                         Gbar_cell_x, Gbar_cell_y)
      !! Per-column semi-implicit solve — the structured `si_column_solve`
      !! ported per-cell (`(nz, ncells)` order).  Predictor: PLAIN old-time
      !! velocity u^n_k = hu_k/max(h_k, H_VANISHED) — NO surface PGF (the
      !! edge-based PGF split: by linearity of the tridiagonal, all surface
      !! PGF applies at edges via the -1 response R_k).  Bottom drag:
      !! Manning quadratic on the bed row (the unstructured barotropic
      !! state carries no linear-drag knob):
      !!   cf = g*n^2 / H_total^{1/3};  lambda = cf*|u_bed| / h_1
      !! folded into the bed diagonal as +dt*lambda.  Tridiagonal, Thomas
      !! factorisation, three back-solves (u^n, v^n, -1), then
      !!   Gbar = sum_k h_k G_k,  Hhat = -sum_k h_k R_k  (>= 0).
      !! R_k is axis-independent (u and v share the matrix).
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nzp1
         !! nz + 1 (the interface-array leading dim, explicit-shape).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! Layer thicknesses [m], bottom-up.
      real(wp), intent(in) :: hu_layer(nz, nct)
         !! Layer x-transports [m^2/s] (predictor tendencies already folded).
      real(wp), intent(in) :: hv_layer(nz, nct)
         !! Layer y-transports [m^2/s].
      real(wp), intent(in) :: nu_iface(nzp1, nct)
         !! Interface viscosity [m^2/s] (vdiff convention).
      real(wp), intent(in) :: manning_n
         !! Manning roughness n [s/m^{1/3}]; <= 0 disables bottom drag.
      real(wp), intent(in) :: dt
         !! Time step [s].
      real(wp), intent(out) :: g_col_u(nz, nct)
         !! Column velocity response to u^n (G_k^u) [m/s].
      real(wp), intent(out) :: g_col_v(nz, nct)
         !! Column velocity response to v^n (G_k^v) [m/s].
      real(wp), intent(out) :: r_col(nz, nct)
         !! Column velocity response to -1 (R_k) [m/s] (<= 0).
      real(wp), intent(out) :: Hhat_cell(nct)
         !! Effective column depth Hhat = -sum_k h_k R_k [m] (>= 0).
      real(wp), intent(out) :: Gbar_cell_x(nct)
         !! Column-integrated predictor transport sum_k h_k G_k^u [m^2/s].
      real(wp), intent(out) :: Gbar_cell_y(nct)
         !! Column-integrated y predictor transport [m^2/s].

      integer :: c, k
      real(wp) :: hk_col(NZ_STACK_MAX)
      real(wp) :: ak_col(NZ_STACK_MAX)
      real(wp) :: bk_col(NZ_STACK_MAX)
      real(wp) :: ck_col(NZ_STACK_MAX)
      real(wp) :: rhs_u(NZ_STACK_MAX)
      real(wp) :: rhs_v(NZ_STACK_MAX)
      real(wp) :: rhs_r(NZ_STACK_MAX)
      real(wp) :: hc, hm, hp, dz_bot, dz_top
      real(wp) :: alph, beta_k, denom_k
      real(wp) :: lambda_drag, h_total, u_bed, v_bed, cf
      real(wp) :: gx_sum, gy_sum, hhat_sum

      !$omp parallel do private(k, hk_col, ak_col, bk_col, ck_col, rhs_u, rhs_v, rhs_r, hc, hm, hp, dz_bot, dz_top, alph, beta_k, denom_k, lambda_drag, h_total, u_bed, v_bed, cf, gx_sum, gy_sum, hhat_sum)
      do c = 1, nco

         ! ---- Floored thicknesses + plain u^n/v^n predictors (no PGF) ----
         do k = 1, nz
            hc = max(h_layer(k, c), H_VANISHED)
            hk_col(k) = hc
            rhs_u(k) = hu_layer(k, c)/hc
            rhs_v(k) = hv_layer(k, c)/hc
            rhs_r(k) = -1.0_wp
         end do

         ! ---- Manning bottom drag on the bed row ----
         lambda_drag = 0.0_wp
         if (manning_n > 0.0_wp) then
            h_total = 0.0_wp
            do k = 1, nz
               h_total = h_total + hk_col(k)
            end do
            cf = GRAVITY*manning_n*manning_n/(h_total**(1.0_wp/3.0_wp))
            u_bed = hu_layer(1, c)/hk_col(1)
            v_bed = hv_layer(1, c)/hk_col(1)
            lambda_drag = cf*sqrt(u_bed*u_bed + v_bed*v_bed)/hk_col(1)
         end if

         ! ---- Tridiagonal build (bottom-up; vanished layers decouple) ----
         hc = hk_col(1)
         if (nz > 1) then
            hp = hk_col(2)
            dz_top = 0.5_wp*(hc + hp)
            alph = dt*nu_iface(2, c)/(hc*dz_top)
            if (h_layer(1, c) <= H_VANISHED .or. h_layer(2, c) <= H_VANISHED) then
               alph = 0.0_wp
            end if
         else
            alph = 0.0_wp
         end if
         ak_col(1) = 0.0_wp
         ck_col(1) = -alph
         bk_col(1) = 1.0_wp + alph + dt*lambda_drag

         do k = 2, nz - 1
            hc = hk_col(k)
            hm = hk_col(k - 1)
            hp = hk_col(k + 1)
            dz_bot = 0.5_wp*(hm + hc)
            dz_top = 0.5_wp*(hc + hp)
            beta_k = dt*nu_iface(k, c)/(hc*dz_bot)
            alph = dt*nu_iface(k + 1, c)/(hc*dz_top)
            if (h_layer(k, c) <= H_VANISHED .or. h_layer(k - 1, c) <= H_VANISHED) then
               beta_k = 0.0_wp
            end if
            if (h_layer(k, c) <= H_VANISHED .or. h_layer(k + 1, c) <= H_VANISHED) then
               alph = 0.0_wp
            end if
            ak_col(k) = -beta_k
            ck_col(k) = -alph
            bk_col(k) = 1.0_wp + alph + beta_k
         end do

         if (nz > 1) then
            hc = hk_col(nz)
            hm = hk_col(nz - 1)
            dz_bot = 0.5_wp*(hm + hc)
            beta_k = dt*nu_iface(nz, c)/(hc*dz_bot)
            if (h_layer(nz, c) <= H_VANISHED .or. h_layer(nz - 1, c) <= H_VANISHED) then
               beta_k = 0.0_wp
            end if
            ak_col(nz) = -beta_k
            ck_col(nz) = 0.0_wp
            bk_col(nz) = 1.0_wp + beta_k
         end if

         ! ---- Thomas factorisation + three-RHS sweep ----
         ck_col(1) = ck_col(1)/bk_col(1)
         do k = 2, nz
            denom_k = bk_col(k) - ak_col(k)*ck_col(k - 1)
            ck_col(k) = ck_col(k)/denom_k
            bk_col(k) = denom_k
         end do
         rhs_u(1) = rhs_u(1)/bk_col(1)
         rhs_v(1) = rhs_v(1)/bk_col(1)
         rhs_r(1) = rhs_r(1)/bk_col(1)
         do k = 2, nz
            rhs_u(k) = (rhs_u(k) - ak_col(k)*rhs_u(k - 1))/bk_col(k)
            rhs_v(k) = (rhs_v(k) - ak_col(k)*rhs_v(k - 1))/bk_col(k)
            rhs_r(k) = (rhs_r(k) - ak_col(k)*rhs_r(k - 1))/bk_col(k)
         end do
         ! Back-substitution (in place: rhs_* become the solutions x_*).
         do k = nz - 1, 1, -1
            rhs_u(k) = rhs_u(k) - ck_col(k)*rhs_u(k + 1)
            rhs_v(k) = rhs_v(k) - ck_col(k)*rhs_v(k + 1)
            rhs_r(k) = rhs_r(k) - ck_col(k)*rhs_r(k + 1)
         end do

         ! ---- Cell integrals + store responses ----
         ! RAW h_layer in the integrals (NOT the floored hk_col): the floor's
         ! only job is arming the tridiagonal 1/h — a vanished layer carries
         ! no transport, so it must contribute 0 to the column conveyance.
         ! The per-layer recovery telescopes with raw h_k*R_k / h_k*G_k, so
         ! floored integrals would make the eta-solve's face conveyance
         ! DISAGREE with the recovered transports by up to nz*H_VANISHED at
         ! dry/thin columns — measured as a ~1e-3 m/step mass<->tracer
         ! closure breach at the wet/dry front (water moving without salt,
         ! ~1.7 PSU freshening over a flood; the Q8 ML+wet/dry lesson).
         ! For any layer with h > H_VANISHED raw == floored, so columns where
         ! EVERY layer exceeds H_VANISHED are bit-identical (a genuinely wet
         ! but sub-H_VANISHED-thin layer is wet yet floored — NOT
         ! bit-identical).
         gx_sum = 0.0_wp
         gy_sum = 0.0_wp
         hhat_sum = 0.0_wp
         do k = 1, nz
            gx_sum = gx_sum + h_layer(k, c)*rhs_u(k)
            gy_sum = gy_sum + h_layer(k, c)*rhs_v(k)
            hhat_sum = hhat_sum - h_layer(k, c)*rhs_r(k)
            g_col_u(k, c) = rhs_u(k)
            g_col_v(k, c) = rhs_v(k)
            r_col(k, c) = rhs_r(k)
         end do
         Gbar_cell_x(c) = gx_sum
         Gbar_cell_y(c) = gy_sum
         Hhat_cell(c) = hhat_sum
      end do
      !$omp end parallel do
   end subroutine si_unstr_column_solve

   subroutine si_unstr_fill_predictor_edges(mesh, nz, nct, nco, &
                                                 Gbar_cell_x, Gbar_cell_y, &
                                                 hu_layer0, hv_layer0, &
                                                 u_n_edge, eta, Hf_edge, &
                                                 eta_ghost_edge, u_star, &
                                                 g_omtheta_dt)
      !! EDGE predictor for the theta-blend RHS — ANCHOR + INCREMENT form:
      !!
      !!   U*_e = U^n_e + 0.5*((Gbar - X0)_c + (Gbar - X0)_nbr)·n_e
      !!               - g(1-theta)*dt*Hf*(deta^n/dn)
      !!
      !! where U^n_e is the PERSISTED per-edge column transport (the 2D
      !! prognostic — full edge-local fidelity, no cell round-trip) and
      !! `Gbar - X0` is the column-integrated IMPLICIT-MIXING + explicit-
      !! tendency INCREMENT (X0 = sum_k hu0_k, the PRE-tendency saved
      !! transport), cell-averaged onto the edge.  This differs from the
      !! structured `si_ml_fill_predictor_faces` (which cell-averages the
      !! WHOLE Gbar — consistent there because the structured 2D predictor
      !! also rebuilds faces from cell momenta) precisely so that the
      !! uniform-density inviscid limit (Gbar == X0) reduces EXACTLY to the
      !! edge-persistent 2D unstructured predictor — the baroclinic-limit
      !! gate.  Mixing / drag / wind / BPG enter through the increment with
      !! the same cell-average quality as the structured path.
      !!
      !! GATED on Hf > 0 (dead edges exactly 0; U^n_e already dead-zeroed,
      !! R3).  Dirichlet boundary edge: ghost increment = 0 (the Dirichlet
      !! DOF has no solved column response and no saved momentum) and the
      !! eta gradient reads the seeded ghost surface.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: Gbar_cell_x(nct), Gbar_cell_y(nct)
      real(wp), intent(in) :: hu_layer0(nz, nct), hv_layer0(nz, nct)
         !! PRE-tendency saved layer transports (the X0 columns).
      real(wp), intent(in) :: u_n_edge(EDGES_PER_TRI, nct)
         !! Persisted U^n edge transports (dead edges zeroed).
      real(wp), intent(in) :: eta(nct)
         !! eta^n at circumcentres.
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(in) :: eta_ghost_edge(EDGES_PER_TRI, nct)
      real(wp), intent(out) :: u_star(EDGES_PER_TRI, nct)
      real(wp), intent(in) :: g_omtheta_dt
         !! g*(1-theta)*dt.

      integer :: c, i, k, kn
      real(wp) :: inv_delta, dgdn, inc_n, enx, eny
      real(wp) :: x0x_c, x0y_c, x0x_n, x0y_n

      !$omp parallel do private(i, k, kn, inv_delta, dgdn, inc_n, enx, eny, x0x_c, x0y_c, x0x_n, x0y_n)
      do c = 1, nco
         ! Own-cell pre-tendency column transport X0 (hoisted out of edges).
         x0x_c = 0.0_wp
         x0y_c = 0.0_wp
         do k = 1, nz
            x0x_c = x0x_c + hu_layer0(k, c)
            x0y_c = x0y_c + hv_layer0(k, c)
         end do
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, c) > 0.0_wp) then
               inv_delta = 1.0_wp/max(mesh%edge_delta(i, c), DELTA_EPS)
               enx = mesh%edge_n1(i, c)
               eny = mesh%edge_n2(i, c)
               kn = mesh%neighbours(i, c)
               if (kn > 0) then
                  x0x_n = 0.0_wp
                  x0y_n = 0.0_wp
                  do k = 1, nz
                     x0x_n = x0x_n + hu_layer0(k, kn)
                     x0y_n = x0y_n + hv_layer0(k, kn)
                  end do
                  inc_n = 0.5_wp*(((Gbar_cell_x(c) - x0x_c) &
                                   + (Gbar_cell_x(kn) - x0x_n))*enx &
                                  + ((Gbar_cell_y(c) - x0y_c) &
                                     + (Gbar_cell_y(kn) - x0y_n))*eny)
                  dgdn = Hf_edge(i, c)*(eta(kn) - eta(c))*inv_delta
               else
                  ! Dirichlet edge: ghost increment 0; gradient vs ghost eta.
                  inc_n = 0.5_wp*((Gbar_cell_x(c) - x0x_c)*enx &
                                  + (Gbar_cell_y(c) - x0y_c)*eny)
                  dgdn = Hf_edge(i, c)*(eta_ghost_edge(i, c) - eta(c))*inv_delta
               end if
               u_star(i, c) = u_n_edge(i, c) + inc_n - g_omtheta_dt*dgdn
            else
               u_star(i, c) = 0.0_wp
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_fill_predictor_edges

   subroutine si_unstr_recover_edge_layers(mesh, nz, nct, nco, h_layer, &
                                                hu_layer0, hv_layer0, dsig, &
                                                bed, g_col_u, g_col_v, &
                                                r_col, eta_n, eta_np1, Hf_edge, &
                                                eta_ghost_edge, dt, theta, &
                                                u_edge_layer)
      !! Recover the per-layer EDGE transports after the eta solve — the
      !! structured `si_recover_face_transports` in the ANCHOR + INCREMENT
      !! form matching `si_unstr_fill_predictor_edges`.  Per edge i of cell
      !! c (gated Hf > 0; dead edges stay exactly 0) and layer k:
      !!
      !!   anchor_k = u_edge_layer(k,i,c)                (the PERSISTED
      !!              per-layer edge transport — read before overwrite;
      !!              own-entry read/write only, race-free)
      !!   inc_k = 0.5*((h_k*G_k - hu0_k)_c·n_e + (h_k*G_k - hu0_k)_nbr·n_e)
      !!   hhf_k = 0.5*(h_k(c)*R_k(c) + h_k(nbr)*R_k(nbr))          (<= 0)
      !!   grad_blend = ((1-theta)*deta^n + theta*deta^{n+1})/delta_e
      !!   u_edge_layer(k,i,c) = anchor_k + inc_k + g*dt*grad_blend*hhf_k
      !!
      !! MASS-CONSISTENCY INVARIANT: sum_k anchor_k = U^n_e (the persisted
      !! column total), sum_k inc_k = the predictor's increment, and
      !! sum_k hhf_k = -Hf, so the column sum telescopes EXACTLY to
      !! U*_e - g*theta*dt*Hf*grad(eta^{n+1}) — the 2D back-substitution.
      !! In the uniform-density inviscid limit inc_k = 0 and
      !! hhf_k = -dsig_k*Hf, so every layer is EXACTLY dsig_k times the 2D
      !! edge update (the baroclinic-limit gate); sheared columns keep
      !! their per-layer structure in the anchor (no cell-average
      !! flattening — the reason for this form over the structured pure
      !! cell-average, whose 2D twin averages cells too).  `h_layer` MUST
      !! be the SAME time-n field that fed the tridiagonal.
      !!
      !! DIRICHLET boundary edge (the structured estuary-mouth lesson,
      !! applied inline): the ghost is a Dirichlet DOF with NO solved
      !! column response and no saved momentum, so ghost increment = 0,
      !! ghost R = -1 exactly, and the ghost per-layer thickness is the
      !! sigma split of the prescribed column h_ghost = max(eta_bc - b_c,
      !! 0): h_k_ghost = dsig_k*h_ghost.  Then sum_k hhf = -0.5*(Hhat_c +
      !! h_ghost) = -Hf exactly — the recovered open-edge transport matches
      !! the Helmholtz flux, so continuity and tracers share the mouth flux
      !! and uniform tracers stay invariant (the freshening-defect guard).
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! Layer thicknesses [m] at time n (the tridiagonal's field).
      real(wp), intent(in) :: hu_layer0(nz, nct), hv_layer0(nz, nct)
         !! PRE-tendency saved layer transports (the increment reference).
      real(wp), intent(in) :: dsig(nz)
         !! Sigma layer fractions (sum = 1).
      real(wp), intent(in) :: bed(nct)
         !! Bed elevation [m] (the Dirichlet ghost column depth).
      real(wp), intent(in) :: g_col_u(nz, nct), g_col_v(nz, nct)
      real(wp), intent(in) :: r_col(nz, nct)
      real(wp), intent(in) :: eta_n(nct)
         !! Old-time surface eta^n at circumcentres.
      real(wp), intent(in) :: eta_np1(nct)
         !! New-time surface eta^{n+1} (the CG solution).
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(in) :: eta_ghost_edge(EDGES_PER_TRI, nct)
         !! Seeded Dirichlet ghost surface (frozen over the step — used for
         !! BOTH the eta^n and eta^{n+1} gradients, as on the structured
         !! path where the ghost persists across the interior-only CG).
      real(wp), intent(in) :: dt, theta
      real(wp), intent(inout) :: u_edge_layer(nz, EDGES_PER_TRI, nct)
         !! Per-layer edge transports [m^2/s], + outward normal: the
         !! persisted anchor on entry, the recovered U^{n+1} on exit.

      integer :: c, i, k, kn
      real(wp) :: inv_delta, grad_blend, inc_k, hhf, enx, eny
      real(wp) :: h_ghost, hk_ghost

      !$omp parallel do private(i, k, kn, inv_delta, grad_blend, inc_k, hhf, enx, eny, h_ghost, hk_ghost)
      do c = 1, nco
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, c) > 0.0_wp) then
               inv_delta = 1.0_wp/max(mesh%edge_delta(i, c), DELTA_EPS)
               enx = mesh%edge_n1(i, c)
               eny = mesh%edge_n2(i, c)
               kn = mesh%neighbours(i, c)
               if (kn > 0) then
                  grad_blend = ((1.0_wp - theta)*(eta_n(kn) - eta_n(c)) &
                                + theta*(eta_np1(kn) - eta_np1(c)))*inv_delta
                  do k = 1, nz
                     inc_k = 0.5_wp*((h_layer(k, c)*g_col_u(k, c) - hu_layer0(k, c) &
                                      + h_layer(k, kn)*g_col_u(k, kn) - hu_layer0(k, kn))*enx &
                                     + (h_layer(k, c)*g_col_v(k, c) - hv_layer0(k, c) &
                                        + h_layer(k, kn)*g_col_v(k, kn) - hv_layer0(k, kn))*eny)
                     hhf = 0.5_wp*(h_layer(k, c)*r_col(k, c) &
                                   + h_layer(k, kn)*r_col(k, kn))
                     u_edge_layer(k, i, c) = u_edge_layer(k, i, c) + inc_k &
                                             + GRAVITY*dt*grad_blend*hhf
                  end do
               else
                  ! Dirichlet edge: frozen ghost eta for both time levels;
                  ! ghost increment = 0, ghost R = -1, h_k = dsig_k*h_ghost.
                  grad_blend = ((1.0_wp - theta)*(eta_ghost_edge(i, c) - eta_n(c)) &
                                + theta*(eta_ghost_edge(i, c) - eta_np1(c)))*inv_delta
                  h_ghost = max(eta_ghost_edge(i, c) - bed(c), 0.0_wp)
                  do k = 1, nz
                     hk_ghost = dsig(k)*h_ghost
                     inc_k = 0.5_wp*((h_layer(k, c)*g_col_u(k, c) - hu_layer0(k, c))*enx &
                                     + (h_layer(k, c)*g_col_v(k, c) - hv_layer0(k, c))*eny)
                     hhf = 0.5_wp*(h_layer(k, c)*r_col(k, c) &
                                   + hk_ghost*(-1.0_wp))
                     u_edge_layer(k, i, c) = u_edge_layer(k, i, c) + inc_k &
                                             + GRAVITY*dt*grad_blend*hhf
                  end do
               end if
            else
               do k = 1, nz
                  u_edge_layer(k, i, c) = 0.0_wp
               end do
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_recover_edge_layers

   subroutine si_unstr_fill_mass_flux_edges(nz, nct, nco, u_edge_layer, &
                                                 u_n_edge, Hf_edge, theta, &
                                                 mass_flux_edge, u_edge)
      !! The blended column edge flux the eta update used:
      !!   F_e = theta*sum_k u_edge_layer_k + (1-theta)*U^n_e
      !! (gated Hf > 0; dead edges exactly 0) — feeds the dsig tracer split
      !! so continuity and tracers share edge fluxes exactly.  Also stores
      !! the NEW column-total edge transport U^{n+1}_e = sum_k u_edge_layer
      !! into `u_edge` (the 2D prognostic the next step's (1-theta) side
      !! reads — the telescoped 2D back-substitution).
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: u_edge_layer(nz, EDGES_PER_TRI, nct)
      real(wp), intent(in) :: u_n_edge(EDGES_PER_TRI, nct)
         !! Saved U^n edge transports (dead edges already zeroed, R3).
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(in) :: theta
      real(wp), intent(out) :: mass_flux_edge(EDGES_PER_TRI, nct)
         !! Blended column edge flux [m^2/s].
      real(wp), intent(out) :: u_edge(EDGES_PER_TRI, nct)
         !! New column-total edge transport [m^2/s].

      integer :: c, i, k
      real(wp) :: fsum

      !$omp parallel do private(i, k, fsum)
      do c = 1, nco
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, c) > 0.0_wp) then
               fsum = 0.0_wp
               do k = 1, nz
                  fsum = fsum + u_edge_layer(k, i, c)
               end do
               u_edge(i, c) = fsum
               mass_flux_edge(i, c) = theta*fsum + (1.0_wp - theta)*u_n_edge(i, c)
            else
               u_edge(i, c) = 0.0_wp
               mass_flux_edge(i, c) = 0.0_wp
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_fill_mass_flux_edges

   subroutine si_unstr_advect_tracer_dsig(mesh, nz, nct, nco, n_bc, dsig, &
                                               h_layer, mass_flux_edge, &
                                               tr_bc_ext, hTr_old, hTr, dt)
      !! Depth-proportional, race-free, first-order-upwind tracer advection
      !! for the sigma SI path — the structured `si_ml_advect_tracer_dsig`,
      !! edge-wise (cell-centric double-visit).  The per-layer edge flux is
      !! EXACTLY `dsig_k * mass_flux_e * len_e * T_upwind` — the only split
      !! consistent with the dsig-slaved thickness update, so a uniform
      !! tracer stays uniform PER LAYER to CG tolerance and total tracer is
      !! conserved to round-off (interior edge fluxes are bit-exactly
      !! antisymmetric between the two visits; dead edges carry 0).
      !!
      !! Upwind values read the PRE-STEP snapshot `hTr_old` + the pre-update
      !! `h_layer` ONLY (no read-neighbour/write-own hazard — the explicit
      !! path's in-place update leaks conservation on sequential backends at
      !! SI step sizes).  Outflow (flux > 0, leaving the cell) upwinds the
      !! OWN cell; inflow through a Dirichlet boundary edge carries the
      !! per-group external tracer `tr_bc_ext(group)` (the clamped-inflow
      !! convention of the explicit unstructured path).
      !!
      !! FLAT-IMPL for the registry outer shim: `hTr_old`/`hTr` are plain
      !! explicit-shape arrays — the step module loops `tracers(:)` on the
      !! host and passes `tracers(it)%hTr` etc. (array-of-derived-types
      !! device rule).
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      integer, intent(in) :: n_bc
         !! BC group table size.
      real(wp), intent(in) :: dsig(nz)
         !! Sigma layer fractions.
      real(wp), intent(in) :: h_layer(nz, nct)
         !! PRE-update layer thicknesses (time n) — upwind T denominators.
      real(wp), intent(in) :: mass_flux_edge(EDGES_PER_TRI, nct)
         !! Blended column edge flux [m^2/s], + outward.
      real(wp), intent(in) :: tr_bc_ext(n_bc)
         !! Per-group external tracer value for Dirichlet inflow.
      real(wp), intent(in) :: hTr_old(nz, nct)
         !! Pre-step tracer snapshot (upwind source; never written here).
      real(wp), intent(inout) :: hTr(nz, nct)
         !! Tracer field updated in place (each cell writes only itself).
      real(wp), intent(in) :: dt

      integer :: c, i, k, kn, g
      real(wp) :: flux_c, tr_up, h_up, h_dn, div_sum, inv_area

      ! WET-SIDE FALLBACK at dry edges (the Q8 ML+wet/dry front lesson): a
      ! freshly-flooded edge can carry a transient REVERSED theta-blend flux
      ! whose upwind layer was EMPTY at time n (the just-wetted cell sloshing
      ! back).  A bare tr_up = 0 there moves real water with no tracer —
      ! measured ~1.7 PSU freshening over a flood tide.  When the upwind
      ! layer is empty the water physically originates from the SAME-step
      ! inflow off the wet side, so the DOWNWIND side's concentration is the
      ! consistent value.  Both double-visits of the edge then read the SAME
      ! (wet-side) tr_up — antisymmetry and conservation hold exactly, a
      ! uniform tracer stays uniform, and any resolved upwind layer
      ! (h > DRY_TOLERANCE) is bit-identical to plain upwinding.
      !$omp parallel do private(i, k, kn, g, flux_c, tr_up, h_up, h_dn, div_sum, inv_area)
      do c = 1, nco
         inv_area = 1.0_wp/mesh%area(c)
         do k = 1, nz
            div_sum = 0.0_wp
            do i = 1, EDGES_PER_TRI
               flux_c = dsig(k)*mass_flux_edge(i, c)*mesh%edge_len(i, c)
               if (flux_c > 0.0_wp) then
                  ! Outflow: upwind the OWN cell; empty own => the downwind
                  ! (neighbour / Dirichlet reservoir) concentration.
                  h_up = h_layer(k, c)
                  if (h_up > DRY_TOLERANCE) then
                     tr_up = hTr_old(k, c)/h_up
                  else
                     kn = mesh%neighbours(i, c)
                     if (kn > 0) then
                        h_dn = h_layer(k, kn)
                        tr_up = merge(hTr_old(k, kn)/h_dn, 0.0_wp, &
                                      h_dn > DRY_TOLERANCE)
                     else
                        g = mesh%bc_group(i, c)
                        if (g >= 1 .and. g <= n_bc) then
                           tr_up = tr_bc_ext(g)
                        else
                           tr_up = 0.0_wp
                        end if
                     end if
                  end if
               else if (flux_c < 0.0_wp) then
                  kn = mesh%neighbours(i, c)
                  if (kn > 0) then
                     ! Inflow: upwind the neighbour; empty neighbour => the
                     ! downwind (own) concentration.
                     h_up = h_layer(k, kn)
                     if (h_up > DRY_TOLERANCE) then
                        tr_up = hTr_old(k, kn)/h_up
                     else
                        h_dn = h_layer(k, c)
                        tr_up = merge(hTr_old(k, c)/h_dn, 0.0_wp, &
                                      h_dn > DRY_TOLERANCE)
                     end if
                  else
                     ! Dirichlet inflow: the per-group external tracer.
                     g = mesh%bc_group(i, c)
                     if (g >= 1 .and. g <= n_bc) then
                        tr_up = tr_bc_ext(g)
                     else
                        tr_up = 0.0_wp
                     end if
                  end if
               else
                  tr_up = 0.0_wp
               end if
               div_sum = div_sum + flux_c*tr_up
            end do
            hTr(k, c) = hTr_old(k, c) - dt*div_sum*inv_area
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_advect_tracer_dsig

   subroutine si_unstr_edge_layers_to_cells(mesh, nz, nct, nco, &
                                                 u_edge_layer, hu_layer, hv_layer)
      !! Per-layer least-squares A-grid cell-momentum reconstruction from
      !! the three per-layer edge-normal transports — the 2D
      !! `si_unstr_reconstruct_cell_momenta`, per layer (the documented
      !! reconstruction compromise; restart round-trips through the cell
      !! layer momenta).
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: u_edge_layer(nz, EDGES_PER_TRI, nct)
      real(wp), intent(inout) :: hu_layer(nz, nct), hv_layer(nz, nct)
         !! Cell layer momenta [m^2/s]; owned cells overwritten.

      integer :: c, i, k
      real(wp) :: m11, m12, m22, rx, ry, det, wlen, enx, eny, ue

      !$omp parallel do private(i, k, m11, m12, m22, rx, ry, det, wlen, enx, eny, ue)
      do c = 1, nco
         ! The 2x2 normal matrix is layer-independent (mesh only).
         m11 = 0.0_wp
         m12 = 0.0_wp
         m22 = 0.0_wp
         do i = 1, EDGES_PER_TRI
            wlen = mesh%edge_len(i, c)
            enx = mesh%edge_n1(i, c)
            eny = mesh%edge_n2(i, c)
            m11 = m11 + wlen*enx*enx
            m12 = m12 + wlen*enx*eny
            m22 = m22 + wlen*eny*eny
         end do
         det = m11*m22 - m12*m12
         do k = 1, nz
            rx = 0.0_wp
            ry = 0.0_wp
            do i = 1, EDGES_PER_TRI
               wlen = mesh%edge_len(i, c)
               enx = mesh%edge_n1(i, c)
               eny = mesh%edge_n2(i, c)
               ue = u_edge_layer(k, i, c)
               rx = rx + wlen*ue*enx
               ry = ry + wlen*ue*eny
            end do
            ! det ~ 0 only for a degenerate cell (collinear edge normals) —
            ! excluded by the fail-loud mesh quality gate at enable_si; the
            ! zero fallback is 1/0 armour for a diagnostic-only quantity
            ! (mass/eta close on the edges regardless).
            if (abs(det) > tiny(1.0_wp)) then
               hu_layer(k, c) = (m22*rx - m12*ry)/det
               hv_layer(k, c) = (m11*ry - m12*rx)/det
            else
               hu_layer(k, c) = 0.0_wp
               hv_layer(k, c) = 0.0_wp
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_edge_layers_to_cells

   subroutine si_unstr_update_thickness_ml(nz, nct, nco, eta, bed, dsig, &
                                                stage, height, h_layer)
      !! h^{n+1} = max(eta^{n+1} - b, 0) (exactly the edge-flux continuity
      !! update) + the sigma slaving h_layer_k = dsig_k*h, and the stage
      !! refresh stage = eta.
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta(nct)
      real(wp), intent(in) :: bed(nct)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(inout) :: stage(nct), height(nct)
      real(wp), intent(inout) :: h_layer(nz, nct)

      integer :: c, k

      !$omp parallel do private(k)
      do c = 1, nco
         stage(c) = eta(c)
         height(c) = max(eta(c) - bed(c), 0.0_wp)
         do k = 1, nz
            h_layer(k, c) = dsig(k)*height(c)
         end do
      end do
      !$omp end parallel do
   end subroutine si_unstr_update_thickness_ml

   subroutine si_unstr_apply_wind_top_layer(nz, nct, nco, h_layer, &
                                                 hu_layer, hv_layer, &
                                                 wind_tau_x, wind_tau_y, dt)
      !! Surface wind stress into the TOP layer (k = nz) transport:
      !! dt*tau/RHO_WATER — the explicit forcing that rides the pure-u^n
      !! predictor into G_k (by linearity R_k / Hhat / the Helmholtz
      !! operator are untouched; zero wind => bit-identical).  The
      !! unstructured wind is the uniform scalar pair on the barotropic
      !! state (`wind_tau_x/y`), so this is a trivially cheap kernel.
      !! Gated on a wet surface layer.
      integer, intent(in) :: nz
         !! Number of layers.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: h_layer(nz, nct)
      real(wp), intent(inout) :: hu_layer(nz, nct), hv_layer(nz, nct)
      real(wp), intent(in) :: wind_tau_x, wind_tau_y
         !! Uniform surface wind stress [Pa].
      real(wp), intent(in) :: dt

      integer :: c

      !$omp parallel do
      do c = 1, nco
         if (h_layer(nz, c) > H_VANISHED) then
            hu_layer(nz, c) = hu_layer(nz, c) + dt*wind_tau_x/RHO_WATER
            hv_layer(nz, c) = hv_layer(nz, c) + dt*wind_tau_y/RHO_WATER
         end if
      end do
      !$omp end parallel do
   end subroutine si_unstr_apply_wind_top_layer

end module rki_si_unstr_column
