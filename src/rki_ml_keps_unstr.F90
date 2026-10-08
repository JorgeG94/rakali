!! Two-equation k-epsilon vertical-mixing closure (unstructured multilayer).
module rki_ml_keps_unstr
   !! Unstructured (triangular) port of the structured `rki_ml_keps`: the same
   !! clean-room k-epsilon closure (Launder-Spalding; Rodi 1987; Warner 2005)
   !! per water column, cell-centric on the unstructured C-grid.
   !!
   !!   dk/dt   = d/dz(nu_t/sig_k dk/dz) + P + B - eps
   !!   deps/dt = d/dz(nu_t/sig_e deps/dz) + (eps/k)(ce1 P + ce3 B - ce2 eps)
   !!   nu_t = c_mu k^2/eps,  kappa_t = c_mu' k^2/eps
   !!
   !! Staggering (bottom-up, ROMS): velocity/density at layer centres
   !! (1..nz); k/eps/nu_t/kappa_t at interfaces (1..nz+1, 1=bed, nz+1=surface).
   !! Per cell: implicit (Thomas) vertical diffusion + semi-implicit (Patankar)
   !! dissipation sinks for positivity.  Stability functions come from the
   !! shared `rki_keps_stability` (constant / Galperin / Canuto) as a pure
   !! `!$acc routine seq` (cross-module — a cheap device call, not inlined).
   !!
   !! The kernel is a flat explicit-shape `_impl` (arrays as `(nz, ncell)`)
   !! invoked through a thin state-passed shim — keeps the device loop free of
   !! derived-type dereferences and is directly unit-testable like the
   !! structured `ml_keps_update`.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, KEPS_STAB_CONSTANT
   use rki_state_unstr, only: unstr_state_t
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR, BPG_NZ_MAX, &
                                  VANISHING_LAYER_TOL_UNSTR
   use rki_keps_stability, only: keps_stability_funcs
   implicit none
   private

   public :: ml_keps_update_unstr
   public :: ml_keps_update_unstr_impl   !! public for unit testing
   public :: ml_keps_apply_momentum_unstr
   public :: ml_keps_apply_momentum_unstr_impl   !! public for unit testing
   public :: ml_keps_apply_diffusion_tracer_unstr
   public :: ml_keps_apply_diffusion_tracer_unstr_impl   !! public for unit testing
   public :: ml_keps_friction_velocities_unstr

   ! Standard k-epsilon constants (mirror rki_ml_keps; physical, never change).
   real(wp), parameter :: CMU = 0.09_wp
   real(wp), parameter :: CE1 = 1.44_wp, CE2 = 1.92_wp, CE3 = -0.4_wp
   real(wp), parameter :: SIG_K = 1.0_wp, SIG_E = 1.3_wp
   real(wp), parameter :: K_MIN = 1.0e-8_wp, EPS_MIN = 1.0e-12_wp
   real(wp), parameter :: NUT_MIN = 1.0e-6_wp, NUT_MAX = 10.0_wp
   real(wp), parameter :: KAPPA_EFF = 0.4327_wp
   real(wp), parameter :: Z0_BED = 0.01_wp, Z0_SURF = 0.01_wp

contains

   subroutine ml_keps_update_unstr(state, dt, use_walls)
      !! State-passed shim: advance k/eps one step per cell and recover
      !! nu_t/kappa_t.  `use_walls` imposes the Dirichlet log-law wall functions
      !! (from keps_ustar_s/b); .false. uses the zero-gradient closure (the
      !! standalone unit tests).  Extracts the flat arrays and calls the impl.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt
      logical, intent(in) :: use_walls

      if (.not. state%use_multilayer) return
      if (state%multilayer%nz_ml < 2) return
      if (.not. allocated(state%multilayer%keps_tke)) return

      call ml_keps_update_unstr_impl(state%multilayer%nz_ml, state%ncells_owned, &
                                     state%multilayer%h_layer, state%multilayer%hu_layer, &
                                     state%multilayer%hv_layer, state%multilayer%rho_layer, &
                                     state%multilayer%keps_tke, state%multilayer%keps_eps, &
                                     state%multilayer%keps_nut, state%multilayer%keps_kt, &
                                     state%multilayer%keps_ustar_s, state%multilayer%keps_ustar_b, &
                                     state%barotropic%height_cv, state%multilayer%rho_0, &
                                     dt, use_walls, state%multilayer%keps_stability)
   end subroutine ml_keps_update_unstr

   pure subroutine ml_keps_update_unstr_impl(nz, n, h_layer, hu_layer, hv_layer, rho_layer, &
                                             tke, eps, nut, kt, ustar_s, ustar_b, height_cv, &
                                             rho_0, dt, use_walls, scheme)
      !! Per-cell k/eps advance on flat explicit-shape arrays.  One
      !! `do concurrent` over cells (thread-per-cell — full occupancy for the
      !! short vertical dimension); sequential vertical sweep per cell with
      !! BPG_NZ_MAX stack workspace, same pattern as the unstructured PP81/KPP.
      integer, intent(in) :: nz, n, scheme
      real(wp), intent(in) :: h_layer(nz, n), hu_layer(nz, n)
      real(wp), intent(in) :: hv_layer(nz, n), rho_layer(nz, n)
      real(wp), intent(inout) :: tke(nz + 1, n), eps(nz + 1, n)
      real(wp), intent(inout) :: nut(nz + 1, n), kt(nz + 1, n)
      real(wp), intent(in) :: ustar_s(n), ustar_b(n), height_cv(n)
      real(wp), intent(in) :: rho_0, dt
      logical, intent(in) :: use_walls

      integer :: cell, m, ni
      real(wp) :: u_dn, u_up, v_dn, v_up, dudz, dvdz, drho
      real(wp) :: prod, buoy, ktx, dd, du, sink, src, eok
      real(wp) :: k_wb, k_ws, e_wb, e_ws, cmu_m, cmup_m
      real(wp) :: S2(BPG_NZ_MAX), N2(BPG_NZ_MAX), dz_if(BPG_NZ_MAX)
      real(wp) :: a(BPG_NZ_MAX), b(BPG_NZ_MAX), c(BPG_NZ_MAX)
      real(wp) :: rhs(BPG_NZ_MAX), cp(BPG_NZ_MAX)

      ni = nz + 1

      do concurrent(cell=1:n) &
         local(m, u_dn, u_up, v_dn, v_up, dudz, dvdz, drho, &
               prod, buoy, ktx, dd, du, sink, src, eok, &
               k_wb, k_ws, e_wb, e_ws, cmu_m, cmup_m, &
               S2, N2, dz_if, a, b, c, rhs, cp)
         if (height_cv(cell) >= ML_MIN_DEPTH_UNSTR) then

            ! ---- interface metrics, shear S^2 and buoyancy N^2 ----
            dz_if(1) = 0.5_wp*h_layer(1, cell)
            dz_if(ni) = 0.5_wp*h_layer(nz, cell)
            S2(1) = 0.0_wp
            N2(1) = 0.0_wp
            S2(ni) = 0.0_wp
            N2(ni) = 0.0_wp
            do m = 2, nz
               dz_if(m) = 0.5_wp*(h_layer(m - 1, cell) + h_layer(m, cell))
               u_dn = hu_layer(m - 1, cell)/max(h_layer(m - 1, cell), DRY_TOLERANCE)
               u_up = hu_layer(m, cell)/max(h_layer(m, cell), DRY_TOLERANCE)
               v_dn = hv_layer(m - 1, cell)/max(h_layer(m - 1, cell), DRY_TOLERANCE)
               v_up = hv_layer(m, cell)/max(h_layer(m, cell), DRY_TOLERANCE)
               dudz = (u_up - u_dn)/dz_if(m)
               dvdz = (v_up - v_dn)/dz_if(m)
               S2(m) = dudz*dudz + dvdz*dvdz
               drho = rho_layer(m, cell) - rho_layer(m - 1, cell)   ! above - below
               N2(m) = -GRAVITY*drho/rho_0/dz_if(m)
            end do

            ! ---- log-law wall values (Dirichlet) when use_walls ----
            k_wb = K_MIN
            e_wb = EPS_MIN
            k_ws = K_MIN
            e_ws = EPS_MIN
            if (use_walls) then
               k_wb = max(ustar_b(cell)*ustar_b(cell)/sqrt(CMU), K_MIN)
               e_wb = max(ustar_b(cell)**3/(KAPPA_EFF*(h_layer(1, cell) + Z0_BED)), EPS_MIN)
               k_ws = max(ustar_s(cell)*ustar_s(cell)/sqrt(CMU), K_MIN)
               e_ws = max(ustar_s(cell)**3/(KAPPA_EFF*(h_layer(nz, cell) + Z0_SURF)), EPS_MIN)
            end if

            ! ---- TKE: implicit diffusion + production - dissipation ----
            do m = 2, nz
               if (use_walls .and. (m == 2 .or. m == nz)) then
                  a(m) = 0.0_wp
                  c(m) = 0.0_wp
                  b(m) = 1.0_wp
                  if (m == nz) then
                     rhs(m) = k_ws
                  else
                     rhs(m) = k_wb
                  end if
               else
                  ktx = kt(m, cell)
                  prod = nut(m, cell)*S2(m)
                  buoy = -ktx*N2(m)
                  dd = 0.5_wp*(nut(m - 1, cell) + nut(m, cell))/SIG_K
                  du = 0.5_wp*(nut(m, cell) + nut(m + 1, cell))/SIG_K
                  a(m) = -dt/dz_if(m)*dd/h_layer(m - 1, cell)
                  c(m) = -dt/dz_if(m)*du/h_layer(m, cell)
                  sink = eps(m, cell)/max(tke(m, cell), K_MIN)
                  b(m) = 1.0_wp - a(m) - c(m) + dt*sink
                  src = max(prod + buoy, 0.0_wp)
                  rhs(m) = tke(m, cell) + dt*src
               end if
            end do
            call keps_tridiag_solve_u(a, b, c, rhs, cp, ni)
            do m = 1, ni
               tke(m, cell) = max(rhs(m), K_MIN)
            end do

            ! ---- eps: implicit diffusion + (eps/k)(ce1 P + ce3 B - ce2 eps) ----
            do m = 2, nz
               if (use_walls .and. (m == 2 .or. m == nz)) then
                  a(m) = 0.0_wp
                  c(m) = 0.0_wp
                  b(m) = 1.0_wp
                  if (m == nz) then
                     rhs(m) = e_ws
                  else
                     rhs(m) = e_wb
                  end if
               else
                  ktx = kt(m, cell)
                  prod = nut(m, cell)*S2(m)
                  buoy = -ktx*N2(m)
                  dd = 0.5_wp*(nut(m - 1, cell) + nut(m, cell))/SIG_E
                  du = 0.5_wp*(nut(m, cell) + nut(m + 1, cell))/SIG_E
                  a(m) = -dt/dz_if(m)*dd/h_layer(m - 1, cell)
                  c(m) = -dt/dz_if(m)*du/h_layer(m, cell)
                  eok = eps(m, cell)/max(tke(m, cell), K_MIN)        ! tke is NEW
                  sink = CE2*eok
                  b(m) = 1.0_wp - a(m) - c(m) + dt*sink
                  src = max(eok*(CE1*prod + CE3*buoy), 0.0_wp)
                  rhs(m) = eps(m, cell) + dt*src
               end if
            end do
            call keps_tridiag_solve_u(a, b, c, rhs, cp, ni)
            do m = 1, ni
               eps(m, cell) = max(rhs(m), EPS_MIN)
            end do

            ! ---- recover nu_t, kappa_t ----
            do m = 1, ni
               call keps_stability_funcs(scheme, S2(m), N2(m), tke(m, cell), &
                                         eps(m, cell), cmu_m, cmup_m)
               nut(m, cell) = min(max(cmu_m*tke(m, cell)*tke(m, cell) &
                                      /max(eps(m, cell), EPS_MIN), NUT_MIN), NUT_MAX)
               kt(m, cell) = min(max(cmup_m*tke(m, cell)*tke(m, cell) &
                                     /max(eps(m, cell), EPS_MIN), NUT_MIN), NUT_MAX)
            end do
         end if
      end do
   end subroutine ml_keps_update_unstr_impl

   subroutine ml_keps_apply_momentum_unstr(state, dt)
      !! State-passed shim: apply the keps eddy viscosity to layer momentum by
      !! an implicit (Thomas) vertical-diffusion sweep.  No-flux at bed/surface
      !! (wind + bottom drag are separate kernels), conserves the column total.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt

      if (.not. state%use_multilayer) return
      if (state%multilayer%nz_ml < 2) return
      if (.not. allocated(state%multilayer%keps_nut)) return

      call ml_keps_apply_momentum_unstr_impl(state%multilayer%nz_ml, state%ncells_owned, &
                                             state%multilayer%h_layer, state%multilayer%hu_layer, &
                                             state%multilayer%hv_layer, state%multilayer%keps_nut, &
                                             state%barotropic%height_cv, dt)
   end subroutine ml_keps_apply_momentum_unstr

   pure subroutine ml_keps_apply_momentum_unstr_impl(nz, n, h_layer, hu_layer, hv_layer, &
                                                     nut, height_cv, dt)
      !! Per-cell implicit vertical diffusion of momentum on flat explicit-shape
      !! arrays.  Implicit (not the explicit PP81 clamp) -- the keps nu_t can
      !! reach NUT_MAX.  One factorisation, two RHS (u, v).
      integer, intent(in) :: nz, n
      real(wp), intent(in) :: h_layer(nz, n)
      real(wp), intent(inout) :: hu_layer(nz, n), hv_layer(nz, n)
      real(wp), intent(in) :: nut(nz + 1, n), height_cv(n)
      real(wp), intent(in) :: dt

      integer :: cell, k
      real(wp) :: h_k, h_dn, h_up, dz_dn, dz_up, cd, cu
      real(wp) :: a(BPG_NZ_MAX), b(BPG_NZ_MAX), c(BPG_NZ_MAX)
      real(wp) :: ru(BPG_NZ_MAX), rv(BPG_NZ_MAX), cp(BPG_NZ_MAX)

      do concurrent(cell=1:n) &
         local(k, h_k, h_dn, h_up, dz_dn, dz_up, cd, cu, a, b, c, ru, rv, cp)
         if (height_cv(cell) >= ML_MIN_DEPTH_UNSTR) then
            do k = 1, nz
               h_k = h_layer(k, cell)
               if (h_k > VANISHING_LAYER_TOL_UNSTR) then
                  cd = 0.0_wp
                  if (k > 1) then
                     h_dn = h_layer(k - 1, cell)
                     dz_dn = 0.5_wp*(h_dn + h_k)
                     if (h_dn > VANISHING_LAYER_TOL_UNSTR) cd = dt*nut(k, cell)/(h_k*dz_dn)
                  end if
                  cu = 0.0_wp
                  if (k < nz) then
                     h_up = h_layer(k + 1, cell)
                     dz_up = 0.5_wp*(h_k + h_up)
                     if (h_up > VANISHING_LAYER_TOL_UNSTR) cu = dt*nut(k + 1, cell)/(h_k*dz_up)
                  end if
                  a(k) = -cd
                  c(k) = -cu
                  b(k) = 1.0_wp + cd + cu
                  ru(k) = hu_layer(k, cell)/h_k
                  rv(k) = hv_layer(k, cell)/h_k
               else
                  a(k) = 0.0_wp
                  c(k) = 0.0_wp
                  b(k) = 1.0_wp
                  if (h_k > 0.0_wp) then
                     ru(k) = hu_layer(k, cell)/h_k
                     rv(k) = hv_layer(k, cell)/h_k
                  else
                     ru(k) = 0.0_wp
                     rv(k) = 0.0_wp
                  end if
               end if
            end do
            call keps_tridiag_solve_noflux_u(a, b, c, ru, cp, nz)
            call keps_tridiag_solve_noflux_u(a, b, c, rv, cp, nz)
            do k = 1, nz
               hu_layer(k, cell) = ru(k)*h_layer(k, cell)
               hv_layer(k, cell) = rv(k)*h_layer(k, cell)
            end do
         end if
      end do
   end subroutine ml_keps_apply_momentum_unstr_impl

   subroutine ml_keps_apply_diffusion_tracer_unstr(state, hTr, kappa_bg, dt)
      !! State-passed shim: implicit vertical diffusion of one tracer with the
      !! keps diffusivity kt + a background floor.  No-flux top/bottom (surface
      !! tracer fluxes are separate), conserves the column tracer content.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout) :: hTr(:, :)
      real(wp), intent(in) :: kappa_bg, dt

      if (.not. state%use_multilayer) return
      if (state%multilayer%nz_ml < 2) return
      if (.not. allocated(state%multilayer%keps_kt)) return

      call ml_keps_apply_diffusion_tracer_unstr_impl(state%multilayer%nz_ml, state%ncells_owned, &
                                                     state%multilayer%h_layer, hTr, &
                                                     state%multilayer%keps_kt, kappa_bg, &
                                                     state%barotropic%height_cv, dt)
   end subroutine ml_keps_apply_diffusion_tracer_unstr

   pure subroutine ml_keps_apply_diffusion_tracer_unstr_impl(nz, n, h_layer, hTr, kt, kappa_bg, &
                                                             height_cv, dt)
      !! Per-cell implicit tracer diffusion on flat explicit-shape arrays.  Same
      !! tridiagonal as the momentum apply but with the scalar diffusivity
      !! kt + kappa_bg and a single RHS (the tracer concentration hTr/h).
      integer, intent(in) :: nz, n
      real(wp), intent(in) :: h_layer(nz, n)
      real(wp), intent(inout) :: hTr(nz, n)
      real(wp), intent(in) :: kt(nz + 1, n), height_cv(n)
      real(wp), intent(in) :: kappa_bg, dt

      integer :: cell, k
      real(wp) :: h_k, h_dn, h_up, dz_dn, dz_up, cd, cu
      real(wp) :: a(BPG_NZ_MAX), b(BPG_NZ_MAX), c(BPG_NZ_MAX)
      real(wp) :: r(BPG_NZ_MAX), cp(BPG_NZ_MAX)

      do concurrent(cell=1:n) &
         local(k, h_k, h_dn, h_up, dz_dn, dz_up, cd, cu, a, b, c, r, cp)
         if (height_cv(cell) >= ML_MIN_DEPTH_UNSTR) then
            do k = 1, nz
               h_k = h_layer(k, cell)
               if (h_k > VANISHING_LAYER_TOL_UNSTR) then
                  cd = 0.0_wp
                  if (k > 1) then
                     h_dn = h_layer(k - 1, cell)
                     dz_dn = 0.5_wp*(h_dn + h_k)
                     if (h_dn > VANISHING_LAYER_TOL_UNSTR) cd = dt*(kt(k, cell) + kappa_bg)/(h_k*dz_dn)
                  end if
                  cu = 0.0_wp
                  if (k < nz) then
                     h_up = h_layer(k + 1, cell)
                     dz_up = 0.5_wp*(h_k + h_up)
                     if (h_up > VANISHING_LAYER_TOL_UNSTR) cu = dt*(kt(k + 1, cell) + kappa_bg)/(h_k*dz_up)
                  end if
                  a(k) = -cd
                  c(k) = -cu
                  b(k) = 1.0_wp + cd + cu
                  r(k) = hTr(k, cell)/h_k
               else
                  a(k) = 0.0_wp
                  c(k) = 0.0_wp
                  b(k) = 1.0_wp
                  if (h_k > 0.0_wp) then
                     r(k) = hTr(k, cell)/h_k
                  else
                     r(k) = 0.0_wp
                  end if
               end if
            end do
            call keps_tridiag_solve_noflux_u(a, b, c, r, cp, nz)
            do k = 1, nz
               hTr(k, cell) = r(k)*h_layer(k, cell)
            end do
         end if
      end do
   end subroutine ml_keps_apply_diffusion_tracer_unstr_impl

   pure subroutine ml_keps_friction_velocities_unstr(state)
      !! Per-cell surface + bed friction velocities feeding the keps wall
      !! functions.  Surface: u_*s = sqrt(|tau_wind|/rho_0) from the (uniform)
      !! unstructured wind stress.  Bed: quadratic Manning u_*b = sqrt(C_d)*|u_b|,
      !! C_d = g n^2 / H^(1/3), from the bottom (k=1) layer speed.  Writes
      !! keps_ustar_s/b on the multilayer state.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, m, nz, n
      real(wp) :: tau_x, tau_y, rho_0, mann, us_s, h_col, h_b, u_b, v_b, speed_b, cf

      if (.not. state%use_multilayer) return
      if (.not. allocated(state%multilayer%keps_ustar_s)) return
      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      tau_x = state%barotropic%wind_tau_x
      tau_y = state%barotropic%wind_tau_y
      rho_0 = state%multilayer%rho_0
      mann = state%barotropic%manning_n
      us_s = sqrt(sqrt(tau_x*tau_x + tau_y*tau_y)/rho_0)   ! uniform wind -> per-cell constant

      do concurrent(cell=1:n) local(m, h_col, h_b, u_b, v_b, speed_b, cf)
         state%multilayer%keps_ustar_s(cell) = us_s
         h_col = 0.0_wp
         do m = 1, nz
            h_col = h_col + state%multilayer%h_layer(m, cell)
         end do
         h_b = state%multilayer%h_layer(1, cell)
         if (h_b > DRY_TOLERANCE .and. h_col > ML_MIN_DEPTH_UNSTR) then
            u_b = state%multilayer%hu_layer(1, cell)/h_b
            v_b = state%multilayer%hv_layer(1, cell)/h_b
            speed_b = sqrt(u_b*u_b + v_b*v_b)
            cf = GRAVITY*mann*mann/(h_col**(1.0_wp/3.0_wp))
            state%multilayer%keps_ustar_b(cell) = sqrt(cf)*speed_b
         else
            state%multilayer%keps_ustar_b(cell) = 0.0_wp
         end if
      end do
   end subroutine ml_keps_friction_velocities_unstr

   pure subroutine keps_tridiag_solve_u(a, b, c, rhs, cp, ni)
      !$acc routine seq
      !! In-place Thomas solve over interfaces 1..ni with zero-gradient
      !! (no-flux) boundary rows.  `rhs` returns the solution; `cp` scratch.
      !! (Duplicate of the structured keps_tridiag_solve; tiny standard
      !! algorithm — can be consolidated into a shared module later.)
      integer, intent(in) :: ni
      real(wp), intent(inout) :: a(ni), b(ni), c(ni), rhs(ni), cp(ni)
      integer :: m
      real(wp) :: denom

      a(1) = 0.0_wp
      b(1) = 1.0_wp
      c(1) = -1.0_wp
      rhs(1) = 0.0_wp
      a(ni) = -1.0_wp
      b(ni) = 1.0_wp
      c(ni) = 0.0_wp
      rhs(ni) = 0.0_wp

      cp(1) = c(1)/b(1)
      rhs(1) = rhs(1)/b(1)
      do m = 2, ni
         denom = b(m) - a(m)*cp(m - 1)
         cp(m) = c(m)/denom
         rhs(m) = (rhs(m) - a(m)*rhs(m - 1))/denom
      end do
      do m = ni - 1, 1, -1
         rhs(m) = rhs(m) - cp(m)*rhs(m + 1)
      end do
   end subroutine keps_tridiag_solve_u

   pure subroutine keps_tridiag_solve_noflux_u(a, b, c, rhs, cp, nn)
      !$acc routine seq
      !! In-place Thomas solve over rows 1..nn with the boundary closure already
      !! baked into a/b/c by the caller (no row override, unlike
      !! `keps_tridiag_solve_u`).  a/b/c read-only so one factorisation serves
      !! multiple RHS.  Duplicate of the structured keps_tridiag_solve_noflux,
      !! kept same-module so it inlines in the per-cell apply loop.
      integer, intent(in) :: nn
      real(wp), intent(in) :: a(nn), b(nn), c(nn)
      real(wp), intent(inout) :: rhs(nn), cp(nn)
      integer :: m
      real(wp) :: denom

      cp(1) = c(1)/b(1)
      rhs(1) = rhs(1)/b(1)
      do m = 2, nn
         denom = b(m) - a(m)*cp(m - 1)
         cp(m) = c(m)/denom
         rhs(m) = (rhs(m) - a(m)*rhs(m - 1))/denom
      end do
      do m = nn - 1, 1, -1
         rhs(m) = rhs(m) - cp(m)*rhs(m + 1)
      end do
   end subroutine keps_tridiag_solve_noflux_u

end module rki_ml_keps_unstr
