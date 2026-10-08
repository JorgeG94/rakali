!! Two-equation k-epsilon vertical-mixing closure (structured multilayer).
module rki_ml_keps
   !! Standard two-equation k-epsilon model (Launder & Spalding 1974; constants
   !! per Rodi 1987 / Warner et al. 2005):
   !!   dk/dt   = d/dz(nu_t/sig_k dk/dz) + P + B - eps
   !!   deps/dt = d/dz(nu_t/sig_e deps/dz) + (eps/k)(ce1 P + ce3 B - ce2 eps)
   !!   nu_t = c_mu k^2/eps,  kappa_t = nu_t/Pr_t,  P = nu_t S^2,  B = -kappa_t N^2
   !!
   !! Staggering (bottom-up): velocity/density at layer centres (1..nz);
   !! k, eps, nu_t, kappa_t at interfaces (1..nz+1, 1=bed, nz+1=surface).
   !! Diffusion implicit (Thomas), dissipation sinks semi-implicit (Patankar)
   !! for positivity. One `do concurrent` over columns; sequential vertical
   !! sweep with NZ_STACK_MAX stack workspace.
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, &
                            KEPS_STAB_CONSTANT
#else
   use rki_constants, only: NZ_STACK_MAX, wp, GRAVITY, DRY_TOLERANCE, &
                            KEPS_STAB_CONSTANT
#endif
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH, VANISHING_LAYER_TOL
   use rki_keps_stability, only: keps_stability_funcs
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ml_keps_update
   public :: ml_keps_apply_momentum
   public :: ml_keps_apply_diffusion_tracer
   public :: ml_keps_friction_velocities

   ! Standard k-epsilon constants.
   real(wp), parameter :: CMU = 0.09_wp
   real(wp), parameter :: CE1 = 1.44_wp, CE2 = 1.92_wp, CE3 = -0.4_wp
   real(wp), parameter :: SIG_K = 1.0_wp, SIG_E = 1.3_wp
   real(wp), parameter :: PR_T = 1.0_wp
   real(wp), parameter :: K_MIN = 1.0e-8_wp, EPS_MIN = 1.0e-12_wp
   real(wp), parameter :: NUT_MIN = 1.0e-6_wp, NUT_MAX = 10.0_wp
   ! Wall-function constants.  KAPPA_EFF is the effective von Karman constant
   ! implied by the k-eps constants, sqrt(SIG_E*(CE2-CE1)*sqrt(CMU)) = 0.4327
   ! (hardcoded — sqrt() is not a standard constant-expression intrinsic).
   ! Z0_* are log-law roughness lengths (m) at the bed / surface.
   real(wp), parameter :: KAPPA_EFF = 0.4327_wp
   real(wp), parameter :: Z0_BED = 0.01_wp, Z0_SURF = 0.01_wp

contains

   pure subroutine ml_keps_update(h_layer, hu_layer, hv_layer, rho_layer, &
                                  tke, eps, nut, kt, rho_0, dt, grid, nz, &
                                  ustar_b, ustar_s, scheme)
      !! Advance k and eps one step per column and recover nu_t / kappa_t.
      !! `tke`, `eps`, `nut`, `kt` are interface-shaped (nz+1); `tke`/`eps`/
      !! `nut` carry state across steps (the diffusion reads the old `nut`).
      !!
      !! Boundary closure: if the optional friction velocities `ustar_b`
      !! (bed) and `ustar_s` (surface) are present, Dirichlet log-law wall
      !! functions are imposed at the first interior interfaces (interface 2
      !! above the bed, interface nz below the surface):
      !!   k = u_*^2 / sqrt(c_mu),   eps = u_*^3 / (kappa_eff*(z + z0))
      !! with z the interface distance from the wall (= the adjacent layer
      !! thickness).  Absent both, the boundary interfaces fall back to the
      !! zero-gradient (no-flux) closure (used by the standalone unit tests).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: tke(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(inout) :: eps(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(inout) :: nut(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(inout) :: kt(grid%nx_total, grid%ny_total, nz + 1)
         !! Scalar diffusivity, carried across steps: the buoyancy term reads
         !! the old value, so a stratification-aware scheme stays consistent.
      real(wp), intent(in) :: rho_0, dt
      real(wp), intent(in), optional :: ustar_b(grid%nx_total, grid%ny_total)
         !! Bed friction velocity sqrt(|tau_b|/rho_0) (m/s), per column.
      real(wp), intent(in), optional :: ustar_s(grid%nx_total, grid%ny_total)
         !! Surface friction velocity sqrt(|tau_wind|/rho_0) (m/s), per column.
      integer, intent(in), optional :: scheme
         !! Stability-function scheme (KEPS_STAB_*): absent / CONSTANT -> the
         !! standard constant c_mu = 0.09 (bit-identical to the original);
         !! GALPERIN -> stratification-aware quasi-equilibrium.

      integer :: i, j, m, nx, ny, nghost, ni, scheme_eff
      real(wp) :: h_col, u_dn, u_up, v_dn, v_up, dudz, dvdz, drho
      real(wp) :: prod, buoy, ktx, dd, du, denom, sink, src, eok
      real(wp) :: k_wb, k_ws, e_wb, e_ws, cmu_m, cmup_m
      logical :: walls
      real(wp) :: S2(NZ_STACK_MAX), N2(NZ_STACK_MAX), dz_if(NZ_STACK_MAX)
      real(wp) :: a(NZ_STACK_MAX), b(NZ_STACK_MAX), c(NZ_STACK_MAX)
      real(wp) :: rhs(NZ_STACK_MAX), cp(NZ_STACK_MAX)

      if (nz < 2) return
      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      ni = nz + 1
      walls = present(ustar_b) .and. present(ustar_s)
      scheme_eff = KEPS_STAB_CONSTANT
      if (present(scheme)) scheme_eff = scheme

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(m, h_col, u_dn, u_up, v_dn, v_up, dudz, dvdz, drho, &
               prod, buoy, ktx, dd, du, denom, sink, src, eok, &
               k_wb, k_ws, e_wb, e_ws, cmu_m, cmup_m, &
               S2, N2, dz_if, a, b, c, rhs, cp)

         h_col = 0.0_wp
         do m = 1, nz
            h_col = h_col + h_layer(i, j, m)
         end do
         if (h_col > ML_MIN_DEPTH) then

            ! ---- interface metrics, shear S^2 and buoyancy N^2 ----
            dz_if(1) = 0.5_wp*h_layer(i, j, 1)
            dz_if(ni) = 0.5_wp*h_layer(i, j, nz)
            S2(1) = 0.0_wp
            N2(1) = 0.0_wp
            S2(ni) = 0.0_wp
            N2(ni) = 0.0_wp
            do m = 2, nz
               dz_if(m) = 0.5_wp*(h_layer(i, j, m - 1) + h_layer(i, j, m))
               u_dn = hu_layer(i, j, m - 1)/max(h_layer(i, j, m - 1), DRY_TOLERANCE)
               u_up = hu_layer(i, j, m)/max(h_layer(i, j, m), DRY_TOLERANCE)
               v_dn = hv_layer(i, j, m - 1)/max(h_layer(i, j, m - 1), DRY_TOLERANCE)
               v_up = hv_layer(i, j, m)/max(h_layer(i, j, m), DRY_TOLERANCE)
               dudz = (u_up - u_dn)/dz_if(m)
               dvdz = (v_up - v_dn)/dz_if(m)
               S2(m) = dudz*dudz + dvdz*dvdz
               drho = rho_layer(i, j, m) - rho_layer(i, j, m - 1)   ! above - below
               N2(m) = -GRAVITY*drho/rho_0/dz_if(m)
            end do

            ! ---- log-law wall values (Dirichlet) when friction velocities
            !      are supplied; z = adjacent-layer thickness ----
            if (walls) then
               k_wb = max(ustar_b(i, j)*ustar_b(i, j)/sqrt(CMU), K_MIN)
               e_wb = max(ustar_b(i, j)**3/(KAPPA_EFF*(h_layer(i, j, 1) + Z0_BED)), EPS_MIN)
               k_ws = max(ustar_s(i, j)*ustar_s(i, j)/sqrt(CMU), K_MIN)
               e_ws = max(ustar_s(i, j)**3/(KAPPA_EFF*(h_layer(i, j, nz) + Z0_SURF)), EPS_MIN)
            end if

            ! ---- TKE: implicit diffusion + production - dissipation ----
            ! At the near-wall interfaces (2 = bed, nz = surface) the wall
            ! functions pin k by Dirichlet; interfaces 1 / nz+1 then copy
            ! the wall value via the zero-gradient rows in the tridiag solve.
            do m = 2, nz
               if (walls .and. (m == 2 .or. m == nz)) then
                  a(m) = 0.0_wp
                  c(m) = 0.0_wp
                  b(m) = 1.0_wp
                  if (m == nz) then
                     rhs(m) = k_ws
                  else
                     rhs(m) = k_wb
                  end if
               else
                  ktx = kt(i, j, m)
                  prod = nut(i, j, m)*S2(m)
                  buoy = -ktx*N2(m)
                  dd = 0.5_wp*(nut(i, j, m - 1) + nut(i, j, m))/SIG_K       ! layer m-1
                  du = 0.5_wp*(nut(i, j, m) + nut(i, j, m + 1))/SIG_K       ! layer m
                  a(m) = -dt/dz_if(m)*dd/h_layer(i, j, m - 1)
                  c(m) = -dt/dz_if(m)*du/h_layer(i, j, m)
                  sink = eps(i, j, m)/max(tke(i, j, m), K_MIN)
                  b(m) = 1.0_wp - a(m) - c(m) + dt*sink
                  src = max(prod + buoy, 0.0_wp)
                  rhs(m) = tke(i, j, m) + dt*src
               end if
            end do
            call keps_tridiag_solve(a, b, c, rhs, cp, ni)
            do m = 1, ni
               tke(i, j, m) = max(rhs(m), K_MIN)
            end do

            ! ---- eps: implicit diffusion + (eps/k)(ce1 P + ce3 B - ce2 eps) ----
            do m = 2, nz
               if (walls .and. (m == 2 .or. m == nz)) then
                  a(m) = 0.0_wp
                  c(m) = 0.0_wp
                  b(m) = 1.0_wp
                  if (m == nz) then
                     rhs(m) = e_ws
                  else
                     rhs(m) = e_wb
                  end if
               else
                  ktx = kt(i, j, m)
                  prod = nut(i, j, m)*S2(m)
                  buoy = -ktx*N2(m)
                  dd = 0.5_wp*(nut(i, j, m - 1) + nut(i, j, m))/SIG_E
                  du = 0.5_wp*(nut(i, j, m) + nut(i, j, m + 1))/SIG_E
                  a(m) = -dt/dz_if(m)*dd/h_layer(i, j, m - 1)
                  c(m) = -dt/dz_if(m)*du/h_layer(i, j, m)
                  eok = eps(i, j, m)/max(tke(i, j, m), K_MIN)              ! tke is NEW
                  sink = CE2*eok
                  b(m) = 1.0_wp - a(m) - c(m) + dt*sink
                  src = max(eok*(CE1*prod + CE3*buoy), 0.0_wp)
                  rhs(m) = eps(i, j, m) + dt*src
               end if
            end do
            call keps_tridiag_solve(a, b, c, rhs, cp, ni)
            do m = 1, ni
               eps(i, j, m) = max(rhs(m), EPS_MIN)
            end do

            ! ---- recover nu_t, kappa_t (constant or stratification-aware) ----
            do m = 1, ni
               call keps_stability_funcs(scheme_eff, S2(m), N2(m), tke(i, j, m), &
                                         eps(i, j, m), cmu_m, cmup_m)
               nut(i, j, m) = min(max(cmu_m*tke(i, j, m)*tke(i, j, m) &
                                      /max(eps(i, j, m), EPS_MIN), NUT_MIN), NUT_MAX)
               kt(i, j, m) = min(max(cmup_m*tke(i, j, m)*tke(i, j, m) &
                                     /max(eps(i, j, m), EPS_MIN), NUT_MIN), NUT_MAX)
            end do
         end if
      end do

   end subroutine ml_keps_update

   pure subroutine keps_tridiag_solve(a, b, c, rhs, cp, ni)
      !$acc routine seq
      !! In-place Thomas solve over interfaces 1..ni with zero-gradient
      !! (no-flux) boundary rows.  `rhs` returns the solution; `cp` scratch.
      integer, intent(in) :: ni
      real(wp), intent(inout) :: a(ni), b(ni), c(ni), rhs(ni), cp(ni)
      integer :: m
      real(wp) :: denom

      ! Zero-gradient closure: var(1)=var(2), var(ni)=var(ni-1).
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
   end subroutine keps_tridiag_solve

   pure subroutine ml_keps_apply_momentum(h_layer, hu_layer, hv_layer, nut, &
                                          grid, nz, dt)
      !! Apply k-epsilon eddy viscosity `nut` (interface-shaped nz+1, from
      !! `ml_keps_update`) to layer momentum via an implicit (Thomas) vertical-
      !! diffusion sweep. No-flux at bed/surface (wind stress and bottom drag
      !! are separate kernels); conserves column total sum_k h_k u_k exactly.
      !! Implicit (not explicit PP81) because nu_t can reach NUT_MAX. One
      !! `do concurrent` over columns; shared matrix, two RHS (u, v).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: nut(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_col, h_k, h_dn, h_up, dz_dn, dz_up, cd, cu
      real(wp) :: a(NZ_STACK_MAX), b(NZ_STACK_MAX), c(NZ_STACK_MAX)
      real(wp) :: ru(NZ_STACK_MAX), rv(NZ_STACK_MAX), cp(NZ_STACK_MAX)

      if (nz < 2) return
      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, h_col, h_k, h_dn, h_up, dz_dn, dz_up, cd, cu, &
               a, b, c, ru, rv, cp)

         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + h_layer(i, j, k)
         end do
         if (h_col > ML_MIN_DEPTH) then

            ! Build the tridiagonal.  Row k couples to k-1 across interface
            ! k and to k+1 across interface k+1; cd / cu are the implicit
            ! coupling coefficients dt*nu_t/(h_k * dz_interface).  A
            ! vanishing layer on either side drops that coupling (no-flux),
            ! so dividing by a near-zero thickness can't blow up (same guard
            ! philosophy as ml_vertical_viscosity).
            do k = 1, nz
               h_k = h_layer(i, j, k)
               if (h_k > VANISHING_LAYER_TOL) then
                  cd = 0.0_wp
                  if (k > 1) then
                     h_dn = h_layer(i, j, k - 1)
                     dz_dn = 0.5_wp*(h_dn + h_k)
                     if (h_dn > VANISHING_LAYER_TOL) then
                        cd = dt*nut(i, j, k)/(h_k*dz_dn)
                     end if
                  end if
                  cu = 0.0_wp
                  if (k < nz) then
                     h_up = h_layer(i, j, k + 1)
                     dz_up = 0.5_wp*(h_k + h_up)
                     if (h_up > VANISHING_LAYER_TOL) then
                        cu = dt*nut(i, j, k + 1)/(h_k*dz_up)
                     end if
                  end if
                  a(k) = -cd
                  c(k) = -cu
                  b(k) = 1.0_wp + cd + cu
                  ru(k) = hu_layer(i, j, k)/h_k
                  rv(k) = hv_layer(i, j, k)/h_k
               else
                  ! Vanishing / dry layer: decouple, leave momentum unchanged.
                  a(k) = 0.0_wp
                  c(k) = 0.0_wp
                  b(k) = 1.0_wp
                  if (h_k > 0.0_wp) then
                     ru(k) = hu_layer(i, j, k)/h_k
                     rv(k) = hv_layer(i, j, k)/h_k
                  else
                     ru(k) = 0.0_wp
                     rv(k) = 0.0_wp
                  end if
               end if
            end do

            call keps_tridiag_solve_noflux(a, b, c, ru, cp, nz)
            call keps_tridiag_solve_noflux(a, b, c, rv, cp, nz)

            do k = 1, nz
               hu_layer(i, j, k) = ru(k)*h_layer(i, j, k)
               hv_layer(i, j, k) = rv(k)*h_layer(i, j, k)
            end do
         end if
      end do

   end subroutine ml_keps_apply_momentum

   pure subroutine ml_keps_apply_diffusion_tracer(h_layer, hTr, kt, kappa_bg, &
                                                  grid, nz, dt)
      !! Implicit (Thomas) vertical diffusion of one tracer using k-eps
      !! diffusivity `kt` (interface-shaped) plus background floor `kappa_bg`.
      !! No-flux at bed/surface (surface tracer fluxes are a separate kernel);
      !! conserves column tracer content sum_k hTr_k exactly. Implicit (not
      !! explicit PP81) because kt can reach NUT_MAX. Pass `tracers(it)%hTr`.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hTr(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: kt(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: kappa_bg, dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_col, h_k, h_dn, h_up, dz_dn, dz_up, cd, cu
      real(wp) :: a(NZ_STACK_MAX), b(NZ_STACK_MAX), c(NZ_STACK_MAX)
      real(wp) :: r(NZ_STACK_MAX), cp(NZ_STACK_MAX)

      if (nz < 2) return
      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, h_col, h_k, h_dn, h_up, dz_dn, dz_up, cd, cu, a, b, c, r, cp)

         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + h_layer(i, j, k)
         end do
         if (h_col > ML_MIN_DEPTH) then
            ! Same tridiagonal as ml_keps_apply_momentum, but with the
            ! scalar diffusivity kt + kappa_bg and a single RHS (the tracer
            ! concentration hTr/h).  Vanishing layers drop the coupling.
            do k = 1, nz
               h_k = h_layer(i, j, k)
               if (h_k > VANISHING_LAYER_TOL) then
                  cd = 0.0_wp
                  if (k > 1) then
                     h_dn = h_layer(i, j, k - 1)
                     dz_dn = 0.5_wp*(h_dn + h_k)
                     if (h_dn > VANISHING_LAYER_TOL) then
                        cd = dt*(kt(i, j, k) + kappa_bg)/(h_k*dz_dn)
                     end if
                  end if
                  cu = 0.0_wp
                  if (k < nz) then
                     h_up = h_layer(i, j, k + 1)
                     dz_up = 0.5_wp*(h_k + h_up)
                     if (h_up > VANISHING_LAYER_TOL) then
                        cu = dt*(kt(i, j, k + 1) + kappa_bg)/(h_k*dz_up)
                     end if
                  end if
                  a(k) = -cd
                  c(k) = -cu
                  b(k) = 1.0_wp + cd + cu
                  r(k) = hTr(i, j, k)/h_k
               else
                  a(k) = 0.0_wp
                  c(k) = 0.0_wp
                  b(k) = 1.0_wp
                  if (h_k > 0.0_wp) then
                     r(k) = hTr(i, j, k)/h_k
                  else
                     r(k) = 0.0_wp
                  end if
               end if
            end do

            call keps_tridiag_solve_noflux(a, b, c, r, cp, nz)

            do k = 1, nz
               hTr(i, j, k) = r(k)*h_layer(i, j, k)
            end do
         end if
      end do

   end subroutine ml_keps_apply_diffusion_tracer

   pure subroutine ml_keps_friction_velocities(wind_x, wind_y, h_layer, &
                                               hu_layer, hv_layer, manning_n, &
                                               use_linear_drag, drag_r, rho_0, &
                                               ustar_s, ustar_b, grid, nz)
      !! Per-column surface and bed friction velocities feeding the k-eps
      !! wall functions (`ml_keps_update`).
      !!   surface: u_*s = sqrt(|tau_wind| / rho_0), tau_wind = (wind_x,
      !!            wind_y) the wind-stress components (Pa).
      !!   bed:     from the active bottom-drag law — quadratic Manning
      !!            u_*b = sqrt(C_d)*|u_b|, C_d = g n^2 / H^(1/3); or linear
      !!            (Rayleigh) u_*b = sqrt(r * h_b * |u_b|).
      !! Bottom velocity is the k=1 (bed) layer; H is the column depth.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: manning_n, drag_r, rho_0
      logical, intent(in) :: use_linear_drag
      real(wp), intent(out) :: ustar_s(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: ustar_b(grid%nx_total, grid%ny_total)

      integer :: i, j, m, nx, ny, nghost
      real(wp) :: h_col, h_b, u_b, v_b, speed_b, cf

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(m, h_col, h_b, u_b, v_b, speed_b, cf)

         ! surface: wind-stress magnitude -> u_*s = sqrt(|tau|/rho_0)
         ustar_s(i, j) = sqrt(sqrt(wind_x(i, j)*wind_x(i, j) &
                                   + wind_y(i, j)*wind_y(i, j))/rho_0)

         ! bed: bottom-layer (k=1) speed + column depth
         h_col = 0.0_wp
         do m = 1, nz
            h_col = h_col + h_layer(i, j, m)
         end do
         h_b = h_layer(i, j, 1)
         if (h_b > DRY_TOLERANCE .and. h_col > ML_MIN_DEPTH) then
            u_b = hu_layer(i, j, 1)/h_b
            v_b = hv_layer(i, j, 1)/h_b
            speed_b = sqrt(u_b*u_b + v_b*v_b)
            if (use_linear_drag) then
               ustar_b(i, j) = sqrt(drag_r*h_b*speed_b)
            else
               cf = GRAVITY*manning_n*manning_n/(h_col**(1.0_wp/3.0_wp))
               ustar_b(i, j) = sqrt(cf)*speed_b
            end if
         else
            ustar_b(i, j) = 0.0_wp
         end if
      end do

   end subroutine ml_keps_friction_velocities

   pure subroutine keps_tridiag_solve_noflux(a, b, c, rhs, cp, n)
      !$acc routine seq
      !! In-place Thomas solve over rows 1..n with the boundary closure
      !! already baked into a/b/c by the caller (no row override, unlike
      !! `keps_tridiag_solve`).  `rhs` returns the solution; `cp` scratch.
      !! a/b/c are read-only so one factorisation can serve multiple RHS.
      integer, intent(in) :: n
      real(wp), intent(in) :: a(n), b(n), c(n)
      real(wp), intent(inout) :: rhs(n), cp(n)
      integer :: m
      real(wp) :: denom

      cp(1) = c(1)/b(1)
      rhs(1) = rhs(1)/b(1)
      do m = 2, n
         denom = b(m) - a(m)*cp(m - 1)
         cp(m) = c(m)/denom
         rhs(m) = (rhs(m) - a(m)*rhs(m - 1))/denom
      end do
      do m = n - 1, 1, -1
         rhs(m) = rhs(m) - cp(m)*rhs(m + 1)
      end do
   end subroutine keps_tridiag_solve_noflux

end module rki_ml_keps
