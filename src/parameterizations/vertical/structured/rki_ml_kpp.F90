!! KPP (Large, McWilliams & Doney 1994) surface-boundary-layer mixing
!! for the structured coastal multilayer path (A-grid).
module rki_ml_kpp
   !! KPP surface-boundary-layer mixing for the coastal A-grid. Unlike the
   !! ocean path (implicit Thomas solve), the coastal path has no implicit
   !! vertical-diffusion solver, so the KPP overlay is fused with the explicit
   !! CFL-clamped diffusive-flux apply:
   !!   * `ml_kpp_boundary_layer` — per column u_*, surface buoyancy flux B_0,
   !!     and bulk-Ri boundary-layer depth h_b (LMD94 eq-23 V_t^2 term);
   !!     `bl_depth` persists across steps to seed V_t^2 (self-bootstrapping).
   !!   * `ml_kpp_viscosity` — momentum: nu_t = max(PP81 interior, KPP overlay
   !!     h_b·w_s·G(sigma)), explicit flux to hu/hv.
   !!   * `ml_kpp_diffusion_tracers` — S+T: same overlay + non-local
   !!     counter-gradient transport under destabilizing forcing.
   !!
   !! Overlay inactive (w_s=0) ⇒ falls through to bit-identical PP81. `vmix_use_kpp`
   !! off (default) leaves all results unchanged. Bottom-up layers (k=1 bed,
   !! k=nz surface); interface "k" sits between layer k and k+1.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH, VANISHING_LAYER_TOL
   implicit none
   private

   public :: ml_kpp_boundary_layer
   public :: ml_kpp_viscosity
   public :: ml_kpp_diffusion_tracers
   public :: kpp_shape

   real(wp), parameter :: EPS = 1.0e-20_wp
      !! Shear² floor in the interior gradient-Ri (matches PP81 kernels).
   real(wp), parameter :: SHEAR2_FLOOR = 1.0e-10_wp
      !! Velocity-difference floor in the bulk-Ri denominator (matches
      !! `rki_ocean_vmix`); stops Ri_bulk blowing up in quiescent water.

contains

   pure function kpp_shape(sigma) result(g)
      !! KPP cubic shape function G(sigma) = sigma·(1-sigma)² (LMD94), valid for
      !! 0 <= sigma <= 1. Exposed only for unit tests; hot kernels inline the
      !! one-line form to avoid a device cross-call.
      real(wp), intent(in) :: sigma
      real(wp) :: g
      g = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
   end function kpp_shape

   subroutine ml_kpp_boundary_layer(h_layer, hu_layer, hv_layer, rho_layer, &
                                         hS_layer, rho_0, alpha_T, beta_S, cp, &
                                         wind_x, wind_y, heat_flux, evap_minus_precip, &
                                         ri_crit, c_vt2, &
                                         u_star, buoy_flux, bl_depth, &
                                         grid, nz)
      !! Pass 1: per-column u_*, B_0, and bulk-Ri boundary-layer depth h_b.
      !!   u_* = sqrt(|tau|/rho_0), tau = surface wind stress (Pa).
      !!   B_0 = g·(alpha_T·F_T - beta_S·F_S), kinematic fluxes
      !!         F_T = Q_heat/(rho_0·cp), F_S = (E-P)·S_top; B_0 < 0 destabilizing.
      !!   h_b : top-down bulk-Ri sweep; crossing of `ri_crit` interpolated.
      !!         V_t²(d) = (c_vt2/ri_crit)·d·N(d)·w_s uses lagged `bl_depth`
      !!         (zero on first call ⇒ self-bootstrapping). No crossing ⇒ h_b = H.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0, alpha_T, beta_S, cp
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: heat_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: evap_minus_precip(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: ri_crit, c_vt2
      real(wp), intent(inout) :: u_star(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: buoy_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: bl_depth(grid%nx_total, grid%ny_total)

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: h_col, h_top, s_top, f_t_kin, f_s_kin, b_0, us
      real(wp) :: u_ref, v_ref, b_ref, d_centre_ref
      real(wp) :: h_b_lagged, wstar3_lagged, w_s_col
      real(wp) :: u_k, v_k, b_k, d_centre_k, delta_d, delta_b
      real(wp) :: n_brunt, vt2, shear2, ri_bulk, ri_prev, d_prev
      real(wp) :: denom, frac, h_b, h_k, h_k1
      logical :: crossing_found

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (nz < 2) return

      !$omp target teams distribute parallel do collapse(2) private(k, kk, h_col, h_top, s_top, f_t_kin, f_s_kin, b_0, us, u_ref, v_ref, b_ref, d_centre_ref, h_b_lagged, wstar3_lagged, w_s_col, u_k, v_k, b_k, d_centre_k, delta_d, delta_b, n_brunt, vt2, shear2, ri_bulk, ri_prev, d_prev, denom, frac, h_b, h_k, h_k1, crossing_found)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col <= ML_MIN_DEPTH) then
            u_star(i, j) = 0.0_wp
            buoy_flux(i, j) = 0.0_wp
            bl_depth(i, j) = 0.0_wp
            cycle
         end if

         ! Surface friction velocity (A-grid wind stress is at the cell
         ! centre, so no face averaging).
         us = sqrt(sqrt(wind_x(i, j)*wind_x(i, j) + wind_y(i, j)*wind_y(i, j))/rho_0)

         ! Surface buoyancy flux.
         h_top = h_layer(i, j, nz)
         if (h_top > DRY_TOLERANCE) then
            s_top = hS_layer(i, j, nz)/h_top
            u_ref = hu_layer(i, j, nz)/h_top
            v_ref = hv_layer(i, j, nz)/h_top
         else
            s_top = 0.0_wp
            u_ref = 0.0_wp
            v_ref = 0.0_wp
         end if
         f_t_kin = heat_flux(i, j)/(rho_0*cp)
         f_s_kin = evap_minus_precip(i, j)*s_top
         ! Buoyancy flux from the *dimensional* coastal EOS derivatives:
         ! b = -g*rho/rho_0, drho = -alpha_T*dT + beta_S*dS (kg/m^3 per
         ! degC / PSU), so B_0 = (g/rho_0)*(alpha_T*F_T - beta_S*F_S).
         ! (The ocean path folds the 1/rho_0 into its expansion
         ! coefficients; coastal keeps it explicit.)
         b_0 = (GRAVITY/rho_0)*(alpha_T*f_t_kin - beta_S*f_s_kin)

         u_star(i, j) = us
         buoy_flux(i, j) = b_0

         ! Reference buoyancy at the surface layer centre.
         b_ref = -GRAVITY*rho_layer(i, j, nz)/rho_0
         d_centre_ref = 0.5_wp*h_top

         ! Lagged convective scale for V_t² (self-bootstraps from 0).
         h_b_lagged = bl_depth(i, j)
         wstar3_lagged = max(0.0_wp, -b_0)*h_b_lagged
         w_s_col = sqrt(us*us + wstar3_lagged**(2.0_wp/3.0_wp))

         ri_prev = 0.0_wp
         d_prev = d_centre_ref
         h_b = 0.0_wp
         crossing_found = .false.
         d_centre_k = d_centre_ref
         do k = nz - 1, 1, -1
            if (.not. crossing_found) then
               h_k = h_layer(i, j, k)
               h_k1 = h_layer(i, j, k + 1)
               d_centre_k = d_centre_k + 0.5_wp*(h_k1 + h_k)
               if (h_k > DRY_TOLERANCE) then
                  u_k = hu_layer(i, j, k)/h_k
                  v_k = hv_layer(i, j, k)/h_k
               else
                  u_k = 0.0_wp
                  v_k = 0.0_wp
               end if
               b_k = -GRAVITY*rho_layer(i, j, k)/rho_0

               delta_d = d_centre_k - d_centre_ref
               delta_b = b_ref - b_k
               n_brunt = sqrt(max(0.0_wp, delta_b/delta_d))
               vt2 = (c_vt2/ri_crit)*d_centre_k*n_brunt*w_s_col
               shear2 = max((u_k - u_ref)*(u_k - u_ref) + &
                            (v_k - v_ref)*(v_k - v_ref) + vt2, SHEAR2_FLOOR)
               ri_bulk = delta_b*delta_d/shear2

               if (ri_bulk >= ri_crit) then
                  denom = ri_bulk - ri_prev
                  if (abs(denom) > 1.0e-12_wp) then
                     frac = (ri_crit - ri_prev)/denom
                  else
                     frac = 0.5_wp
                  end if
                  h_b = d_prev + frac*(d_centre_k - d_prev)
                  crossing_found = .true.
               end if

               ri_prev = ri_bulk
               d_prev = d_centre_k
            end if
         end do

         if (.not. crossing_found) h_b = h_col
         bl_depth(i, j) = h_b
      end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine ml_kpp_boundary_layer

   subroutine ml_kpp_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                    rho_0, nu_max, nu_bg, ri_alpha, &
                                    u_star, buoy_flux, bl_depth, &
                                    grid, nz, dt)
      !! Pass 2 (momentum): PP81 interior viscosity with the KPP
      !! shape-function overlay max'd in, then the explicit CFL-clamped
      !! diffusive flux applied to hu/hv.  Identical to
      !! `ml_vertical_viscosity` everywhere the overlay is inactive
      !! (w_s = 0 or interface below the BL).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0, nu_max, nu_bg, ri_alpha
      real(wp), intent(in) :: u_star(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: buoy_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: bl_depth(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, h_col
      real(wp) :: du, dv, drho, Ri, denom, nu_t, flux_u, flux_v
      real(wp) :: h_b, b_0, us, wstar3, w_star, w_s
      real(wp) :: d_below, d_iface, sigma, g_shape, nu_kpp

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (nz < 2) return

      !$omp target teams distribute parallel do collapse(2) private(k, kk, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, du, dv, drho, Ri, denom, nu_t, flux_u, flux_v, h_col, h_b, b_0, us, wstar3, w_star, w_s, d_below, d_iface, sigma, g_shape, nu_kpp)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then
            h_b = bl_depth(i, j)
            b_0 = buoy_flux(i, j)
            us = u_star(i, j)
            wstar3 = max(0.0_wp, -b_0)*h_b
            w_star = wstar3**(1.0_wp/3.0_wp)
            w_s = sqrt(us*us + w_star*w_star)

            d_below = 0.0_wp
            do k = 1, nz - 1
               h_k = h_layer(i, j, k)
               h_k1 = h_layer(i, j, k + 1)
               h_avg = 0.5_wp*(h_k + h_k1)
               d_below = d_below + h_k
               d_iface = h_col - d_below

               if (h_avg > VANISHING_LAYER_TOL .and. &
                   h_k > VANISHING_LAYER_TOL .and. &
                   h_k1 > VANISHING_LAYER_TOL) then
                  if (h_k > DRY_TOLERANCE) then
                     u_k = hu_layer(i, j, k)/h_k
                     v_k = hv_layer(i, j, k)/h_k
                  else
                     u_k = 0.0_wp
                     v_k = 0.0_wp
                  end if
                  if (h_k1 > DRY_TOLERANCE) then
                     u_k1 = hu_layer(i, j, k + 1)/h_k1
                     v_k1 = hv_layer(i, j, k + 1)/h_k1
                  else
                     u_k1 = 0.0_wp
                     v_k1 = 0.0_wp
                  end if

                  du = u_k1 - u_k
                  dv = v_k1 - v_k
                  drho = rho_layer(i, j, k) - rho_layer(i, j, k + 1)

                  Ri = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + EPS)
                  if (Ri < 0.0_wp) Ri = 0.0_wp

                  denom = 1.0_wp + ri_alpha*Ri
                  nu_t = nu_max/(denom*denom) + nu_bg

                  ! KPP overlay inside the boundary layer.
                  if (w_s > 0.0_wp .and. h_b > 0.0_wp .and. d_iface < h_b) then
                     sigma = d_iface/h_b
                     g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
                     nu_kpp = h_b*w_s*g_shape
                     if (nu_kpp > nu_t) nu_t = nu_kpp
                  end if

                  ! Explicit-flux stability: nu_t*dt/h_avg² < 0.5.
                  if (nu_t*dt > 0.5_wp*h_avg*h_avg) then
                     nu_t = 0.5_wp*h_avg*h_avg/dt
                  end if

                  flux_u = nu_t*du/h_avg*dt
                  flux_v = nu_t*dv/h_avg*dt

                  hu_layer(i, j, k) = hu_layer(i, j, k) + flux_u
                  hu_layer(i, j, k + 1) = hu_layer(i, j, k + 1) - flux_u
                  hv_layer(i, j, k) = hv_layer(i, j, k) + flux_v
                  hv_layer(i, j, k + 1) = hv_layer(i, j, k + 1) - flux_v
               end if
            end do
         end if
      end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine ml_kpp_viscosity

   subroutine ml_kpp_diffusion_tracers(h_layer, hu_layer, hv_layer, &
                                            hS_layer, hT_layer, rho_layer, rho_0, &
                                            nu_max, kappa_S_bg, kappa_T_bg, ri_alpha, &
                                            cs_nonlocal, heat_flux, evap_minus_precip, &
                                            cp, u_star, buoy_flux, bl_depth, &
                                            grid, nz, dt)
      !! Pass 2 (tracers): PP81 cubic interior diffusivity with the KPP
      !! overlay max'd in (same overlay as momentum), explicit flux on
      !! salinity + temperature, plus the non-local counter-gradient
      !! transport under destabilizing forcing.  Fuses the two
      !! `ml_vertical_diffusion_tracer` calls; bit-identical to them
      !! where the overlay + non-local terms vanish.
      !!
      !! Non-local: gamma(iface) = cs·G(sigma)·F_kin is a downward flux;
      !! its divergence redistributes the surface flux through the BL.
      !! Applied in a second column pass so the diffusive flux above is
      !! computed from un-redistributed tracer (matches the ocean order
      !! of vdiff-then-nonlocal).  gamma = 0 at the bed and free
      !! surface ⇒ the divergence conserves the column integral.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0, nu_max, kappa_S_bg, kappa_T_bg, ri_alpha
      real(wp), intent(in) :: cs_nonlocal, cp
      real(wp), intent(in) :: heat_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: evap_minus_precip(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: u_star(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: buoy_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: bl_depth(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, h_col
      real(wp) :: s_k, s_k1, t_k, t_k1, du, dv, drho, Ri, denom
      real(wp) :: kappa_S, kappa_T, flux_S, flux_T
      real(wp) :: h_b, b_0, us, wstar3, w_star, w_s
      real(wp) :: h_top, s_top, f_t_kin, f_s_kin
      real(wp) :: d_below, d_iface, sigma, g_shape, nu_kpp, gamma_T, gamma_S
      logical :: destabilizing

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (nz < 2) return

      !$omp target teams distribute parallel do collapse(2) private(k, kk, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, h_col, s_k, s_k1, t_k, t_k1, du, dv, drho, Ri, denom, kappa_S, kappa_T, flux_S, flux_T, h_b, b_0, us, wstar3, w_star, w_s, h_top, s_top, f_t_kin, f_s_kin, d_below, d_iface, sigma, g_shape, nu_kpp, gamma_T, gamma_S, destabilizing)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then
            h_b = bl_depth(i, j)
            b_0 = buoy_flux(i, j)
            us = u_star(i, j)
            wstar3 = max(0.0_wp, -b_0)*h_b
            w_star = wstar3**(1.0_wp/3.0_wp)
            w_s = sqrt(us*us + w_star*w_star)
            destabilizing = (b_0 < 0.0_wp)

            ! ---- diffusive flux (PP81 cubic + KPP overlay) ----
            d_below = 0.0_wp
            do k = 1, nz - 1
               h_k = h_layer(i, j, k)
               h_k1 = h_layer(i, j, k + 1)
               h_avg = 0.5_wp*(h_k + h_k1)
               d_below = d_below + h_k
               d_iface = h_col - d_below

               if (h_avg > VANISHING_LAYER_TOL .and. &
                   h_k > VANISHING_LAYER_TOL .and. &
                   h_k1 > VANISHING_LAYER_TOL) then
                  if (h_k > DRY_TOLERANCE) then
                     u_k = hu_layer(i, j, k)/h_k
                     v_k = hv_layer(i, j, k)/h_k
                     s_k = hS_layer(i, j, k)/h_k
                     t_k = hT_layer(i, j, k)/h_k
                  else
                     u_k = 0.0_wp
                     v_k = 0.0_wp
                     s_k = 0.0_wp
                     t_k = 0.0_wp
                  end if
                  if (h_k1 > DRY_TOLERANCE) then
                     u_k1 = hu_layer(i, j, k + 1)/h_k1
                     v_k1 = hv_layer(i, j, k + 1)/h_k1
                     s_k1 = hS_layer(i, j, k + 1)/h_k1
                     t_k1 = hT_layer(i, j, k + 1)/h_k1
                  else
                     u_k1 = 0.0_wp
                     v_k1 = 0.0_wp
                     s_k1 = 0.0_wp
                     t_k1 = 0.0_wp
                  end if

                  du = u_k1 - u_k
                  dv = v_k1 - v_k
                  drho = rho_layer(i, j, k) - rho_layer(i, j, k + 1)
                  Ri = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + EPS)
                  if (Ri < 0.0_wp) Ri = 0.0_wp

                  denom = 1.0_wp + ri_alpha*Ri
                  kappa_S = nu_max/(denom*denom*denom) + kappa_S_bg
                  kappa_T = nu_max/(denom*denom*denom) + kappa_T_bg

                  if (w_s > 0.0_wp .and. h_b > 0.0_wp .and. d_iface < h_b) then
                     sigma = d_iface/h_b
                     g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
                     nu_kpp = h_b*w_s*g_shape
                     if (nu_kpp > kappa_S) kappa_S = nu_kpp
                     if (nu_kpp > kappa_T) kappa_T = nu_kpp
                  end if

                  if (kappa_S*dt > 0.5_wp*h_avg*h_avg) kappa_S = 0.5_wp*h_avg*h_avg/dt
                  if (kappa_T*dt > 0.5_wp*h_avg*h_avg) kappa_T = 0.5_wp*h_avg*h_avg/dt

                  flux_S = kappa_S*(s_k1 - s_k)/h_avg*dt
                  flux_T = kappa_T*(t_k1 - t_k)/h_avg*dt

                  hS_layer(i, j, k) = hS_layer(i, j, k) + flux_S
                  hS_layer(i, j, k + 1) = hS_layer(i, j, k + 1) - flux_S
                  hT_layer(i, j, k) = hT_layer(i, j, k) + flux_T
                  hT_layer(i, j, k + 1) = hT_layer(i, j, k + 1) - flux_T
               end if
            end do

            ! ---- non-local counter-gradient transport ----
            if (destabilizing .and. w_s > 0.0_wp .and. h_b > 0.0_wp) then
               h_top = h_layer(i, j, nz)
               if (h_top > DRY_TOLERANCE) then
                  s_top = hS_layer(i, j, nz)/h_top
               else
                  s_top = 0.0_wp
               end if
               f_t_kin = heat_flux(i, j)/(rho_0*cp)
               f_s_kin = evap_minus_precip(i, j)*s_top

               d_below = 0.0_wp
               do k = 1, nz - 1
                  h_k = h_layer(i, j, k)
                  h_k1 = h_layer(i, j, k + 1)
                  h_avg = 0.5_wp*(h_k + h_k1)
                  d_below = d_below + h_k
                  d_iface = h_col - d_below

                  if (h_avg > VANISHING_LAYER_TOL .and. &
                      h_k > VANISHING_LAYER_TOL .and. &
                      h_k1 > VANISHING_LAYER_TOL .and. d_iface < h_b) then
                     sigma = d_iface/h_b
                     g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
                     gamma_T = cs_nonlocal*g_shape*f_t_kin
                     gamma_S = cs_nonlocal*g_shape*f_s_kin
                     hT_layer(i, j, k) = hT_layer(i, j, k) + dt*gamma_T
                     hT_layer(i, j, k + 1) = hT_layer(i, j, k + 1) - dt*gamma_T
                     hS_layer(i, j, k) = hS_layer(i, j, k) + dt*gamma_S
                     hS_layer(i, j, k + 1) = hS_layer(i, j, k + 1) - dt*gamma_S
                  end if
               end do
            end if
         end if
      end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine ml_kpp_diffusion_tracers

end module rki_ml_kpp
