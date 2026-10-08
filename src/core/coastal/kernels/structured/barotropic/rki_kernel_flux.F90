!! HLL / HLLC Riemann solver with hydrostatic reconstruction and minmod limiter
module rki_kernel_flux
   !! Net flux divergence for the shallow water equations via HLL or HLLC
   !! Riemann solver with hydrostatic reconstruction (Audusse et al. 2004)
   !! and minmod slope-limited reconstruction. HLL/HLLC differ only in the
   !! transverse momentum flux. Free surface eta = h + B is reconstructed
   !! (not h) for the well-balanced property.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, THIN_LAYER_THRESHOLD
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: compute_flux_hll, compute_flux_hllc

contains

   subroutine compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                                    mass_flux_x, mass_flux_y, grid)
      !! Net flux divergence at each interior cell (HLL Riemann solver).
      !! Also writes the per-face HLL h-flux at the right face
      !! `mass_flux_x(i, j)` and top face `mass_flux_y(i, j)`, reused by
      !! the multilayer tracer kernels for consistency-with-continuity (CWC).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total), &
                              hu(grid%nx_total, grid%ny_total), &
                              hv(grid%nx_total, grid%ny_total), &
                              b(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: flux_h(grid%nx_total, grid%ny_total), &
                               flux_hu(grid%nx_total, grid%ny_total), &
                               flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: mass_flux_x(grid%nx_total, grid%ny_total), &
                               mass_flux_y(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: fh, fhu, fhv, fmxL, fmxR, fmyB, fmyT
      logical, parameter :: use_hllc = .false.

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      !$omp target teams distribute parallel do collapse(2) private(fh, fhu, fhv, fmxL, fmxR, fmyB, fmyT)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         call flux_cell(h, hu, hv, b, nx, ny, i, j, dx, dy, use_hllc, &
                        fh, fhu, fhv, fmxL, fmxR, fmyB, fmyT)
         flux_h(i, j) = fh
         flux_hu(i, j) = fhu
         flux_hv(i, j) = fhv
         ! mass_flux_x(i, j) = right face of cell i (face i+1/2); the
         ! leftmost interior cell also writes its left face into (i-1, j).
         mass_flux_x(i, j) = fmxR
         mass_flux_y(i, j) = fmyT
         if (i == nghost + 1) mass_flux_x(i - 1, j) = fmxL
         if (j == nghost + 1) mass_flux_y(i, j - 1) = fmyB
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine compute_flux_hll

   subroutine compute_flux_hllc(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                                     mass_flux_x, mass_flux_y, grid)
      !! Net flux divergence at each interior cell (HLLC Riemann solver).
      !! h/hu flux and `mass_flux_x/y` match the HLL version; only the
      !! transverse momentum flux differs.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total), &
                              hu(grid%nx_total, grid%ny_total), &
                              hv(grid%nx_total, grid%ny_total), &
                              b(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: flux_h(grid%nx_total, grid%ny_total), &
                               flux_hu(grid%nx_total, grid%ny_total), &
                               flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: mass_flux_x(grid%nx_total, grid%ny_total), &
                               mass_flux_y(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: fh, fhu, fhv, fmxL, fmxR, fmyB, fmyT
      logical, parameter :: use_hllc = .true.

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      !$omp target teams distribute parallel do collapse(2) private(fh, fhu, fhv, fmxL, fmxR, fmyB, fmyT)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         call flux_cell(h, hu, hv, b, nx, ny, i, j, dx, dy, use_hllc, &
                        fh, fhu, fhv, fmxL, fmxR, fmyB, fmyT)
         flux_h(i, j) = fh
         flux_hu(i, j) = fhu
         flux_hv(i, j) = fhv
         mass_flux_x(i, j) = fmxR
         mass_flux_y(i, j) = fmyT
         if (i == nghost + 1) mass_flux_x(i - 1, j) = fmxL
         if (j == nghost + 1) mass_flux_y(i, j - 1) = fmyB
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine compute_flux_hllc

   pure subroutine flux_cell(h, hu, hv, b, nx, ny, i, j, dx, dy, use_hllc, &
                             out_fh, out_fhu, out_fhv, &
                             out_fmass_xL, out_fmass_xR, &
                             out_fmass_yB, out_fmass_yT)
      !! Net flux divergence for a single cell (i, j). Also returns the
      !! four per-face HLL h-fluxes for the multilayer tracer kernels' CWC.
      integer, intent(in) :: nx, ny, i, j
      real(wp), intent(in) :: h(nx, ny), hu(nx, ny), hv(nx, ny), b(nx, ny)
      real(wp), intent(in) :: dx, dy
      logical, intent(in) :: use_hllc
         !! Selects the Riemann solver. Passed as a literal constant from
         !! the two entry points, so the in-kernel branch constant-folds.
      real(wp), intent(out) :: out_fh, out_fhu, out_fhv
      real(wp), intent(out) :: out_fmass_xL
         !! HLL h-flux at left face (i-1/2).
      real(wp), intent(out) :: out_fmass_xR
         !! HLL h-flux at right face (i+1/2).
      real(wp), intent(out) :: out_fmass_yB
         !! HLL h-flux at bottom face (j-1/2).
      real(wp), intent(out) :: out_fmass_yT
         !! HLL h-flux at top face (j+1/2).

      real(wp) :: eta_c, eta_im1, eta_ip1, eta_im2, eta_ip2
      real(wp) :: eta_jm1, eta_jp1, eta_jm2, eta_jp2
      real(wp) :: seta_x, shu_x, shv_x, seta_y, shu_y, shv_y
      real(wp) :: seta_x_nb, shu_x_nb, shv_x_nb
      real(wp) :: seta_y_nb, shu_y_nb, shv_y_nb
      real(wp) :: eta_L, eta_R, hu_L, hv_L, hu_R, hv_R
      real(wp) :: h_L, h_R
      real(wp) :: b_max, h_hyd_l, h_hyd_r
      real(wp) :: u_l, v_l, u_r, v_r
      real(wp) :: fx_h_l, fx_hu_l, fx_hv_l
      real(wp) :: fx_h_r, fx_hu_r, fx_hv_r
      real(wp) :: fy_h_b, fy_hu_b, fy_hv_b
      real(wp) :: fy_h_t, fy_hu_t, fy_hv_t
      real(wp) :: h_hyd_left, h_pre_left
      real(wp) :: h_hyd_right, h_pre_right
      real(wp) :: h_hyd_bottom, h_pre_bottom
      real(wp) :: h_hyd_top, h_pre_top
      real(wp) :: half_g

      half_g = 0.5_wp*GRAVITY

      ! Precompute free surface eta = h + B
      eta_c = h(i, j) + b(i, j)
      eta_im1 = h(i - 1, j) + b(i - 1, j)
      eta_ip1 = h(i + 1, j) + b(i + 1, j)
      eta_im2 = h(i - 2, j) + b(i - 2, j)
      eta_ip2 = h(i + 2, j) + b(i + 2, j)
      eta_jm1 = h(i, j - 1) + b(i, j - 1)
      eta_jp1 = h(i, j + 1) + b(i, j + 1)
      eta_jm2 = h(i, j - 2) + b(i, j - 2)
      eta_jp2 = h(i, j + 2) + b(i, j + 2)

      ! Fall back to first-order (zero slopes) near wet/dry fronts.
      if (h(i, j) < THIN_LAYER_THRESHOLD &
          .or. h(i - 1, j) < DRY_TOLERANCE &
          .or. h(i + 1, j) < DRY_TOLERANCE) then
         seta_x = 0.0_wp
         shu_x = 0.0_wp
         shv_x = 0.0_wp
      else
         seta_x = minmod(eta_c - eta_im1, eta_ip1 - eta_c)
         shu_x = minmod(hu(i, j) - hu(i - 1, j), hu(i + 1, j) - hu(i, j))
         shv_x = minmod(hv(i, j) - hv(i - 1, j), hv(i + 1, j) - hv(i, j))
      end if

      if (h(i, j) < THIN_LAYER_THRESHOLD &
          .or. h(i, j - 1) < DRY_TOLERANCE &
          .or. h(i, j + 1) < DRY_TOLERANCE) then
         seta_y = 0.0_wp
         shu_y = 0.0_wp
         shv_y = 0.0_wp
      else
         seta_y = minmod(eta_c - eta_jm1, eta_jp1 - eta_c)
         shu_y = minmod(hu(i, j) - hu(i, j - 1), hu(i, j + 1) - hu(i, j))
         shv_y = minmod(hv(i, j) - hv(i, j - 1), hv(i, j + 1) - hv(i, j))
      end if

      ! ============================================================
      ! x-direction: LEFT interface (i-1/2)
      ! ============================================================
      eta_R = eta_c - 0.5_wp*seta_x
      hu_R = hu(i, j) - 0.5_wp*shu_x
      hv_R = hv(i, j) - 0.5_wp*shv_x
      h_R = max(0.0_wp, eta_R - b(i, j))
      h_pre_left = h_R

      if (h(i - 1, j) < THIN_LAYER_THRESHOLD &
          .or. h(i - 2, j) < DRY_TOLERANCE &
          .or. h(i, j) < DRY_TOLERANCE) then
         seta_x_nb = 0.0_wp
         shu_x_nb = 0.0_wp
         shv_x_nb = 0.0_wp
      else
         seta_x_nb = minmod(eta_im1 - eta_im2, eta_c - eta_im1)
         shu_x_nb = minmod(hu(i - 1, j) - hu(i - 2, j), hu(i, j) - hu(i - 1, j))
         shv_x_nb = minmod(hv(i - 1, j) - hv(i - 2, j), hv(i, j) - hv(i - 1, j))
      end if
      eta_L = eta_im1 + 0.5_wp*seta_x_nb
      hu_L = hu(i - 1, j) + 0.5_wp*shu_x_nb
      hv_L = hv(i - 1, j) + 0.5_wp*shv_x_nb
      h_L = max(0.0_wp, eta_L - b(i - 1, j))

      b_max = max(b(i - 1, j), b(i, j))
      h_hyd_l = max(0.0_wp, h_L + b(i - 1, j) - b_max)
      h_hyd_r = max(0.0_wp, h_R + b(i, j) - b_max)
      h_hyd_left = h_hyd_r

      call extract_velocity(h_L, hu_L, hv_L, u_l, v_l)
      call extract_velocity(h_R, hu_R, hv_R, u_r, v_r)

      if (use_hllc) then
         call hllc_flux_x(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                          h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                          fx_h_l, fx_hu_l, fx_hv_l)
      else
         call hll_flux_x(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                         h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                         fx_h_l, fx_hu_l, fx_hv_l)
      end if

      ! ============================================================
      ! x-direction: RIGHT interface (i+1/2)
      ! ============================================================
      eta_L = eta_c + 0.5_wp*seta_x
      hu_L = hu(i, j) + 0.5_wp*shu_x
      hv_L = hv(i, j) + 0.5_wp*shv_x
      h_L = max(0.0_wp, eta_L - b(i, j))
      h_pre_right = h_L

      if (h(i + 1, j) < THIN_LAYER_THRESHOLD &
          .or. h(i, j) < DRY_TOLERANCE &
          .or. h(i + 2, j) < DRY_TOLERANCE) then
         seta_x_nb = 0.0_wp
         shu_x_nb = 0.0_wp
         shv_x_nb = 0.0_wp
      else
         seta_x_nb = minmod(eta_ip1 - eta_c, eta_ip2 - eta_ip1)
         shu_x_nb = minmod(hu(i + 1, j) - hu(i, j), hu(i + 2, j) - hu(i + 1, j))
         shv_x_nb = minmod(hv(i + 1, j) - hv(i, j), hv(i + 2, j) - hv(i + 1, j))
      end if
      eta_R = eta_ip1 - 0.5_wp*seta_x_nb
      hu_R = hu(i + 1, j) - 0.5_wp*shu_x_nb
      hv_R = hv(i + 1, j) - 0.5_wp*shv_x_nb
      h_R = max(0.0_wp, eta_R - b(i + 1, j))

      b_max = max(b(i, j), b(i + 1, j))
      h_hyd_l = max(0.0_wp, h_L + b(i, j) - b_max)
      h_hyd_r = max(0.0_wp, h_R + b(i + 1, j) - b_max)
      h_hyd_right = h_hyd_l

      call extract_velocity(h_L, hu_L, hv_L, u_l, v_l)
      call extract_velocity(h_R, hu_R, hv_R, u_r, v_r)

      if (use_hllc) then
         call hllc_flux_x(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                          h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                          fx_h_r, fx_hu_r, fx_hv_r)
      else
         call hll_flux_x(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                         h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                         fx_h_r, fx_hu_r, fx_hv_r)
      end if

      ! ============================================================
      ! y-direction: BOTTOM interface (j-1/2)
      ! ============================================================
      eta_R = eta_c - 0.5_wp*seta_y
      hu_R = hu(i, j) - 0.5_wp*shu_y
      hv_R = hv(i, j) - 0.5_wp*shv_y
      h_R = max(0.0_wp, eta_R - b(i, j))
      h_pre_bottom = h_R

      if (h(i, j - 1) < THIN_LAYER_THRESHOLD &
          .or. h(i, j - 2) < DRY_TOLERANCE &
          .or. h(i, j) < DRY_TOLERANCE) then
         seta_y_nb = 0.0_wp
         shu_y_nb = 0.0_wp
         shv_y_nb = 0.0_wp
      else
         seta_y_nb = minmod(eta_jm1 - eta_jm2, eta_c - eta_jm1)
         shu_y_nb = minmod(hu(i, j - 1) - hu(i, j - 2), hu(i, j) - hu(i, j - 1))
         shv_y_nb = minmod(hv(i, j - 1) - hv(i, j - 2), hv(i, j) - hv(i, j - 1))
      end if
      eta_L = eta_jm1 + 0.5_wp*seta_y_nb
      hu_L = hu(i, j - 1) + 0.5_wp*shu_y_nb
      hv_L = hv(i, j - 1) + 0.5_wp*shv_y_nb
      h_L = max(0.0_wp, eta_L - b(i, j - 1))

      b_max = max(b(i, j - 1), b(i, j))
      h_hyd_l = max(0.0_wp, h_L + b(i, j - 1) - b_max)
      h_hyd_r = max(0.0_wp, h_R + b(i, j) - b_max)
      h_hyd_bottom = h_hyd_r

      call extract_velocity(h_L, hu_L, hv_L, u_l, v_l)
      call extract_velocity(h_R, hu_R, hv_R, u_r, v_r)

      if (use_hllc) then
         call hllc_flux_y(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                          h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                          fy_h_b, fy_hu_b, fy_hv_b)
      else
         call hll_flux_y(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                         h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                         fy_h_b, fy_hu_b, fy_hv_b)
      end if

      ! ============================================================
      ! y-direction: TOP interface (j+1/2)
      ! ============================================================
      eta_L = eta_c + 0.5_wp*seta_y
      hu_L = hu(i, j) + 0.5_wp*shu_y
      hv_L = hv(i, j) + 0.5_wp*shv_y
      h_L = max(0.0_wp, eta_L - b(i, j))
      h_pre_top = h_L

      if (h(i, j + 1) < THIN_LAYER_THRESHOLD &
          .or. h(i, j) < DRY_TOLERANCE &
          .or. h(i, j + 2) < DRY_TOLERANCE) then
         seta_y_nb = 0.0_wp
         shu_y_nb = 0.0_wp
         shv_y_nb = 0.0_wp
      else
         seta_y_nb = minmod(eta_jp1 - eta_c, eta_jp2 - eta_jp1)
         shu_y_nb = minmod(hu(i, j + 1) - hu(i, j), hu(i, j + 2) - hu(i, j + 1))
         shv_y_nb = minmod(hv(i, j + 1) - hv(i, j), hv(i, j + 2) - hv(i, j + 1))
      end if
      eta_R = eta_jp1 - 0.5_wp*seta_y_nb
      hu_R = hu(i, j + 1) - 0.5_wp*shu_y_nb
      hv_R = hv(i, j + 1) - 0.5_wp*shv_y_nb
      h_R = max(0.0_wp, eta_R - b(i, j + 1))

      b_max = max(b(i, j), b(i, j + 1))
      h_hyd_l = max(0.0_wp, h_L + b(i, j) - b_max)
      h_hyd_r = max(0.0_wp, h_R + b(i, j + 1) - b_max)
      h_hyd_top = h_hyd_l

      call extract_velocity(h_L, hu_L, hv_L, u_l, v_l)
      call extract_velocity(h_R, hu_R, hv_R, u_r, v_r)

      if (use_hllc) then
         call hllc_flux_y(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                          h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                          fy_h_t, fy_hu_t, fy_hv_t)
      else
         call hll_flux_y(h_hyd_l, h_hyd_l*u_l, h_hyd_l*v_l, &
                         h_hyd_r, h_hyd_r*u_r, h_hyd_r*v_r, &
                         fy_h_t, fy_hu_t, fy_hv_t)
      end if

      ! ============================================================
      ! Net flux divergence with well-balanced correction
      ! ============================================================
      out_fh = -(fx_h_r - fx_h_l)/dx &
               - (fy_h_t - fy_h_b)/dy

      out_fhu = -(fx_hu_r - fx_hu_l)/dx &
                - (fy_hu_t - fy_hu_b)/dy &
                - half_g*(h_hyd_left**2 - h_pre_left**2 &
                          - h_hyd_right**2 + h_pre_right**2)/dx

      out_fhv = -(fx_hv_r - fx_hv_l)/dx &
                - (fy_hv_t - fy_hv_b)/dy &
                - half_g*(h_hyd_bottom**2 - h_pre_bottom**2 &
                          - h_hyd_top**2 + h_pre_top**2)/dy

      ! Per-face HLL mass fluxes (each interior face computed twice but
      ! symmetric, so values match exactly).
      out_fmass_xL = fx_h_l
      out_fmass_xR = fx_h_r
      out_fmass_yB = fy_h_b
      out_fmass_yT = fy_h_t

   end subroutine flux_cell

   pure function minmod(a, b) result(m)
      !! Minmod slope limiter
      real(wp), intent(in) :: a, b
      real(wp) :: m

      if (a*b <= 0.0_wp) then
         m = 0.0_wp
      else if (abs(a) < abs(b)) then
         m = a
      else
         m = b
      end if

   end function minmod

   pure subroutine extract_velocity(h_val, hu_val, hv_val, u_val, v_val)
      !! Extract velocity from conserved variables, handling dry cells
      real(wp), intent(in) :: h_val, hu_val, hv_val
      real(wp), intent(out) :: u_val, v_val

      if (h_val > DRY_TOLERANCE) then
         u_val = hu_val/h_val
         v_val = hv_val/h_val
      else
         u_val = 0.0_wp
         v_val = 0.0_wp
      end if

   end subroutine extract_velocity

   pure subroutine hll_flux_x(h_l, hu_l, hv_l, h_r, hu_r, hv_r, &
                              f_h, f_hu, f_hv)
      !! HLL numerical flux in x-direction between left and right states
      real(wp), intent(in) :: h_l, hu_l, hv_l
      real(wp), intent(in) :: h_r, hu_r, hv_r
      real(wp), intent(out) :: f_h, f_hu, f_hv

      real(wp) :: u_l, u_r, v_l, v_r, a_l, a_r
      real(wp) :: s_l, s_r, denom
      real(wp) :: fl_h, fl_hu, fl_hv
      real(wp) :: fr_h, fr_hu, fr_hv

      if (h_l > DRY_TOLERANCE) then
         u_l = hu_l/h_l
         v_l = hv_l/h_l
      else
         u_l = 0.0_wp
         v_l = 0.0_wp
      end if

      if (h_r > DRY_TOLERANCE) then
         u_r = hu_r/h_r
         v_r = hv_r/h_r
      else
         u_r = 0.0_wp
         v_r = 0.0_wp
      end if

      a_l = sqrt(GRAVITY*max(h_l, 0.0_wp))
      a_r = sqrt(GRAVITY*max(h_r, 0.0_wp))
      ! Wave-speed estimate.  Davis (min/max of the one-sided speeds) is right
      ! for a genuine two-sided WET Riemann problem, but wrong at a wet/dry
      ! front and must be replaced there.
      !
      ! With the right state dry, `a_R = 0` and `u_R = 0`, so Davis returns
      ! `s_r = u_L + a_L`.  The exact dry-front speed is `u_L + 2 a_L` — the
      ! Ritter (1892) result — so Davis underestimates the right-going wave by
      ! up to a factor of two and the front cannot advance at the right rate.
      ! The dry-state speeds below are the standard fix (Toro, "Riemann
      ! Solvers and Numerical Methods for Fluid Dynamics", 3rd ed. S10.4;
      ! Einfeldt et al. 1991).
      !
      ! Measured consequence before this branch existed: the Ritter dry-bed
      ! tip advanced at 1.49 c0 rather than 2 c0, the leading tongue below
      ! ~1 % of the upstream depth was absent, and on a sloping beach a thin
      ! film would not recede at all (validation phases 4 and 5).
      if (h_r <= DRY_TOLERANCE) then
         ! Dry on the RIGHT: the right-going wave IS the dry front.
         s_l = u_l - a_l
         s_r = u_l + 2.0_wp*a_l
      else if (h_l <= DRY_TOLERANCE) then
         ! Dry on the LEFT: mirror image.
         s_l = u_r - 2.0_wp*a_r
         s_r = u_r + a_r
      else
         s_l = min(u_l - a_l, u_r - a_r)
         s_r = max(u_l + a_l, u_r + a_r)
      end if

      fl_h = hu_l
      fl_hu = hu_l*u_l + 0.5_wp*GRAVITY*h_l*h_l
      fl_hv = hu_l*v_l

      fr_h = hu_r
      fr_hu = hu_r*u_r + 0.5_wp*GRAVITY*h_r*h_r
      fr_hv = hu_r*v_r

      if (s_l >= 0.0_wp) then
         f_h = fl_h
         f_hu = fl_hu
         f_hv = fl_hv
      else if (s_r <= 0.0_wp) then
         f_h = fr_h
         f_hu = fr_hu
         f_hv = fr_hv
      else
         denom = 1.0_wp/(s_r - s_l)
         f_h = (s_r*fl_h - s_l*fr_h &
                + s_l*s_r*(h_r - h_l))*denom
         f_hu = (s_r*fl_hu - s_l*fr_hu &
                 + s_l*s_r*(hu_r - hu_l))*denom
         f_hv = (s_r*fl_hv - s_l*fr_hv &
                 + s_l*s_r*(hv_r - hv_l))*denom
      end if

   end subroutine hll_flux_x

   pure subroutine hll_flux_y(h_b, hu_b, hv_b, h_t, hu_t, hv_t, &
                              g_h, g_hu, g_hv)
      !! HLL numerical flux in y-direction between bottom and top states
      real(wp), intent(in) :: h_b, hu_b, hv_b
      real(wp), intent(in) :: h_t, hu_t, hv_t
      real(wp), intent(out) :: g_h, g_hu, g_hv

      real(wp) :: u_b, u_t, v_b, v_t, a_b, a_t
      real(wp) :: s_l, s_r, denom
      real(wp) :: gl_h, gl_hu, gl_hv
      real(wp) :: gr_h, gr_hu, gr_hv

      if (h_b > DRY_TOLERANCE) then
         u_b = hu_b/h_b
         v_b = hv_b/h_b
      else
         u_b = 0.0_wp
         v_b = 0.0_wp
      end if

      if (h_t > DRY_TOLERANCE) then
         u_t = hu_t/h_t
         v_t = hv_t/h_t
      else
         u_t = 0.0_wp
         v_t = 0.0_wp
      end if

      a_b = sqrt(GRAVITY*max(h_b, 0.0_wp))
      a_t = sqrt(GRAVITY*max(h_t, 0.0_wp))
      ! Wave-speed estimate.  Davis (min/max of the one-sided speeds) is right
      ! for a genuine two-sided WET Riemann problem, but wrong at a wet/dry
      ! front and must be replaced there.
      !
      ! With the right state dry, `a_R = 0` and `u_R = 0`, so Davis returns
      ! `s_r = u_L + a_L`.  The exact dry-front speed is `u_L + 2 a_L` — the
      ! Ritter (1892) result — so Davis underestimates the right-going wave by
      ! up to a factor of two and the front cannot advance at the right rate.
      ! The dry-state speeds below are the standard fix (Toro, "Riemann
      ! Solvers and Numerical Methods for Fluid Dynamics", 3rd ed. S10.4;
      ! Einfeldt et al. 1991).
      !
      ! Measured consequence before this branch existed: the Ritter dry-bed
      ! tip advanced at 1.49 c0 rather than 2 c0, the leading tongue below
      ! ~1 % of the upstream depth was absent, and on a sloping beach a thin
      ! film would not recede at all (validation phases 4 and 5).
      if (h_t <= DRY_TOLERANCE) then
         ! Dry on the RIGHT: the right-going wave IS the dry front.
         s_l = v_b - a_b
         s_r = v_b + 2.0_wp*a_b
      else if (h_b <= DRY_TOLERANCE) then
         ! Dry on the LEFT: mirror image.
         s_l = v_t - 2.0_wp*a_t
         s_r = v_t + a_t
      else
         s_l = min(v_b - a_b, v_t - a_t)
         s_r = max(v_b + a_b, v_t + a_t)
      end if

      gl_h = hv_b
      gl_hu = hu_b*v_b
      gl_hv = hv_b*v_b + 0.5_wp*GRAVITY*h_b*h_b

      gr_h = hv_t
      gr_hu = hu_t*v_t
      gr_hv = hv_t*v_t + 0.5_wp*GRAVITY*h_t*h_t

      if (s_l >= 0.0_wp) then
         g_h = gl_h
         g_hu = gl_hu
         g_hv = gl_hv
      else if (s_r <= 0.0_wp) then
         g_h = gr_h
         g_hu = gr_hu
         g_hv = gr_hv
      else
         denom = 1.0_wp/(s_r - s_l)
         g_h = (s_r*gl_h - s_l*gr_h &
                + s_l*s_r*(h_t - h_b))*denom
         g_hu = (s_r*gl_hu - s_l*gr_hu &
                 + s_l*s_r*(hu_t - hu_b))*denom
         g_hv = (s_r*gl_hv - s_l*gr_hv &
                 + s_l*s_r*(hv_t - hv_b))*denom
      end if

   end subroutine hll_flux_y

   pure subroutine hllc_flux_x(h_l, hu_l, hv_l, h_r, hu_r, hv_r, &
                               f_h, f_hu, f_hv)
      !! HLLC numerical flux in x-direction. Transverse momentum hv is
      !! upwinded across the contact wave: f_hv = f_h * v_upwind
      !! (v_l if s_star >= 0, else v_r); h and hu match HLL.
      real(wp), intent(in) :: h_l, hu_l, hv_l
      real(wp), intent(in) :: h_r, hu_r, hv_r
      real(wp), intent(out) :: f_h, f_hu, f_hv

      real(wp) :: u_l, u_r, v_l, v_r, a_l, a_r
      real(wp) :: s_l, s_r, denom, denom_star, s_star
      real(wp) :: fl_h, fl_hu, fl_hv
      real(wp) :: fr_h, fr_hu, fr_hv

      if (h_l > DRY_TOLERANCE) then
         u_l = hu_l/h_l
         v_l = hv_l/h_l
      else
         u_l = 0.0_wp
         v_l = 0.0_wp
      end if

      if (h_r > DRY_TOLERANCE) then
         u_r = hu_r/h_r
         v_r = hv_r/h_r
      else
         u_r = 0.0_wp
         v_r = 0.0_wp
      end if

      a_l = sqrt(GRAVITY*max(h_l, 0.0_wp))
      a_r = sqrt(GRAVITY*max(h_r, 0.0_wp))
      ! Wave-speed estimate.  Davis (min/max of the one-sided speeds) is right
      ! for a genuine two-sided WET Riemann problem, but wrong at a wet/dry
      ! front and must be replaced there.
      !
      ! With the right state dry, `a_R = 0` and `u_R = 0`, so Davis returns
      ! `s_r = u_L + a_L`.  The exact dry-front speed is `u_L + 2 a_L` — the
      ! Ritter (1892) result — so Davis underestimates the right-going wave by
      ! up to a factor of two and the front cannot advance at the right rate.
      ! The dry-state speeds below are the standard fix (Toro, "Riemann
      ! Solvers and Numerical Methods for Fluid Dynamics", 3rd ed. S10.4;
      ! Einfeldt et al. 1991).
      !
      ! Measured consequence before this branch existed: the Ritter dry-bed
      ! tip advanced at 1.49 c0 rather than 2 c0, the leading tongue below
      ! ~1 % of the upstream depth was absent, and on a sloping beach a thin
      ! film would not recede at all (validation phases 4 and 5).
      if (h_r <= DRY_TOLERANCE) then
         ! Dry on the RIGHT: the right-going wave IS the dry front.
         s_l = u_l - a_l
         s_r = u_l + 2.0_wp*a_l
      else if (h_l <= DRY_TOLERANCE) then
         ! Dry on the LEFT: mirror image.
         s_l = u_r - 2.0_wp*a_r
         s_r = u_r + a_r
      else
         s_l = min(u_l - a_l, u_r - a_r)
         s_r = max(u_l + a_l, u_r + a_r)
      end if

      fl_h = hu_l
      fl_hu = hu_l*u_l + 0.5_wp*GRAVITY*h_l*h_l
      fl_hv = hu_l*v_l

      fr_h = hu_r
      fr_hu = hu_r*u_r + 0.5_wp*GRAVITY*h_r*h_r
      fr_hv = hu_r*v_r

      if (s_l >= 0.0_wp) then
         f_h = fl_h
         f_hu = fl_hu
         f_hv = fl_hv
      else if (s_r <= 0.0_wp) then
         f_h = fr_h
         f_hu = fr_hu
         f_hv = fr_hv
      else
         denom = 1.0_wp/(s_r - s_l)
         f_h = (s_r*fl_h - s_l*fr_h &
                + s_l*s_r*(h_r - h_l))*denom
         f_hu = (s_r*fl_hu - s_l*fr_hu &
                 + s_l*s_r*(hu_r - hu_l))*denom

         ! Contact-wave speed (Toro). Fall back to the arithmetic mean
         ! when the denominator vanishes.
         denom_star = h_r*(u_r - s_r) - h_l*(u_l - s_l)
         if (abs(denom_star) > DRY_TOLERANCE) then
            s_star = (s_l*h_r*(u_r - s_r) - s_r*h_l*(u_l - s_l))/denom_star
         else
            s_star = 0.5_wp*(u_l + u_r)
         end if

         if (s_star >= 0.0_wp) then
            f_hv = f_h*v_l
         else
            f_hv = f_h*v_r
         end if
      end if

   end subroutine hllc_flux_x

   pure subroutine hllc_flux_y(h_b, hu_b, hv_b, h_t, hu_t, hv_t, &
                               g_h, g_hu, g_hv)
      !! HLLC numerical flux in y-direction. Mirrors hllc_flux_x with
      !! (u, hu) and (v, hv) swapped; transverse hu is upwinded.
      real(wp), intent(in) :: h_b, hu_b, hv_b
      real(wp), intent(in) :: h_t, hu_t, hv_t
      real(wp), intent(out) :: g_h, g_hu, g_hv

      real(wp) :: u_b, u_t, v_b, v_t, a_b, a_t
      real(wp) :: s_l, s_r, denom, denom_star, s_star
      real(wp) :: gl_h, gl_hu, gl_hv
      real(wp) :: gr_h, gr_hu, gr_hv

      if (h_b > DRY_TOLERANCE) then
         u_b = hu_b/h_b
         v_b = hv_b/h_b
      else
         u_b = 0.0_wp
         v_b = 0.0_wp
      end if

      if (h_t > DRY_TOLERANCE) then
         u_t = hu_t/h_t
         v_t = hv_t/h_t
      else
         u_t = 0.0_wp
         v_t = 0.0_wp
      end if

      a_b = sqrt(GRAVITY*max(h_b, 0.0_wp))
      a_t = sqrt(GRAVITY*max(h_t, 0.0_wp))
      ! Wave-speed estimate.  Davis (min/max of the one-sided speeds) is right
      ! for a genuine two-sided WET Riemann problem, but wrong at a wet/dry
      ! front and must be replaced there.
      !
      ! With the right state dry, `a_R = 0` and `u_R = 0`, so Davis returns
      ! `s_r = u_L + a_L`.  The exact dry-front speed is `u_L + 2 a_L` — the
      ! Ritter (1892) result — so Davis underestimates the right-going wave by
      ! up to a factor of two and the front cannot advance at the right rate.
      ! The dry-state speeds below are the standard fix (Toro, "Riemann
      ! Solvers and Numerical Methods for Fluid Dynamics", 3rd ed. S10.4;
      ! Einfeldt et al. 1991).
      !
      ! Measured consequence before this branch existed: the Ritter dry-bed
      ! tip advanced at 1.49 c0 rather than 2 c0, the leading tongue below
      ! ~1 % of the upstream depth was absent, and on a sloping beach a thin
      ! film would not recede at all (validation phases 4 and 5).
      if (h_t <= DRY_TOLERANCE) then
         ! Dry on the RIGHT: the right-going wave IS the dry front.
         s_l = v_b - a_b
         s_r = v_b + 2.0_wp*a_b
      else if (h_b <= DRY_TOLERANCE) then
         ! Dry on the LEFT: mirror image.
         s_l = v_t - 2.0_wp*a_t
         s_r = v_t + a_t
      else
         s_l = min(v_b - a_b, v_t - a_t)
         s_r = max(v_b + a_b, v_t + a_t)
      end if

      gl_h = hv_b
      gl_hu = hu_b*v_b
      gl_hv = hv_b*v_b + 0.5_wp*GRAVITY*h_b*h_b

      gr_h = hv_t
      gr_hu = hu_t*v_t
      gr_hv = hv_t*v_t + 0.5_wp*GRAVITY*h_t*h_t

      if (s_l >= 0.0_wp) then
         g_h = gl_h
         g_hu = gl_hu
         g_hv = gl_hv
      else if (s_r <= 0.0_wp) then
         g_h = gr_h
         g_hu = gr_hu
         g_hv = gr_hv
      else
         denom = 1.0_wp/(s_r - s_l)
         g_h = (s_r*gl_h - s_l*gr_h &
                + s_l*s_r*(h_t - h_b))*denom
         g_hv = (s_r*gl_hv - s_l*gr_hv &
                 + s_l*s_r*(hv_t - hv_b))*denom

         denom_star = h_t*(v_t - s_r) - h_b*(v_b - s_l)
         if (abs(denom_star) > DRY_TOLERANCE) then
            s_star = (s_l*h_t*(v_t - s_r) - s_r*h_b*(v_b - s_l))/denom_star
         else
            s_star = 0.5_wp*(v_b + v_t)
         end if

         if (s_star >= 0.0_wp) then
            g_hu = g_h*u_b
         else
            g_hu = g_h*u_t
         end if
      end if

   end subroutine hllc_flux_y

end module rki_kernel_flux
