!! Conservative state update kernels
module rki_kernel_update
   !! Advances the conserved variables: forward Euler update (each RK
   !! stage) and the SSP-RK2 averaging step.
   !!
   !! Thin-cell treatment (all three kernels, in order):
   !!  1. negative depths are clamped to zero (momentum zeroed);
   !!  2. genuinely dry cells (`h < DRY_TOLERANCE`) have momentum zeroed;
   !!  3. thin-but-wet cells (`h < eps = DESING_EPS`) get the
   !!     Kurganov & Petrova (2007) desingularized velocity
   !!     `u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4))` with momentum
   !!     reconstructed as `hu := h*u`.  For `h >= eps` the formula is the
   !!     identity `u = hu/h`, so well-resolved cells are left bit-exact;
   !!  4. a last-resort Froude backstop `|u| <= FROUDE_BACKSTOP*sqrt(g h)`
   !!     arrests the thin-film velocity runaway that develops at wet/dry
   !!     fronts on sloping beds (measured: u -> 4e3 m/s at h ~ 3.5 mm and
   !!     a ~2000x dt collapse without it).  It is far above any depth
   !!     ratio this scheme can resolve, so it never shapes the physics —
   !!     see the constant's docstring for the measurements.
   !!
   !! `eps` is the same DESING_EPS used by `compute_timestep`, so the
   !! update and the CFL estimator agree about what "thin" means.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, DESING_EPS, &
                            FROUDE_BACKSTOP
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: compute_update
   public :: compute_rk2_average
   public :: compute_ssp_combine

contains

   pure subroutine compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, &
                                  grid, dt)
      !! Forward Euler update U = U + dt*RHS (each SSP-RK2 stage). Also
      !! enforces non-negative depth, zeroes momentum in dry cells, and
      !! desingularizes thin-cell velocity (KP07) with a Froude backstop.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: flux_h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: flux_hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt
         !! Timestep (s)

      integer :: i, j, nx, ny, nghost
      real(wp) :: eps4, h4, u_val, v_val, speed_sq, fr_max_sq, fr_scale

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      eps4 = DESING_EPS**4

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(h4, u_val, v_val, speed_sq, fr_max_sq, fr_scale)
         h(i, j) = h(i, j) + dt*flux_h(i, j)
         hu(i, j) = hu(i, j) + dt*flux_hu(i, j)
         hv(i, j) = hv(i, j) + dt*flux_hv(i, j)

         ! Non-negative depth + dry-cell momentum zeroing; clamp h only
         ! for negative undershoots (preserves mass).  Thin-but-wet cells
         ! keep a desingularized (KP07) momentum instead of being zeroed.
         if (h(i, j) < 0.0_wp) then
            h(i, j) = 0.0_wp
            hu(i, j) = 0.0_wp
            hv(i, j) = 0.0_wp
         else if (h(i, j) < DRY_TOLERANCE) then
            hu(i, j) = 0.0_wp
            hv(i, j) = 0.0_wp
         else
            h4 = h(i, j)**4
            if (h4 < eps4) then
               ! u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4)); hu := h*u.
               ! Identity for h >= eps, so applied only below eps.
               u_val = sqrt(2.0_wp)*h(i, j)*hu(i, j)/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h(i, j)*hv(i, j)/sqrt(h4 + eps4)
               hu(i, j) = h(i, j)*u_val
               hv(i, j) = h(i, j)*v_val
            end if
            ! Froude backstop: blowup arrester only (triggers on a TRUE
            ! comparison, so NaN passes through and fails loudly).
            speed_sq = (hu(i, j)**2 + hv(i, j)**2)/h(i, j)**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h(i, j)
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu(i, j) = hu(i, j)*fr_scale
               hv(i, j) = hv(i, j)*fr_scale
            end if
         end if
      end do

   end subroutine compute_update

   pure subroutine compute_rk2_average(h, hu, hv, h0, hu0, hv0, grid)
      !! SSP-RK2 averaging step: U^{n+1} = 0.5*(U^n + U**), averaging the
      !! original state (h0, hu0, hv0) with the twice-advanced (h, hu, hv).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: h0(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu0(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv0(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny, nghost
      real(wp) :: eps4, h4, u_val, v_val, speed_sq, fr_max_sq, fr_scale

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      eps4 = DESING_EPS**4

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(h4, u_val, v_val, speed_sq, fr_max_sq, fr_scale)
         h(i, j) = 0.5_wp*(h0(i, j) + h(i, j))
         hu(i, j) = 0.5_wp*(hu0(i, j) + hu(i, j))
         hv(i, j) = 0.5_wp*(hv0(i, j) + hv(i, j))

         ! Same dry clamp + KP07 desingularization + backstop as
         ! compute_update.
         if (h(i, j) < 0.0_wp) then
            h(i, j) = 0.0_wp
            hu(i, j) = 0.0_wp
            hv(i, j) = 0.0_wp
         else if (h(i, j) < DRY_TOLERANCE) then
            hu(i, j) = 0.0_wp
            hv(i, j) = 0.0_wp
         else
            h4 = h(i, j)**4
            if (h4 < eps4) then
               u_val = sqrt(2.0_wp)*h(i, j)*hu(i, j)/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h(i, j)*hv(i, j)/sqrt(h4 + eps4)
               hu(i, j) = h(i, j)*u_val
               hv(i, j) = h(i, j)*v_val
            end if
            speed_sq = (hu(i, j)**2 + hv(i, j)**2)/h(i, j)**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h(i, j)
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu(i, j) = hu(i, j)*fr_scale
               hv(i, j) = hv(i, j)*fr_scale
            end if
         end if
      end do

   end subroutine compute_rk2_average

   pure subroutine compute_ssp_combine(h, hu, hv, h0, hu0, hv0, a_n, grid)
      !! Generalized SSP Shu-Osher combine: U = a_n*U^n + (1 - a_n)*U,
      !! where U is the post-Euler-stage iterate and (h0,hu0,hv0) = U^n.
      !! Reduces exactly to compute_rk2_average for a_n = 0.5.  Applies the
      !! same positivity + thin-cell treatment so each stage stays
      !! wet/dry-safe.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: h0(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu0(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv0(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: a_n
         !! Blend weight on the start-of-step state U^n.

      integer :: i, j, nx, ny, nghost
      real(wp) :: a_prev, eps4, h4, u_val, v_val, speed_sq, fr_max_sq, fr_scale

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      a_prev = 1.0_wp - a_n
      eps4 = DESING_EPS**4

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(h4, u_val, v_val, speed_sq, fr_max_sq, fr_scale)
         h(i, j) = a_n*h0(i, j) + a_prev*h(i, j)
         hu(i, j) = a_n*hu0(i, j) + a_prev*hu(i, j)
         hv(i, j) = a_n*hv0(i, j) + a_prev*hv(i, j)

         if (h(i, j) < 0.0_wp) then
            h(i, j) = 0.0_wp
            hu(i, j) = 0.0_wp
            hv(i, j) = 0.0_wp
         else if (h(i, j) < DRY_TOLERANCE) then
            hu(i, j) = 0.0_wp
            hv(i, j) = 0.0_wp
         else
            h4 = h(i, j)**4
            if (h4 < eps4) then
               u_val = sqrt(2.0_wp)*h(i, j)*hu(i, j)/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h(i, j)*hv(i, j)/sqrt(h4 + eps4)
               hu(i, j) = h(i, j)*u_val
               hv(i, j) = h(i, j)*v_val
            end if
            speed_sq = (hu(i, j)**2 + hv(i, j)**2)/h(i, j)**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h(i, j)
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu(i, j) = hu(i, j)*fr_scale
               hv(i, j) = hv(i, j)*fr_scale
            end if
         end if
      end do

   end subroutine compute_ssp_combine

end module rki_kernel_update
