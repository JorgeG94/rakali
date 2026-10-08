!! Per-layer Manning bottom friction and per-layer Coriolis rotation.
module rki_ml_friction
   !! Manning quadratic bottom friction on the bottom layer (k=1) and per-layer
   !! Coriolis rotation of every layer's momentum (`ml_apply_coriolis`).
   !!
   !! Per-layer (vs depth-averaged) friction/Coriolis damps and rotates the
   !! baroclinic shear (gravity currents, dense overflows) using each layer's
   !! own velocity. Replaces the depth-averaged `apply_coriolis` on the
   !! multilayer path; the barotropic-only path keeps the depth-averaged version.
   !!
   !! Layers are bottom-up (k=1 bed, k=nz surface): bed friction acts at k=1.
   !! Semi-implicit friction: hu_new = hu_old / (1 + dt*C_d*|u_b|/h_k),
   !! C_d = g*n^2/H^{1/3}.
   use rki_constants, only: wp, GRAVITY, THIN_LAYER_THRESHOLD
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: VANISHING_LAYER_TOL
   implicit none
   private

   public :: ml_apply_bottom_friction
   public :: ml_apply_coriolis

contains

   subroutine ml_apply_bottom_friction(h_layer, hu_layer, hv_layer, &
                                            manning_n, use_linear_drag, drag_r, &
                                            grid, nz, dt)
      !! Apply bottom friction to the lowest active layer (k=1, bottom-up) using
      !! that layer's velocity. Two semi-implicit drag laws:
      !!   quadratic Manning (default):  hu_new = hu / (1 + dt*C_d*|u_b|/h_k),
      !!     C_d = g*n^2/H^{1/3}, H = total column depth.
      !!   linear Rayleigh (`use_linear_drag`):  hu_new = hu / (1 + dt*r).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: manning_n, dt
      logical, intent(in) :: use_linear_drag
         !! Use linear (Rayleigh) drag with rate `drag_r` instead of Manning.
      real(wp), intent(in) :: drag_r
         !! Linear bottom-drag rate (1/s); only read when use_linear_drag.

      integer :: i, j, k_bot, kk, nx, ny, nghost
      real(wp) :: h_k, h_total, u_b, v_b, speed_b, cf, denom

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (use_linear_drag) then
         if (drag_r <= 0.0_wp) return
      else
         if (manning_n <= 0.0_wp) return
      end if

      !$omp parallel do collapse(2) private(k_bot, kk, h_k, h_total, u_b, v_b, speed_b, cf, denom)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost

         ! Find the lowest active layer (k_bot).  For VCOORD_SIGMA /
         ! ZSIGMA / ZSTAR (lite), k=1 is always active (the deepest
         ! non-vanishing layer = k=1).  For VCOORD_ZSTAR_FULL the
         ! bottom layer can be parked at h_min when the column is
         ! shallower than the reference; in that case skip up to the
         ! lowest layer with thickness above the vanishing threshold.
         k_bot = 1
         do kk = 1, nz
            if (h_layer(i, j, kk) > VANISHING_LAYER_TOL) then
               k_bot = kk
               exit
            end if
         end do

         h_k = h_layer(i, j, k_bot)
         if (h_k < THIN_LAYER_THRESHOLD) cycle

         ! Total column depth for Manning coefficient
         h_total = 0.0_wp
         do kk = 1, nz
            h_total = h_total + h_layer(i, j, kk)
         end do
         if (h_total < THIN_LAYER_THRESHOLD) cycle

         u_b = hu_layer(i, j, k_bot)/h_k
         v_b = hv_layer(i, j, k_bot)/h_k
         speed_b = sqrt(u_b*u_b + v_b*v_b)

         ! Semi-implicit divide: linear (constant rate) or Manning (speed-
         ! and depth-dependent C_d).
         if (use_linear_drag) then
            denom = 1.0_wp + dt*drag_r
         else
            cf = GRAVITY*manning_n*manning_n/(h_total**(1.0_wp/3.0_wp))
            denom = 1.0_wp + dt*cf*speed_b/h_k
         end if
         hu_layer(i, j, k_bot) = hu_layer(i, j, k_bot)/denom
         hv_layer(i, j, k_bot) = hv_layer(i, j, k_bot)/denom
      end do
      end do
      !$omp end parallel do

   end subroutine ml_apply_bottom_friction

   subroutine ml_apply_coriolis(h_layer, hu_layer, hv_layer, grid, nz, dt, f)
      !! Per-layer Coriolis rotation: rotate every layer's `(hu_k, hv_k)`
      !! through angle `f*dt` via the closed-form (energy-conserving, det=1)
      !! rotation matrix
      !!   hu_new =  hu_old*cos(f*dt) + hv_old*sin(f*dt)
      !!   hv_new = -hu_old*sin(f*dt) + hv_old*cos(f*dt)
      !! `f` is a scalar (latitude-only). Skips layers below the dry-cell
      !! threshold so already-zeroed momentum is untouched.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dt, f

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: cosf, sinf, hu_old, hv_old

      if (f == 0.0_wp) return

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      cosf = cos(f*dt)
      sinf = sin(f*dt)

      !$omp parallel do collapse(2) private(k, hu_old, hv_old)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         do k = 1, nz
            if (h_layer(i, j, k) > THIN_LAYER_THRESHOLD) then
               hu_old = hu_layer(i, j, k)
               hv_old = hv_layer(i, j, k)
               hu_layer(i, j, k) = hu_old*cosf + hv_old*sinf
               hv_layer(i, j, k) = -hu_old*sinf + hv_old*cosf
            end if
         end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_apply_coriolis

end module rki_ml_friction
