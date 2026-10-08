!! CFL-based adaptive timestep computation kernel
module rki_kernel_timestep
   !! Maximum stable timestep from the CFL condition:
   !! dt = cfl * min(dx, dy) / max_wavespeed,
   !! max_wavespeed = max(|u| + sqrt(g*h), |v| + sqrt(g*h)).
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, DESING_EPS
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: compute_timestep
   public :: compute_timestep_baroclinic

contains

   subroutine compute_timestep(h, hu, hv, grid, cfl, dt)
      !! Maximum stable timestep: reduces over interior cells for the max
      !! wave speed, then dt = cfl * min(dx,dy) / max_wavespeed.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: cfl
         !! CFL number
      real(wp), intent(out) :: dt
         !! Computed timestep (s)

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: max_speed, speed_x, speed_y, u_abs, v_abs, wave_speed
      real(wp) :: h4, eps4, denom, eps

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      ! Desingularization parameter (Kurganov & Petrova, 2007): the SAME
      ! eps as the update kernels' thin-cell desingularization, so the dt
      ! bound and the update agree about what "thin" means.  (The old
      ! `max(eps, 0.01*min(dx,dy))` grid term is gone: on coarse grids it
      ! made this estimator blind to real velocities in cells shallower
      ! than 1% of the cell width — a latent CFL bug.)  The denominator
      ! here omits the KP07 `max(h^4, eps^4)`, so the speed estimated
      ! below is >= the desingularized velocity the update kernels produce
      ! for every h (by up to sqrt(2) in deep water) — a strict, mildly
      ! conservative upper bound on the CFL.
      eps = DESING_EPS
      eps4 = eps**4
      max_speed = 0.0_wp

      !$omp parallel do collapse(2) private(wave_speed, h4, denom, u_abs, v_abs, speed_x, speed_y)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         if (h(i, j) > DRY_TOLERANCE) then
            wave_speed = sqrt(GRAVITY*h(i, j))

            ! Desingularized velocity (Kurganov & Petrova, 2007):
            !   u = sqrt(2)*h*hu / sqrt(h^4 + eps^4)
            h4 = h(i, j)**4
            denom = sqrt(h4 + eps4)
            u_abs = abs(hu(i, j))*sqrt(2.0_wp)*h(i, j)/denom
            v_abs = abs(hv(i, j))*sqrt(2.0_wp)*h(i, j)/denom

            speed_x = u_abs + wave_speed
            speed_y = v_abs + wave_speed
            max_speed = max(max_speed, speed_x/dx, speed_y/dy)
         end if
      end do
      end do
      !$omp end parallel do

      if (max_speed > 0.0_wp) then
         dt = cfl/max_speed
      else
         dt = cfl*min(dx, dy)/sqrt(GRAVITY*DRY_TOLERANCE)
      end if

   end subroutine compute_timestep

   subroutine compute_timestep_baroclinic(h, u_layer, v_layer, dsig, &
                                          grid, nz, cfl, dt)
      !! Baroclinic timestep from the internal-wave CFL:
      !! max(|u_layer_k| + c_int)/dx, c_int = sqrt(g*H/nz) (uniform sigma).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: cfl
      real(wp), intent(out) :: dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: max_speed, c_int, speed_x, speed_y

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      max_speed = 0.0_wp

      !$omp parallel do collapse(2) private(c_int, speed_x, speed_y, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         if (h(i, j) > DRY_TOLERANCE) then
            c_int = sqrt(GRAVITY*h(i, j)/real(nz, wp))
            do k = 1, nz
               speed_x = (abs(u_layer(i, j, k)) + c_int)/dx
               speed_y = (abs(v_layer(i, j, k)) + c_int)/dy
               max_speed = max(max_speed, speed_x, speed_y)
            end do
         end if
      end do
      end do
      !$omp end parallel do

      if (max_speed > 0.0_wp) then
         dt = cfl/max_speed
      else
         dt = cfl*min(dx, dy)/sqrt(GRAVITY*DRY_TOLERANCE/real(nz, wp))
      end if

   end subroutine compute_timestep_baroclinic

end module rki_kernel_timestep
