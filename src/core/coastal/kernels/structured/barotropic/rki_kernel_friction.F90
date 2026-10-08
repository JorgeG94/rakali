!! Bottom friction, wind stress, atmospheric pressure, and precipitation kernel
module rki_kernel_friction
   !! Source terms for the momentum and continuity equations:
   !! Friction -c_f*u*sqrt(u^2+v^2), c_f = g*n^2/h^(1/3) (semi-implicit IMEX);
   !! wind tau_s/rho; pressure -(h/rho)*grad(p_atm); rain dh/dt += P (all explicit).
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, RHO_WATER
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: apply_friction
   public :: apply_coriolis

contains

   subroutine apply_friction(h, hu, hv, grid, dt, &
                                  manning_n, wind_x, wind_y, &
                                  atm_pressure, precipitation, &
                                  forcing_active, include_wind)
      !! Apply bottom friction, wind stress, atmospheric pressure gradient,
      !! and precipitation to (h, hu, hv).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt
         !! Timestep (s)
      real(wp), intent(in) :: manning_n
         !! Manning roughness coefficient (s/m^{1/3})
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: atm_pressure(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: precipitation(grid%nx_total, grid%ny_total)
      logical, intent(in), optional :: forcing_active
         !! When .false., skip wind/pressure/precipitation (default: .true.)
      logical, intent(in), optional :: include_wind
         !! When .false., skip the wind stress (it is applied to the surface
         !! layer instead — wind_surface_layer).  Default: .true.

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: u_val, v_val, speed, cf, denom, h_val
      real(wp) :: dpdx, dpdy
      logical :: has_friction, has_forcing, apply_wind

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      has_friction = (manning_n > 0.0_wp)
      has_forcing = .true.
      if (present(forcing_active)) has_forcing = forcing_active
      apply_wind = .true.
      if (present(include_wind)) apply_wind = include_wind

      ! Early exit: no friction and no atmospheric forcing
      if (.not. has_friction .and. .not. has_forcing) return

      if (.not. has_forcing) then
         ! Fast path: friction only, no forcing array reads
         !$omp parallel do collapse(2) private(u_val, v_val, speed, cf, denom, h_val)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            h_val = h(i, j)
            if (h_val > DRY_TOLERANCE) then
               u_val = hu(i, j)/h_val
               v_val = hv(i, j)/h_val
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*manning_n*manning_n/(h_val**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h_val
               hu(i, j) = hu(i, j)/denom
               hv(i, j) = hv(i, j)/denom
            end if
         end do
         end do
         !$omp end parallel do
      else
         ! Full path: friction + wind + pressure + precipitation
         !$omp parallel do collapse(2) private(u_val, v_val, speed, cf, denom, h_val, dpdx, dpdy)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            h_val = h(i, j)

            if (h_val > DRY_TOLERANCE) then
               u_val = hu(i, j)/h_val
               v_val = hv(i, j)/h_val
               speed = sqrt(u_val*u_val + v_val*v_val)

               if (manning_n > 0.0_wp) then
                  cf = GRAVITY*manning_n*manning_n/(h_val**(1.0_wp/3.0_wp))

                  ! Semi-implicit: divide by (1 + dt*cf*|u|/h)
                  denom = 1.0_wp + dt*cf*speed/h_val
                  hu(i, j) = hu(i, j)/denom
                  hv(i, j) = hv(i, j)/denom
               end if

               ! Wind stress tau/rho; skipped when wind goes to the surface layer.
               if (apply_wind) then
                  hu(i, j) = hu(i, j) + dt*wind_x(i, j)/RHO_WATER
                  hv(i, j) = hv(i, j) + dt*wind_y(i, j)/RHO_WATER
               end if

               ! Atmospheric pressure gradient -(h/rho)*dp/dx, central diff.
               dpdx = (atm_pressure(i + 1, j) - atm_pressure(i - 1, j))/(2.0_wp*dx)
               dpdy = (atm_pressure(i, j + 1) - atm_pressure(i, j - 1))/(2.0_wp*dy)
               hu(i, j) = hu(i, j) - dt*h_val*dpdx/RHO_WATER
               hv(i, j) = hv(i, j) - dt*h_val*dpdy/RHO_WATER
            end if

            ! Precipitation: mass source (applied to all cells including dry)
            h(i, j) = h(i, j) + dt*precipitation(i, j)
         end do
         end do
         !$omp end parallel do
      end if

   end subroutine apply_friction

   subroutine apply_coriolis(h, hu, hv, grid, dt, f)
      !! Apply Coriolis via the exact (energy-conserving) rotation:
      !!   hu_new =  hu*cos(f*dt) + hv*sin(f*dt)
      !!   hv_new = -hu*sin(f*dt) + hv*cos(f*dt)
      type(hgrid_t), intent(in) :: grid
      ! assumed-shape-ok: Coriolis kernel is on the coastal barotropic path
      ! (not the ocean dycore hot path); h/hu/hv grid sizes are read from
      ! grid%nx_total/ny_total, not from the array descriptors.
      real(wp), intent(in) :: h(:, :)
      real(wp), intent(inout) :: hu(:, :), hv(:, :)  ! assumed-shape-ok: coastal barotropic; grid dims from grid%nx_total
      real(wp), intent(in) :: dt, f

      integer :: i, j, nx, ny, nghost
      real(wp) :: cosf, sinf, hu_old, hv_old

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (abs(f) < 1.0e-20_wp) return

      cosf = cos(f*dt)
      sinf = sin(f*dt)

      !$omp parallel do collapse(2) private(hu_old, hv_old)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         if (h(i, j) > DRY_TOLERANCE) then
            hu_old = hu(i, j)
            hv_old = hv(i, j)
            hu(i, j) = hu_old*cosf + hv_old*sinf
            hv(i, j) = -hu_old*sinf + hv_old*cosf
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine apply_coriolis

end module rki_kernel_friction
