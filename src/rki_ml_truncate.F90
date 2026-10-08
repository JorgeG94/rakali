!! CFL-based velocity truncation for layer momentum (MOM6-style limiter).
module rki_ml_truncate
   !! Per-cell-per-layer limiter (not a closure) that clips runaway velocities
   !! from any upstream operator (BPG, vertical exchange, viscosity) to
   !! 0.9*cfl_trunc*dx/dt so a few bad cells don't tank the global timestep.
   use rki_constants, only: wp, THIN_LAYER_THRESHOLD
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: ml_truncate_velocity

contains

   pure subroutine ml_truncate_velocity(h_layer, hu_layer, hv_layer, &
                                        grid, nz, dt, &
                                        cfl_trunc)
      !! Clip any layer velocity with CFL > cfl_trunc to 0.9*cfl_trunc*dx/dt.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: cfl_trunc

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: u_abs, v_abs, h_k, vel_max_x, vel_max_y

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      vel_max_x = 0.9_wp*cfl_trunc*dx/max(dt, 1.0e-10_wp)
      vel_max_y = 0.9_wp*cfl_trunc*dy/max(dt, 1.0e-10_wp)

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, u_abs, v_abs, h_k)

         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k < THIN_LAYER_THRESHOLD) cycle

            u_abs = abs(hu_layer(i, j, k))/h_k
            v_abs = abs(hv_layer(i, j, k))/h_k

            if (u_abs > vel_max_x) then
               hu_layer(i, j, k) = sign(vel_max_x*h_k, hu_layer(i, j, k))
            end if
            if (v_abs > vel_max_y) then
               hv_layer(i, j, k) = sign(vel_max_y*h_k, hv_layer(i, j, k))
            end if
         end do
      end do

   end subroutine ml_truncate_velocity

end module rki_ml_truncate
