!! Sponge layer relaxation kernel
module rki_kernel_sponge
   !! Volumetric damping in interior cells near sponge boundaries, with a
   !! quadratic profile alpha(d) = strength*((width - d)/width)^2 where d
   !! is the distance in cells from the nearest sponge boundary.
   use rki_constants, only: wp
   use rki_boundary_types, only: bc_params_t, BC_SPONGE
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: apply_sponge

contains

   subroutine apply_sponge(h, hu, hv, grid, dt, bc)
      !! Apply sponge layer relaxation to h, hu, hv
      !!
      !! Cells within sponge_width of a sponge boundary are relaxed toward
      !! h_ref with hu_ref = hv_ref = 0. No-op when sponge_width <= 0.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt
         !! Timestep (s)
      type(bc_params_t), intent(in) :: bc
         !! Boundary condition parameters (sponge config)

      integer :: i, j, i_first, i_last, j_first, j_last
      integer :: nx, ny, nghost
      integer :: sw, d_west, d_east, d_south, d_north, d_min
      real(wp) :: alpha, ratio, h_ref, strength
      logical :: has_west, has_east, has_south, has_north

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      sw = bc%sponge_width
      if (sw <= 0) return

      strength = bc%sponge_strength
      h_ref = bc%sponge_h_ref
      if (strength <= 0.0_wp) return

      i_first = nghost + 1
      i_last = nx - nghost
      j_first = nghost + 1
      j_last = ny - nghost

      has_west = (bc%west == BC_SPONGE)
      has_east = (bc%east == BC_SPONGE)
      has_south = (bc%south == BC_SPONGE)
      has_north = (bc%north == BC_SPONGE)

      !$omp target teams distribute parallel do collapse(2) private(d_west, d_east, d_south, d_north, d_min, alpha, ratio)
      do j = j_first, j_last
      do i = i_first, i_last
         ! Compute distance from each sponge boundary (in cells)
         d_min = sw + 1  ! Start beyond sponge width (no damping)

         if (has_west) then
            d_west = i - i_first
            if (d_west < sw) d_min = min(d_min, d_west)
         end if
         if (has_east) then
            d_east = i_last - i
            if (d_east < sw) d_min = min(d_min, d_east)
         end if
         if (has_south) then
            d_south = j - j_first
            if (d_south < sw) d_min = min(d_min, d_south)
         end if
         if (has_north) then
            d_north = j_last - j
            if (d_north < sw) d_min = min(d_min, d_north)
         end if

         if (d_min < sw) then
            ratio = real(sw - d_min, wp)/real(sw, wp)
            alpha = strength*ratio*ratio
            h(i, j) = h(i, j) - dt*alpha*(h(i, j) - h_ref)
            hu(i, j) = hu(i, j) - dt*alpha*hu(i, j)
            hv(i, j) = hv(i, j) - dt*alpha*hv(i, j)
         end if
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine apply_sponge

end module rki_kernel_sponge
