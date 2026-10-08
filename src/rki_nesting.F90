!! One-way nesting from a coarser parent domain
module rki_nesting
   !! Reads solution (h, hu, hv) from a parent domain NetCDF output file
   !! and fills ghost cells on BC_NESTED edges using Flather radiation with
   !! spatially-varying external values. Parent file: dims (x, y, time),
   !! variables h/hu/hv(x,y,time) + coordinate variables. Interpolation is
   !! linear along each edge in space and between bracketing snapshots in time.
   use rki_constants, only: wp, GRAVITY
   use rki_config, only: config_t
   use rki_state, only: state_t
   use rki_boundary_types, only: bc_params_t, BC_NESTED
   use rki_io_netcdf, only: nc_check, nc_open_read, nc_close, &
                            nc_get_dim_len, nc_get_varid
   use netcdf, only: nf90_get_var
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: nesting_t
   public :: nesting_init
   public :: nesting_init_from_config
   public :: nesting_update
   public :: nesting_destroy

   type :: nesting_t
      !! Container for parent domain data and interpolation state

      ! NetCDF file handle
      integer :: ncid = -1
         !! NetCDF file ID (kept open for lazy reads)
      logical :: active = .false.
         !! Whether nesting is enabled

      ! Parent grid coordinates
      integer :: p_nx = 0
         !! Number of x-points in parent file
      integer :: p_ny = 0
         !! Number of y-points in parent file
      integer :: p_nt = 0
         !! Number of time steps in parent file
      real(wp), allocatable :: p_x(:)
         !! Parent x-coordinates (m)
      real(wp), allocatable :: p_y(:)
         !! Parent y-coordinates (m)
      real(wp), allocatable :: p_time(:)
         !! Parent time coordinates (s)

      ! Variable IDs
      integer :: varid_h = -1
      integer :: varid_hu = -1
      integer :: varid_hv = -1

      ! Two time-bracketing snapshots on parent grid
      integer :: idx_lo = 0
      integer :: idx_hi = 0
      real(wp), allocatable :: h_lo(:, :), h_hi(:, :)
      real(wp), allocatable :: hu_lo(:, :), hu_hi(:, :)
      real(wp), allocatable :: hv_lo(:, :), hv_hi(:, :)

      ! Child-to-parent coordinate offset
      real(wp) :: x_origin = 0.0_wp
         !! Child grid x-origin in parent coordinates
      real(wp) :: y_origin = 0.0_wp
         !! Child grid y-origin in parent coordinates

      ! Precomputed interpolation indices and weights for each boundary edge
      ! West/East edges: indexed by child j (1:ny_total)
      integer, allocatable :: we_ix0(:)
         !! Parent x-index for west/east edge interpolation
      real(wp), allocatable :: we_wx(:)
         !! Parent x-weight for west/east edge
      integer, allocatable :: we_iy0(:)
         !! Lower parent y-index per child j
      real(wp), allocatable :: we_wy(:)
         !! Parent y-weight per child j

      ! South/North edges: indexed by child i (1:nx_total)
      integer, allocatable :: sn_iy0(:)
         !! Parent y-index for south/north edge interpolation
      real(wp), allocatable :: sn_wy(:)
         !! Parent y-weight for south/north edge
      integer, allocatable :: sn_ix0(:)
         !! Lower parent x-index per child i
      real(wp), allocatable :: sn_wx(:)
         !! Parent x-weight per child i

      ! Parent x-indices for west and east edges
      integer :: ix_west = 1
      real(wp) :: wx_west = 0.0_wp
      integer :: ix_east = 1
      real(wp) :: wx_east = 0.0_wp

      ! Parent y-indices for south and north edges
      integer :: iy_south = 1
      real(wp) :: wy_south = 0.0_wp
      integer :: iy_north = 1
      real(wp) :: wy_north = 0.0_wp

   end type nesting_t

contains

   subroutine nesting_init(nest, filename, state, x_origin, y_origin)
      !! Initialise nesting from a parent domain NetCDF file.
      !! Public only for the unit-test suite.
      type(nesting_t), intent(inout) :: nest
      character(len=*), intent(in) :: filename
      type(state_t), intent(in) :: state
      real(wp), intent(in), optional :: x_origin
      real(wp), intent(in), optional :: y_origin

      integer :: ncid

      if (len_trim(filename) == 0) return

      nest%x_origin = 0.0_wp
      nest%y_origin = 0.0_wp
      if (present(x_origin)) nest%x_origin = x_origin
      if (present(y_origin)) nest%y_origin = y_origin

      call logger%info("Loading nesting data from: "//trim(filename))

      call nc_open_read(filename, ncid)
      nest%ncid = ncid

      ! Read parent grid
      call read_parent_grid(nest, ncid)

      ! Find variable IDs
      call nc_get_varid(ncid, "h", nest%varid_h)
      call nc_get_varid(ncid, "hu", nest%varid_hu)
      call nc_get_varid(ncid, "hv", nest%varid_hv)

      ! Allocate snapshot arrays
      call allocate_parent_snapshots(nest)

      ! Precompute interpolation weights for boundary edges
      call compute_edge_weights(nest, state)

      ! Load first two snapshots
      if (nest%p_nt >= 2) then
         nest%idx_lo = 1
         nest%idx_hi = 2
         call load_parent_snapshot(nest, ncid, 1, &
                                   nest%h_lo, nest%hu_lo, nest%hv_lo)
         call load_parent_snapshot(nest, ncid, 2, &
                                   nest%h_hi, nest%hu_hi, nest%hv_hi)
      else if (nest%p_nt == 1) then
         nest%idx_lo = 1
         nest%idx_hi = 1
         call load_parent_snapshot(nest, ncid, 1, &
                                   nest%h_lo, nest%hu_lo, nest%hv_lo)
         nest%h_hi = nest%h_lo
         nest%hu_hi = nest%hu_lo
         nest%hv_hi = nest%hv_lo
      end if

      nest%active = .true.

      call logger%info("Parent grid: "//to_string(nest%p_nx)//" x "// &
                       to_string(nest%p_ny)//", "// &
                       to_string(nest%p_nt)//" time steps")
      call logger%info("Parent time range: "//to_string(nest%p_time(1))// &
                       " to "//to_string(nest%p_time(nest%p_nt))//" s")

   end subroutine nesting_init

   subroutine nesting_init_from_config(nest, cfg, state)
      !! Initialise nesting from config parameters
      type(nesting_t), intent(inout) :: nest
      type(config_t), intent(in) :: cfg
      type(state_t), intent(in) :: state

      call nesting_init(nest, trim(cfg%nesting_file), state, &
                        cfg%nesting_x_origin, cfg%nesting_y_origin)

   end subroutine nesting_init_from_config

   subroutine nesting_update(nest, state, sim_time)
      !! Fill ghost cells on nested edges by interpolating from the parent,
      !! using Flather radiation with spatially-varying external values.
      type(nesting_t), intent(inout) :: nest
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: sim_time

      real(wp) :: alpha

      if (.not. nest%active) return

      ! Check if any edge is nested
      if (state%bc%west /= BC_NESTED .and. state%bc%east /= BC_NESTED .and. &
          state%bc%south /= BC_NESTED .and. state%bc%north /= BC_NESTED) return

      ! Advance time bracket
      call advance_parent_bracket(nest, sim_time)

      ! Temporal interpolation weight
      if (nest%idx_lo == nest%idx_hi) then
         alpha = 0.0_wp
      else
         alpha = (sim_time - nest%p_time(nest%idx_lo))/ &
                 (nest%p_time(nest%idx_hi) - nest%p_time(nest%idx_lo))
         alpha = max(0.0_wp, min(1.0_wp, alpha))
      end if

      ! Apply nested BCs per edge
      if (state%bc%west == BC_NESTED) then
         call apply_nested_west(nest, state, alpha)
      end if
      if (state%bc%east == BC_NESTED) then
         call apply_nested_east(nest, state, alpha)
      end if
      if (state%bc%south == BC_NESTED) then
         call apply_nested_south(nest, state, alpha)
      end if
      if (state%bc%north == BC_NESTED) then
         call apply_nested_north(nest, state, alpha)
      end if

   end subroutine nesting_update

   subroutine nesting_destroy(nest)
      !! Clean up nesting state and close file
      type(nesting_t), intent(inout) :: nest

      if (nest%ncid >= 0) then
         call nc_close(nest%ncid)
         nest%ncid = -1
      end if

      if (allocated(nest%p_x)) deallocate (nest%p_x)
      if (allocated(nest%p_y)) deallocate (nest%p_y)
      if (allocated(nest%p_time)) deallocate (nest%p_time)
      if (allocated(nest%h_lo)) deallocate (nest%h_lo)
      if (allocated(nest%h_hi)) deallocate (nest%h_hi)
      if (allocated(nest%hu_lo)) deallocate (nest%hu_lo)
      if (allocated(nest%hu_hi)) deallocate (nest%hu_hi)
      if (allocated(nest%hv_lo)) deallocate (nest%hv_lo)
      if (allocated(nest%hv_hi)) deallocate (nest%hv_hi)
      if (allocated(nest%we_iy0)) deallocate (nest%we_iy0)
      if (allocated(nest%we_wy)) deallocate (nest%we_wy)
      if (allocated(nest%sn_ix0)) deallocate (nest%sn_ix0)
      if (allocated(nest%sn_wx)) deallocate (nest%sn_wx)

      nest%active = .false.

   end subroutine nesting_destroy

   ! ================================================================
   ! Private helpers
   ! ================================================================

   subroutine read_parent_grid(nest, ncid)
      type(nesting_t), intent(inout) :: nest
      integer, intent(in) :: ncid

      integer :: varid

      call nc_get_dim_len(ncid, "x", nest%p_nx)
      call nc_get_dim_len(ncid, "y", nest%p_ny)
      call nc_get_dim_len(ncid, "time", nest%p_nt)

      allocate (nest%p_x(nest%p_nx))
      allocate (nest%p_y(nest%p_ny))
      allocate (nest%p_time(nest%p_nt))

      call nc_get_varid(ncid, "x", varid)
      call nc_check(nf90_get_var(ncid, varid, nest%p_x), "reading parent x")

      call nc_get_varid(ncid, "y", varid)
      call nc_check(nf90_get_var(ncid, varid, nest%p_y), "reading parent y")

      call nc_get_varid(ncid, "time", varid)
      call nc_check(nf90_get_var(ncid, varid, nest%p_time), "reading parent time")

   end subroutine read_parent_grid

   subroutine allocate_parent_snapshots(nest)
      type(nesting_t), intent(inout) :: nest

      integer :: nx, ny

      nx = nest%p_nx
      ny = nest%p_ny

      allocate (nest%h_lo(nx, ny), source=0.0_wp)
      allocate (nest%h_hi(nx, ny), source=0.0_wp)
      allocate (nest%hu_lo(nx, ny), source=0.0_wp)
      allocate (nest%hu_hi(nx, ny), source=0.0_wp)
      allocate (nest%hv_lo(nx, ny), source=0.0_wp)
      allocate (nest%hv_hi(nx, ny), source=0.0_wp)

   end subroutine allocate_parent_snapshots

   subroutine compute_edge_weights(nest, state)
      !! Precompute interpolation weights mapping child edge cells to parent grid
      type(nesting_t), intent(inout) :: nest
      type(state_t), intent(in) :: state

      integer :: i, j, nx, ny, ng
      real(wp) :: x_child, y_child

      nx = state%grid%nx_total
      ny = state%grid%ny_total
      ng = state%grid%nghost

      ! West/East edges: need y-interpolation along j, and x-position for the edge
      allocate (nest%we_iy0(ny), source=1)
      allocate (nest%we_wy(ny), source=0.0_wp)

      do j = 1, ny
         y_child = nest%y_origin + (real(j - ng, wp) - 0.5_wp)*state%grid%dy
         call find_bracket_1d(nest%p_y, nest%p_ny, y_child, &
                              nest%we_iy0(j), nest%we_wy(j))
      end do

      ! x-position for west edge (first physical cell centre)
      x_child = nest%x_origin + 0.5_wp*state%grid%dx
      call find_bracket_1d(nest%p_x, nest%p_nx, x_child, &
                           nest%ix_west, nest%wx_west)

      ! x-position for east edge (last physical cell centre)
      x_child = nest%x_origin + (real(state%grid%nx_phys, wp) - 0.5_wp)*state%grid%dx
      call find_bracket_1d(nest%p_x, nest%p_nx, x_child, &
                           nest%ix_east, nest%wx_east)

      ! South/North edges: need x-interpolation along i, and y-position for the edge
      allocate (nest%sn_ix0(nx), source=1)
      allocate (nest%sn_wx(nx), source=0.0_wp)

      do i = 1, nx
         x_child = nest%x_origin + (real(i - ng, wp) - 0.5_wp)*state%grid%dx
         call find_bracket_1d(nest%p_x, nest%p_nx, x_child, &
                              nest%sn_ix0(i), nest%sn_wx(i))
      end do

      ! y-position for south edge
      y_child = nest%y_origin + 0.5_wp*state%grid%dy
      call find_bracket_1d(nest%p_y, nest%p_ny, y_child, &
                           nest%iy_south, nest%wy_south)

      ! y-position for north edge
      y_child = nest%y_origin + (real(state%grid%ny_phys, wp) - 0.5_wp)*state%grid%dy
      call find_bracket_1d(nest%p_y, nest%p_ny, y_child, &
                           nest%iy_north, nest%wy_north)

   end subroutine compute_edge_weights

   pure subroutine find_bracket_1d(coords, n, val, idx, frac)
      !! Find bracketing interval: coords(idx) <= val < coords(idx+1)
      real(wp), intent(in) :: coords(:)
      integer, intent(in) :: n
      real(wp), intent(in) :: val
      integer, intent(out) :: idx
      real(wp), intent(out) :: frac

      integer :: k

      if (n <= 1) then
         idx = 1
         frac = 0.0_wp
         return
      end if

      if (val <= coords(1)) then
         idx = 1
         frac = 0.0_wp
         return
      end if

      if (val >= coords(n)) then
         idx = n
         frac = 0.0_wp
         return
      end if

      do k = 1, n - 1
         if (val >= coords(k) .and. val < coords(k + 1)) then
            idx = k
            frac = (val - coords(k))/(coords(k + 1) - coords(k))
            return
         end if
      end do

      idx = n - 1
      frac = 1.0_wp

   end subroutine find_bracket_1d

   subroutine load_parent_snapshot(nest, ncid, tidx, h, hu, hv)
      type(nesting_t), intent(in) :: nest
      integer, intent(in) :: ncid, tidx
      real(wp), intent(out) :: h(:, :), hu(:, :), hv(:, :)

      integer :: start3(3), count3(3)

      start3 = [1, 1, tidx]
      count3 = [nest%p_nx, nest%p_ny, 1]

      call nc_check(nf90_get_var(ncid, nest%varid_h, h, &
                                 start=start3, count=count3), "reading parent h")
      call nc_check(nf90_get_var(ncid, nest%varid_hu, hu, &
                                 start=start3, count=count3), "reading parent hu")
      call nc_check(nf90_get_var(ncid, nest%varid_hv, hv, &
                                 start=start3, count=count3), "reading parent hv")

   end subroutine load_parent_snapshot

   subroutine advance_parent_bracket(nest, sim_time)
      type(nesting_t), intent(inout) :: nest
      real(wp), intent(in) :: sim_time

      do while (nest%idx_hi < nest%p_nt .and. &
                sim_time >= nest%p_time(nest%idx_hi))
         nest%idx_lo = nest%idx_hi
         nest%h_lo = nest%h_hi
         nest%hu_lo = nest%hu_hi
         nest%hv_lo = nest%hv_hi

         nest%idx_hi = nest%idx_hi + 1
         call load_parent_snapshot(nest, nest%ncid, nest%idx_hi, &
                                   nest%h_hi, nest%hu_hi, nest%hv_hi)
      end do

   end subroutine advance_parent_bracket

   pure function interp_parent(fld_lo, fld_hi, ix, iy, wx, wy, alpha, pnx, pny) result(val)
      !! Bilinear spatial + linear temporal interpolation of parent field
      real(wp), intent(in) :: fld_lo(:, :), fld_hi(:, :)
      integer, intent(in) :: ix, iy, pnx, pny
      real(wp), intent(in) :: wx, wy, alpha
      real(wp) :: val

      integer :: ix1, iy1
      real(wp) :: v_lo, v_hi

      ix1 = min(ix + 1, pnx)
      iy1 = min(iy + 1, pny)

      v_lo = (1.0_wp - wx)*(1.0_wp - wy)*fld_lo(ix, iy) + &
             wx*(1.0_wp - wy)*fld_lo(ix1, iy) + &
             (1.0_wp - wx)*wy*fld_lo(ix, iy1) + &
             wx*wy*fld_lo(ix1, iy1)

      v_hi = (1.0_wp - wx)*(1.0_wp - wy)*fld_hi(ix, iy) + &
             wx*(1.0_wp - wy)*fld_hi(ix1, iy) + &
             (1.0_wp - wx)*wy*fld_hi(ix, iy1) + &
             wx*wy*fld_hi(ix1, iy1)

      val = (1.0_wp - alpha)*v_lo + alpha*v_hi

   end function interp_parent

   subroutine apply_nested_west(nest, state, alpha)
      !! Fill west ghost cells using Flather with parent values: eta_ext
      !! prescribed from parent, h_ghost = max(0, eta_ext - b_ghost),
      !! hu_ghost = h_ghost * (u_int - c*(eta_int - eta_ext)).
      type(nesting_t), intent(in) :: nest
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: alpha

      integer :: j, k, ng, ny, i_first
      integer :: iy0, pnx, pny
      real(wp) :: wy_j, h_ext, hu_ext, hv_ext, eta_ext
      real(wp) :: h_int, eta_int, speed

      ng = state%grid%nghost
      ny = state%grid%ny_total
      i_first = ng + 1
      pnx = nest%p_nx
      pny = nest%p_ny

      do j = 1, ny
         iy0 = nest%we_iy0(j)
         wy_j = nest%we_wy(j)

         ! Interpolate parent values at west edge
         h_ext = interp_parent(nest%h_lo, nest%h_hi, &
                               nest%ix_west, iy0, nest%wx_west, wy_j, &
                               alpha, pnx, pny)
         hu_ext = interp_parent(nest%hu_lo, nest%hu_hi, &
                                nest%ix_west, iy0, nest%wx_west, wy_j, &
                                alpha, pnx, pny)
         hv_ext = interp_parent(nest%hv_lo, nest%hv_hi, &
                                nest%ix_west, iy0, nest%wx_west, wy_j, &
                                alpha, pnx, pny)

         eta_ext = h_ext + state%barotropic%b(i_first, j)

         ! Interior values
         h_int = state%barotropic%h(i_first, j)
         eta_int = h_int + state%barotropic%b(i_first, j)
         speed = sqrt(GRAVITY*max(h_int, 1.0e-6_wp))

         ! Fill ghost cells with Flather BC
         do k = 1, ng
            state%barotropic%h(i_first - k, j) = max(0.0_wp, &
                                                     eta_ext - state%barotropic%b(i_first - k, j))
            state%barotropic%hu(i_first - k, j) = state%barotropic%h(i_first - k, j)* &
                                                  (hu_ext/max(h_ext, 1.0e-6_wp) - &
                                                   speed*(eta_int - eta_ext)/max(h_int, 1.0e-6_wp))
            state%barotropic%hv(i_first - k, j) = hv_ext
         end do
      end do

   end subroutine apply_nested_west

   subroutine apply_nested_east(nest, state, alpha)
      !! Fill east ghost cells using Flather with parent values
      type(nesting_t), intent(in) :: nest
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: alpha

      integer :: j, k, ng, ny, nx, i_last
      integer :: iy0, pnx, pny
      real(wp) :: wy_j, h_ext, hu_ext, hv_ext, eta_ext
      real(wp) :: h_int, eta_int, speed

      ng = state%grid%nghost
      nx = state%grid%nx_total
      ny = state%grid%ny_total
      i_last = nx - ng
      pnx = nest%p_nx
      pny = nest%p_ny

      do j = 1, ny
         iy0 = nest%we_iy0(j)
         wy_j = nest%we_wy(j)

         h_ext = interp_parent(nest%h_lo, nest%h_hi, &
                               nest%ix_east, iy0, nest%wx_east, wy_j, &
                               alpha, pnx, pny)
         hu_ext = interp_parent(nest%hu_lo, nest%hu_hi, &
                                nest%ix_east, iy0, nest%wx_east, wy_j, &
                                alpha, pnx, pny)
         hv_ext = interp_parent(nest%hv_lo, nest%hv_hi, &
                                nest%ix_east, iy0, nest%wx_east, wy_j, &
                                alpha, pnx, pny)

         eta_ext = h_ext + state%barotropic%b(i_last, j)

         h_int = state%barotropic%h(i_last, j)
         eta_int = h_int + state%barotropic%b(i_last, j)
         speed = sqrt(GRAVITY*max(h_int, 1.0e-6_wp))

         do k = 1, ng
            state%barotropic%h(i_last + k, j) = max(0.0_wp, &
                                                    eta_ext - state%barotropic%b(i_last + k, j))
            state%barotropic%hu(i_last + k, j) = state%barotropic%h(i_last + k, j)* &
                                                 (hu_ext/max(h_ext, 1.0e-6_wp) + &
                                                  speed*(eta_int - eta_ext)/max(h_int, 1.0e-6_wp))
            state%barotropic%hv(i_last + k, j) = hv_ext
         end do
      end do

   end subroutine apply_nested_east

   subroutine apply_nested_south(nest, state, alpha)
      !! Fill south ghost cells using Flather with parent values
      type(nesting_t), intent(in) :: nest
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: alpha

      integer :: i, k, ng, nx, j_first
      integer :: ix0, pnx, pny
      real(wp) :: wx_i, h_ext, hu_ext, hv_ext, eta_ext
      real(wp) :: h_int, eta_int, speed

      ng = state%grid%nghost
      nx = state%grid%nx_total
      j_first = ng + 1
      pnx = nest%p_nx
      pny = nest%p_ny

      do i = 1, nx
         ix0 = nest%sn_ix0(i)
         wx_i = nest%sn_wx(i)

         h_ext = interp_parent(nest%h_lo, nest%h_hi, &
                               ix0, nest%iy_south, wx_i, nest%wy_south, &
                               alpha, pnx, pny)
         hu_ext = interp_parent(nest%hu_lo, nest%hu_hi, &
                                ix0, nest%iy_south, wx_i, nest%wy_south, &
                                alpha, pnx, pny)
         hv_ext = interp_parent(nest%hv_lo, nest%hv_hi, &
                                ix0, nest%iy_south, wx_i, nest%wy_south, &
                                alpha, pnx, pny)

         eta_ext = h_ext + state%barotropic%b(i, j_first)

         h_int = state%barotropic%h(i, j_first)
         eta_int = h_int + state%barotropic%b(i, j_first)
         speed = sqrt(GRAVITY*max(h_int, 1.0e-6_wp))

         do k = 1, ng
            state%barotropic%h(i, j_first - k) = max(0.0_wp, &
                                                     eta_ext - state%barotropic%b(i, j_first - k))
            state%barotropic%hu(i, j_first - k) = hu_ext
            state%barotropic%hv(i, j_first - k) = state%barotropic%h(i, j_first - k)* &
                                                  (hv_ext/max(h_ext, 1.0e-6_wp) - &
                                                   speed*(eta_int - eta_ext)/max(h_int, 1.0e-6_wp))
         end do
      end do

   end subroutine apply_nested_south

   subroutine apply_nested_north(nest, state, alpha)
      !! Fill north ghost cells using Flather with parent values
      type(nesting_t), intent(in) :: nest
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: alpha

      integer :: i, k, ng, nx, ny, j_last
      integer :: ix0, pnx, pny
      real(wp) :: wx_i, h_ext, hu_ext, hv_ext, eta_ext
      real(wp) :: h_int, eta_int, speed

      ng = state%grid%nghost
      nx = state%grid%nx_total
      ny = state%grid%ny_total
      j_last = ny - ng
      pnx = nest%p_nx
      pny = nest%p_ny

      do i = 1, nx
         ix0 = nest%sn_ix0(i)
         wx_i = nest%sn_wx(i)

         h_ext = interp_parent(nest%h_lo, nest%h_hi, &
                               ix0, nest%iy_north, wx_i, nest%wy_north, &
                               alpha, pnx, pny)
         hu_ext = interp_parent(nest%hu_lo, nest%hu_hi, &
                                ix0, nest%iy_north, wx_i, nest%wy_north, &
                                alpha, pnx, pny)
         hv_ext = interp_parent(nest%hv_lo, nest%hv_hi, &
                                ix0, nest%iy_north, wx_i, nest%wy_north, &
                                alpha, pnx, pny)

         eta_ext = h_ext + state%barotropic%b(i, j_last)

         h_int = state%barotropic%h(i, j_last)
         eta_int = h_int + state%barotropic%b(i, j_last)
         speed = sqrt(GRAVITY*max(h_int, 1.0e-6_wp))

         do k = 1, ng
            state%barotropic%h(i, j_last + k) = max(0.0_wp, &
                                                    eta_ext - state%barotropic%b(i, j_last + k))
            state%barotropic%hu(i, j_last + k) = hu_ext
            state%barotropic%hv(i, j_last + k) = state%barotropic%h(i, j_last + k)* &
                                                 (hv_ext/max(h_ext, 1.0e-6_wp) + &
                                                  speed*(eta_int - eta_ext)/max(h_int, 1.0e-6_wp))
         end do
      end do

   end subroutine apply_nested_north

end module rki_nesting
