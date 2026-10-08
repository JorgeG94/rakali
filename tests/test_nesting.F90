!! Unit tests for the one-way nesting module
module test_nesting
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_boundary, only: BC_NESTED, BC_WALL, bc_type_from_string
   use rki_nesting, only: nesting_t, nesting_init, nesting_update, nesting_destroy
   use netcdf, only: nf90_create, nf90_def_dim, nf90_def_var, nf90_enddef, &
                     nf90_put_var, nf90_close, nf90_clobber, nf90_netcdf4, nf90_double
   implicit none
   private

   public :: collect_nesting_tests

contains

   subroutine collect_nesting_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_destroy_lifecycle", test_init_destroy), &
                  new_unittest("no_file_noop", test_no_file), &
                  new_unittest("bc_nested_string", test_bc_string), &
                  new_unittest("uniform_parent_west", test_uniform_west), &
                  new_unittest("uniform_parent_all_edges", test_uniform_all), &
                  new_unittest("time_interpolation", test_time_interp) &
                  ]
   end subroutine collect_nesting_tests

   ! ================================================================
   ! Helper: create a parent domain NetCDF file
   ! ================================================================

   subroutine create_parent_file(filename, p_nx, p_ny, p_nt, &
                                 x_vals, y_vals, t_vals, &
                                 h_vals, hu_vals, hv_vals)
      !! Create a test parent solution NetCDF file
      character(len=*), intent(in) :: filename
      integer, intent(in) :: p_nx, p_ny, p_nt
      real(wp), intent(in) :: x_vals(p_nx), y_vals(p_ny), t_vals(p_nt)
      real(wp), intent(in) :: h_vals(p_nx, p_ny, p_nt)
      real(wp), intent(in) :: hu_vals(p_nx, p_ny, p_nt)
      real(wp), intent(in) :: hv_vals(p_nx, p_ny, p_nt)

      integer :: ncid, dim_x, dim_y, dim_t
      integer :: var_x, var_y, var_t, var_h, var_hu, var_hv
      integer :: status, dimids(3)

      status = nf90_create(trim(filename), ior(nf90_clobber, nf90_netcdf4), ncid)

      status = nf90_def_dim(ncid, "x", p_nx, dim_x)
      status = nf90_def_dim(ncid, "y", p_ny, dim_y)
      status = nf90_def_dim(ncid, "time", p_nt, dim_t)

      status = nf90_def_var(ncid, "x", nf90_double, [dim_x], var_x)
      status = nf90_def_var(ncid, "y", nf90_double, [dim_y], var_y)
      status = nf90_def_var(ncid, "time", nf90_double, [dim_t], var_t)

      dimids = [dim_x, dim_y, dim_t]
      status = nf90_def_var(ncid, "h", nf90_double, dimids, var_h)
      status = nf90_def_var(ncid, "hu", nf90_double, dimids, var_hu)
      status = nf90_def_var(ncid, "hv", nf90_double, dimids, var_hv)

      status = nf90_enddef(ncid)

      status = nf90_put_var(ncid, var_x, x_vals)
      status = nf90_put_var(ncid, var_y, y_vals)
      status = nf90_put_var(ncid, var_t, t_vals)
      status = nf90_put_var(ncid, var_h, h_vals)
      status = nf90_put_var(ncid, var_hu, hu_vals)
      status = nf90_put_var(ncid, var_hv, hv_vals)

      status = nf90_close(ncid)

   end subroutine create_parent_file

   ! ================================================================
   ! Tests
   ! ================================================================

   subroutine test_init_destroy(error)
      !! Nesting init/destroy should not crash with valid file
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(nesting_t) :: nest
      character(len=*), parameter :: FNAME = "test_nest_lifecycle.nc"
      integer, parameter :: PNX = 3, PNY = 3, PNT = 2
      real(wp) :: px(PNX), py(PNY), pt(PNT)
      real(wp) :: h(PNX, PNY, PNT), hu(PNX, PNY, PNT), hv(PNX, PNY, PNT)

      px = [0.0_wp, 50.0_wp, 100.0_wp]
      py = [0.0_wp, 50.0_wp, 100.0_wp]
      pt = [0.0_wp, 3600.0_wp]
      h = 2.0_wp; hu = 0.0_wp; hv = 0.0_wp

      call create_parent_file(FNAME, PNX, PNY, PNT, px, py, pt, h, hu, hv)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call nesting_init(nest, FNAME, state)

      call check(error, nest%active, "Nesting should be active after init")
      if (allocated(error)) return

      call check(error, nest%p_nx == PNX, "p_nx should match file")
      if (allocated(error)) return

      call nesting_destroy(nest)
      call check(error,.not. nest%active, "Nesting should be inactive after destroy")

      call state%destroy()

   end subroutine test_init_destroy

   subroutine test_no_file(error)
      !! Empty nesting file path should leave nesting inactive
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(nesting_t) :: nest

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      call nesting_init(nest, "", state)

      call check(error,.not. nest%active, &
                 "Nesting should be inactive with empty filename")

      call nesting_destroy(nest)
      call state%destroy()

   end subroutine test_no_file

   subroutine test_bc_string(error)
      !! "nested" string should map to BC_NESTED
      type(error_type), allocatable, intent(out) :: error

      call check(error, bc_type_from_string("nested") == BC_NESTED, &
                 "'nested' should map to BC_NESTED")

   end subroutine test_bc_string

   subroutine test_uniform_west(error)
      !! Uniform parent h=5 on west nested edge should fill ghost cells
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(nesting_t) :: nest
      character(len=*), parameter :: FNAME = "test_nest_west.nc"
      integer, parameter :: PNX = 3, PNY = 3, PNT = 1
      real(wp) :: px(PNX), py(PNY), pt(PNT)
      real(wp) :: h_p(PNX, PNY, PNT), hu_p(PNX, PNY, PNT), hv_p(PNX, PNY, PNT)
      real(wp), parameter :: H_EXT = 5.0_wp
      integer :: ng, j_mid

      ! Parent grid covers child domain
      px = [0.0_wp, 50.0_wp, 100.0_wp]
      py = [0.0_wp, 50.0_wp, 100.0_wp]
      pt = [0.0_wp]
      h_p = H_EXT; hu_p = 0.0_wp; hv_p = 0.0_wp

      call create_parent_file(FNAME, PNX, PNY, PNT, px, py, pt, &
                              h_p, hu_p, hv_p)

      ! Child grid: 10x10, dx=dy=10m, domain 0-100m
      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      state%bc%west = BC_NESTED

      ! Set child interior to a different depth
      state%barotropic%h = 2.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp
      state%barotropic%b = 0.0_wp

      call nesting_init(nest, FNAME, state)
      call nesting_update(nest, state, 0.0_wp)

      ng = state%grid%nghost
      j_mid = ng + state%grid%ny_phys/2

      ! Ghost cells on west edge should have h = H_EXT (flat bottom, b=0)
      ! With Flather and zero velocity, interior h=2:
      ! eta_ext = H_EXT + b = 5, eta_int = h_int + b = 2
      ! h_ghost = max(0, eta_ext - b_ghost) = 5
      call check(error, abs(state%barotropic%h(ng, j_mid) - H_EXT) < 1.0e-6_wp, &
                 "West ghost cell h should reflect parent h")
      if (allocated(error)) return

      ! Flather momentum check — NOT just "positive depth". Earlier
      ! versions of this test only asserted h_ghost > 0, which is
      ! trivially true for any nonzero parent h. The whole point of
      ! Flather is the momentum radiation correction:
      !
      !   hu_ghost = h_ghost * (u_ext - c * (eta_int - eta_ext) / h_int)
      !
      ! With parent u_ext = 0, h_int = 2, eta_int = 2, eta_ext = 5,
      ! c = sqrt(g * h_int) ~ 4.43, h_ghost = 5, the formula gives
      ! hu_ghost ~ 5 * (0 - 4.43 * (2 - 5) / 2) = +33.2 m^2/s. A bug
      ! that zeroed the momentum branch (or flipped the sign of the
      ! radiation correction) would fail here.
      block
         real(wp) :: h_int, h_ghost, eta_int, eta_ext
         real(wp) :: speed, hu_expected
         real(wp), parameter :: H_INT_VAL = 2.0_wp
         h_int = H_INT_VAL
         h_ghost = state%barotropic%h(ng, j_mid)
         eta_int = h_int + state%barotropic%b(ng + 1, j_mid)
         eta_ext = H_EXT + state%barotropic%b(ng + 1, j_mid)
         speed = sqrt(GRAVITY*max(h_int, 1.0e-6_wp))
         hu_expected = h_ghost*(0.0_wp &
                                - speed*(eta_int - eta_ext) &
                                /max(h_int, 1.0e-6_wp))
         call check(error, abs(state%barotropic%hu(ng, j_mid) - hu_expected) < 1.0e-6_wp, &
                    "West Flather hu_ghost must match radiation formula")
         if (allocated(error)) return
         call check(error, state%barotropic%hu(ng, j_mid) > 0.0_wp, &
                    "West Flather: higher parent eta drives inflow (+hu)")
      end block

      call nesting_destroy(nest)
      call state%destroy()

   end subroutine test_uniform_west

   subroutine test_uniform_all(error)
      !! All four edges nested with uniform parent: ghost cells filled
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(nesting_t) :: nest
      character(len=*), parameter :: FNAME = "test_nest_all.nc"
      integer, parameter :: PNX = 3, PNY = 3, PNT = 1
      real(wp) :: px(PNX), py(PNY), pt(PNT)
      real(wp) :: h_p(PNX, PNY, PNT), hu_p(PNX, PNY, PNT), hv_p(PNX, PNY, PNT)
      real(wp), parameter :: H_EXT = 3.0_wp
      integer :: ng, nx, ny, i_mid, j_mid

      px = [0.0_wp, 50.0_wp, 100.0_wp]
      py = [0.0_wp, 50.0_wp, 100.0_wp]
      pt = [0.0_wp]
      h_p = H_EXT; hu_p = 0.0_wp; hv_p = 0.0_wp

      call create_parent_file(FNAME, PNX, PNY, PNT, px, py, pt, &
                              h_p, hu_p, hv_p)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      state%bc%west = BC_NESTED
      state%bc%east = BC_NESTED
      state%bc%south = BC_NESTED
      state%bc%north = BC_NESTED
      state%barotropic%h = H_EXT  ! Match parent so Flather gives zero velocity
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp
      state%barotropic%b = 0.0_wp

      call nesting_init(nest, FNAME, state)
      call nesting_update(nest, state, 0.0_wp)

      ng = state%grid%nghost
      nx = state%grid%nx_total
      ny = state%grid%ny_total
      i_mid = ng + state%grid%nx_phys/2
      j_mid = ng + state%grid%ny_phys/2

      ! West ghost
      call check(error, abs(state%barotropic%h(ng, j_mid) - H_EXT) < 1.0e-6_wp, &
                 "West ghost h should match parent")
      if (allocated(error)) return

      ! East ghost
      call check(error, abs(state%barotropic%h(nx - ng + 1, j_mid) - H_EXT) < 1.0e-6_wp, &
                 "East ghost h should match parent")
      if (allocated(error)) return

      ! South ghost
      call check(error, abs(state%barotropic%h(i_mid, ng) - H_EXT) < 1.0e-6_wp, &
                 "South ghost h should match parent")
      if (allocated(error)) return

      ! North ghost
      call check(error, abs(state%barotropic%h(i_mid, ny - ng + 1) - H_EXT) < 1.0e-6_wp, &
                 "North ghost h should match parent")

      call nesting_destroy(nest)
      call state%destroy()

   end subroutine test_uniform_all

   subroutine test_time_interp(error)
      !! Time interpolation between two parent snapshots
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      type(nesting_t) :: nest
      character(len=*), parameter :: FNAME = "test_nest_time.nc"
      integer, parameter :: PNX = 3, PNY = 3, PNT = 2
      real(wp) :: px(PNX), py(PNY), pt(PNT)
      real(wp) :: h_p(PNX, PNY, PNT), hu_p(PNX, PNY, PNT), hv_p(PNX, PNY, PNT)
      integer :: ng, j_mid
      real(wp) :: expected_h

      px = [0.0_wp, 50.0_wp, 100.0_wp]
      py = [0.0_wp, 50.0_wp, 100.0_wp]
      pt = [0.0_wp, 100.0_wp]
      ! h = 2.0 at t=0, h = 4.0 at t=100
      h_p(:, :, 1) = 2.0_wp
      h_p(:, :, 2) = 4.0_wp
      hu_p = 0.0_wp; hv_p = 0.0_wp

      call create_parent_file(FNAME, PNX, PNY, PNT, px, py, pt, &
                              h_p, hu_p, hv_p)

      call state%init(10, 10, 2, 10.0_wp, 10.0_wp, 0.45_wp)
      state%bc%west = BC_NESTED
      ! Interior h deliberately distinct from any expected ghost value
      ! at the test times, so a nesting_update that silently did
      ! nothing (or hit a wrong code path) would leave the ghost at
      ! 1.5 and be caught by every check below.
      state%barotropic%h = 1.5_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp
      state%barotropic%b = 0.0_wp

      call nesting_init(nest, FNAME, state)

      ! At t=50 (halfway), expect h=3.0 from parent (= (2.0 + 4.0)/2)
      call nesting_update(nest, state, 50.0_wp)

      ng = state%grid%nghost
      j_mid = ng + state%grid%ny_phys/2
      expected_h = 3.0_wp

      ! Ghost h should be the interpolated parent value. With the
      ! interior initialised to 1.5 a missing ghost write would leave
      ! 1.5 here (fail by 1.5), so this is a real interpolation test,
      ! not a "IC matches target by coincidence" test.
      call check(error, abs(state%barotropic%h(ng, j_mid) - expected_h) < 1.0e-6_wp, &
                 "Time-interpolated parent h should be 3.0 at t=50s")
      if (allocated(error)) return

      ! At t=25 (quarter), expect h=2.5 (linear between 2.0 and 4.0)
      call nesting_update(nest, state, 25.0_wp)
      expected_h = 2.5_wp
      call check(error, abs(state%barotropic%h(ng, j_mid) - expected_h) < 1.0e-6_wp, &
                 "Time-interpolated parent h should be 2.5 at t=25s")
      if (allocated(error)) return

      ! At t=75, expect h=3.5
      call nesting_update(nest, state, 75.0_wp)
      expected_h = 3.5_wp
      call check(error, abs(state%barotropic%h(ng, j_mid) - expected_h) < 1.0e-6_wp, &
                 "Time-interpolated parent h should be 3.5 at t=75s")

      call nesting_destroy(nest)
      call state%destroy()

   end subroutine test_time_interp

end module test_nesting
