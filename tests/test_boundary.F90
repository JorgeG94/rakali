!! Unit tests for boundary condition application
module test_boundary
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, TWO_PI, GRAVITY
   use rki_boundary, only: apply_boundary_conditions, bc_params_t, &
                           BC_WALL, BC_OPEN, BC_TIDAL, BC_INFLOW, &
                           BC_DISCHARGE, BC_CLAMPED, BC_SPONGE, BC_CHAPMAN
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_boundary_tests

contains

   subroutine collect_boundary_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("wall_negates_normal_velocity", test_wall_normal), &
                  new_unittest("wall_mirrors_tangential_velocity", test_wall_tangential), &
                  new_unittest("wall_mirrors_eta_flat_bottom", test_wall_eta_flat), &
                  new_unittest("wall_mirrors_eta_sloped_bottom", test_wall_eta_sloped), &
                  new_unittest("wall_symmetry", test_wall_symmetry), &
                  new_unittest("wall_no_negative_depth", test_wall_no_negative_h), &
                  new_unittest("open_zero_gradient_h", test_open_zero_gradient), &
                  new_unittest("tidal_prescribes_eta", test_tidal_eta), &
                  new_unittest("tidal_zero_at_half_period", test_tidal_zero_at_half_period), &
                  new_unittest("tidal_multi_constituent", test_tidal_multi_constituent), &
                  new_unittest("inflow_prescribes_velocity", test_inflow_velocity), &
                  new_unittest("inflow_extrapolates_depth", test_inflow_depth), &
                  new_unittest("discharge_prescribes_flux", test_discharge_flux), &
                  new_unittest("clamped_prescribes_eta", test_clamped_eta), &
                  new_unittest("sponge_falls_through_to_wall", test_sponge_wall_ghost), &
                  new_unittest("chapman_zero_dt_extrapolates", test_chapman_zero_dt), &
                  new_unittest("wall_normal_high_edges", test_wall_normal_high_edges), &
                  new_unittest("discharge_normal_all_edges", test_discharge_normal_all_edges), &
                  new_unittest("open_radiation_outward_all_edges", test_open_radiation_outward_all_edges), &
                  new_unittest("inflow_normal_all_edges", test_inflow_normal_all_edges), &
                  new_unittest("chapman_dt_positive_two_edges", test_chapman_dt_positive) &
                  ]
   end subroutine collect_boundary_tests

   subroutine apply_bc_wrapper(h, hu, hv, b, nx, ny, ng, bc)
      !! Helper to wrap OpenACC/CUDA data management
      integer, intent(in) :: nx, ny, ng
      real(wp), intent(inout) :: h(nx, ny), hu(nx, ny), hv(nx, ny)
      real(wp), intent(inout) :: b(nx, ny)
         !! INOUT since wall/sponge edges now mirror the bed into the ghosts.
      type(bc_params_t), intent(inout) :: bc

      type(hgrid_t) :: grid

      call grid%init(nx - 2*ng, ny - 2*ng, ng, 1.0_wp, 1.0_wp)

      !$omp target data map(tofrom:  h, hu, hv, b)
      call apply_boundary_conditions(h, hu, hv, b, grid, bc)
      !$omp end target data
   end subroutine apply_bc_wrapper

   subroutine test_wall_normal(error)
      !! Normal velocity component should be negated at wall boundaries
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      hu(NG + 1, 4) = -5.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hu(NG, 4) - 5.0_wp) < 1.0e-14_wp, &
                 "West ghost hu should negate interior hu")
      if (allocated(error)) return

      hu = 0.0_wp
      hv(4, NG + 1) = -3.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hv(4, NG) - 3.0_wp) < 1.0e-14_wp, &
                 "South ghost hv should negate interior hv")
   end subroutine test_wall_normal

   subroutine test_wall_tangential(error)
      !! Tangential velocity component should be mirrored (same sign) at walls
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      hv(NG + 1, 4) = 7.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hv(NG, 4) - 7.0_wp) < 1.0e-14_wp, &
                 "West ghost hv (tangential) should mirror interior hv")
      if (allocated(error)) return

      hu = 0.0_wp
      hv = 0.0_wp
      hu(4, NG + 1) = 4.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hu(4, NG) - 4.0_wp) < 1.0e-14_wp, &
                 "South ghost hu (tangential) should mirror interior hu")
   end subroutine test_wall_tangential

   subroutine test_wall_eta_flat(error)
      !! On flat bottom, ghost h should equal interior h (eta mirror = h mirror)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - 2.0_wp) < 1.0e-14_wp, &
                 "West ghost h should equal interior h on flat bottom")
      if (allocated(error)) return
      call check(error, abs(h(NX - NG + 1, 4) - 2.0_wp) < 1.0e-14_wp, &
                 "East ghost h should equal interior h on flat bottom")
      if (allocated(error)) return
      call check(error, abs(h(4, NG) - 2.0_wp) < 1.0e-14_wp, &
                 "South ghost h should equal interior h on flat bottom")
      if (allocated(error)) return
      call check(error, abs(h(4, NY - NG + 1) - 2.0_wp) < 1.0e-14_wp, &
                 "North ghost h should equal interior h on flat bottom")
   end subroutine test_wall_eta_flat

   subroutine test_wall_eta_sloped(error)
      !! On sloped bottom, ghost h should mirror eta = h + B, not h directly
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: eta_interior, h_ghost_expected
      integer :: i
      type(bc_params_t) :: bc

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp

      do i = 1, NX
         b(i, :) = real(i, wp)*0.1_wp
      end do

      do i = NG + 1, NX - NG
         h(i, :) = 2.0_wp - b(i, 1)
      end do

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      eta_interior = 2.0_wp
      h_ghost_expected = max(0.0_wp, eta_interior - b(NG, 1))

      call check(error, abs(h(NG, 4) - h_ghost_expected) < 1.0e-12_wp, &
                 "Ghost h should mirror eta on sloped bottom")
   end subroutine test_wall_eta_sloped

   subroutine test_wall_symmetry(error)
      !! Uniform quiescent state should remain symmetric after BC application
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 10, NY = 10, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      integer :: k
      type(bc_params_t) :: bc

      h = 1.5_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      do k = 1, NG
         call check(error, abs(h(k, 5) - 1.5_wp) < 1.0e-14_wp, &
                    "West ghost h should match interior for uniform state")
         if (allocated(error)) return
         call check(error, abs(h(NX - k + 1, 5) - 1.5_wp) < 1.0e-14_wp, &
                    "East ghost h should match interior for uniform state")
         if (allocated(error)) return
      end do
   end subroutine test_wall_symmetry

   subroutine test_wall_no_negative_h(error)
      !! Ghost cell depth should never be negative
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      integer :: i
      type(bc_params_t) :: bc

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      do i = 1, NG
         b(i, :) = 10.0_wp
      end do
      h(NG + 1, :) = 0.5_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, h(NG, 4) >= 0.0_wp, &
                 "Ghost cell depth should not be negative")
      if (allocated(error)) return
      call check(error, h(1, 4) >= 0.0_wp, &
                 "Outer ghost cell depth should not be negative")
   end subroutine test_wall_no_negative_h

   subroutine test_open_zero_gradient(error)
      !! Open boundary should copy interior h to ghost cells (Flather)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      h(NG + 1, :) = 3.0_wp

      bc%west = BC_OPEN
      bc%east = BC_OPEN
      bc%south = BC_OPEN
      bc%north = BC_OPEN

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - 3.0_wp) < 1.0e-12_wp, &
                 "Open west ghost h should equal first interior h")
      if (allocated(error)) return
      call check(error, abs(h(1, 4) - 3.0_wp) < 1.0e-12_wp, &
                 "Open west outer ghost h should equal first interior h")
   end subroutine test_open_zero_gradient

   subroutine test_tidal_eta(error)
      !! Tidal boundary should prescribe eta = A*sin(2*pi*t/T) in ghost cells
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc
      real(wp) :: expected_eta, expected_h

      h = 5.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = -3.0_wp

      bc%west = BC_TIDAL
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%tidal_amplitude = 1.0_wp
      bc%tidal_period = 100.0_wp
      bc%time = 25.0_wp

      expected_eta = 1.0_wp*sin(TWO_PI*25.0_wp/100.0_wp)
      expected_h = max(0.0_wp, expected_eta - b(NG, 4))

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - expected_h) < 1.0e-12_wp, &
                 "Tidal west ghost h should match prescribed eta")
   end subroutine test_tidal_eta

   subroutine test_tidal_zero_at_half_period(error)
      !! At half period, tidal elevation should be ~0 (sin(pi) = 0)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc
      real(wp) :: expected_h

      h = 5.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = -3.0_wp

      bc%west = BC_TIDAL
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%tidal_amplitude = 1.0_wp
      bc%tidal_period = 100.0_wp
      bc%time = 50.0_wp

      expected_h = max(0.0_wp, 0.0_wp - b(NG, 4))

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - expected_h) < 1.0e-10_wp, &
                 "Tidal ghost h at half period should be ~reference depth")
   end subroutine test_tidal_zero_at_half_period

   subroutine test_tidal_multi_constituent(error)
      !! Multi-constituent tidal BC should prescribe eta = sum A_n cos(w_n t - phi_n)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc
      real(wp) :: expected_eta, expected_h, t_sim

      h = 5.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = -3.0_wp

      bc%west = BC_TIDAL
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL

      ! Two constituents: M2-like and S2-like
      bc%n_tidal_constituents = 2
      bc%tidal_amp(1) = 0.5_wp
      bc%tidal_omega(1) = 1.405189e-4_wp
      bc%tidal_phase(1) = -0.9_wp
      bc%tidal_amp(2) = 0.12_wp
      bc%tidal_omega(2) = 1.454441e-4_wp
      bc%tidal_phase(2) = -0.8_wp

      t_sim = 20000.0_wp
      bc%time = t_sim

      expected_eta = 0.5_wp*cos(1.405189e-4_wp*t_sim - (-0.9_wp)) &
                     + 0.12_wp*cos(1.454441e-4_wp*t_sim - (-0.8_wp))
      expected_h = max(0.0_wp, expected_eta - b(NG, 4))

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - expected_h) < 1.0e-10_wp, &
                 "Multi-constituent tidal ghost h should match superposition")
   end subroutine test_tidal_multi_constituent

   subroutine test_inflow_velocity(error)
      !! Inflow BC should prescribe hu = h * inflow_velocity in ghost cells
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc
      real(wp) :: expected_hu

      h = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      bc%west = BC_INFLOW
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%inflow_velocity = 3.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      ! Ghost h should be extrapolated from interior (= 2.0)
      ! Ghost hu should be h * inflow_velocity = 2.0 * 3.0 = 6.0
      expected_hu = h(NG, 4)*3.0_wp
      call check(error, abs(hu(NG, 4) - expected_hu) < 1.0e-12_wp, &
                 "Inflow west ghost hu should be h * inflow_velocity")
      if (allocated(error)) return

      ! East inflow: sign flips (negative hu = into domain from east)
      bc%west = BC_WALL
      bc%east = BC_INFLOW
      h = 2.0_wp
      hu = 0.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      expected_hu = -h(NX - NG + 1, 4)*3.0_wp
      call check(error, abs(hu(NX - NG + 1, 4) - expected_hu) < 1.0e-12_wp, &
                 "Inflow east ghost hu should be -h * inflow_velocity")
   end subroutine test_inflow_velocity

   subroutine test_inflow_depth(error)
      !! Inflow BC should extrapolate h from interior
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      ! Set a non-uniform interior depth
      h(NG + 1, :) = 4.0_wp

      bc%west = BC_INFLOW
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%inflow_velocity = 1.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      ! Ghost h should equal first interior h (extrapolation)
      call check(error, abs(h(NG, 4) - 4.0_wp) < 1.0e-12_wp, &
                 "Inflow ghost h should be extrapolated from interior")
      if (allocated(error)) return
      call check(error, abs(h(1, 4) - 4.0_wp) < 1.0e-12_wp, &
                 "Inflow outer ghost h should be extrapolated from interior")
   end subroutine test_inflow_depth

   subroutine test_discharge_flux(error)
      !! Discharge BC should prescribe hu directly as the flux value
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      bc%west = BC_DISCHARGE
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%discharge_flux = 5.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      ! Ghost hu should be exactly the discharge flux value
      call check(error, abs(hu(NG, 4) - 5.0_wp) < 1.0e-12_wp, &
                 "Discharge west ghost hu should be discharge_flux")
      if (allocated(error)) return

      ! East discharge: sign flips
      bc%west = BC_WALL
      bc%east = BC_DISCHARGE
      h = 2.0_wp
      hu = 0.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hu(NX - NG + 1, 4) - (-5.0_wp)) < 1.0e-12_wp, &
                 "Discharge east ghost hu should be -discharge_flux")
   end subroutine test_discharge_flux

   subroutine test_clamped_eta(error)
      !! Clamped BC should prescribe h = max(0, eta_clamped - b)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc
      real(wp) :: expected_h

      h = 5.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = -3.0_wp

      bc%west = BC_CLAMPED
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%clamped_elevation = 1.5_wp

      expected_h = max(0.0_wp, 1.5_wp - b(NG, 4))  ! = 1.5 - (-3.0) = 4.5

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - expected_h) < 1.0e-12_wp, &
                 "Clamped west ghost h should match prescribed elevation")
      if (allocated(error)) return

      ! Velocity should be extrapolated (zero-gradient)
      call check(error, abs(hu(NG, 4) - hu(NG + 1, 4)) < 1.0e-12_wp, &
                 "Clamped west ghost hu should be zero-gradient from interior")
   end subroutine test_clamped_eta

   subroutine test_sponge_wall_ghost(error)
      !! Sponge BC ghost cells should behave like wall BC
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      hu(NG + 1, 4) = -5.0_wp

      bc%west = BC_SPONGE
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%sponge_width = 5
      bc%sponge_strength = 1.0_wp
      bc%sponge_h_ref = 1.0_wp

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      ! Ghost cells should mirror like wall: hu negated
      call check(error, abs(hu(NG, 4) - 5.0_wp) < 1.0e-14_wp, &
                 "Sponge west ghost hu should negate like wall BC")
   end subroutine test_sponge_wall_ghost

   subroutine test_chapman_zero_dt(error)
      !! Chapman BC with dt=0 should extrapolate eta from interior
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      bc%west = BC_CHAPMAN
      bc%east = BC_WALL
      bc%south = BC_WALL
      bc%north = BC_WALL
      bc%dt = 0.0_wp
      bc%dx = 1.0_wp
      bc%dy = 1.0_wp

      ! Allocate eta_old arrays
      allocate (bc%eta_old_west(NY), source=0.0_wp)

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      ! With dt=0, eta_ghost_new = eta_int = h_int + b_int = 2.0
      ! h_ghost = max(0, 2.0 - 0.0) = 2.0
      call check(error, abs(h(NG, 4) - 2.0_wp) < 1.0e-12_wp, &
                 "Chapman with dt=0 should extrapolate interior eta to ghost")

      deallocate (bc%eta_old_west)
   end subroutine test_chapman_zero_dt

   ! -----------------------------------------------------------------
   ! High-edge + y-edge coverage (added to guard the per-edge
   ! decomposition refactor of apply_boundary_conditions).  The pre-
   ! existing tests assert almost exclusively on the WEST edge, leaving
   ! the east/south/north index arithmetic, the hu<->hv normal/tangential
   ! swap (normal = hu on W/E, hv on S/N), and the outward-radiation
   ! sign untested.  These close those gaps with exact values derived
   ! from the kernel formulas.
   ! -----------------------------------------------------------------

   subroutine test_wall_normal_high_edges(error)
      !! Wall reflection of the NORMAL momentum on the EAST (hu) and
      !! NORTH (hv) edges — the high-index edges the existing wall test
      !! never asserts.  Covers the i_last / j_last mirror arithmetic.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 8, NG = 2
      integer, parameter :: I_LAST = NX - NG, J_LAST = NY - NG
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 1.0_wp; hu = 0.0_wp; hv = 0.0_wp; b = 0.0_wp
      hu(I_LAST, 4) = 5.0_wp     ! interior normal velocity at the east face
      hv(4, J_LAST) = 4.0_wp     ! interior normal velocity at the north face

      bc%west = BC_WALL; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_WALL

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hu(I_LAST + 1, 4) - (-5.0_wp)) < 1.0e-14_wp, &
                 "East wall ghost hu should negate interior hu")
      if (allocated(error)) return
      call check(error, abs(hv(4, J_LAST + 1) - (-4.0_wp)) < 1.0e-14_wp, &
                 "North wall ghost hv should negate interior hv")
   end subroutine test_wall_normal_high_edges

   subroutine test_discharge_normal_all_edges(error)
      !! BC_DISCHARGE prescribes the NORMAL momentum (hu on W/E, hv on
      !! S/N) to +flux on the low edges (W/S) and -flux on the high
      !! edges (E/N).  This is the most direct check of the hu<->hv
      !! normal-field selection on every edge — the part a per-edge
      !! refactor is most likely to get wrong.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 8, NG = 2
      integer, parameter :: I_LAST = NX - NG, J_LAST = NY - NG
      real(wp), parameter :: DF = 3.0_wp
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = 2.0_wp; hu = 0.0_wp; hv = 0.0_wp; b = 0.0_wp
      bc%west = BC_DISCHARGE; bc%east = BC_DISCHARGE
      bc%south = BC_DISCHARGE; bc%north = BC_DISCHARGE
      bc%discharge_flux = DF

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(hu(NG, 4) - DF) < 1.0e-14_wp, &
                 "West discharge: hu ghost = +flux")
      if (allocated(error)) return
      call check(error, abs(hu(I_LAST + 1, 4) - (-DF)) < 1.0e-14_wp, &
                 "East discharge: hu ghost = -flux")
      if (allocated(error)) return
      call check(error, abs(hv(4, NG) - DF) < 1.0e-14_wp, &
                 "South discharge: hv ghost = +flux (normal momentum is hv)")
      if (allocated(error)) return
      call check(error, abs(hv(4, J_LAST + 1) - (-DF)) < 1.0e-14_wp, &
                 "North discharge: hv ghost = -flux (normal momentum is hv)")
   end subroutine test_discharge_normal_all_edges

   subroutine test_open_radiation_outward_all_edges(error)
      !! Flather (BC_OPEN) at rest (u=v=0, flat bottom, uniform depth H)
      !! radiates the surface outward: the ghost normal momentum is
      !! -/+ sqrt(g*H)*H, signed to point OUT of the domain on every
      !! edge (W:-x, E:+x, S:-y, N:+y).  Exact value from the kernel;
      !! covers per-edge radiation sign + normal-field selection.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 8, NG = 2
      integer, parameter :: I_LAST = NX - NG, J_LAST = NY - NG
      real(wp), parameter :: H0 = 2.0_wp
      real(wp) :: rad
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = H0; hu = 0.0_wp; hv = 0.0_wp; b = 0.0_wp
      bc%west = BC_OPEN; bc%east = BC_OPEN
      bc%south = BC_OPEN; bc%north = BC_OPEN

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      rad = sqrt(GRAVITY*H0)*H0     ! magnitude of the radiated normal momentum
      call check(error, abs(hu(NG, 4) - (-rad)) < 1.0e-10_wp, &
                 "West open: hu ghost radiates outward (-x) = -sqrt(gH)H")
      if (allocated(error)) return
      call check(error, abs(hu(I_LAST + 1, 4) - rad) < 1.0e-10_wp, &
                 "East open: hu ghost radiates outward (+x) = +sqrt(gH)H")
      if (allocated(error)) return
      call check(error, abs(hv(4, NG) - (-rad)) < 1.0e-10_wp, &
                 "South open: hv ghost radiates outward (-y), normal = hv")
      if (allocated(error)) return
      call check(error, abs(hv(4, J_LAST + 1) - rad) < 1.0e-10_wp, &
                 "North open: hv ghost radiates outward (+y), normal = hv")
      if (allocated(error)) return
      call check(error, abs(h(NG, 4) - H0) < 1.0e-12_wp, &
                 "Open ghost h equals interior depth")
   end subroutine test_open_radiation_outward_all_edges

   subroutine test_inflow_normal_all_edges(error)
      !! BC_INFLOW sets the NORMAL momentum to +/- h_ghost*vel with
      !! h_ghost extrapolated from the interior.  Uniform depth H, so
      !! the normal momentum is +/- H*vel; +on W/S, -on E/N.  Covers
      !! the normal-field selection + sign on all four edges.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 8, NG = 2
      integer, parameter :: I_LAST = NX - NG, J_LAST = NY - NG
      real(wp), parameter :: H0 = 2.0_wp, V = 0.5_wp
      real(wp) :: mag
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = H0; hu = 0.0_wp; hv = 0.0_wp; b = 0.0_wp
      bc%west = BC_INFLOW; bc%east = BC_INFLOW
      bc%south = BC_INFLOW; bc%north = BC_INFLOW
      bc%inflow_velocity = V

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      mag = H0*V
      call check(error, abs(hu(NG, 4) - mag) < 1.0e-12_wp, &
                 "West inflow: hu ghost = +h*vel")
      if (allocated(error)) return
      call check(error, abs(hu(I_LAST + 1, 4) - (-mag)) < 1.0e-12_wp, &
                 "East inflow: hu ghost = -h*vel")
      if (allocated(error)) return
      call check(error, abs(hv(4, NG) - mag) < 1.0e-12_wp, &
                 "South inflow: hv ghost = +h*vel (normal momentum is hv)")
      if (allocated(error)) return
      call check(error, abs(hv(4, J_LAST + 1) - (-mag)) < 1.0e-12_wp, &
                 "North inflow: hv ghost = -h*vel (normal momentum is hv)")
   end subroutine test_inflow_normal_all_edges

   subroutine test_chapman_dt_positive(error)
      !! Chapman radiation with dt>0 on WEST and NORTH — the radiation
      !! branch the existing dt=0 test never enters.  Flat bottom +
      !! uniform interior depth makes the along-normal eta gradient
      !! zero, so the phase speed clamps to zero and the ghost surface
      !! holds the persistent eta_old.  With b=0 that means
      !! h_ghost = eta_old.  Also confirms each edge reads its OWN
      !! eta_old array (west <- eta_old_west, north <- eta_old_north);
      !! a swap there would surface the wrong seed.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NX = 8, NY = 8, NG = 2
      integer, parameter :: J_LAST = NY - NG
      real(wp), parameter :: H0 = 2.0_wp, E0_W = 1.5_wp, E0_N = 0.75_wp
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      type(bc_params_t) :: bc

      h = H0; hu = 0.0_wp; hv = 0.0_wp; b = 0.0_wp
      bc%west = BC_CHAPMAN; bc%east = BC_WALL
      bc%south = BC_WALL; bc%north = BC_CHAPMAN
      bc%dt = 1.0_wp; bc%dx = 1.0_wp; bc%dy = 1.0_wp

      allocate (bc%eta_old_west(NY), source=E0_W)
      allocate (bc%eta_old_north(NX), source=E0_N)

      call apply_bc_wrapper(h, hu, hv, b, NX, NY, NG, bc)

      call check(error, abs(h(NG, 4) - E0_W) < 1.0e-12_wp, &
                 "West Chapman (dt>0): zero gradient => ghost h holds eta_old_west")
      if (allocated(error)) return
      call check(error, abs(h(4, J_LAST + 1) - E0_N) < 1.0e-12_wp, &
                 "North Chapman (dt>0): zero gradient => ghost h holds eta_old_north")

      deallocate (bc%eta_old_west)
      deallocate (bc%eta_old_north)
   end subroutine test_chapman_dt_positive

end module test_boundary
