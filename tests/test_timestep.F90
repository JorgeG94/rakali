!! Unit tests for the CFL timestep computation
module test_timestep
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, DESING_EPS
   use rki_kernel_timestep, only: compute_timestep
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_timestep_tests

   integer, parameter :: NX = 14, NY = 14, NG = 2
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp
   type(hgrid_t), save :: grid
   logical, save :: grid_ready = .false.

contains

   subroutine ensure_grid()
      if (grid_ready) return
      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
      grid_ready = .true.
   end subroutine ensure_grid

   subroutine collect_timestep_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      call ensure_grid()

      testsuite = [ &
                  new_unittest("uniform_depth", test_uniform_depth), &
                  new_unittest("dry_cells_ignored", test_dry_cells_ignored), &
                  new_unittest("thin_layer_bounded", test_thin_layer_bounded), &
                  new_unittest("high_velocity_cell", test_high_velocity_cell), &
                  new_unittest("all_dry", test_all_dry) &
                  ]
   end subroutine collect_timestep_tests

   subroutine test_uniform_depth(error)
      !! Uniform h=1, u=v=0: dt should be cfl*dx/sqrt(g*h)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: dt, expected_dt
      real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      expected_dt = CFL*DX/sqrt(GRAVITY*1.0_wp)
      call check(error, abs(dt - expected_dt) < 1.0e-12_wp, &
                 "dt should match CFL*dx/sqrt(g*h)")
   end subroutine test_uniform_depth

   subroutine test_dry_cells_ignored(error)
      !! Cells below DRY_TOLERANCE should not affect dt
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: dt1, dt2
      real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp

      ! Case 1: all wet at h=1
      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt1)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! Case 2: same but add dry cells with huge momentum
      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      h(5, 5) = 0.5_wp*DRY_TOLERANCE
      hu(5, 5) = 1000.0_wp
      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt2)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, abs(dt1 - dt2) < 1.0e-12_wp, &
                 "Dry cells with huge momentum should not change dt")
   end subroutine test_dry_cells_ignored

   subroutine test_thin_layer_bounded(error)
      !! Thin-layer cells should have bounded velocity via desingularization
      !! and not dominate the CFL computation
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: dt
      real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp
      real(wp) :: dt_quiescent

      ! Quiescent baseline at h=1
      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt_quiescent)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! Sub-DESING_EPS cell with momentum that would give u = hu/h =
      ! 10 m/s.  The estimator's desingularized velocity rolls off below
      ! eps (~3.4 m/s here), so the thin cell must not dominate the CFL.
      ! (Cells ABOVE eps with real velocity now DO constrain dt — that is
      ! the fix, not a regression: the update kernels no longer zero
      ! their momentum, so the estimator must count them.)
      h(5, 5) = 0.5_wp*DESING_EPS
      hu(5, 5) = 0.5_wp*DESING_EPS*10.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! With desingularization, the thin cell should not dominate over
      ! the deep cells — dt should not drop catastrophically
      call check(error, dt > 0.5_wp*dt_quiescent, &
                 "Desingularized thin-layer cell should not halve dt")
   end subroutine test_thin_layer_bounded

   subroutine test_high_velocity_cell(error)
      !! A single fast cell should drive dt down
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: dt, dt_quiescent
      real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp

      ! dt with no velocity
      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt_quiescent)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! One cell with high velocity
      hu(7, 7) = 10.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      ! dt should be smaller than quiescent case
      call check(error, dt < dt_quiescent, &
                 "Fast cell should reduce dt below quiescent value")
      if (allocated(error)) return
      call check(error, dt > 0.0_wp, "dt should be positive")
   end subroutine test_high_velocity_cell

   subroutine test_all_dry(error)
      !! All dry grid should give a fallback dt, not crash
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: dt
      real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp

      ! [acc->omp-cpu] dropped device directive (host no-op): data
      call compute_timestep(h, hu, hv, grid, CFL, dt)
      ! [acc->omp-cpu] dropped device directive (host no-op): end data

      call check(error, dt > 0.0_wp, "All-dry grid should give positive fallback dt")
      call check(error, dt < 1.0e10_wp, "Fallback dt should be finite")
   end subroutine test_all_dry

end module test_timestep
