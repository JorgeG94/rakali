!! Unit tests for barotropic-baroclinic mode splitting
module test_mode_split
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_kernel_timestep, only: compute_timestep, compute_timestep_baroclinic
   use rki_kernel_pressure, only: compute_nh_baro_tendency, &
                                  apply_baro_nh_tendency, &
                                  reconcile_barotropic, &
                                  compute_depth_averaged_momentum
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_init, only: init_dam_break, compute_total_mass
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_mode_split_tests

   integer, parameter :: NX = 12, NY = 12, NG = 2, NZ = 3
   real(wp), parameter :: DX = 10.0_wp, DY = 10.0_wp
   type(hgrid_t), save :: grid
   logical, save :: grid_ready = .false.

contains

   subroutine ensure_grid()
      if (grid_ready) return
      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
      grid_ready = .true.
   end subroutine ensure_grid

   subroutine collect_mode_split_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      call ensure_grid()

      testsuite = [ &
                  new_unittest("baroclinic_cfl_slower", test_baroclinic_cfl), &
                  new_unittest("nh_baro_tendency_known_gradient", test_nh_baro_tendency), &
                  new_unittest("apply_baro_tendency", test_apply_baro_tendency), &
                  new_unittest("reconcile_consistent", test_reconcile_consistent), &
                  new_unittest("split_mass_conservation", test_split_mass_conservation), &
                  new_unittest("split_ratio_1_matches_unsplit", test_split_ratio_1_matches_unsplit), &
                  new_unittest("baroclinic_cfl_nz1", test_baroclinic_cfl_nz1) &
                  ]
   end subroutine collect_mode_split_tests

   subroutine test_baroclinic_cfl(error)
      !! Baroclinic CFL should give a larger timestep than barotropic CFL
      !! because internal wave speed c_int = sqrt(g*H/nz) < sqrt(g*H)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: dt_fast, dt_slow
      real(wp), parameter :: CFL = 0.45_wp
      integer :: k

      h = 10.0_wp
      hu = 10.0_wp   ! h*u = 10, so u = 1 m/s
      hv = 0.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      do k = 1, NZ
         u_layer(:, :, k) = 1.0_wp
         v_layer(:, :, k) = 0.0_wp
      end do

      !$acc data copyin( h, hu, hv, u_layer, v_layer, dsig)
      call compute_timestep(h, hu, hv, grid, CFL, dt_fast)
      call compute_timestep_baroclinic(h, u_layer, v_layer, dsig, &
                                       grid, NZ, CFL, dt_slow)
      !$acc end data

      ! dt_slow should be larger because c_int < c_barotropic
      call check(error, dt_slow > dt_fast, &
                 "Baroclinic timestep should be larger than barotropic")
      if (allocated(error)) return

      ! The ratio should be approximately sqrt(nz) for uniform velocity
      ! (c_baro = |u| + sqrt(g*H), c_int = |u| + sqrt(g*H/nz))
      ! Not exactly sqrt(nz) because velocity contributes, but should be > 1
      call check(error, dt_slow/dt_fast > 1.2_wp, &
                 "dt_slow/dt_fast should be significantly > 1 for nz=3")
   end subroutine test_baroclinic_cfl

   subroutine test_nh_baro_tendency(error)
      !! Known linear q_nh gradient should produce correct tendency
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), q_nh(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: tend_x(NX, NY), tend_y(NX, NY)
      integer :: i, j, k
      real(wp) :: expected_tend_x, dqdx_per_layer

      h = 1.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      q_nh = 0.0_wp

      ! Set up linear gradient in x: q_nh(i,j,k) = 2.0 * (i - NG - 0.5) * DX
      ! so dq/dx via central diff = (q(i+1) - q(i-1)) / (2*DX) = 2.0 for all k
      do k = 1, NZ
         do j = 1, NY
            do i = 1, NX
               q_nh(i, j, k) = 2.0_wp*real(i, wp)*DX
            end do
         end do
      end do

      ! Expected: tend_x = -sum(dsig_k * dq/dx) = -sum(1/nz * 2.0) = -2.0
      expected_tend_x = -2.0_wp

      !$acc data copyin( h, q_nh, dsig) copy( tend_x, tend_y)
      call compute_nh_baro_tendency(h, q_nh, dsig, tend_x, tend_y, &
                                    grid, NZ)
      !$acc end data

      call check(error, abs(tend_x(5, 5) - expected_tend_x) < 1.0e-10_wp, &
                 "NH baro tendency x should match analytical gradient")
      if (allocated(error)) return
      ! y gradient is zero -> tend_y should be zero
      call check(error, abs(tend_y(5, 5)) < 1.0e-10_wp, &
                 "NH baro tendency y should be zero for x-only gradient")
   end subroutine test_nh_baro_tendency

   subroutine test_apply_baro_tendency(error)
      !! apply_baro_nh_tendency adds dt * h * tend to hu, hv
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: hu(NX, NY), hv(NX, NY), h(NX, NY)
      real(wp) :: tend_x(NX, NY), tend_y(NX, NY)
      real(wp), parameter :: DT = 0.1_wp
      real(wp) :: expected_hu

      hu = 5.0_wp
      hv = 3.0_wp
      h = 2.0_wp
      tend_x = 1.0_wp
      tend_y = -0.5_wp

      ! Expected: hu += dt * h * tend_x = 5 + 0.1 * 2 * 1 = 5.2
      expected_hu = 5.0_wp + DT*2.0_wp*1.0_wp

      !$acc data copy( hu, hv) copyin( h, tend_x, tend_y)
      call apply_baro_nh_tendency(hu, hv, h, tend_x, tend_y, &
                                  grid, DT)
      !$acc end data

      call check(error, abs(hu(5, 5) - expected_hu) < 1.0e-12_wp, &
                 "hu += dt*h*tend_x")
      if (allocated(error)) return
      call check(error, abs(hv(5, 5) - (3.0_wp + DT*2.0_wp*(-0.5_wp))) < 1.0e-12_wp, &
                 "hv += dt*h*tend_y")
   end subroutine test_apply_baro_tendency

   subroutine test_reconcile_consistent(error)
      !! After reconciliation, hu == h * sum(dsig_k * u_layer_k)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: u_layer(NX, NY, NZ), v_layer(NX, NY, NZ)
      real(wp) :: dsig(NZ)
      real(wp) :: expected_hu
      integer :: k

      h = 3.0_wp
      dsig = 1.0_wp/real(NZ, wp)
      ! Set different velocities per layer
      do k = 1, NZ
         u_layer(:, :, k) = real(k, wp)  ! 1, 2, 3
         v_layer(:, :, k) = -real(k, wp)
      end do
      hu = 999.0_wp  ! will be overwritten
      hv = 999.0_wp

      ! Expected: hu = h * sum(dsig * u) = 3 * (1/3*1 + 1/3*2 + 1/3*3) = 3 * 2 = 6
      expected_hu = 3.0_wp*2.0_wp

      !$acc data copyin( h, u_layer, v_layer, dsig) copy( hu, hv)
      call reconcile_barotropic(h, hu, hv, u_layer, v_layer, dsig, &
                                grid, NZ)
      !$acc end data

      call check(error, abs(hu(5, 5) - expected_hu) < 1.0e-12_wp, &
                 "hu = h * sum(dsig*u) = 6 after reconciliation")
      if (allocated(error)) return
      call check(error, abs(hv(5, 5) - (-expected_hu)) < 1.0e-12_wp, &
                 "hv = h * sum(dsig*v) = -6 after reconciliation")
   end subroutine test_reconcile_consistent

   subroutine test_split_mass_conservation(error)
      !! Run 10 steps with split_ratio=2 on dam break IC with NH, verify mass conservation
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp) :: mass_before, mass_after
      integer :: n

      ! Set up NH state with split_ratio=2
      s%use_nonhydrostatic = .true.
      s%nh%nz = 2
      s%nh%split_ratio = 2
      s%nh%split_averaging = .true.
      call s%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      call init_dam_break(s)
      mass_before = compute_total_mass(s)

      call solver_enter_data(s)
      do n = 1, 10
         call solver_step(s)
      end do
      call solver_exit_data(s)

      mass_after = compute_total_mass(s)

      call check(error, abs(mass_after - mass_before)/mass_before < 1.0e-10_wp, &
                 "Mass should be conserved with mode splitting (split_ratio=2)")

      call s%destroy()
   end subroutine test_split_mass_conservation

   subroutine test_split_ratio_1_matches_unsplit(error)
      !! Mass conservation on both the unsplit (`split_ratio = 1`) and
      !! the auto-split (`split_ratio = 0`) NH paths, independently.
      !! Earlier versions of this test compared the two h fields with a
      !! `max_diff < 0.5 m` tolerance and a comment admitting the
      !! comparison was basically nonsense — the two paths integrate
      !! different discrete equations, so pointwise agreement was never
      !! the right thing to check. Conservation on each path IS the
      !! right invariant, and tight (closed basin + well-balanced
      !! barotropic step keeps mass to machine precision).
      !!
      !! Both runs also must stay physically bounded (h in the dam-break
      !! range, no NaN) so a blow-up in the split path can't slip past
      !! a pure conservation check.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s_split, s_unsplit
      real(wp) :: mass_before_unsplit, mass_after_unsplit
      real(wp) :: mass_before_split, mass_after_split
      integer :: n

      ! Unsplit path: split_ratio = 1 collapses mode splitting back to
      ! a single slow step (no fast substepping).
      s_unsplit%use_nonhydrostatic = .true.
      s_unsplit%nh%nz = 2
      s_unsplit%nh%split_ratio = 1
      call s_unsplit%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call init_dam_break(s_unsplit)

      ! Auto-split path: split_ratio = 0 lets the solver pick
      ! `ceiling(dt_slow/dt_fast)` substeps per slow step; for nz = 2
      ! that is at least 1 and usually > 1.
      s_split%use_nonhydrostatic = .true.
      s_split%nh%nz = 2
      s_split%nh%split_ratio = 0
      s_split%nh%split_averaging = .true.
      call s_split%init(20, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call init_dam_break(s_split)

      mass_before_unsplit = compute_total_mass(s_unsplit)
      mass_before_split = compute_total_mass(s_split)

      call solver_enter_data(s_unsplit)
      call solver_enter_data(s_split)

      do n = 1, 5
         call solver_step(s_unsplit)
         call solver_step(s_split)
      end do

      call solver_exit_data(s_unsplit)
      call solver_exit_data(s_split)

      mass_after_unsplit = compute_total_mass(s_unsplit)
      mass_after_split = compute_total_mass(s_split)

      call check(error, &
                 abs(mass_after_unsplit - mass_before_unsplit)/mass_before_unsplit &
                 < 1.0e-10_wp, &
                 "unsplit NH path must conserve mass to 1e-10")
      if (allocated(error)) return
      call check(error, &
                 abs(mass_after_split - mass_before_split)/mass_before_split &
                 < 1.0e-10_wp, &
                 "auto-split NH path must conserve mass to 1e-10")
      if (allocated(error)) return

      ! Physical-bound check — rules out NaN/Inf and gross blow-up
      ! without the pointwise-match trap.
      call check(error, &
                 maxval(s_unsplit%barotropic%h) < 3.0_wp .and. minval(s_unsplit%barotropic%h) >= 0.0_wp, &
                 "unsplit h stays in physical range after 5 NH steps")
      if (allocated(error)) return
      call check(error, &
                 maxval(s_split%barotropic%h) < 3.0_wp .and. minval(s_split%barotropic%h) >= 0.0_wp, &
                 "auto-split h stays in physical range after 5 NH steps")

      call s_unsplit%destroy()
      call s_split%destroy()
   end subroutine test_split_ratio_1_matches_unsplit

   subroutine test_baroclinic_cfl_nz1(error)
      !! When nz=1, baroclinic CFL should equal barotropic CFL
      !! because c_int = sqrt(g*H/1) = sqrt(g*H) = c_baro
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: u_layer(NX, NY, 1), v_layer(NX, NY, 1)
      real(wp) :: dsig(1)
      real(wp) :: dt_fast, dt_slow
      real(wp), parameter :: CFL = 0.45_wp

      h = 10.0_wp
      hu = 10.0_wp
      hv = 0.0_wp
      dsig(1) = 1.0_wp
      u_layer(:, :, 1) = 1.0_wp
      v_layer(:, :, 1) = 0.0_wp

      !$acc data copyin( h, hu, hv, u_layer, v_layer, dsig)
      call compute_timestep(h, hu, hv, grid, CFL, dt_fast)
      call compute_timestep_baroclinic(h, u_layer, v_layer, dsig, &
                                       grid, 1, CFL, dt_slow)
      !$acc end data

      ! With nz=1, wave speed is the same (sqrt(g*H)), but baroclinic uses
      ! raw |u_layer| while barotropic uses desingularized velocity, so
      ! dt_slow >= dt_fast. The ratio should be close to 1 (same wave speed).
      call check(error, dt_slow >= dt_fast, &
                 "dt_slow >= dt_fast when nz=1")
      if (allocated(error)) return
      call check(error, dt_slow/dt_fast < 2.0_wp, &
                 "dt_slow/dt_fast < 2 when nz=1 (same wave speed, similar velocity)")
   end subroutine test_baroclinic_cfl_nz1

end module test_mode_split
