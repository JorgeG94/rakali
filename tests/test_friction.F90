!! Unit tests for the bottom friction, wind stress, pressure, and precipitation kernel
module test_friction
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, RHO_WATER
   use rki_kernel_friction, only: apply_friction, apply_coriolis
   use rki_ml_friction, only: ml_apply_coriolis, ml_apply_bottom_friction
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_friction_tests

   integer, parameter :: NX = 8, NY = 8, NG = 2
   real(wp), parameter :: DX_TEST = 1.0_wp, DY_TEST = 1.0_wp

   type(hgrid_t), save :: grid
   logical, save :: grid_ready = .false.

contains

   subroutine ensure_grid()
      if (grid_ready) return
      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX_TEST, DY_TEST)
      grid_ready = .true.
   end subroutine ensure_grid

   subroutine collect_friction_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      call ensure_grid()

      testsuite = [ &
                  new_unittest("no_friction_no_wind_unchanged", test_no_friction_no_wind), &
                  new_unittest("manning_reduces_momentum", test_manning_reduces), &
                  new_unittest("ml_linear_bottom_drag_decay", test_ml_linear_bottom_drag), &
                  new_unittest("manning_semi_implicit_stable", test_manning_stable), &
                  new_unittest("wind_adds_momentum", test_wind_adds_momentum), &
                  new_unittest("dry_cells_unaffected", test_dry_unaffected), &
                  new_unittest("friction_plus_wind", test_friction_plus_wind), &
                  new_unittest("friction_symmetric", test_friction_symmetric), &
                  new_unittest("precipitation_adds_mass", test_precipitation_adds_mass), &
                  new_unittest("pressure_gradient_accelerates", test_pressure_gradient), &
                  new_unittest("coriolis_zero_f_unchanged", test_coriolis_zero_f), &
                  new_unittest("coriolis_rotates_momentum", test_coriolis_rotates), &
                  new_unittest("coriolis_conserves_energy", test_coriolis_energy), &
                  new_unittest("ml_coriolis_uniform_matches_depth_averaged", &
                               test_ml_coriolis_uniform_matches_depth_averaged), &
                  new_unittest("ml_coriolis_baroclinic_per_layer_energy", &
                               test_ml_coriolis_baroclinic_energy) &
                  ]
   end subroutine collect_friction_tests

   subroutine make_uniform_forcing(wx_val, wy_val, wx, wy, p_atm, precip)
      !! Helper: create uniform forcing arrays
      real(wp), intent(in) :: wx_val, wy_val
      real(wp), intent(out) :: wx(NX, NY), wy(NX, NY)
      real(wp), intent(out) :: p_atm(NX, NY), precip(NX, NY)

      wx = wx_val
      wy = wy_val
      p_atm = 0.0_wp
      precip = 0.0_wp
   end subroutine make_uniform_forcing

   subroutine test_no_friction_no_wind(error)
      !! Zero Manning and zero wind should leave momentum unchanged
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)

      h = 1.0_wp
      hu = 5.0_wp
      hv = 3.0_wp
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)

      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu, hv)
      call apply_friction(h, hu, hv, grid, 0.1_wp, 0.0_wp, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      call check(error, abs(hu(4, 4) - 5.0_wp) < 1.0e-14_wp, &
                 "hu unchanged with zero friction and zero wind")
      if (allocated(error)) return
      call check(error, abs(hv(4, 4) - 3.0_wp) < 1.0e-14_wp, &
                 "hv unchanged with zero friction and zero wind")
   end subroutine test_no_friction_no_wind

   subroutine test_manning_reduces(error)
      !! Analytic semi-implicit Manning decay:
      !!    hu^{n+1} = hu^n / (1 + dt * c_f * |u| / h)
      !!    c_f     = g * n^2 / h^{1/3}
      !! so the ratio (hu_new / hu_old) must equal 1 / denom exactly.
      !! Inputs: h=1, hu=5, hv=3, n=0.025, dt=0.1 => |u| = sqrt(34).
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: H_VAL = 1.0_wp, HU_VAL = 5.0_wp, HV_VAL = 3.0_wp
      real(wp), parameter :: N_MAN = 0.025_wp, DT_VAL = 0.1_wp
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)
      real(wp) :: speed, cf, denom, hu_expected, hv_expected

      h = H_VAL
      hu = HU_VAL
      hv = HV_VAL
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)

      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu, hv)
      call apply_friction(h, hu, hv, grid, DT_VAL, N_MAN, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      speed = sqrt((HU_VAL/H_VAL)**2 + (HV_VAL/H_VAL)**2)
      cf = GRAVITY*N_MAN*N_MAN/(H_VAL**(1.0_wp/3.0_wp))
      denom = 1.0_wp + DT_VAL*cf*speed/H_VAL
      hu_expected = HU_VAL/denom
      hv_expected = HV_VAL/denom

      call check(error, abs(hu(4, 4) - hu_expected) < 1.0e-12_wp, &
                 "hu must match analytic semi-implicit Manning decay")
      if (allocated(error)) return
      call check(error, abs(hv(4, 4) - hv_expected) < 1.0e-12_wp, &
                 "hv must match analytic semi-implicit Manning decay")
      if (allocated(error)) return
      ! hu/hv must decay by the same factor: friction acts isotropically.
      call check(error, abs(hu(4, 4)/HU_VAL - hv(4, 4)/HV_VAL) < 1.0e-12_wp, &
                 "hu and hv must decay by the same factor")
   end subroutine test_manning_reduces

   subroutine test_ml_linear_bottom_drag(error)
      !! Linear (Rayleigh) bottom drag on the multilayer bottom layer: the
      !! bottom layer decays by the exact semi-implicit factor 1/(1 + dt*r),
      !! independent of speed, and the interior layers are untouched.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 3
      real(wp), parameter :: HU_VAL = 2.0_wp, HV_VAL = -1.5_wp
      real(wp), parameter :: R = 1.0e-3_wp, DT = 50.0_wp
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: denom

      h_layer = 1.0_wp
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp
      hu_layer(:, :, 1) = HU_VAL    ! bottom layer (k=1, ROMS bottom-up)
      hv_layer(:, :, 1) = HV_VAL

      !$omp target data map(to:  h_layer) map(tofrom:  hu_layer, hv_layer)
      call ml_apply_bottom_friction(h_layer, hu_layer, hv_layer, &
                                    0.0_wp, .true., R, grid, NZ, DT)
      !$omp end target data

      denom = 1.0_wp + DT*R
      call check(error, abs(hu_layer(4, 4, 1) - HU_VAL/denom) < 1.0e-12_wp, &
                 "linear drag: bottom hu = hu / (1 + dt*r)")
      if (allocated(error)) return
      call check(error, abs(hv_layer(4, 4, 1) - HV_VAL/denom) < 1.0e-12_wp, &
                 "linear drag: bottom hv = hv / (1 + dt*r)")
      if (allocated(error)) return
      ! Drag is bottom-only: surface layer must be untouched.
      call check(error, abs(hu_layer(4, 4, NZ)) < 1.0e-14_wp, &
                 "linear drag: surface layer unchanged")
   end subroutine test_ml_linear_bottom_drag

   subroutine test_manning_stable(error)
      !! Semi-implicit friction should not overshoot even with large dt or thin layer
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)

      h = 0.01_wp
      hu = 1.0_wp
      hv = 0.5_wp
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)

      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu, hv)
      call apply_friction(h, hu, hv, grid, 100.0_wp, 0.05_wp, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      call check(error, hu(4, 4) >= 0.0_wp, &
                 "Semi-implicit friction should not overshoot (hu)")
      if (allocated(error)) return
      call check(error, hv(4, 4) >= 0.0_wp, &
                 "Semi-implicit friction should not overshoot (hv)")
      if (allocated(error)) return
      call check(error, hu(4, 4) < 1.0_wp, "Friction should reduce hu")
   end subroutine test_manning_stable

   subroutine test_wind_adds_momentum(error)
      !! Wind stress should add momentum proportional to dt * tau / rho
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)
      real(wp), parameter :: DT = 0.5_wp
      real(wp), parameter :: WIND_X = 0.1_wp, WIND_Y = -0.05_wp
      real(wp) :: expected_hu, expected_hv

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      call make_uniform_forcing(WIND_X, WIND_Y, wx, wy, p_atm, precip)

      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu, hv)
      call apply_friction(h, hu, hv, grid, DT, 0.0_wp, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      expected_hu = DT*WIND_X/RHO_WATER
      expected_hv = DT*WIND_Y/RHO_WATER

      call check(error, abs(hu(4, 4) - expected_hu) < 1.0e-14_wp, &
                 "hu should equal dt * wind_x / rho")
      if (allocated(error)) return
      call check(error, abs(hv(4, 4) - expected_hv) < 1.0e-14_wp, &
                 "hv should equal dt * wind_y / rho")
   end subroutine test_wind_adds_momentum

   subroutine test_dry_unaffected(error)
      !! Dry cells should not be modified by friction or wind
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)

      h = 0.0_wp
      hu = 100.0_wp
      hv = 50.0_wp
      call make_uniform_forcing(1.0_wp, 1.0_wp, wx, wy, p_atm, precip)

      !$omp target data map(tofrom:  h) map(tofrom:  hu, hv) map(to:  wx, wy, p_atm, precip)
      call apply_friction(h, hu, hv, grid, 1.0_wp, 0.025_wp, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      call check(error, abs(hu(4, 4) - 100.0_wp) < 1.0e-14_wp, &
                 "Dry cell hu should be unchanged by friction")
      if (allocated(error)) return
      call check(error, abs(hv(4, 4) - 50.0_wp) < 1.0e-14_wp, &
                 "Dry cell hv should be unchanged by friction")
   end subroutine test_dry_unaffected

   subroutine test_friction_plus_wind(error)
      !! Friction reduces momentum, then wind adds -- verify combined effect
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)
      real(wp) :: hu_no_wind, hu_with_wind
      real(wp), parameter :: DT = 0.1_wp, MANNING = 0.025_wp

      h = 1.0_wp

      ! First: friction only
      hu = 5.0_wp
      hv = 0.0_wp
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)
      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu, hv)
      call apply_friction(h, hu, hv, grid, DT, MANNING, &
                          wx, wy, p_atm, precip)
      !$omp end target data
      hu_no_wind = hu(4, 4)

      ! Second: friction + positive wind
      hu = 5.0_wp
      hv = 0.0_wp
      call make_uniform_forcing(10.0_wp, 0.0_wp, wx, wy, p_atm, precip)
      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu, hv)
      call apply_friction(h, hu, hv, grid, DT, MANNING, &
                          wx, wy, p_atm, precip)
      !$omp end target data
      hu_with_wind = hu(4, 4)

      call check(error, hu_with_wind > hu_no_wind, &
                 "Positive wind should increase hu compared to friction-only")
   end subroutine test_friction_plus_wind

   subroutine test_friction_symmetric(error)
      !! Friction should be direction-symmetric: same magnitude reduction for +/- velocity
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu_pos(NX, NY), hv(NX, NY)
      real(wp) :: hu_neg(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)
      real(wp), parameter :: DT = 0.1_wp, MANNING = 0.03_wp

      h = 1.0_wp
      hv = 0.0_wp
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)

      hu_pos = 5.0_wp
      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu_pos, hv)
      call apply_friction(h, hu_pos, hv, grid, DT, MANNING, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      hv = 0.0_wp
      hu_neg = -5.0_wp
      !$omp target data map(to:  h, wx, wy, p_atm, precip) map(tofrom:  hu_neg, hv)
      call apply_friction(h, hu_neg, hv, grid, DT, MANNING, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      call check(error, abs(hu_pos(4, 4) + hu_neg(4, 4)) < 1.0e-12_wp, &
                 "Friction on +hu and -hu should give equal and opposite results")
   end subroutine test_friction_symmetric

   subroutine test_precipitation_adds_mass(error)
      !! Precipitation should increase water depth by dt * P
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)
      real(wp), parameter :: DT = 1.0_wp
      real(wp), parameter :: RAIN_RATE = 0.001_wp  ! 1 mm/s
      real(wp) :: h_before

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      h_before = 1.0_wp
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)
      precip = RAIN_RATE

      !$omp target data map(tofrom:  h) map(tofrom:  hu, hv) map(to:  wx, wy, p_atm, precip)
      call apply_friction(h, hu, hv, grid, DT, 0.0_wp, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      call check(error, abs(h(4, 4) - (h_before + DT*RAIN_RATE)) < 1.0e-14_wp, &
                 "Precipitation should add dt * P to water depth")
   end subroutine test_precipitation_adds_mass

   subroutine test_pressure_gradient(error)
      !! Atmospheric pressure gradient should accelerate flow
      !! High pressure on left, low on right -> flow in +x direction
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), p_atm(NX, NY), precip(NX, NY)
      real(wp), parameter :: DT = 1.0_wp
      integer :: i

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      call make_uniform_forcing(0.0_wp, 0.0_wp, wx, wy, p_atm, precip)

      ! Linear pressure gradient in x: p = 101325 - 100*x
      ! dp/dx = -100 Pa/m -> force in +x direction (-(h/rho)*dp/dx > 0)
      do i = 1, NX
         p_atm(i, :) = 101325.0_wp - 100.0_wp*real(i, wp)*DX_TEST
      end do

      !$omp target data map(tofrom:  h) map(tofrom:  hu, hv) map(to:  wx, wy, p_atm, precip)
      call apply_friction(h, hu, hv, grid, DT, 0.0_wp, &
                          wx, wy, p_atm, precip)
      !$omp end target data

      ! dp/dx = -100, so -(h/rho)*dp/dx = -(1.0/1025)*(-100) > 0
      call check(error, hu(4, 4) > 0.0_wp, &
                 "Negative pressure gradient should accelerate flow in +x")
      if (allocated(error)) return

      ! Expected: hu += dt * (-h/rho) * dp/dx = 1.0 * (-1.0/1025.0) * (-100) = 0.0976
      call check(error, abs(hu(4, 4) - DT*1.0_wp*100.0_wp/RHO_WATER) < 1.0e-10_wp, &
                 "Pressure gradient acceleration should match -(h/rho)*dp/dx")
      if (allocated(error)) return

      ! No acceleration in y (uniform p in y)
      call check(error, abs(hv(4, 4)) < 1.0e-14_wp, &
                 "No y-acceleration with x-only pressure gradient")
   end subroutine test_pressure_gradient

   subroutine test_coriolis_zero_f(error)
      !! Zero Coriolis parameter should leave momentum unchanged
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)

      h = 1.0_wp
      hu = 5.0_wp
      hv = 3.0_wp

      !$omp target data map(to:  h) map(tofrom:  hu, hv)
      call apply_coriolis(h, hu, hv, grid, 0.1_wp, 0.0_wp)
      !$omp end target data

      call check(error, abs(hu(4, 4) - 5.0_wp) < 1.0e-14_wp, &
                 "Zero f should not change hu")
      if (allocated(error)) return

      call check(error, abs(hv(4, 4) - 3.0_wp) < 1.0e-14_wp, &
                 "Zero f should not change hv")
   end subroutine test_coriolis_zero_f

   subroutine test_coriolis_rotates(error)
      !! Coriolis should rotate momentum: positive f with (hu,0) should
      !! produce positive hu and negative hv (northern hemisphere deflection)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp), parameter :: F_COR = 1.0e-4_wp  ! ~mid-latitude
      real(wp), parameter :: DT = 100.0_wp

      h = 1.0_wp
      hu = 1.0_wp
      hv = 0.0_wp

      !$omp target data map(to:  h) map(tofrom:  hu, hv)
      call apply_coriolis(h, hu, hv, grid, DT, F_COR)
      !$omp end target data

      ! Positive f with initial u>0 should deflect flow to the right (hv < 0)
      call check(error, hv(4, 4) < 0.0_wp, &
                 "Coriolis with f>0 and u>0 should produce hv<0 (rightward deflection)")
      if (allocated(error)) return

      ! hu should decrease slightly (energy transferred to hv)
      call check(error, hu(4, 4) < 1.0_wp, &
                 "Coriolis rotation should reduce hu when hv starts at zero")
   end subroutine test_coriolis_rotates

   subroutine test_coriolis_energy(error)
      !! Semi-implicit Coriolis should conserve kinetic energy exactly:
      !! hu^2 + hv^2 should be unchanged
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp), parameter :: F_COR = 1.0e-4_wp
      real(wp), parameter :: DT = 500.0_wp
      real(wp) :: ke_before, ke_after

      h = 1.0_wp
      hu = 3.0_wp
      hv = 4.0_wp
      ke_before = hu(4, 4)**2 + hv(4, 4)**2

      !$omp target data map(to:  h) map(tofrom:  hu, hv)
      call apply_coriolis(h, hu, hv, grid, DT, F_COR)
      !$omp end target data

      ke_after = hu(4, 4)**2 + hv(4, 4)**2

      call check(error, abs(ke_after - ke_before)/ke_before < 1.0e-12_wp, &
                 "Semi-implicit Coriolis should conserve hu^2 + hv^2 exactly")
   end subroutine test_coriolis_energy

   subroutine test_ml_coriolis_uniform_matches_depth_averaged(error)
      !! For a column with the SAME velocity in every layer, the per-layer
      !! kernel must produce a result whose depth-averaged sum equals what
      !! the existing depth-averaged apply_coriolis produces from the same
      !! IC.  This is the regression check that today's barotropic Coriolis
      !! is preserved as a special case of the new per-layer kernel.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 4
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp) :: h2d(NX, NY), hu2d(NX, NY), hv2d(NX, NY)
      real(wp), parameter :: F_COR = 1.0e-4_wp
      real(wp), parameter :: DT = 500.0_wp
      real(wp) :: hu_sum, hv_sum
      integer :: k

      ! Uniform sigma: each layer is H/NZ thick, same velocity in each layer.
      h_layer = 1.0_wp/real(NZ, wp)   ! layer thickness
      hu_layer = 0.75_wp              ! per-layer hu (depth-integrated -> 3.0 across NZ=4)
      hv_layer = 1.0_wp               ! per-layer hv (depth-integrated -> 4.0)

      ! Depth-averaged equivalent: h = sum h_layer = 1, hu = sum hu_layer = 3, hv = 4
      h2d = 1.0_wp
      hu2d = 3.0_wp
      hv2d = 4.0_wp

      !$omp target data map(to:  h_layer) map(tofrom:  hu_layer, hv_layer)
      call ml_apply_coriolis(h_layer, hu_layer, hv_layer, grid, NZ, DT, F_COR)
      !$omp end target data

      !$omp target data map(to:  h2d) map(tofrom:  hu2d, hv2d)
      call apply_coriolis(h2d, hu2d, hv2d, grid, DT, F_COR)
      !$omp end target data

      ! Sum the rotated layer momenta and compare to the depth-averaged result
      hu_sum = 0.0_wp
      hv_sum = 0.0_wp
      do k = 1, NZ
         hu_sum = hu_sum + hu_layer(4, 4, k)
         hv_sum = hv_sum + hv_layer(4, 4, k)
      end do

      call check(error, abs(hu_sum - hu2d(4, 4)) < 1.0e-12_wp, &
                 "Per-layer rotation summed to depth-averaged matches barotropic rotation (hu)")
      if (allocated(error)) return
      call check(error, abs(hv_sum - hv2d(4, 4)) < 1.0e-12_wp, &
                 "Per-layer rotation summed to depth-averaged matches barotropic rotation (hv)")
   end subroutine test_ml_coriolis_uniform_matches_depth_averaged

   subroutine test_ml_coriolis_baroclinic_energy(error)
      !! With baroclinic shear (different velocity in each layer), each
      !! layer's kinetic energy hu_k^2 + hv_k^2 must be conserved exactly
      !! by the closed-form rotation.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 4
      real(wp) :: h_layer(NX, NY, NZ), hu_layer(NX, NY, NZ), hv_layer(NX, NY, NZ)
      real(wp), parameter :: F_COR = 1.0e-4_wp
      real(wp), parameter :: DT = 500.0_wp
      real(wp) :: ke_before(NZ), ke_after(NZ)
      integer :: k

      h_layer = 0.5_wp

      ! Baroclinic shear: layer k=1 (bottom) at rest, k=NZ (surface) strongest
      do k = 1, NZ
         hu_layer(:, :, k) = real(k - 1, wp)
         hv_layer(:, :, k) = real(NZ - k, wp)
         ke_before(k) = hu_layer(4, 4, k)**2 + hv_layer(4, 4, k)**2
      end do

      !$omp target data map(to:  h_layer) map(tofrom:  hu_layer, hv_layer)
      call ml_apply_coriolis(h_layer, hu_layer, hv_layer, grid, NZ, DT, F_COR)
      !$omp end target data

      do k = 1, NZ
         ke_after(k) = hu_layer(4, 4, k)**2 + hv_layer(4, 4, k)**2
         if (ke_before(k) > 0.0_wp) then
            call check(error, abs(ke_after(k) - ke_before(k))/ke_before(k) < 1.0e-12_wp, &
                       "Per-layer Coriolis conserves kinetic energy in every layer")
         else
            call check(error, abs(ke_after(k) - ke_before(k)) < 1.0e-12_wp, &
                       "Per-layer Coriolis preserves zero KE in resting layers")
         end if
         if (allocated(error)) return
      end do
   end subroutine test_ml_coriolis_baroclinic_energy

end module test_friction
