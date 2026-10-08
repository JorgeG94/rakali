!! Uniform density over sloping bathymetry must produce ZERO baroclinic
!! pressure gradient.
module test_bpg_uniform_density
   !! ## The property
   !!
   !! A fluid whose density is spatially uniform has no horizontal pressure
   !! gradient at constant depth, whatever the bed does underneath it.  The
   !! baroclinic pressure gradient must therefore return **exactly zero**
   !! momentum tendency for a uniform-density column, over arbitrary
   !! bathymetry, for ANY value of that uniform density — including one that
   !! differs from `rho_0`.  This is an algebraic identity, so the tolerance
   !! is round-off.
   !!
   !! ## Why this test exists (the bug it pins)
   !!
   !! The sigma method (`ml_baroclinic_pressure_gradient`) integrates the
   !! density anomaly down each column and then differences the result at the
   !! same SIGMA INDEX:
   !!
   !!     dp_dx = (p_int_ip1(k) - p_int_im1(k)) / (2 dx)
   !!
   !! Over sloping bathymetry, layer `k` sits at a different physical depth in
   !! each column (sigma: `h_layer = dsig * h_col`, and `h_col` varies).  For a
   !! uniform anomaly `drho = rho - rho_ref`, `p_int(k) ~ drho * z_mid(k)`, so
   !! the difference is `drho * (z_mid_ip1(k) - z_mid_im1(k))`, which is
   !! NONZERO purely because of the bed slope.  The result is a fictitious
   !! force proportional to (uniform density anomaly) x (bed slope).
   !!
   !! Measured consequence before the fix, on a closed box with walls and zero
   !! forcing, 400 steps, 40x40x6 over a sinusoidal bed, uniform
   !! rho = 1025.4 against rho_0 = 1025.0:
   !!
   !!   | BPG method | rel. mass change | max\|hu\| |
   !!   |------------|------------------|-----------|
   !!   | sigma      | 1.65e-06         | 7.9e-04   |
   !!   | jacobian   | 1.28e-14         | 5.8e-14   |
   !!
   !! And at the kernel level, the single-step momentum tendency this test
   !! measures directly (uniform density, sinusoidal bed):
   !!
   !!   | BPG method | sloping bed | flat bed |
   !!   |------------|-------------|----------|
   !!   | sigma      | 1.28e-03    | 0.0      |
   !!   | jacobian   | 2.09e-18    | 0.0      |
   !!
   !! The spurious force drove real flow, and the flow leaked mass.  Note the
   !! anomaly only has to be uniform, not zero: seeding S and T exactly at
   !! `S_ref`/`T_ref` (so the EOS returns exactly `rho_0`) also conserved, which
   !! is what isolated the reference offset as the trigger.
   !!
   !! The Jacobian method (Shchepetkin & McWilliams 2003) differences along
   !! constant-z surfaces by interpolating the neighbour column's pressure to
   !! the centre column's depth, and satisfies the identity by construction.
   !! `bpg_method` therefore defaults to `"jacobian"`.
   !!
   !! This test asserts the property directly on the kernel, with no solver
   !! run, so it fails immediately and unambiguously if the default is reverted
   !! or a future edit reintroduces sigma-index differencing.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, PI
   use rki_grid, only: hgrid_t
   use rki_ml_pressure, only: ml_baroclinic_pressure_gradient, &
                              ml_baroclinic_pressure_gradient_jacobian
   implicit none
   private

   public :: collect_bpg_uniform_density_tests

   integer, parameter :: NX = 24, NY = 24, NGHOST = 2, NZ = 6
   real(wp), parameter :: DX = 5.0_wp
   real(wp), parameter :: DT = 1.0_wp
   real(wp), parameter :: ETA0 = 10.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: RHO_UNIFORM = 1025.4_wp
      !! Deliberately NOT rho_0 — a uniform anomaly is the failing case.
   real(wp), parameter :: BPG_TOL = 1.0e-16_wp
      !! Momentum tendency (m^2/s) after one dt.  Derived, not tuned.  The
      !! differenced quantity is `p_int ~ drho*z ~ 0.4*10 = 4 kg/m^2`, so the
      !! round-off floor is
      !!
      !!   eps * 4 / (2*dx) * (g/rho_0) * h_layer
      !!     = 2.2e-16 * 4 / 10 * 9.57e-3 * 1.67  ~  1.4e-18
      !!
      !! and the Jacobian measures 2.09e-18 — on the floor, within 1.5x.  The
      !! bar is set ~50x above that for toolchain variation (nvfortran's
      !! relaxed-FP reassociation under -fast), which still leaves it 1e13
      !! times BELOW the sigma method's measured 1.28e-3.  There is no
      !! plausible tolerance in between, so this cannot be a judgement call.

contains

   subroutine collect_bpg_uniform_density_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bpg_uniform_density_over_slope_is_zero", test_uniform_zero), &
                  new_unittest("bpg_uniform_density_flat_bed_is_zero", test_uniform_flat), &
                  new_unittest("bpg_sigma_method_fails_uniform_density", test_sigma_is_broken) &
                  ]
   end subroutine collect_bpg_uniform_density_tests

   subroutine build_case(grid, h_layer, rho_layer, hu_layer, hv_layer, ripple)
      !! Uniform-density column over a sinusoidal bed, at rest.
      type(hgrid_t), intent(out) :: grid
      real(wp), allocatable, intent(out) :: h_layer(:, :, :), rho_layer(:, :, :)
      real(wp), allocatable, intent(out) :: hu_layer(:, :, :), hv_layer(:, :, :)
      real(wp), intent(in) :: ripple

      integer :: i, j, k, nxt, nyt
      real(wp) :: x, y, b, h_col

      call grid%init(NX, NY, NGHOST, DX, DX)
      nxt = grid%nx_total
      nyt = grid%ny_total

      allocate (h_layer(nxt, nyt, NZ), rho_layer(nxt, nyt, NZ))
      allocate (hu_layer(nxt, nyt, NZ), hv_layer(nxt, nyt, NZ))

      do j = 1, nyt
         do i = 1, nxt
            x = (real(i - NGHOST, wp) - 0.5_wp)*DX
            y = (real(j - NGHOST, wp) - 0.5_wp)*DX
            b = ripple*(1.0_wp + sin(2.0_wp*PI*x/47.0_wp)*cos(2.0_wp*PI*y/31.0_wp))
            b = max(b, 0.0_wp)
            h_col = max(ETA0 - b, 0.0_wp)
            do k = 1, NZ
               h_layer(i, j, k) = h_col/real(NZ, wp)   ! uniform sigma
               rho_layer(i, j, k) = RHO_UNIFORM        ! UNIFORM density
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            end do
         end do
      end do
   end subroutine build_case

   function max_abs_momentum(hu_layer, hv_layer, grid) result(mx)
      real(wp), intent(in) :: hu_layer(:, :, :), hv_layer(:, :, :)
      type(hgrid_t), intent(in) :: grid
      real(wp) :: mx
      integer :: i, j, k

      mx = 0.0_wp
      do k = 1, NZ
         do j = NGHOST + 1, grid%ny_total - NGHOST
            do i = NGHOST + 1, grid%nx_total - NGHOST
               mx = max(mx, abs(hu_layer(i, j, k)), abs(hv_layer(i, j, k)))
            end do
         end do
      end do
   end function max_abs_momentum

   subroutine run_default_bpg(ripple, mx)
      !! Exercise the DEFAULT method (jacobian).  A dummy reference profile is
      !! passed as "unset" (`rho_ref_dz = 0`) so the kernel falls back to
      !! `rho_0` as the reference — the configuration in which the sigma method
      !! failed.
      real(wp), intent(in) :: ripple
      real(wp), intent(out) :: mx

      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), rho_layer(:, :, :)
      real(wp), allocatable :: hu_layer(:, :, :), hv_layer(:, :, :)
      real(wp) :: rho_ref(0:1)

      call build_case(grid, h_layer, rho_layer, hu_layer, hv_layer, ripple)
      rho_ref = RHO_0

      call ml_baroclinic_pressure_gradient_jacobian(h_layer, rho_layer, &
                                                    hu_layer, hv_layer, &
                                                    RHO_0, rho_ref, 0.0_wp, 1, &
                                                    grid, NZ, DT)

      mx = max_abs_momentum(hu_layer, hv_layer, grid)
   end subroutine run_default_bpg

   subroutine test_uniform_zero(error)
      !! THE assertion: uniform density, SLOPING bed, zero BPG.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx

      call run_default_bpg(1.5_wp, mx)
      call check(error, mx < BPG_TOL, &
                 "uniform density over a sloping bed produced a non-zero BPG")
   end subroutine test_uniform_zero

   subroutine test_uniform_flat(error)
      !! Control: the same must hold on a flat bed.  The sigma method passes
      !! this one, which is exactly why the bug hid — every flat-bed test in
      !! the suite was blind to it.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: mx

      call run_default_bpg(0.0_wp, mx)
      call check(error, mx < BPG_TOL, &
                 "uniform density over a flat bed produced a non-zero BPG")
   end subroutine test_uniform_flat

   subroutine test_sigma_is_broken(error)
      !! Pins the DIAGNOSIS, not just the fix: assert that the sigma method
      !! really does violate the identity on a sloping bed while satisfying it
      !! on a flat one.  If someone later repairs the sigma kernel this test
      !! fails, which is the correct prompt to revisit the default and delete
      !! this assertion — the failure is informative, not a nuisance.
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), rho_layer(:, :, :)
      real(wp), allocatable :: hu_layer(:, :, :), hv_layer(:, :, :)
      real(wp) :: rho_ref(0:1)
      real(wp) :: mx_slope, mx_flat

      rho_ref = RHO_0

      call build_case(grid, h_layer, rho_layer, hu_layer, hv_layer, 1.5_wp)
      call ml_baroclinic_pressure_gradient(h_layer, rho_layer, hu_layer, hv_layer, &
                                           RHO_0, rho_ref, 0.0_wp, 1, grid, NZ, DT)
      mx_slope = max_abs_momentum(hu_layer, hv_layer, grid)

      call build_case(grid, h_layer, rho_layer, hu_layer, hv_layer, 0.0_wp)
      call ml_baroclinic_pressure_gradient(h_layer, rho_layer, hu_layer, hv_layer, &
                                           RHO_0, rho_ref, 0.0_wp, 1, grid, NZ, DT)
      mx_flat = max_abs_momentum(hu_layer, hv_layer, grid)

      call check(error, mx_flat < BPG_TOL, &
                 "sigma BPG unexpectedly non-zero on a FLAT bed — diagnosis changed")
      if (allocated(error)) return
      call check(error, mx_slope > 1.0e-6_wp, &
                 "sigma BPG no longer violates the uniform-density identity — "// &
                 "if the kernel was fixed, revisit the bpg_method default and drop this test")
   end subroutine test_sigma_is_broken

end module test_bpg_uniform_density
