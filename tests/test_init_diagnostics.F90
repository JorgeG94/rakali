!! Unit tests for the diagnostic sums in rki_init
module test_init_diagnostics
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_init, only: compute_total_mass, compute_layer_sums
   implicit none
   private

   public :: collect_init_diagnostics_tests

   integer, parameter :: NX = 6, NY = 4, NG = 2, NZ = 3
   real(wp), parameter :: DX = 2.0_wp, DY = 3.0_wp

contains

   subroutine collect_init_diagnostics_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("total_mass_constant_h", test_total_mass_constant), &
                  new_unittest("total_mass_excludes_ghost", test_total_mass_ghost), &
                  new_unittest("total_mass_skips_garbage", test_total_mass_garbage), &
                  new_unittest("layer_sums_constant", test_layer_sums_constant), &
                  new_unittest("layer_sums_no_multilayer", test_layer_sums_disabled) &
                  ]
   end subroutine collect_init_diagnostics_tests

   subroutine test_total_mass_constant(error)
      !! Uniform h over interior cells: mass = h * nx_phys * ny_phys * dx * dy.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp), parameter :: H_VAL = 2.5_wp
      real(wp) :: mass, expected

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)
      s%barotropic%h = H_VAL

      mass = compute_total_mass(s)
      expected = H_VAL*real(NX, wp)*real(NY, wp)*DX*DY

      call check(error, abs(mass - expected) < 1.0e-12_wp*expected, &
                 "mass should equal h * nx_phys * ny_phys * dx * dy")

      call s%destroy()
   end subroutine test_total_mass_constant

   subroutine test_total_mass_ghost(error)
      !! Ghost cells must be excluded from the sum even if they hold
      !! large values. Sets interior=1.0 and the ghost band to a poison
      !! value; if the loop bounds were wrong the result would be huge.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp) :: mass, expected
      integer :: i, j

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      s%barotropic%h = 999.0_wp
      do j = NG + 1, s%grid%ny_total - NG
         do i = NG + 1, s%grid%nx_total - NG
            s%barotropic%h(i, j) = 1.0_wp
         end do
      end do

      mass = compute_total_mass(s)
      expected = real(NX, wp)*real(NY, wp)*DX*DY

      call check(error, abs(mass - expected) < 1.0e-12_wp*expected, &
                 "mass should sum only interior cells")

      call s%destroy()
   end subroutine test_total_mass_ghost

   subroutine test_total_mass_garbage(error)
      !! The NaN / dry-cell guard (h <= 0 .or. h >= 1e6) must drop
      !! poisoned cells from the sum.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp) :: mass, expected
      integer :: i_lo, j_lo

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)
      s%barotropic%h = 1.0_wp

      i_lo = NG + 1
      j_lo = NG + 1
      s%barotropic%h(i_lo, j_lo) = -5.0_wp           ! negative — skipped
      s%barotropic%h(i_lo + 1, j_lo) = 1.0e10_wp     ! unphysical — skipped

      mass = compute_total_mass(s)
      ! 2 of NX*NY interior cells are poisoned and dropped
      expected = (real(NX*NY, wp) - 2.0_wp)*DX*DY

      call check(error, abs(mass - expected) < 1.0e-12_wp*expected, &
                 "mass should skip negative and unphysical h values")

      call s%destroy()
   end subroutine test_total_mass_garbage

   subroutine test_layer_sums_constant(error)
      !! Uniform layer fields: each sum = value * nx_phys * ny_phys * nz_ml.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp), parameter :: H_LAYER = 2.0_wp
      real(wp), parameter :: HS_LAYER = 70.0_wp   ! h*S with S=35 PSU
      real(wp), parameter :: HT_LAYER = 30.0_wp   ! h*T with T=15 degC
      real(wp) :: sum_h, sum_hS, sum_hT, n_cells

      s%use_multilayer = .true.
      s%multilayer%nz_ml = NZ
      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      s%multilayer%h_layer = H_LAYER
      s%multilayer%tracers(s%multilayer%idx_salinity)%hTr = HS_LAYER
      s%multilayer%tracers(s%multilayer%idx_temperature)%hTr = HT_LAYER

      call compute_layer_sums(s, sum_h, sum_hS, sum_hT)
      n_cells = real(NX*NY*NZ, wp)

      call check(error, abs(sum_h - H_LAYER*n_cells) < 1.0e-12_wp*sum_h, &
                 "sum_h should equal h_layer * nx_phys * ny_phys * nz_ml")
      if (allocated(error)) then; call s%destroy(); return; end if
      call check(error, abs(sum_hS - HS_LAYER*n_cells) < 1.0e-12_wp*sum_hS, &
                 "sum_hS should equal hS_layer * nx_phys * ny_phys * nz_ml")
      if (allocated(error)) then; call s%destroy(); return; end if
      call check(error, abs(sum_hT - HT_LAYER*n_cells) < 1.0e-12_wp*sum_hT, &
                 "sum_hT should equal hT_layer * nx_phys * ny_phys * nz_ml")

      call s%destroy()
   end subroutine test_layer_sums_constant

   subroutine test_layer_sums_disabled(error)
      !! With multilayer off, the routine returns zeros without touching
      !! the (unallocated) layer arrays.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp) :: sum_h, sum_hS, sum_hT

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      sum_h = 1.0_wp; sum_hS = 1.0_wp; sum_hT = 1.0_wp
      call compute_layer_sums(s, sum_h, sum_hS, sum_hT)

      call check(error, sum_h == 0.0_wp .and. sum_hS == 0.0_wp .and. &
                 sum_hT == 0.0_wp, &
                 "layer sums should return zero when use_multilayer is false")

      call s%destroy()
   end subroutine test_layer_sums_disabled

end module test_init_diagnostics
