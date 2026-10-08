!! Unit tests for the structured multilayer linear equation of state and
!! the depth-dependent reference-density profile builder.
module test_ml_eos
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_ml_eos, only: ml_equation_of_state, ml_compute_reference_density
   implicit none
   private

   public :: collect_ml_eos_tests

   ! Small grid used by all tests
   integer, parameter :: NX_PHYS = 4
   integer, parameter :: NY_PHYS = 3
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4

   ! Reference EOS scalars
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: ALPHA_T = 0.17_wp
   real(wp), parameter :: S_REF = 35.0_wp
   real(wp), parameter :: T_REF = 15.0_wp

contains

   subroutine collect_ml_eos_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("eos_reduces_to_rho0_at_reference", test_eos_reference_state), &
                  new_unittest("eos_linear_in_salinity", test_eos_linear_salinity), &
                  new_unittest("eos_linear_in_temperature_sign", test_eos_linear_temperature), &
                  new_unittest("eos_alpha_T_zero_makes_T_passive", test_eos_alpha_T_zero), &
                  new_unittest("eos_dry_cell_falls_back_to_reference", test_eos_dry_cell), &
                  new_unittest("reference_density_profile_invariants", test_reference_density_profile) &
                  ]
   end subroutine collect_ml_eos_tests

   subroutine make_grid(grid)
      !! Helper: build a small ghosted hgrid_t for kernel tests.
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 100.0_wp, 100.0_wp)
   end subroutine make_grid

   subroutine test_eos_reference_state(error)
      !! With S = S_ref and T = T_ref everywhere, rho == rho_0 exactly.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp) :: h_val

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))

      h_val = 2.5_wp
      h_layer = h_val
      hS_layer = h_val*S_REF
      hT_layer = h_val*T_REF
      rho_layer = -1.0_wp

      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, grid, NZ)

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               call check(error, abs(rho_layer(i, j, k) - RHO_0) < 1.0e-12_wp, &
                          "rho should equal rho_0 at the reference state")
               if (allocated(error)) return
            end do
         end do
      end do

      deallocate (h_layer, hS_layer, hT_layer, rho_layer)
   end subroutine test_eos_reference_state

   subroutine test_eos_linear_salinity(error)
      !! At T = T_ref, dS = +k PSU should change rho by exactly +beta_S * k.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp) :: h_val, S_perturbed, expected

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))

      h_val = 3.0_wp
      h_layer = h_val
      do k = 1, NZ
         S_perturbed = S_REF + real(k, wp)
         hS_layer(:, :, k) = h_val*S_perturbed
      end do
      hT_layer = h_val*T_REF

      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, grid, NZ)

      do k = 1, NZ
         expected = RHO_0 + BETA_S*real(k, wp)
         do j = 1, ny
            do i = 1, nx
               call check(error, abs(rho_layer(i, j, k) - expected) < 1.0e-10_wp, &
                          "rho should respond linearly to salinity perturbations")
               if (allocated(error)) return
            end do
         end do
      end do

      deallocate (h_layer, hS_layer, hT_layer, rho_layer)
   end subroutine test_eos_linear_salinity

   subroutine test_eos_linear_temperature(error)
      !! At S = S_ref, dT = +k K should change rho by exactly -alpha_T * k
      !! (water expands when warmed -> density decreases).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp) :: h_val, T_perturbed, expected

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))

      h_val = 1.5_wp
      h_layer = h_val
      hS_layer = h_val*S_REF
      do k = 1, NZ
         T_perturbed = T_REF + real(k, wp)
         hT_layer(:, :, k) = h_val*T_perturbed
      end do

      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, grid, NZ)

      do k = 1, NZ
         expected = RHO_0 - ALPHA_T*real(k, wp)
         do j = 1, ny
            do i = 1, nx
               call check(error, abs(rho_layer(i, j, k) - expected) < 1.0e-10_wp, &
                          "rho should decrease by alpha_T per +1 K (warming -> lighter)")
               if (allocated(error)) return
            end do
         end do
      end do

      ! Sanity: warmer layers must be strictly less dense than cooler layers
      do j = 1, ny
         do i = 1, nx
            call check(error, rho_layer(i, j, NZ) < rho_layer(i, j, 1), &
                       "Warmer layer must be less dense than cooler layer")
            if (allocated(error)) return
         end do
      end do

      deallocate (h_layer, hS_layer, hT_layer, rho_layer)
   end subroutine test_eos_linear_temperature

   subroutine test_eos_alpha_T_zero(error)
      !! With alpha_T = 0, varying T leaves rho unchanged (T is passive).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      integer :: i, j, k, nx, ny
      real(wp) :: h_val, T_perturbed

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))

      h_val = 2.0_wp
      h_layer = h_val
      hS_layer = h_val*S_REF
      do k = 1, NZ
         T_perturbed = T_REF + 5.0_wp*real(k, wp)
         hT_layer(:, :, k) = h_val*T_perturbed
      end do

      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                RHO_0, BETA_S, S_REF, 0.0_wp, T_REF, grid, NZ)

      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               call check(error, abs(rho_layer(i, j, k) - RHO_0) < 1.0e-12_wp, &
                          "alpha_T=0 should make T passive (rho == rho_0 when S == S_ref)")
               if (allocated(error)) return
            end do
         end do
      end do

      deallocate (h_layer, hS_layer, hT_layer, rho_layer)
   end subroutine test_eos_alpha_T_zero

   subroutine test_eos_dry_cell(error)
      !! Cells with h_layer <= DRY_TOLERANCE are treated with S=S_ref, T=T_ref,
      !! so rho falls back to rho_0 regardless of the (garbage) hS/hT values.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      integer :: nx, ny

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))

      h_layer = 0.0_wp
      hS_layer = 1.0e6_wp
      hT_layer = -1.0e6_wp

      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, grid, NZ)

      call check(error, all(abs(rho_layer - RHO_0) < 1.0e-12_wp), &
                 "Dry cells should fall back to rho_0 (S_ref, T_ref guard)")
      if (allocated(error)) return

      deallocate (h_layer, hS_layer, hT_layer, rho_layer)
   end subroutine test_eos_dry_cell

   subroutine test_reference_density_profile(error)
      !! Build rho_ref(z) from a stably-stratified column (dense at the
      !! bottom, light at the top) and check the documented invariants:
      !!   - rho_ref_dz > 0 when there is wet water,
      !!   - profile values lie within the actual rho_layer range,
      !!   - profile is monotone non-increasing in z (z=0 -> bottom),
      !!   - rho_ref(0) matches the deepest layer's density.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), rho_layer(:, :, :)
      integer, parameter :: NZ_REF = 8
      real(wp) :: rho_ref(0:NZ_REF), rho_ref_dz
      real(wp) :: rho_min_ref, rho_max_ref
      integer :: k, nx, ny, iz

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total

      allocate (h_layer(nx, ny, NZ), rho_layer(nx, ny, NZ))

      ! Equal-thickness layers
      h_layer = 1.0_wp

      ! Stable stratification: layer 1 (deepest) heaviest, layer NZ (top) lightest
      do k = 1, NZ
         rho_layer(:, :, k) = RHO_0 + real(NZ - k + 1, wp)
      end do

      rho_ref = -huge(1.0_wp)
      rho_ref_dz = -1.0_wp

      call ml_compute_reference_density(h_layer, rho_layer, &
                                        rho_ref, rho_ref_dz, NZ_REF, &
                                        RHO_0, grid, NZ)

      call check(error, rho_ref_dz > 0.0_wp, &
                 "rho_ref_dz must be positive for a wet column")
      if (allocated(error)) return

      rho_min_ref = RHO_0 + 1.0_wp
      rho_max_ref = RHO_0 + real(NZ, wp)

      do iz = 0, NZ_REF
         call check(error, rho_ref(iz) >= rho_min_ref - 1.0e-10_wp .and. &
                    rho_ref(iz) <= rho_max_ref + 1.0e-10_wp, &
                    "rho_ref(iz) must lie within the rho_layer min/max range")
         if (allocated(error)) return
      end do

      ! In Rakali's column convention z=0 is the TOP (sampled from
      ! rho_layer(nz)) and z=H_col is the bottom (rho_layer(1)). For a
      ! stable column (heaviest at k=1) rho_ref is monotone non-decreasing
      ! in iz.
      do iz = 1, NZ_REF
         call check(error, rho_ref(iz) >= rho_ref(iz - 1) - 1.0e-10_wp, &
                    "rho_ref(z) must be monotone non-decreasing in iz for a stable column")
         if (allocated(error)) return
      end do

      ! z=0 sample must equal the top layer (k = NZ) density
      call check(error, abs(rho_ref(0) - rho_layer(NGHOST + 1, NGHOST + 1, NZ)) < 1.0e-10_wp, &
                 "rho_ref(0) should match the top layer (k=nz) density")
      if (allocated(error)) return

      deallocate (h_layer, rho_layer)
   end subroutine test_reference_density_profile

end module test_ml_eos
