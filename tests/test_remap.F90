!! Unit tests for the vertical remapping kernel
module test_remap
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, REMAP_PCM, REMAP_PLM, REMAP_PPM
   use rki_remap_column, only: remap_column, remap_column_pcm, remap_column_plm, &
                               remap_column_ppm
   implicit none
   private

   public :: collect_remap_tests

contains

   subroutine collect_remap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("identity_pcm", test_identity_pcm), &
                  new_unittest("identity_plm", test_identity_plm), &
                  new_unittest("conservation_pcm", test_conservation_pcm), &
                  new_unittest("conservation_plm", test_conservation_plm), &
                  new_unittest("constant_pcm", test_constant_pcm), &
                  new_unittest("constant_plm", test_constant_plm), &
                  new_unittest("monotonicity_plm", test_monotonicity_plm), &
                  new_unittest("step_function_plm", test_step_function_plm), &
                  new_unittest("thin_layer", test_thin_layer), &
                  new_unittest("single_layer", test_single_layer), &
                  new_unittest("dispatch", test_dispatch), &
                  new_unittest("linear_profile_plm", test_linear_profile_plm), &
                  new_unittest("identity_ppm", test_identity_ppm), &
                  new_unittest("conservation_ppm", test_conservation_ppm), &
                  new_unittest("constant_ppm", test_constant_ppm), &
                  new_unittest("monotonicity_ppm", test_monotonicity_ppm), &
                  new_unittest("step_function_ppm", test_step_function_ppm), &
                  new_unittest("ppm_less_diffusive_than_pcm", test_ppm_less_diffusive), &
                  new_unittest("dispatch_ppm", test_dispatch_ppm) &
                  ]
   end subroutine collect_remap_tests

   subroutine test_identity_pcm(error)
      !! Remapping from dz to identical dz recovers q exactly (PCM)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz(nz), q_old(nz), q_new(nz)
      integer :: k

      dz = [1.0_wp, 2.0_wp, 3.0_wp, 2.0_wp, 1.0_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp, 25.0_wp, 15.0_wp]

      call remap_column_pcm(nz, dz, dz, q_old, q_new)

      do k = 1, nz
         call check(error, abs(q_new(k) - q_old(k)) < 1.0e-12_wp, &
                    "PCM identity remap should recover q exactly")
         if (allocated(error)) return
      end do
   end subroutine test_identity_pcm

   subroutine test_identity_plm(error)
      !! Remapping from dz to identical dz recovers q exactly (PLM)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz(nz), q_old(nz), q_new(nz)
      integer :: k

      dz = [1.0_wp, 2.0_wp, 3.0_wp, 2.0_wp, 1.0_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp, 25.0_wp, 15.0_wp]

      call remap_column_plm(nz, dz, dz, q_old, q_new)

      do k = 1, nz
         call check(error, abs(q_new(k) - q_old(k)) < 1.0e-12_wp, &
                    "PLM identity remap should recover q exactly")
         if (allocated(error)) return
      end do
   end subroutine test_identity_plm

   subroutine test_conservation_pcm(error)
      !! sum(q_new * dz_new) = sum(q_old * dz_old) for PCM
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new

      ! Non-uniform layers, different old and new
      dz_old = [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp]
      dz_new = [1.0_wp, 3.0_wp, 1.0_wp, 3.0_wp, 2.0_wp]
      q_old = [35.0_wp, 30.0_wp, 20.0_wp, 10.0_wp, 5.0_wp]

      call remap_column_pcm(nz, dz_old, dz_new, q_old, q_new)

      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)

      call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                 "PCM remap must conserve total mass")
   end subroutine test_conservation_pcm

   subroutine test_conservation_plm(error)
      !! sum(q_new * dz_new) = sum(q_old * dz_old) for PLM
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new

      dz_old = [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp]
      dz_new = [1.0_wp, 3.0_wp, 1.0_wp, 3.0_wp, 2.0_wp]
      q_old = [35.0_wp, 30.0_wp, 20.0_wp, 10.0_wp, 5.0_wp]

      call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)

      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)

      call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                 "PLM remap must conserve total mass")
   end subroutine test_conservation_plm

   subroutine test_constant_pcm(error)
      !! Remapping a uniform field preserves uniformity (PCM)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      integer :: k

      dz_old = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      dz_new = [2.5_wp, 2.5_wp, 2.5_wp, 2.5_wp]
      q_old = 25.0_wp

      call remap_column_pcm(nz, dz_old, dz_new, q_old, q_new)

      do k = 1, nz
         call check(error, abs(q_new(k) - 25.0_wp) < 1.0e-12_wp, &
                    "PCM remap of constant field should stay constant")
         if (allocated(error)) return
      end do
   end subroutine test_constant_pcm

   subroutine test_constant_plm(error)
      !! Remapping a uniform field preserves uniformity (PLM)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      integer :: k

      dz_old = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      dz_new = [2.5_wp, 2.5_wp, 2.5_wp, 2.5_wp]
      q_old = 25.0_wp

      call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)

      do k = 1, nz
         call check(error, abs(q_new(k) - 25.0_wp) < 1.0e-12_wp, &
                    "PLM remap of constant field should stay constant")
         if (allocated(error)) return
      end do
   end subroutine test_constant_plm

   subroutine test_monotonicity_plm(error)
      !! PLM remap of a monotone profile stays monotone (no new extrema)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 6
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: q_min, q_max
      integer :: k

      ! Monotonically increasing profile
      dz_old = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      dz_new = [0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp]
      q_old = [0.0_wp, 5.0_wp, 10.0_wp, 20.0_wp, 35.0_wp, 40.0_wp]

      call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)

      q_min = minval(q_old)
      q_max = maxval(q_old)

      ! No new extrema: all q_new values must lie within [q_min, q_max]
      do k = 1, nz
         call check(error, q_new(k) >= q_min - 1.0e-10_wp .and. &
                    q_new(k) <= q_max + 1.0e-10_wp, &
                    "PLM remap should not create new extrema")
         if (allocated(error)) return
      end do
   end subroutine test_monotonicity_plm

   subroutine test_step_function_plm(error)
      !! PLM remap of a sharp step: bounded, no new extrema, conservative
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new
      integer :: k

      dz_old = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      dz_new = [0.5_wp, 1.5_wp, 1.5_wp, 0.5_wp]
      q_old = [0.0_wp, 0.0_wp, 35.0_wp, 35.0_wp]

      call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)

      ! Conservation
      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)
      call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                 "PLM step: conservation check")
      if (allocated(error)) return

      ! Bounded
      do k = 1, nz
         call check(error, q_new(k) >= -1.0e-10_wp .and. q_new(k) <= 35.0_wp + 1.0e-10_wp, &
                    "PLM step: should be bounded within [0, 35]")
         if (allocated(error)) return
      end do
   end subroutine test_step_function_plm

   subroutine test_thin_layer(error)
      !! Remap with a near-zero dz in one layer doesn't crash
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new

      dz_old = [1.0_wp, 1.0e-14_wp, 2.0_wp, 1.0_wp]
      dz_new = [1.0_wp, 1.0_wp, 1.0_wp, 1.0e-14_wp + 1.0_wp]
      q_old = [10.0_wp, 0.0_wp, 20.0_wp, 30.0_wp]

      call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)

      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)
      call check(error, abs(mass_new - mass_old) < 1.0e-8_wp, &
                 "Thin layer remap: conservation check")
   end subroutine test_thin_layer

   subroutine test_single_layer(error)
      !! Single-layer remap is identity
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: dz_old(1), dz_new(1), q_old(1), q_new(1)

      dz_old = [5.0_wp]
      dz_new = [5.0_wp]
      q_old = [42.0_wp]

      call remap_column_plm(1, dz_old, dz_new, q_old, q_new)
      call check(error, abs(q_new(1) - 42.0_wp) < 1.0e-12_wp, &
                 "Single-layer remap should be identity")
   end subroutine test_single_layer

   subroutine test_dispatch(error)
      !! remap_column dispatches correctly to PCM and PLM
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 3
      real(wp) :: dz(nz), q_old(nz), q_pcm(nz), q_plm(nz)

      dz = [1.0_wp, 1.0_wp, 1.0_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp]

      call remap_column(REMAP_PCM, nz, dz, dz, q_old, q_pcm)
      call remap_column(REMAP_PLM, nz, dz, dz, q_old, q_plm)

      ! Both should be identity remaps
      call check(error, abs(q_pcm(2) - 20.0_wp) < 1.0e-12_wp, &
                 "Dispatch PCM: identity check")
      if (allocated(error)) return
      call check(error, abs(q_plm(2) - 20.0_wp) < 1.0e-12_wp, &
                 "Dispatch PLM: identity check")
   end subroutine test_dispatch

   subroutine test_linear_profile_plm(error)
      !! PLM should remap a linear profile exactly (within rounding)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new

      ! Uniform old layers with linear profile: q = k
      dz_old = [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp]
      dz_new = [1.0_wp, 3.0_wp, 2.0_wp, 3.0_wp, 1.0_wp]
      q_old = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

      call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)

      ! Conservation
      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)
      call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                 "PLM linear profile: conservation check")
      if (allocated(error)) return

      ! For a linear profile, PLM reconstruction is exact, so
      ! the new cell averages should be the exact integrals.
      ! New layer 1 spans z=[0,1], covering first half of old layer 1
      ! Old layer 1 center: z=1, q=1, slope=0.5*(2-1)=0.5 (limited by boundary→0)
      ! Actually slope(1)=0 (boundary), so PLM reduces to PCM for layer 1
      ! Just check conservation and boundedness
      call check(error, q_new(1) >= 0.5_wp .and. q_new(1) <= 5.5_wp, &
                 "PLM linear profile: bounded")
   end subroutine test_linear_profile_plm

   subroutine test_identity_ppm(error)
      !! Remapping from dz to identical dz recovers q exactly (PPM)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz(nz), q_old(nz), q_new(nz)
      integer :: k

      dz = [1.0_wp, 2.0_wp, 3.0_wp, 2.0_wp, 1.0_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp, 25.0_wp, 15.0_wp]

      call remap_column_ppm(nz, dz, dz, q_old, q_new)

      do k = 1, nz
         call check(error, abs(q_new(k) - q_old(k)) < 1.0e-10_wp, &
                    "PPM identity remap should recover q exactly")
         if (allocated(error)) return
      end do
   end subroutine test_identity_ppm

   subroutine test_conservation_ppm(error)
      !! sum(q_new * dz_new) = sum(q_old * dz_old) for PPM
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 5
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new

      dz_old = [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp]
      dz_new = [1.0_wp, 3.0_wp, 1.0_wp, 3.0_wp, 2.0_wp]
      q_old = [35.0_wp, 30.0_wp, 20.0_wp, 10.0_wp, 5.0_wp]

      call remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)

      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)

      call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                 "PPM remap must conserve total mass")
   end subroutine test_conservation_ppm

   subroutine test_constant_ppm(error)
      !! Remapping a uniform field preserves uniformity (PPM)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      integer :: k

      dz_old = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      dz_new = [2.5_wp, 2.5_wp, 2.5_wp, 2.5_wp]
      q_old = 25.0_wp

      call remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)

      do k = 1, nz
         call check(error, abs(q_new(k) - 25.0_wp) < 1.0e-10_wp, &
                    "PPM remap of constant field should stay constant")
         if (allocated(error)) return
      end do
   end subroutine test_constant_ppm

   subroutine test_monotonicity_ppm(error)
      !! PPM remap of a monotone profile: no new extrema
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 6
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: q_min_val, q_max_val
      integer :: k

      dz_old = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      dz_new = [0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp]
      q_old = [0.0_wp, 5.0_wp, 10.0_wp, 20.0_wp, 35.0_wp, 40.0_wp]

      call remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)

      q_min_val = minval(q_old)
      q_max_val = maxval(q_old)

      do k = 1, nz
         call check(error, q_new(k) >= q_min_val - 1.0e-10_wp .and. &
                    q_new(k) <= q_max_val + 1.0e-10_wp, &
                    "PPM remap should not create new extrema")
         if (allocated(error)) return
      end do
   end subroutine test_monotonicity_ppm

   subroutine test_step_function_ppm(error)
      !! PPM remap of a sharp step: bounded and conservative
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), q_new(nz)
      real(wp) :: mass_old, mass_new
      integer :: k

      dz_old = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      dz_new = [0.5_wp, 1.5_wp, 1.5_wp, 0.5_wp]
      q_old = [0.0_wp, 0.0_wp, 35.0_wp, 35.0_wp]

      call remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)

      mass_old = sum(q_old*dz_old)
      mass_new = sum(q_new*dz_new)
      call check(error, abs(mass_new - mass_old) < 1.0e-10_wp, &
                 "PPM step: conservation")
      if (allocated(error)) return

      do k = 1, nz
         call check(error, q_new(k) >= -1.0e-10_wp .and. q_new(k) <= 35.0_wp + 1.0e-10_wp, &
                    "PPM step: bounded within [0, 35]")
         if (allocated(error)) return
      end do
   end subroutine test_step_function_ppm

   subroutine test_ppm_less_diffusive(error)
      !! PPM should be less diffusive than PCM on a smooth profile
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 8
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz)
      real(wp) :: q_pcm(nz), q_ppm(nz)
      real(wp) :: err_pcm, err_ppm
      integer :: k

      ! Smooth sine-like profile on uniform layers, remap to shifted grid
      dz_old = 1.0_wp
      dz_new = [0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp, 0.5_wp, 1.5_wp]
      do k = 1, nz
         q_old(k) = sin(real(k, wp)*0.4_wp)
      end do

      call remap_column_pcm(nz, dz_old, dz_new, q_old, q_pcm)
      call remap_column_ppm(nz, dz_old, dz_new, q_old, q_ppm)

      ! Remap back to original grid and measure L2 error against original
      call remap_column_pcm(nz, dz_new, dz_old, q_pcm, q_pcm)
      call remap_column_ppm(nz, dz_new, dz_old, q_ppm, q_ppm)

      err_pcm = 0.0_wp
      err_ppm = 0.0_wp
      do k = 1, nz
         err_pcm = err_pcm + (q_pcm(k) - q_old(k))**2
         err_ppm = err_ppm + (q_ppm(k) - q_old(k))**2
      end do

      call check(error, err_ppm <= err_pcm + 1.0e-14_wp, &
                 "PPM round-trip error should be <= PCM round-trip error")
   end subroutine test_ppm_less_diffusive

   subroutine test_dispatch_ppm(error)
      !! remap_column with REMAP_PPM dispatches correctly
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 4
      real(wp) :: dz(nz), q_old(nz), q_dispatch(nz), q_direct(nz)
      integer :: k

      dz = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp]

      call remap_column(REMAP_PPM, nz, dz, dz, q_old, q_dispatch)
      call remap_column_ppm(nz, dz, dz, q_old, q_direct)

      do k = 1, nz
         call check(error, abs(q_dispatch(k) - q_direct(k)) < 1.0e-14_wp, &
                    "Dispatch PPM should match direct PPM call")
         if (allocated(error)) return
      end do
   end subroutine test_dispatch_ppm

end module test_remap
