!! Unit tests for the vertical coordinate type
module test_vcoord
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_SIGMA, VCOORD_ZSIGMA, &
                            VCOORD_ZSTAR, VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            REMAP_PLM, REMAP_PPM, REMAP_PCM
   use rki_vcoord, only: vcoord_t, parse_vcoord_type, parse_remap_method
   implicit none
   private

   public :: collect_vcoord_tests

contains

   subroutine collect_vcoord_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_sigma", test_init_sigma), &
                  new_unittest("init_allocates", test_init_allocates), &
                  new_unittest("needs_remap_sigma", test_needs_remap_sigma), &
                  new_unittest("needs_remap_zsigma", test_needs_remap_zsigma), &
                  new_unittest("needs_remap_zsigma", test_needs_remap_zsigma), &
                  new_unittest("cleanup_idempotent", test_cleanup_idempotent), &
                  new_unittest("cleanup_before_init", test_cleanup_before_init), &
                  new_unittest("reinit", test_reinit), &
                  new_unittest("dsig_sum_invariant", test_dsig_sum_invariant), &
                  new_unittest("parse_vcoord_type", test_parse_vcoord_type), &
                  new_unittest("parse_remap_method", test_parse_remap_method) &
                  ]
   end subroutine collect_vcoord_tests

   subroutine test_init_sigma(error)
      !! vcoord_init with VCOORD_SIGMA allocates dsig_target(nz) = 1/nz
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc
      integer :: nz

      do nz = 1, 20
         call vc%init(nz, VCOORD_SIGMA, REMAP_PLM)
         call check(error, size(vc%dsig_target) == nz, &
                    "dsig_target should have nz elements")
         if (allocated(error)) return
         call check(error, abs(vc%dsig_target(1) - 1.0_wp/real(nz, wp)) < epsilon(1.0_wp), &
                    "dsig_target(1) should be 1/nz")
         if (allocated(error)) return
         call vc%cleanup()
      end do
   end subroutine test_init_sigma

   subroutine test_init_allocates(error)
      !! After init, dsig_target is allocated; after cleanup, it is not
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc

      call vc%init(5, VCOORD_SIGMA, REMAP_PLM)
      call check(error, allocated(vc%dsig_target), &
                 "dsig_target should be allocated after init")
      if (allocated(error)) return
      call check(error, vc%nz == 5, "nz should be 5")
      if (allocated(error)) return

      call vc%cleanup()
      call check(error,.not. allocated(vc%dsig_target), &
                 "dsig_target should be deallocated after cleanup")
      if (allocated(error)) return
      call check(error, vc%nz == 0, "nz should be 0 after cleanup")
   end subroutine test_init_allocates

   subroutine test_needs_remap_sigma(error)
      !! needs_remap returns .false. for VCOORD_SIGMA
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc

      call vc%init(5, VCOORD_SIGMA, REMAP_PLM)
      call check(error,.not. vc%needs_remap(), &
                 "sigma should not need remapping")
      call vc%cleanup()
   end subroutine test_needs_remap_sigma

   subroutine test_needs_remap_zsigma(error)
      !! needs_remap returns .true. for z-sigma (anything not pure sigma)
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc

      call vc%init(5, VCOORD_ZSIGMA, REMAP_PLM)
      call check(error, vc%needs_remap(), &
                 "zsigma should need remapping")
      call vc%cleanup()
   end subroutine test_needs_remap_zsigma

   subroutine test_cleanup_idempotent(error)
      !! Calling cleanup twice does not crash
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc

      call vc%init(5, VCOORD_SIGMA, REMAP_PLM)
      call vc%cleanup()
      call vc%cleanup()
      call check(error,.not. allocated(vc%dsig_target), &
                 "dsig_target should remain deallocated")
   end subroutine test_cleanup_idempotent

   subroutine test_cleanup_before_init(error)
      !! Calling cleanup on a fresh (uninitialised) vcoord_t does not crash
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc

      call vc%cleanup()
      call check(error,.not. allocated(vc%dsig_target), &
                 "dsig_target should not be allocated")
   end subroutine test_cleanup_before_init

   subroutine test_reinit(error)
      !! Calling init a second time with different nz reallocates correctly
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc

      call vc%init(3, VCOORD_SIGMA, REMAP_PLM)
      call check(error, size(vc%dsig_target) == 3, &
                 "first init: size should be 3")
      if (allocated(error)) return

      call vc%init(10, VCOORD_SIGMA, REMAP_PPM)
      call check(error, size(vc%dsig_target) == 10, &
                 "second init: size should be 10")
      if (allocated(error)) return
      call check(error, vc%remap_method == REMAP_PPM, &
                 "remap_method should be PPM after reinit")
      if (allocated(error)) return
      call vc%cleanup()
   end subroutine test_reinit

   subroutine test_dsig_sum_invariant(error)
      !! sum(dsig_target) = 1.0 for sigma coordinate
      type(error_type), allocatable, intent(out) :: error
      type(vcoord_t) :: vc
      integer :: nz
      real(wp) :: s

      do nz = 1, 20
         call vc%init(nz, VCOORD_SIGMA, REMAP_PLM)
         s = sum(vc%dsig_target)
         call check(error, abs(s - 1.0_wp) < 100.0_wp*epsilon(1.0_wp), &
                    "sum(dsig_target) should be 1.0")
         if (allocated(error)) return
         call vc%cleanup()
      end do
   end subroutine test_dsig_sum_invariant

   subroutine test_parse_vcoord_type(error)
      !! parse_vcoord_type returns correct enum values
      type(error_type), allocatable, intent(out) :: error

      call check(error, parse_vcoord_type("sigma") == VCOORD_SIGMA, &
                 "'sigma' should parse to VCOORD_SIGMA")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zsigma") == VCOORD_ZSIGMA, &
                 "'zsigma' should parse to VCOORD_ZSIGMA")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z-sigma") == VCOORD_ZSIGMA, &
                 "'z-sigma' should parse to VCOORD_ZSIGMA")
      if (allocated(error)) return
      ! "scoord" / VCOORD_SCOORD removed (D7); "scoord" now falls back to VCOORD_SIGMA
      call check(error, parse_vcoord_type("scoord") == VCOORD_SIGMA, &
                 "'scoord' should now fall back to VCOORD_SIGMA (D7 removal)")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar") == VCOORD_ZSTAR, &
                 "'zstar' should parse to VCOORD_ZSTAR")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar_full") == VCOORD_ZSTAR_FULL, &
                 "'zstar_full' should parse to VCOORD_ZSTAR_FULL")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z-star-full") == VCOORD_ZSTAR_FULL, &
                 "'z-star-full' should parse to VCOORD_ZSTAR_FULL")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("zstar_sigma") == VCOORD_ZSTAR_SIGMA, &
                 "'zstar_sigma' should parse to VCOORD_ZSTAR_SIGMA")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("z-star-sigma") == VCOORD_ZSTAR_SIGMA, &
                 "'z-star-sigma' should parse to VCOORD_ZSTAR_SIGMA")
      if (allocated(error)) return
      call check(error, parse_vcoord_type("unknown") == VCOORD_SIGMA, &
                 "unknown string should default to VCOORD_SIGMA")
   end subroutine test_parse_vcoord_type

   subroutine test_parse_remap_method(error)
      !! parse_remap_method returns correct enum values
      type(error_type), allocatable, intent(out) :: error

      call check(error, parse_remap_method("pcm") == REMAP_PCM, &
                 "'pcm' should parse to REMAP_PCM")
      if (allocated(error)) return
      call check(error, parse_remap_method("plm") == REMAP_PLM, &
                 "'plm' should parse to REMAP_PLM")
      if (allocated(error)) return
      call check(error, parse_remap_method("ppm") == REMAP_PPM, &
                 "'ppm' should parse to REMAP_PPM")
      if (allocated(error)) return
      call check(error, parse_remap_method("unknown") == REMAP_PLM, &
                 "unknown string should default to REMAP_PLM")
   end subroutine test_parse_remap_method

end module test_vcoord
