!! Unit tests for vcoord_target_dz_column — the per-column thickness
!! generator that the existing `test_vcoord` suite did not exercise.
module test_vcoord_target
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_SIGMA, VCOORD_ZSIGMA, &
                            VCOORD_ZSTAR, VCOORD_ZSTAR_SIGMA
   use rki_vcoord, only: vcoord_target_dz_column
   implicit none
   private

   public :: collect_vcoord_target_tests

contains

   subroutine collect_vcoord_target_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("sigma_uniform", test_sigma_uniform), &
                  new_unittest("sigma_nonuniform", test_sigma_nonuniform), &
                  new_unittest("zsigma_shallow_pure_sigma", test_zsigma_shallow), &
                  new_unittest("zsigma_deep_blend", test_zsigma_deep), &
                  new_unittest("zsigma_zero_blend_width", test_zsigma_zero_blend_width), &
                  new_unittest("zsigma_sum_equals_H", test_zsigma_sum), &
                  new_unittest("zstar_proportional_scaling", test_zstar_proportional_scaling), &
                  new_unittest("zstar_conservation", test_zstar_conservation), &
                  new_unittest("zstar_relative_spacing_preserved", test_zstar_relative_spacing), &
                  new_unittest("zstar_sigma_shallow_pure_sigma", test_zstar_sigma_shallow), &
                  new_unittest("zstar_sigma_deep_pure_zstar", test_zstar_sigma_deep), &
                  new_unittest("zstar_sigma_blend_zone", test_zstar_sigma_blend), &
                  new_unittest("zstar_sigma_sum_equals_H", test_zstar_sigma_sum), &
                  new_unittest("unknown_falls_back_to_sigma", test_default_fallback) &
                  ]
   end subroutine collect_vcoord_target_tests

   subroutine test_sigma_uniform(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 10.0_wp
      integer :: k

      dsig = 1.0_wp/real(NZ, wp)
      do k = 0, NZ
         z_ref(k) = real(k, wp)
      end do

      call vcoord_target_dz_column(VCOORD_SIGMA, NZ, H, dsig, z_ref, &
                                   1000.0_wp, 100.0_wp, dz)

      do k = 1, NZ
         call check(error, abs(dz(k) - H/real(NZ, wp)) < 1.0e-12_wp, &
                    "uniform sigma should produce equal layers")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(dz) - H) < 1.0e-12_wp, "sigma layers sum to H")
   end subroutine test_sigma_uniform

   subroutine test_sigma_nonuniform(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 8.0_wp
      integer :: k

      dsig = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
      do k = 0, NZ
         z_ref(k) = real(k, wp)
      end do

      call vcoord_target_dz_column(VCOORD_SIGMA, NZ, H, dsig, z_ref, &
                                   1000.0_wp, 100.0_wp, dz)

      do k = 1, NZ
         call check(error, abs(dz(k) - dsig(k)*H) < 1.0e-12_wp, &
                    "sigma layer thickness = dsig * H")
         if (allocated(error)) return
      end do
   end subroutine test_sigma_nonuniform

   subroutine test_zsigma_shallow(error)
      !! When H <= depth_transition the zsigma path returns pure sigma layers.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 5.0_wp
      real(wp), parameter :: DEPTH_TRANS = 50.0_wp, BLEND = 25.0_wp
      integer :: k

      dsig = 0.25_wp
      do k = 0, NZ
         z_ref(k) = real(k, wp)*5.0_wp
      end do

      call vcoord_target_dz_column(VCOORD_ZSIGMA, NZ, H, dsig, z_ref, &
                                   DEPTH_TRANS, BLEND, dz)

      do k = 1, NZ
         call check(error, abs(dz(k) - dsig(k)*H) < 1.0e-12_wp, &
                    "shallow zsigma should equal pure sigma")
         if (allocated(error)) return
      end do
   end subroutine test_zsigma_shallow

   subroutine test_zsigma_deep(error)
      !! When H >> depth_transition + blend_width, alpha = 1 and dz comes
      !! from clipped z-levels (with deficit redistributed to bottom).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 1000.0_wp
      real(wp), parameter :: DEPTH_TRANS = 50.0_wp, BLEND = 25.0_wp
      integer :: k

      dsig = 0.25_wp
      ! Three z-level interfaces well within H, fourth is much deeper than H.
      ! Pre-deficit: 10 + 20 + 30 + (very large) = clipped to H beyond bottom.
      z_ref = [0.0_wp, 10.0_wp, 30.0_wp, 60.0_wp, 1.0e8_wp]

      call vcoord_target_dz_column(VCOORD_ZSIGMA, NZ, H, dsig, z_ref, &
                                   DEPTH_TRANS, BLEND, dz)

      ! Sum must still equal H (the deficit redistribution guarantees this).
      call check(error, abs(sum(dz) - H) < 1.0e-9_wp, "deep zsigma layers sum to H")
      if (allocated(error)) return

      ! Output is ROMS-ordered (k=1=bottom, k=nz=surface).  In the alpha=1
      ! limit, the SURFACE layer (k=NZ) is set to z_ref(1) - z_ref(0) = 10.
      call check(error, abs(dz(NZ) - 10.0_wp) < 1.0e-9_wp, &
                 "surface layer (k=NZ) hits z-level thickness in deep limit")
      if (allocated(error)) return

      ! Bottom layer (k=1) absorbs the deficit and is positive
      call check(error, dz(1) > 0.0_wp, "bottom layer (k=1) absorbs deficit")
   end subroutine test_zsigma_deep

   subroutine test_zsigma_zero_blend_width(error)
      !! blend_width = 0 must still be safe — the code uses an unguarded
      !! division otherwise, but the implementation special-cases it.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 3
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 200.0_wp

      dsig = 1.0_wp/real(NZ, wp)
      z_ref = [0.0_wp, 50.0_wp, 100.0_wp, 1.0e8_wp]

      call vcoord_target_dz_column(VCOORD_ZSIGMA, NZ, H, dsig, z_ref, &
                                   100.0_wp, 0.0_wp, dz)

      call check(error, abs(sum(dz) - H) < 1.0e-9_wp, &
                 "zsigma with blend_width=0 sum to H")
   end subroutine test_zsigma_zero_blend_width

   subroutine test_zsigma_sum(error)
      !! Sweep the blend region — sum(dz) must equal H at every depth.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp) :: H_test
      integer :: i

      dsig = 0.25_wp
      z_ref = [0.0_wp, 10.0_wp, 30.0_wp, 60.0_wp, 1.0e8_wp]

      do i = 1, 9
         H_test = real(i, wp)*25.0_wp   ! 25, 50, 75, ..., 225
         call vcoord_target_dz_column(VCOORD_ZSIGMA, NZ, H_test, dsig, z_ref, &
                                      50.0_wp, 25.0_wp, dz)
         call check(error, abs(sum(dz) - H_test) < 1.0e-9_wp, &
                    "sum(dz) must equal H at every depth")
         if (allocated(error)) return
      end do
   end subroutine test_zsigma_sum

   subroutine test_zstar_proportional_scaling(error)
      !! z*-lite: scaling H by an arbitrary factor must scale every dz(k)
      !! by exactly the same factor.  This is the defining property of
      !! z*'s SSH-tracking — relative layer spacing is preserved as the
      !! column depth changes.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz1(NZ), dz2(NZ), z_ref(0:NZ)
      real(wp), parameter :: H1 = 50.0_wp, H2 = 75.0_wp
      real(wp), parameter :: SCALE = H2/H1
      integer :: k

      dsig = 1.0_wp/real(NZ, wp)
      ! Non-uniform z_ref to make the test more discriminating
      z_ref = [0.0_wp, 5.0_wp, 15.0_wp, 35.0_wp, 75.0_wp]

      call vcoord_target_dz_column(VCOORD_ZSTAR, NZ, H1, dsig, z_ref, &
                                   100.0_wp, 50.0_wp, dz1)
      call vcoord_target_dz_column(VCOORD_ZSTAR, NZ, H2, dsig, z_ref, &
                                   100.0_wp, 50.0_wp, dz2)

      do k = 1, NZ
         call check(error, abs(dz2(k) - SCALE*dz1(k)) < 1.0e-12_wp, &
                    "zstar: dz scales linearly with H (SSH-tracking)")
         if (allocated(error)) return
      end do
   end subroutine test_zstar_proportional_scaling

   subroutine test_zstar_conservation(error)
      !! sum(dz) must equal H exactly across a sweep of H values.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ), H_test
      integer :: i

      dsig = 1.0_wp/real(NZ, wp)
      z_ref = [0.0_wp, 2.0_wp, 6.0_wp, 14.0_wp, 30.0_wp, 60.0_wp]

      do i = 1, 9
         H_test = real(i, wp)*15.0_wp   ! 15, 30, ..., 135 m
         call vcoord_target_dz_column(VCOORD_ZSTAR, NZ, H_test, dsig, z_ref, &
                                      100.0_wp, 50.0_wp, dz)
         call check(error, abs(sum(dz) - H_test) < 1.0e-12_wp, &
                    "zstar: sum(dz) must equal H at every depth")
         if (allocated(error)) return
      end do
   end subroutine test_zstar_conservation

   subroutine test_zstar_relative_spacing(error)
      !! Layer-thickness ratios dz(k1)/dz(k2) must equal the corresponding
      !! z_ref-difference ratios, independent of H — this confirms layers
      !! preserve their relative position in the column under SSH change.
      !! Also checks ROMS ordering: dz(1) is the BOTTOM layer (largest
      !! z_ref interval at the deepest end), dz(NZ) is the SURFACE layer
      !! (smallest z_ref interval at the shallow end).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 25.0_wp
      real(wp) :: ratio_dz, ratio_zref
      integer :: k

      dsig = 1.0_wp/real(NZ, wp)
      ! Concentrate layers near the surface: dz_ref(top) < dz_ref(bottom)
      z_ref = [0.0_wp, 2.0_wp, 6.0_wp, 14.0_wp, 30.0_wp]
      ! Reference layer thicknesses (top-down): 2, 4, 8, 16
      ! Under ROMS ordering (bottom-up): dz(1)=16*scale, dz(2)=8*scale,
      ! dz(3)=4*scale, dz(4)=2*scale  with scale = H/z_ref(NZ) = 25/30

      call vcoord_target_dz_column(VCOORD_ZSTAR, NZ, H, dsig, z_ref, &
                                   100.0_wp, 50.0_wp, dz)

      ! Check ROMS ordering: dz(1) > dz(NZ) (bottom layer is thickest here)
      call check(error, dz(1) > dz(NZ), &
                 "zstar ROMS ordering: dz(1)=bottom > dz(NZ)=surface for surface-concentrated z_ref")
      if (allocated(error)) return

      ! Compare ratio of any two layers to the corresponding z_ref ratio
      do k = 1, NZ - 1
         ratio_dz = dz(k)/dz(k + 1)
         ratio_zref = (z_ref(NZ - k + 1) - z_ref(NZ - k)) &
                      /(z_ref(NZ - k) - z_ref(NZ - k - 1))
         call check(error, abs(ratio_dz - ratio_zref) < 1.0e-12_wp, &
                    "zstar layer-thickness ratios match z_ref-difference ratios")
         if (allocated(error)) return
      end do
   end subroutine test_zstar_relative_spacing

   ! ---- VCOORD_ZSTAR_SIGMA hybrid ----

   subroutine test_zstar_sigma_shallow(error)
      !! For H ≤ depth_transition the hybrid reduces to pure sigma:
      !! equal dz across layers when dsig is uniform.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 5.0_wp
      real(wp), parameter :: DEPTH_T = 50.0_wp
      real(wp), parameter :: BLEND = 20.0_wp
      integer :: k

      dsig = 1.0_wp/real(NZ, wp)
      do k = 0, NZ
         z_ref(k) = real(k, wp)*14.0_wp  ! arbitrary nonuniform pattern
      end do

      call vcoord_target_dz_column(VCOORD_ZSTAR_SIGMA, NZ, H, dsig, z_ref, &
                                   DEPTH_T, BLEND, dz)

      do k = 1, NZ
         call check(error, abs(dz(k) - H/real(NZ, wp)) < 1.0e-12_wp, &
                    "shallow hybrid should reduce to uniform sigma")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(dz) - H) < 1.0e-12_wp, "sum(dz) = H (shallow)")
   end subroutine test_zstar_sigma_shallow

   subroutine test_zstar_sigma_deep(error)
      !! For H ≥ depth_transition + blend_width the hybrid reduces to
      !! pure z*-lite: dz(k) = (z_ref(nz-k+1) - z_ref(nz-k)) * H / z_ref(nz).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ), dz_expected(NZ)
      real(wp), parameter :: H = 200.0_wp
      real(wp), parameter :: DEPTH_T = 50.0_wp
      real(wp), parameter :: BLEND = 20.0_wp
      real(wp) :: dz_z
      integer :: k

      dsig = 1.0_wp/real(NZ, wp)
      ! Surface-concentrated reference pattern (z_ref(0)=0, z_ref(NZ)=20)
      z_ref = [0.0_wp, 2.0_wp, 6.0_wp, 12.0_wp, 20.0_wp]

      call vcoord_target_dz_column(VCOORD_ZSTAR_SIGMA, NZ, H, dsig, z_ref, &
                                   DEPTH_T, BLEND, dz)

      do k = 1, NZ
         dz_z = (z_ref(NZ - k + 1) - z_ref(NZ - k))*H/z_ref(NZ)
         dz_expected(k) = dz_z
         call check(error, abs(dz(k) - dz_z) < 1.0e-10_wp, &
                    "deep hybrid should reduce to pure z*-lite")
         if (allocated(error)) return
      end do
      call check(error, abs(sum(dz) - H) < 1.0e-10_wp, "sum(dz) = H (deep)")
   end subroutine test_zstar_sigma_deep

   subroutine test_zstar_sigma_blend(error)
      !! In the blend zone (depth_transition < H < depth_transition+blend_width)
      !! the result should fall strictly between the pure-sigma and
      !! pure-z*-lite values, layer-by-layer.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz_hybrid(NZ), dz_sigma(NZ), dz_zstar(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 60.0_wp     ! between 50 and 70
      real(wp), parameter :: DEPTH_T = 50.0_wp
      real(wp), parameter :: BLEND = 20.0_wp
      integer :: k

      dsig = 1.0_wp/real(NZ, wp)
      z_ref = [0.0_wp, 2.0_wp, 6.0_wp, 12.0_wp, 20.0_wp]

      call vcoord_target_dz_column(VCOORD_SIGMA, NZ, H, dsig, z_ref, &
                                   DEPTH_T, BLEND, dz_sigma)
      call vcoord_target_dz_column(VCOORD_ZSTAR, NZ, H, dsig, z_ref, &
                                   DEPTH_T, BLEND, dz_zstar)
      call vcoord_target_dz_column(VCOORD_ZSTAR_SIGMA, NZ, H, dsig, z_ref, &
                                   DEPTH_T, BLEND, dz_hybrid)

      ! Layer thicknesses where sigma != z*-lite: hybrid is strictly between
      do k = 1, NZ
         if (abs(dz_sigma(k) - dz_zstar(k)) > 1.0e-10_wp) then
            call check(error, &
                       (dz_hybrid(k) - dz_sigma(k))*(dz_hybrid(k) - dz_zstar(k)) < 0.0_wp, &
                       "blend value lies strictly between sigma and z*-lite")
            if (allocated(error)) return
         end if
      end do
   end subroutine test_zstar_sigma_blend

   subroutine test_zstar_sigma_sum(error)
      !! sum(dz) = H across a sweep of H values that cross the blend zone.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 4
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: DEPTH_T = 50.0_wp
      real(wp), parameter :: BLEND = 20.0_wp
      real(wp), parameter :: H_test(6) = [5.0_wp, 30.0_wp, 55.0_wp, 65.0_wp, 100.0_wp, 500.0_wp]
      integer :: i

      dsig = 1.0_wp/real(NZ, wp)
      z_ref = [0.0_wp, 2.0_wp, 6.0_wp, 12.0_wp, 20.0_wp]

      do i = 1, size(H_test)
         call vcoord_target_dz_column(VCOORD_ZSTAR_SIGMA, NZ, H_test(i), &
                                      dsig, z_ref, DEPTH_T, BLEND, dz)
         call check(error, abs(sum(dz) - H_test(i)) < 1.0e-10_wp, &
                    "hybrid: sum(dz) = H across blend zone")
         if (allocated(error)) return
      end do
   end subroutine test_zstar_sigma_sum

   subroutine test_default_fallback(error)
      !! Unknown coord_type falls through to the sigma branch.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 3
      real(wp) :: dsig(NZ), dz(NZ), z_ref(0:NZ)
      real(wp), parameter :: H = 10.0_wp
      integer, parameter :: BOGUS = 9999
      integer :: k

      dsig = [0.2_wp, 0.3_wp, 0.5_wp]
      do k = 0, NZ
         z_ref(k) = real(k, wp)
      end do

      call vcoord_target_dz_column(BOGUS, NZ, H, dsig, z_ref, &
                                   100.0_wp, 25.0_wp, dz)

      do k = 1, NZ
         call check(error, abs(dz(k) - dsig(k)*H) < 1.0e-12_wp, &
                    "unknown coord_type should behave like sigma")
         if (allocated(error)) return
      end do
   end subroutine test_default_fallback

end module test_vcoord_target
