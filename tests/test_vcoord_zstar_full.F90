!! Unit tests for VCOORD_ZSTAR_FULL — `zstar_full_build_column` and
!! `vcoord_target_dz_column_zstar_full`.
module test_vcoord_zstar_full
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_vcoord, only: zstar_full_build_column, &
                         vcoord_target_dz_column_zstar_full, &
                         STRETCH_LOG, STRETCH_UNIFORM
   implicit none
   private

   public :: collect_vcoord_zstar_full_tests

contains

   subroutine collect_vcoord_zstar_full_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("builder_uniform_monotonic", &
                               test_builder_uniform_monotonic), &
                  new_unittest("builder_log_surface_concentration", &
                               test_builder_log_surface), &
                  new_unittest("builder_anchors_bed", test_builder_anchors_bed), &
                  new_unittest("builder_dry_cell_degenerate_nonneg", &
                               test_builder_dry_cell_nonneg), &
                  new_unittest("step_sum_equals_H_no_eta", &
                               test_step_sum_no_eta), &
                  new_unittest("step_positive_eta_in_surface", &
                               test_step_positive_eta), &
                  new_unittest("step_negative_eta_creates_vanishing", &
                               test_step_negative_eta_vanishing), &
                  new_unittest("step_sum_equals_H_under_eta", &
                               test_step_sum_under_eta) &
                  ]
   end subroutine collect_vcoord_zstar_full_tests

   ! ---- zstar_full_build_column ----

   subroutine test_builder_uniform_monotonic(error)
      !! Uniform stretching: z_ref_col(0)=0, z_ref_col(nz)=h_bed,
      !! strictly monotonic, equal increments (within the fine zone).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp) :: z_ref_col(0:NZ)
      real(wp), parameter :: H_BED = 100.0_wp
      integer :: k

      call zstar_full_build_column(NZ, H_BED, 1.0_wp, 0, STRETCH_UNIFORM, &
                                   z_ref_col)

      call check(error, abs(z_ref_col(0)) < 1.0e-12_wp, &
                 "z_ref_col(0) should be 0")
      if (allocated(error)) return
      call check(error, abs(z_ref_col(NZ) - H_BED) < 1.0e-12_wp, &
                 "z_ref_col(nz) should equal h_bed")
      if (allocated(error)) return
      do k = 1, NZ
         call check(error, z_ref_col(k) > z_ref_col(k - 1) - 1.0e-12_wp, &
                    "z_ref_col should be monotonically increasing")
         if (allocated(error)) return
      end do
   end subroutine test_builder_uniform_monotonic

   subroutine test_builder_log_surface(error)
      !! Log stretching: surface layers should be thinner than
      !! mid-column layers.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp) :: z_ref_col(0:NZ)
      real(wp), parameter :: H_BED = 500.0_wp
      real(wp), parameter :: H_SURF_TARGET = 1.0_wp
      real(wp) :: dz_surf, dz_mid

      call zstar_full_build_column(NZ, H_BED, H_SURF_TARGET, 4, &
                                   STRETCH_LOG, z_ref_col)

      dz_surf = z_ref_col(1) - z_ref_col(0)
      dz_mid = z_ref_col(NZ) - z_ref_col(NZ - 1)

      call check(error, dz_surf < dz_mid, &
                 "log stretching should give thinner surface than bottom layers")
      if (allocated(error)) return
      call check(error, abs(z_ref_col(NZ) - H_BED) < 1.0e-9_wp, &
                 "log builder must anchor bottom at h_bed")
   end subroutine test_builder_log_surface

   subroutine test_builder_anchors_bed(error)
      !! For a range of bed depths, the deepest interface should
      !! land exactly at h_bed, the surface at 0.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 8
      real(wp) :: z_ref_col(0:NZ)
      real(wp), parameter :: H_BEDS(4) = [5.0_wp, 50.0_wp, 250.0_wp, 1000.0_wp]
      integer :: i

      do i = 1, size(H_BEDS)
         call zstar_full_build_column(NZ, H_BEDS(i), 0.5_wp, 0, &
                                      STRETCH_LOG, z_ref_col)
         call check(error, abs(z_ref_col(0)) < 1.0e-12_wp, &
                    "z_ref_col(0) = 0 across bed depths")
         if (allocated(error)) return
         call check(error, abs(z_ref_col(NZ) - H_BEDS(i)) < 1.0e-9_wp, &
                    "z_ref_col(NZ) = h_bed across bed depths")
         if (allocated(error)) return
      end do
   end subroutine test_builder_anchors_bed

   subroutine test_builder_dry_cell_nonneg(error)
      !! REGRESSION (Lakes Entrance 2026-05-11): land cells in the
      !! structured solver have bathymetric `b > 0` (above sea level),
      !! so the production wrapper `ml_build_vcoord_zstar_full` clamps
      !! `-b` to a non-negative value before calling this kernel.  But
      !! the kernel itself must ALSO be robust to `h_bed = 0` — the
      !! degenerate / dry-land case — by producing a column of zeros,
      !! NOT negative interfaces.  A previous bug propagated negative
      !! `z_ref_col` values from a `h_bed < 0` input, which then surfaced
      !! as negative `dz_new` in the per-step builder and crashed the
      !! solver via NaN CFL.  This test asserts non-negativity of all
      !! interfaces and monotonicity.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 6
      real(wp) :: z_ref_col(0:NZ)
      integer :: k

      ! h_bed = 0 (dry land sentinel): all interfaces should be 0.
      call zstar_full_build_column(NZ, 0.0_wp, 1.0_wp, 2, STRETCH_UNIFORM, &
                                   z_ref_col)
      do k = 0, NZ
         call check(error, z_ref_col(k) >= 0.0_wp, &
                    "z_ref_col(k) must be >= 0 for h_bed = 0 (dry land)")
         if (allocated(error)) return
      end do
      do k = 1, NZ
         call check(error, z_ref_col(k) >= z_ref_col(k - 1), &
                    "z_ref_col must be monotonic non-decreasing")
         if (allocated(error)) return
      end do

      ! Defensive: even if a caller passes a *negative* h_bed (e.g. by
      ! forgetting to negate `b` on the structured side), the kernel must
      ! still emit non-negative interfaces.  Otherwise the per-step
      ! builder produces negative dz_new which the solver can't handle.
      call zstar_full_build_column(NZ, -20.0_wp, 1.0_wp, 2, STRETCH_UNIFORM, &
                                   z_ref_col)
      do k = 0, NZ
         call check(error, z_ref_col(k) >= 0.0_wp, &
                    "z_ref_col(k) must be >= 0 even for negative h_bed input")
         if (allocated(error)) return
      end do
   end subroutine test_builder_dry_cell_nonneg

   ! ---- vcoord_target_dz_column_zstar_full ----

   subroutine test_step_sum_no_eta(error)
      !! eta = 0 (H = h_bed): each ROMS-ordered dz should equal the
      !! reference (top-down) thickness with no SSH offset.
      !! sum(dz) = H exactly.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: z_ref_col(0:NZ), dz(NZ)
      real(wp), parameter :: H_BED = 50.0_wp
      real(wp), parameter :: H_MIN = 1.0e-4_wp
      integer :: k

      call zstar_full_build_column(NZ, H_BED, 0.0_wp, 0, STRETCH_UNIFORM, &
                                   z_ref_col)

      call vcoord_target_dz_column_zstar_full(NZ, H_BED, z_ref_col, &
                                              H_MIN, dz)

      call check(error, abs(sum(dz) - H_BED) < 1.0e-9_wp, &
                 "sum(dz) = H when eta = 0")
      if (allocated(error)) return
      ! All layers should have equal thickness (uniform reference,
      ! uniform stretching, no SSH offset)
      do k = 1, NZ
         call check(error, abs(dz(k) - H_BED/real(NZ, wp)) < 1.0e-9_wp, &
                    "uniform reference + no eta -> equal dz")
         if (allocated(error)) return
      end do
   end subroutine test_step_sum_no_eta

   subroutine test_step_positive_eta(error)
      !! eta > 0 (H > h_bed): surface layer absorbs the entire offset.
      !! sum(dz) = H exactly; non-surface layers unchanged.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 5
      real(wp) :: z_ref_col(0:NZ), dz(NZ), dz_ref(NZ)
      real(wp), parameter :: H_BED = 50.0_wp
      real(wp), parameter :: ETA = 3.0_wp
      real(wp), parameter :: H_MIN = 1.0e-4_wp
      integer :: k

      call zstar_full_build_column(NZ, H_BED, 0.0_wp, 0, STRETCH_UNIFORM, &
                                   z_ref_col)
      call vcoord_target_dz_column_zstar_full(NZ, H_BED, z_ref_col, H_MIN, dz_ref)
      call vcoord_target_dz_column_zstar_full(NZ, H_BED + ETA, z_ref_col, H_MIN, dz)

      call check(error, abs(sum(dz) - (H_BED + ETA)) < 1.0e-9_wp, &
                 "sum(dz) = H when eta > 0")
      if (allocated(error)) return
      ! Subsurface layers (ROMS k=1..nz-1) should match the reference
      do k = 1, NZ - 1
         call check(error, abs(dz(k) - dz_ref(k)) < 1.0e-9_wp, &
                    "subsurface layers unchanged when eta > 0")
         if (allocated(error)) return
      end do
      ! Surface layer (ROMS k=NZ) absorbs the offset
      call check(error, abs(dz(NZ) - dz_ref(NZ) - ETA) < 1.0e-9_wp, &
                 "surface layer absorbs eta")
   end subroutine test_step_positive_eta

   subroutine test_step_negative_eta_vanishing(error)
      !! Column much shallower than reference: deep layers should
      !! collapse to h_min ("vanishing"); upper layers keep reference
      !! thickness.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp) :: z_ref_col(0:NZ), dz(NZ)
      real(wp), parameter :: H_BED = 100.0_wp
      real(wp), parameter :: H_ACTUAL = 30.0_wp  ! shallow
      real(wp), parameter :: H_MIN = 1.0e-4_wp
      integer :: k, n_vanish

      call zstar_full_build_column(NZ, H_BED, 0.0_wp, 0, STRETCH_UNIFORM, &
                                   z_ref_col)
      call vcoord_target_dz_column_zstar_full(NZ, H_ACTUAL, z_ref_col, &
                                              H_MIN, dz)

      ! ROMS order: dz(1) is bottom.  Bed reference is at 100m,
      ! actual column is 30m, so layers below z=30m (which are the
      ! ROMS-bottom layers) should vanish to h_min.
      n_vanish = 0
      do k = 1, NZ
         if (abs(dz(k) - H_MIN) < 10.0_wp*H_MIN) n_vanish = n_vanish + 1
      end do
      call check(error, n_vanish >= 5, &
                 "expected several vanishing layers when H is 30% of h_bed")
      if (allocated(error)) return

      ! At least one surface-side layer should retain its reference thickness
      call check(error, dz(NZ) > 5.0_wp*H_MIN, &
                 "surface layer should be substantially thicker than h_min")
   end subroutine test_step_negative_eta_vanishing

   subroutine test_step_sum_under_eta(error)
      !! sum(dz) = H exactly even when vanishing layers are present
      !! (mass conservation for downstream remap).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: NZ = 10
      real(wp) :: z_ref_col(0:NZ), dz(NZ)
      real(wp), parameter :: H_BED = 100.0_wp
      real(wp), parameter :: H_MIN = 1.0e-4_wp
      real(wp), parameter :: H_test(4) = [80.0_wp, 50.0_wp, 20.0_wp, 5.0_wp]
      integer :: i

      call zstar_full_build_column(NZ, H_BED, 0.0_wp, 0, STRETCH_UNIFORM, &
                                   z_ref_col)
      do i = 1, size(H_test)
         call vcoord_target_dz_column_zstar_full(NZ, H_test(i), z_ref_col, &
                                                 H_MIN, dz)
         call check(error, abs(sum(dz) - H_test(i)) < 1.0e-9_wp, &
                    "sum(dz) = H under negative eta with vanishing layers")
         if (allocated(error)) return
      end do
   end subroutine test_step_sum_under_eta

end module test_vcoord_zstar_full
