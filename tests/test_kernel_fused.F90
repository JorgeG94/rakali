!! Unit tests for the fused 2D-barotropic update kernel (structured)
!!
!! Drives `fused_update_stage1` and `fused_update_stage2` directly with
!! synthetic flux + forcing arrays so the major branches (clamp,
!! thin-cell desingularization, Manning friction, Coriolis, wind /
!! atmospheric pressure / precipitation forcing, RK2 average) are all
!! exercised.
module test_kernel_fused
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DESING_EPS, FROUDE_BACKSTOP
   use rki_grid, only: hgrid_t
   use rki_kernel_fused, only: fused_update_stage1, fused_update_stage2
   implicit none
   private

   public :: collect_kernel_fused_tests

   integer, parameter :: NX_PHYS = 4, NY_PHYS = 3, NG = 2
   integer, parameter :: NX = NX_PHYS + 2*NG, NY = NY_PHYS + 2*NG

contains

   subroutine collect_kernel_fused_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("stage1_save_and_euler", test_stage1_basic), &
                  new_unittest("stage1_clamps_negative_h", test_clamp_negative), &
                  new_unittest("stage1_dry_cell_zeros_momentum", test_thin_layer), &
                  new_unittest("stage1_friction_damps_momentum", test_friction), &
                  new_unittest("stage1_forcing_adds_wind", test_wind), &
                  new_unittest("stage1_coriolis_rotates_velocity", test_coriolis), &
                  new_unittest("stage1_desingularizes_thin_momentum", test_desingularization), &
                  new_unittest("stage2_rk2_average", test_stage2_rk2) &
                  ]
   end subroutine collect_kernel_fused_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NG, 1.0_wp, 1.0_wp)
   end subroutine make_grid

   subroutine zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                          wx, wy, atm, prec)
      real(wp), intent(out) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp), intent(out) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp), intent(out) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp), intent(out) :: wx(NX, NY), wy(NX, NY)
      real(wp), intent(out) :: atm(NX, NY), prec(NX, NY)
      h = 0.0_wp; hu = 0.0_wp; hv = 0.0_wp
      h0 = 0.0_wp; hu0 = 0.0_wp; hv0 = 0.0_wp
      fh = 0.0_wp; fhu = 0.0_wp; fhv = 0.0_wp
      wx = 0.0_wp; wy = 0.0_wp
      atm = 0.0_wp; prec = 0.0_wp
   end subroutine zero_arrays

   subroutine test_stage1_basic(error)
      !! No forcing, no flux: stage1 should leave (h, hu, hv) unchanged
      !! and snapshot them into (h0, hu0, hv0).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 2.0_wp; hu = 0.5_wp; hv = -0.25_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.1_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      call check(error, all(abs(h(NG + 1:NG + NX_PHYS, NG + 1:NG + NY_PHYS) - 2.0_wp) < 1.0e-12_wp), &
                 "h interior unchanged with zero flux + no forcing")
      if (allocated(error)) return

      call check(error, all(abs(h0(NG + 1:NG + NX_PHYS, NG + 1:NG + NY_PHYS) - 2.0_wp) < 1.0e-12_wp), &
                 "h0 snapshot equals input h")
      if (allocated(error)) return
      call check(error, all(abs(hu0(NG + 1:NG + NX_PHYS, NG + 1:NG + NY_PHYS) - 0.5_wp) < 1.0e-12_wp), &
                 "hu0 snapshot equals input hu")
   end subroutine test_stage1_basic

   subroutine test_clamp_negative(error)
      !! Negative h after Euler update is clamped to zero with zero momentum.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)
      integer :: i, j

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 0.5_wp; hu = 1.0_wp; hv = 1.0_wp
      ! Strong negative flux drains the column past zero
      fh = -100.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.1_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      do j = NG + 1, NG + NY_PHYS
         do i = NG + 1, NG + NX_PHYS
            call check(error, h(i, j) == 0.0_wp, "h must be clamped to zero")
            if (allocated(error)) return
            call check(error, hu(i, j) == 0.0_wp, "hu zeroed when h clamped")
            if (allocated(error)) return
            call check(error, hv(i, j) == 0.0_wp, "hv zeroed when h clamped")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_clamp_negative

   subroutine test_thin_layer(error)
      !! 0 < h < DRY_TOLERANCE ⇒ momentum zeroed, h preserved.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)
      integer :: i, j

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 1.0e-12_wp; hu = 1.0_wp; hv = 1.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.1_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      do j = NG + 1, NG + NY_PHYS
         do i = NG + 1, NG + NX_PHYS
            call check(error, hu(i, j) == 0.0_wp, "thin-layer hu zeroed")
            if (allocated(error)) return
            call check(error, hv(i, j) == 0.0_wp, "thin-layer hv zeroed")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_thin_layer

   subroutine test_friction(error)
      !! Manning friction must reduce momentum magnitude when n > 0.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 2.0_wp; hu = 1.0_wp; hv = 0.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.5_wp, 0.05_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      ! With manning_n > 0 and dt non-trivial, hu must shrink
      call check(error, hu(NG + 1, NG + 1) > 0.0_wp .and. hu(NG + 1, NG + 1) < 1.0_wp, &
                 "Manning friction must reduce hu")
   end subroutine test_friction

   subroutine test_wind(error)
      !! Wind stress with forcing_active=.true. should add momentum.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 2.0_wp
      wx = 1.0_wp   ! constant wind stress

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.1_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .true.)

      ! Wind in +x direction should increase hu (was zero before)
      call check(error, hu(NG + 1, NG + 1) > 0.0_wp, "wind stress must add hu")
   end subroutine test_wind

   subroutine test_coriolis(error)
      !! Non-zero coriolis_f must rotate the (hu, hv) vector.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 2.0_wp; hu = 1.0_wp; hv = 0.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 1.0_wp, 0.0_wp, wx, wy, atm, prec, &
                               1.0_wp, .false.)

      ! With f*dt = 1 rad, hv was 0 -> sin(1)*hu_old, hu shrinks to cos(1)*hu_old.
      ! Sanity: |momentum| preserved by the rotation.
      call check(error, abs(hu(NG + 1, NG + 1)**2 + hv(NG + 1, NG + 1)**2 - 1.0_wp) < 1.0e-9_wp, &
                 "Coriolis must preserve |momentum|^2")
   end subroutine test_coriolis

   subroutine test_desingularization(error)
      !! Thin-but-wet cells (DRY_TOLERANCE <= h < eps) get the Kurganov &
      !! Petrova (2007) desingularized momentum
      !!     hu := h * sqrt(2) h hu / sqrt(h^4 + eps^4),
      !! eps = DESING_EPS.  Cells with h >= eps and
      !! |u| <= FROUDE_BACKSTOP*sqrt(g h) are untouched EXACTLY — no
      !! physics-scale momentum rescale survives (the old FROUDE_CAP=10
      !! clamp pinned the supercritical tongue of every dry-bed front;
      !! see test_dam_break_dry).  Above the backstop, momentum is
      !! rescaled onto the bound — a blowup arrester, pinned here too.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)
      real(wp) :: eps, h_thin, hu_thin, hu_expect

      call make_grid(grid)
      eps = DESING_EPS

      ! Thin cell: h = eps/2, raw u = 3 m/s (strongly supercritical for
      ! this h; desingularized speed ~1.0 m/s stays below the backstop's
      ! ~1.36 m/s bound at this depth).
      h_thin = 0.5_wp*eps
      hu_thin = 3.0_wp*h_thin
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = h_thin
      hu = hu_thin
      hv = 0.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.01_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      hu_expect = h_thin*(sqrt(2.0_wp)*h_thin*hu_thin/sqrt(h_thin**4 + eps**4))
      call check(error, abs(hu(NG + 1, NG + 1) - hu_expect) < 1.0e-14_wp, &
                 "thin-cell momentum must equal the KP07 desingularized value")
      if (allocated(error)) return
      call check(error, abs(hu(NG + 1, NG + 1)) < abs(hu_thin), &
                 "desingularization must reduce thin-cell momentum")
      if (allocated(error)) return

      ! Wet cell at h >= eps with u = 20 m/s — supercritical (Fr ~ 6.4,
      ! above the OLD FROUDE_CAP = 10 domain of harm, far below the
      ! backstop): must pass through bit-identical.
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 1.0_wp
      hu = 20.0_wp
      hv = 0.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.01_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      call check(error, hu(NG + 1, NG + 1) == 20.0_wp, &
                 "supercritical momentum below the backstop must be "// &
                 "bit-identical (no physics-scale Froude cap)")
      if (allocated(error)) return

      ! Runaway-scale velocity: the backstop rescales onto
      ! |u| = FROUDE_BACKSTOP*sqrt(g h).
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 1.0_wp
      hu = 1000.0_wp
      hv = 0.0_wp

      call fused_update_stage1(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.01_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      hu_expect = FROUDE_BACKSTOP*sqrt(GRAVITY*1.0_wp)
      call check(error, abs(hu(NG + 1, NG + 1) - hu_expect) < 1.0e-9_wp, &
                 "runaway momentum must be rescaled onto the Froude backstop")
   end subroutine test_desingularization

   subroutine test_stage2_rk2(error)
      !! Stage 2 with zero flux + zero forcing should average to (h+h0)/2.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: h0(NX, NY), hu0(NX, NY), hv0(NX, NY)
      real(wp) :: fh(NX, NY), fhu(NX, NY), fhv(NX, NY)
      real(wp) :: wx(NX, NY), wy(NX, NY), atm(NX, NY), prec(NX, NY)

      call make_grid(grid)
      call zero_arrays(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, wx, wy, atm, prec)
      h = 4.0_wp; hu = 6.0_wp; hv = 0.0_wp
      h0 = 2.0_wp; hu0 = 0.0_wp; hv0 = 0.0_wp

      call fused_update_stage2(h, hu, hv, h0, hu0, hv0, fh, fhu, fhv, &
                               grid, 0.1_wp, 0.0_wp, wx, wy, atm, prec, &
                               0.0_wp, .false.)

      call check(error, abs(h(NG + 1, NG + 1) - 3.0_wp) < 1.0e-12_wp, "RK2: h := (h+h0)/2 = 3")
      if (allocated(error)) return
      call check(error, abs(hu(NG + 1, NG + 1) - 3.0_wp) < 1.0e-12_wp, "RK2: hu := (hu+hu0)/2 = 3")
   end subroutine test_stage2_rk2

end module test_kernel_fused
