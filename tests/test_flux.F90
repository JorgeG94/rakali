!! Unit tests for the HLL / HLLC flux computation kernels
module test_flux
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_kernel_flux, only: compute_flux_hll, compute_flux_hllc
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_flux_tests

contains

   subroutine collect_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("uniform_state_zero_flux", test_uniform_zero_flux), &
                  new_unittest("symmetric_dam_break", test_symmetric_dam_break), &
                  new_unittest("all_dry_zero_flux", test_all_dry_zero_flux), &
                  new_unittest("lake_at_rest_zero_flux", test_lake_at_rest), &
                  new_unittest("x_y_symmetry", test_xy_symmetry), &
                  new_unittest("flux_bounded", test_flux_bounded), &
                  new_unittest("hllc_lake_at_rest", test_hllc_lake_at_rest), &
                  new_unittest("hllc_matches_hll_no_shear", test_hllc_matches_hll_no_shear), &
                  new_unittest("hllc_shear_jump", test_hllc_shear_jump) &
                  ]
   end subroutine collect_flux_tests

   subroutine test_uniform_zero_flux(error)
      !! Uniform quiescent state should produce zero flux divergence
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 12, NY = 12, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: max_flux
      integer :: i, j

      h = 1.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                            mass_flux_x, mass_flux_y, grid)
      !$acc end data

      max_flux = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            max_flux = max(max_flux, abs(flux_h(i, j)))
            max_flux = max(max_flux, abs(flux_hu(i, j)))
            max_flux = max(max_flux, abs(flux_hv(i, j)))
         end do
      end do

      call check(error, max_flux < 1.0e-12_wp, &
                 "Uniform quiescent state should give zero flux")
   end subroutine test_uniform_zero_flux

   subroutine test_symmetric_dam_break(error)
      !! Classic dam break with h_L=2, h_R=1, u=0 at the i_mid / i_mid+1
      !! interface (flat bed, dx=dy=1). Minmod slopes vanish at the dam
      !! (one side of the minmod pair is zero, the other non-zero, so
      !! the limiter returns 0), so the reconstructed states are the
      !! cell-centre values and the HLL mass flux is analytic:
      !!     a_L = sqrt(2g), a_R = sqrt(g), s_L = -sqrt(2g), s_R = sqrt(2g)
      !!     F_h = s_L*s_R*(h_R - h_L)/(s_R - s_L) = sqrt(g/2)
      !! Net flux divergence:
      !!     flux_h(i_mid)   = -sqrt(g/2)   (cell loses mass)
      !!     flux_h(i_mid+1) = +sqrt(g/2)   (cell gains mass)
      !! and conservation demands flux_h(i_mid) + flux_h(i_mid+1) = 0.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 16, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: expected_mass_flux
      integer :: i, j, i_mid

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      i_mid = NX/2

      do j = 1, NY
         do i = 1, NX
            if (i <= i_mid) then
               h(i, j) = 2.0_wp
            else
               h(i, j) = 1.0_wp
            end if
         end do
      end do

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                            mass_flux_x, mass_flux_y, grid)
      !$acc end data

      expected_mass_flux = sqrt(0.5_wp*GRAVITY)

      ! flux_hv must be zero everywhere (no y-momentum in a 1-D x-problem).
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            call check(error, abs(flux_hv(i, j)) < 1.0e-12_wp, &
                       "flux_hv should be zero for x-direction dam break")
            if (allocated(error)) return
         end do
      end do

      ! Every interior j-row must see the analytic HLL mass flux at the
      ! dam cells and exact anti-symmetric conservation across the dam.
      do j = NG + 1, NY - NG
         call check(error, abs(flux_h(i_mid, j) + expected_mass_flux) < 1.0e-12_wp, &
                    "flux_h(i_mid) should equal -sqrt(g/2)")
         if (allocated(error)) return
         call check(error, abs(flux_h(i_mid + 1, j) - expected_mass_flux) < 1.0e-12_wp, &
                    "flux_h(i_mid+1) should equal +sqrt(g/2)")
         if (allocated(error)) return
         call check(error, abs(flux_h(i_mid, j) + flux_h(i_mid + 1, j)) < 1.0e-13_wp, &
                    "mass must be conserved across the dam interface")
         if (allocated(error)) return
      end do

      ! All other interior columns are uniform and must give zero flux.
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            if (i == i_mid .or. i == i_mid + 1) cycle
            call check(error, abs(flux_h(i, j)) < 1.0e-12_wp, &
                       "flux_h should be zero away from the dam")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_symmetric_dam_break

   subroutine test_all_dry_zero_flux(error)
      !! All-dry domain should produce zero flux
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 10, NY = 10, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: max_flux
      integer :: i, j

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                            mass_flux_x, mass_flux_y, grid)
      !$acc end data

      max_flux = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            max_flux = max(max_flux, abs(flux_h(i, j)))
            max_flux = max(max_flux, abs(flux_hu(i, j)))
            max_flux = max(max_flux, abs(flux_hv(i, j)))
         end do
      end do

      call check(error, max_flux < 1.0e-15_wp, &
                 "All-dry domain should give zero flux")
   end subroutine test_all_dry_zero_flux

   subroutine test_lake_at_rest(error)
      !! Lake at rest over variable bathymetry should give zero flux (C-property)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: max_flux, eta_0
      integer :: i, j

      eta_0 = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp

      ! Parabolic bump bathymetry
      do j = 1, NY
         do i = 1, NX
            b(i, j) = 0.5_wp*exp(-0.5_wp*real((i - NX/2)**2 + (j - NY/2)**2, wp))
            h(i, j) = max(0.0_wp, eta_0 - b(i, j))
         end do
      end do

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                            mass_flux_x, mass_flux_y, grid)
      !$acc end data

      max_flux = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            max_flux = max(max_flux, abs(flux_h(i, j)))
            max_flux = max(max_flux, abs(flux_hu(i, j)))
            max_flux = max(max_flux, abs(flux_hv(i, j)))
         end do
      end do

      call check(error, max_flux < 1.0e-10_wp, &
                 "Lake at rest should give near-zero flux (well-balanced)")
   end subroutine test_lake_at_rest

   subroutine test_xy_symmetry(error)
      !! A dam break in y should give the same flux magnitude as one in x
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 12, NG = 2
      real(wp) :: hx(N, N), hux(N, N), hvx(N, N), bx(N, N)
      real(wp) :: fhx(N, N), fhux(N, N), fhvx(N, N)
      real(wp) :: fmxx(N, N), fmxy(N, N)
      real(wp) :: hy(N, N), huy(N, N), hvy(N, N), by(N, N)
      real(wp) :: fhy(N, N), fhuy(N, N), fhvy(N, N)
      real(wp) :: fmyx(N, N), fmyy(N, N)
      type(hgrid_t) :: grid
      integer :: i, j

      call grid%init(N - 2*NG, N - 2*NG, NG, 1.0_wp, 1.0_wp)

      ! x-direction dam break
      hx = 1.0_wp
      hux = 0.0_wp
      hvx = 0.0_wp
      bx = 0.0_wp
      do j = 1, N
         do i = 1, N/2
            hx(i, j) = 2.0_wp
         end do
      end do

      ! y-direction dam break (transposed)
      hy = 1.0_wp
      huy = 0.0_wp
      hvy = 0.0_wp
      by = 0.0_wp
      do j = 1, N/2
         do i = 1, N
            hy(i, j) = 2.0_wp
         end do
      end do

      !$acc data copyin( hx, hux, hvx, bx) &
      !$acc&     copyout( fhx, fhux, fhvx, fmxx, fmxy)
      call compute_flux_hll(hx, hux, hvx, bx, fhx, fhux, fhvx, &
                            fmxx, fmxy, grid)
      !$acc end data
      !$acc data copyin( hy, huy, hvy, by) &
      !$acc&     copyout( fhy, fhuy, fhvy, fmyx, fmyy)
      call compute_flux_hll(hy, huy, hvy, by, fhy, fhuy, fhvy, &
                            fmyx, fmyy, grid)
      !$acc end data

      ! flux_h at the interface should be the same magnitude
      ! x-break at i=N/2, j=N/2; y-break at i=N/2, j=N/2
      call check(error, abs(abs(fhx(N/2, N/2)) - abs(fhy(N/2, N/2))) < 1.0e-12_wp, &
                 "x and y dam break flux_h magnitude should match")
   end subroutine test_xy_symmetry

   subroutine test_flux_bounded(error)
      !! A spatially uniform state (h=1, hu=5, hv=3 everywhere, flat bed)
      !! must give exact zero flux divergence: consistent Riemann flux
      !! means left and right interfaces see identical states, minmod
      !! slopes vanish, and the hydrostatic correction vanishes too.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 10, NY = 10, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      integer :: i, j

      h = 1.0_wp
      hu = 5.0_wp
      hv = 3.0_wp
      b = 0.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                            mass_flux_x, mass_flux_y, grid)
      !$acc end data

      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            call check(error, abs(flux_h(i, j)) < 1.0e-12_wp, &
                       "flux_h must vanish for uniform moving state")
            if (allocated(error)) return
            call check(error, abs(flux_hu(i, j)) < 1.0e-11_wp, &
                       "flux_hu must vanish for uniform moving state")
            if (allocated(error)) return
            call check(error, abs(flux_hv(i, j)) < 1.0e-11_wp, &
                       "flux_hv must vanish for uniform moving state")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_flux_bounded

   subroutine test_hllc_lake_at_rest(error)
      !! HLLC must satisfy the lake-at-rest C-property over variable
      !! bathymetry (well-balanced is a property of the hydrostatic
      !! reconstruction, not the Riemann solver, so HLLC inherits it).
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: max_flux, eta_0
      integer :: i, j

      eta_0 = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            b(i, j) = 0.5_wp*exp(-0.5_wp*real((i - NX/2)**2 + (j - NY/2)**2, wp))
            h(i, j) = max(0.0_wp, eta_0 - b(i, j))
         end do
      end do

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hllc(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                             mass_flux_x, mass_flux_y, grid)
      !$acc end data

      max_flux = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            max_flux = max(max_flux, abs(flux_h(i, j)))
            max_flux = max(max_flux, abs(flux_hu(i, j)))
            max_flux = max(max_flux, abs(flux_hv(i, j)))
         end do
      end do

      call check(error, max_flux < 1.0e-10_wp, &
                 "HLLC lake at rest must give near-zero flux (well-balanced)")
   end subroutine test_hllc_lake_at_rest

   subroutine test_hllc_matches_hll_no_shear(error)
      !! Without transverse shear (hv == 0 everywhere) HLLC and HLL
      !! must agree to round-off — the contact-wave contribution to
      !! hv is f_h * v, and v == 0 collapses it to the HLL average.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 16, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: fh_a(NX, NY), fhu_a(NX, NY), fhv_a(NX, NY)
      real(wp) :: fh_b(NX, NY), fhu_b(NX, NY), fhv_b(NX, NY)
      real(wp) :: mxa(NX, NY), mya(NX, NY), mxb(NX, NY), myb(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: maxdiff
      integer :: i, j, i_mid

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      i_mid = NX/2
      do j = 1, NY
         do i = 1, NX
            if (i <= i_mid) then
               h(i, j) = 2.0_wp
            else
               h(i, j) = 1.0_wp
            end if
         end do
      end do

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( fh_a, fhu_a, fhv_a, mxa, mya)
      call compute_flux_hll(h, hu, hv, b, fh_a, fhu_a, fhv_a, mxa, mya, grid)
      !$acc end data
      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( fh_b, fhu_b, fhv_b, mxb, myb)
      call compute_flux_hllc(h, hu, hv, b, fh_b, fhu_b, fhv_b, mxb, myb, grid)
      !$acc end data

      maxdiff = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            maxdiff = max(maxdiff, abs(fh_a(i, j) - fh_b(i, j)))
            maxdiff = max(maxdiff, abs(fhu_a(i, j) - fhu_b(i, j)))
            maxdiff = max(maxdiff, abs(fhv_a(i, j) - fhv_b(i, j)))
         end do
      end do

      call check(error, maxdiff < 1.0e-13_wp, &
                 "HLLC must match HLL in the absence of transverse shear")
   end subroutine test_hllc_matches_hll_no_shear

   subroutine test_hllc_shear_jump(error)
      !! Shear test: identical h, hu but a transverse-velocity jump in
      !! hv across the dam.  HLL averages the jump (smearing it); HLLC
      !! upwinds it across the contact, so flux_hv at the dam differs
      !! between the two solvers.  The total mass flux flux_h must
      !! agree (transverse momentum doesn't enter the h equation).
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 16, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: fh_a(NX, NY), fhu_a(NX, NY), fhv_a(NX, NY)
      real(wp) :: fh_b(NX, NY), fhu_b(NX, NY), fhv_b(NX, NY)
      real(wp) :: mxa(NX, NY), mya(NX, NY), mxb(NX, NY), myb(NX, NY)
      type(hgrid_t) :: grid
      real(wp) :: max_h_diff, max_hv_diff
      integer :: i, j, i_mid

      ! Subcritical uniform-h flow with a transverse-velocity jump.
      h = 1.0_wp
      hu = 0.5_wp
      hv = 0.0_wp
      b = 0.0_wp
      i_mid = NX/2
      do j = 1, NY
         do i = 1, NX
            if (i <= i_mid) then
               hv(i, j) = 1.0_wp
            else
               hv(i, j) = -1.0_wp
            end if
         end do
      end do

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( fh_a, fhu_a, fhv_a, mxa, mya)
      call compute_flux_hll(h, hu, hv, b, fh_a, fhu_a, fhv_a, mxa, mya, grid)
      !$acc end data
      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( fh_b, fhu_b, fhv_b, mxb, myb)
      call compute_flux_hllc(h, hu, hv, b, fh_b, fhu_b, fhv_b, mxb, myb, grid)
      !$acc end data

      ! flux_h must be identical (transverse momentum doesn't enter
      ! the h equation), but flux_hv must differ at the dam.
      max_h_diff = 0.0_wp
      max_hv_diff = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            max_h_diff = max(max_h_diff, abs(fh_a(i, j) - fh_b(i, j)))
            max_hv_diff = max(max_hv_diff, abs(fhv_a(i, j) - fhv_b(i, j)))
         end do
      end do

      call check(error, max_h_diff < 1.0e-12_wp, &
                 "flux_h must match between HLL and HLLC under shear")
      if (allocated(error)) return
      call check(error, max_hv_diff > 1.0e-3_wp, &
                 "flux_hv must differ between HLL and HLLC under shear")
   end subroutine test_hllc_shear_jump

end module test_flux
