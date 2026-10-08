!! Unit tests for the coastal A-grid KPP overlay (`rki_ml_kpp`).
!!
!! KPP runs on top of the PP81 interior closure and overlays a
!! shape-function-shaped diffusivity inside the surface boundary
!! layer, plus a non-local counter-gradient tracer transport.  Because
!! the coastal path applies an explicit diffusive flux (no implicit
!! solver), the KPP kernels fuse the overlay with the apply step; the
!! tests therefore observe the *effect* on hu/hv/hS/hT and the stored
!! boundary-layer depth rather than a standalone kv field.
!!
!! Coverage (the five KPP pieces called out in the v0.1.0 TODO):
!!   * interior   — KPP with no forcing reproduces PP81 bit-for-bit
!!   * BL-depth   — bounded under stratification, fills an unstratified
!!                  column, deepens with the V_t² term
!!   * shape fn   — analytic G(sigma); overlay lifts surface mixing
!!   * convective — surface cooling activates the overlay with no wind
!!   * non-local  — destabilizing flux redistributes a uniform tracer
!!                  conservatively; stabilizing flux is a no-op
module test_coastal_kpp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, CP_WATER
   use rki_grid, only: hgrid_t
   use rki_ml_kpp, only: ml_kpp_boundary_layer, ml_kpp_viscosity, &
                         ml_kpp_diffusion_tracers, kpp_shape
   use rki_ml_vertical_viscosity, only: ml_vertical_viscosity
   use rki_ml_tracers, only: ml_vertical_diffusion_tracer
   implicit none
   private

   public :: collect_coastal_kpp_tests

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: RHO0 = 1000.0_wp
   real(wp), parameter :: ALPHA_T = 0.17_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: NU_MAX = 0.01_wp
   real(wp), parameter :: NU_BG = 1.0e-4_wp
   real(wp), parameter :: KAPPA_BG = 1.0e-5_wp
   real(wp), parameter :: RI_ALPHA = 5.0_wp
   real(wp), parameter :: RI_CRIT = 0.3_wp
   real(wp), parameter :: CS_NONLOCAL = 6.3_wp
   real(wp), parameter :: C_VT2 = 1.8_wp

contains

   subroutine collect_coastal_kpp_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("kpp_shape_analytic", test_shape_analytic), &
                  new_unittest("kpp_no_forcing_equals_pp81", test_no_forcing_equals_pp81), &
                  new_unittest("kpp_stratified_bl_in_column", test_strat_bl_depth), &
                  new_unittest("kpp_unstratified_fills_column", test_unstrat_fills_column), &
                  new_unittest("kpp_overlay_lifts_surface_mixing", test_overlay_lifts_mixing), &
                  new_unittest("kpp_convective_activates_no_wind", test_convective_no_wind), &
                  new_unittest("kpp_vt2_deepens_bl", test_vt2_deepens_bl), &
                  new_unittest("kpp_vt2_no_forcing_invariant", test_vt2_no_forcing_invariant), &
                  new_unittest("kpp_nonlocal_redistributes_conserves", test_nonlocal_redistributes), &
                  new_unittest("kpp_stabilizing_flux_no_nonlocal", test_stabilizing_no_nonlocal) &
                  ]
   end subroutine collect_coastal_kpp_tests

   ! -----------------------------------------------------------------
   ! Helpers
   ! -----------------------------------------------------------------

   subroutine run_kpp(grid, nz, dt, h, hu, hv, rho, hS, hT, &
                      wind_x, wind_y, heat_flux, emp, c_vt2_in, &
                      u_star, buoy_flux, bl_depth)
      !! Run the three KPP kernels (BL-depth sweep, momentum, tracers)
      !! with explicit GPU data management (inert on host builds).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hS(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hT(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: heat_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: emp(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: c_vt2_in
      real(wp), intent(inout) :: u_star(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: buoy_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: bl_depth(grid%nx_total, grid%ny_total)

      !$acc enter data copyin(h, hu, hv, rho, hS, hT, wind_x, wind_y, &
      !$acc&                  heat_flux, emp, u_star, buoy_flux, bl_depth)
      call ml_kpp_boundary_layer(h, hu, hv, rho, hS, RHO0, ALPHA_T, BETA_S, &
                                 CP_WATER, wind_x, wind_y, heat_flux, emp, &
                                 RI_CRIT, c_vt2_in, u_star, buoy_flux, bl_depth, &
                                 grid, nz)
      call ml_kpp_viscosity(h, hu, hv, rho, RHO0, NU_MAX, NU_BG, RI_ALPHA, &
                            u_star, buoy_flux, bl_depth, grid, nz, dt)
      call ml_kpp_diffusion_tracers(h, hu, hv, hS, hT, rho, RHO0, NU_MAX, &
                                    KAPPA_BG, KAPPA_BG, RI_ALPHA, CS_NONLOCAL, &
                                    heat_flux, emp, CP_WATER, u_star, buoy_flux, &
                                    bl_depth, grid, nz, dt)
      !$acc update self(hu, hv, hS, hT, u_star, buoy_flux, bl_depth)
      !$acc exit data delete(h, hu, hv, rho, hS, hT, wind_x, wind_y, &
      !$acc&                 heat_flux, emp, u_star, buoy_flux, bl_depth)
   end subroutine run_kpp

   subroutine run_pp81(grid, nz, dt, h, hu, hv, rho, hS, hT)
      !! Run the legacy PP81 kernels (momentum + S + T) the KPP path
      !! replaces, with matching GPU data management.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hS(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hT(grid%nx_total, grid%ny_total, nz)

      !$acc enter data copyin(h, hu, hv, rho, hS, hT)
      call ml_vertical_viscosity(h, hu, hv, rho, RHO0, NU_MAX, NU_BG, RI_ALPHA, &
                                 grid, nz, dt)
      call ml_vertical_diffusion_tracer(h, hu, hv, hS, rho, RHO0, NU_MAX, &
                                        KAPPA_BG, RI_ALPHA, grid, nz, dt)
      call ml_vertical_diffusion_tracer(h, hu, hv, hT, rho, RHO0, NU_MAX, &
                                        KAPPA_BG, RI_ALPHA, grid, nz, dt)
      !$acc update self(hu, hv, hS, hT)
      !$acc exit data delete(h, hu, hv, rho, hS, hT)
   end subroutine run_pp81

   integer function iprobe(grid) result(ip)
      type(hgrid_t), intent(in) :: grid
      ip = grid%nx_total/2
   end function iprobe

   integer function jprobe(grid) result(jp)
      type(hgrid_t), intent(in) :: grid
      jp = grid%ny_total/2
   end function jprobe

   ! -----------------------------------------------------------------
   ! Cases
   ! -----------------------------------------------------------------

   subroutine test_shape_analytic(error)
      !! The cubic shape function G(sigma) = sigma(1-sigma)² vanishes at
      !! the surface (sigma=0) and BL base (sigma=1) and peaks at
      !! sigma=1/3 with G=4/27.
      type(error_type), allocatable, intent(out) :: error
      real(wp), parameter :: TOL = 1.0e-12_wp
      checks: block
         call check(error, abs(kpp_shape(0.0_wp)) < TOL, "G(0) /= 0")
         if (allocated(error)) exit checks
         call check(error, abs(kpp_shape(1.0_wp)) < TOL, "G(1) /= 0")
         if (allocated(error)) exit checks
         call check(error, abs(kpp_shape(1.0_wp/3.0_wp) - 4.0_wp/27.0_wp) < TOL, &
                    "G(1/3) /= 4/27 (peak)")
         if (allocated(error)) exit checks
         call check(error, kpp_shape(0.5_wp) < kpp_shape(1.0_wp/3.0_wp), &
                    "G not peaked near sigma=1/3")
      end block checks
   end subroutine test_shape_analytic

   subroutine test_no_forcing_equals_pp81(error)
      !! KPP with zero wind + zero surface flux ⇒ u_* = 0, B_0 = 0 ⇒
      !! w_s = 0 ⇒ overlay inactive and non-local skipped at every
      !! interface, so the momentum + tracer fluxes must reproduce the
      !! PP81 kernels bit-for-bit on a sheared, stratified column.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 5
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 1.0_wp, TOL = 1.0e-12_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :)
      real(wp), allocatable :: hu_a(:, :, :), hv_a(:, :, :), hS_a(:, :, :), hT_a(:, :, :)
      real(wp), allocatable :: hu_b(:, :, :), hv_b(:, :, :), hS_b(:, :, :), hT_b(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny, k
      checks: block
         call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (h(nx, ny, NZ), rho(nx, ny, NZ))
         allocate (hu_a(nx, ny, NZ), hv_a(nx, ny, NZ), hS_a(nx, ny, NZ), hT_a(nx, ny, NZ))
         allocate (hu_b(nx, ny, NZ), hv_b(nx, ny, NZ), hS_b(nx, ny, NZ), hT_b(nx, ny, NZ))
         allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
         allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_depth(nx, ny))

         h = H_LAYER
         do k = 1, NZ
            rho(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)        ! densest at bed
            hu_a(:, :, k) = H_LAYER*0.1_wp*real(k - 1, wp)           ! sheared
            hv_a(:, :, k) = 0.0_wp
            hS_a(:, :, k) = H_LAYER*(35.0_wp - 0.1_wp*real(k - 1, wp))
            hT_a(:, :, k) = H_LAYER*(15.0_wp + 0.2_wp*real(k - 1, wp))
         end do
         hu_b = hu_a
         hv_b = hv_a
         hS_b = hS_a
         hT_b = hT_a
         wind_x = 0.0_wp
         wind_y = 0.0_wp
         heat_flux = 0.0_wp
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_depth = 0.0_wp

         call run_pp81(grid, NZ, DT, h, hu_a, hv_a, rho, hS_a, hT_a)
         call run_kpp(grid, NZ, DT, h, hu_b, hv_b, rho, hS_b, hT_b, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

         call check(error, maxval(abs(hu_b - hu_a)) < TOL, "KPP hu /= PP81 hu (no forcing)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(hv_b - hv_a)) < TOL, "KPP hv /= PP81 hv (no forcing)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(hS_b - hS_a)) < TOL, "KPP hS /= PP81 hS (no forcing)")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(hT_b - hT_a)) < TOL, "KPP hT /= PP81 hT (no forcing)")
      end block checks
      if (allocated(h)) deallocate (h, rho, hu_a, hv_a, hS_a, hT_a, hu_b, hv_b, hS_b, hT_b, &
                                    wind_x, wind_y, heat_flux, emp, u_star, buoy_flux, bl_depth)
   end subroutine test_no_forcing_equals_pp81

   subroutine test_strat_bl_depth(error)
      !! Mixed-layer-over-thermocline + wind: bulk-Ri stays near zero in
      !! the well-mixed top, jumps at the thermocline.  BL must extend
      !! below the surface layer and stop short of the bed.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 8, NZ_ML = 3
      real(wp), parameter :: H_LAYER = 25.0_wp, DT = 1.0_wp, TAU = 0.1_wp, U_ML = 0.5_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :), hu(:, :, :), hv(:, :, :)
      real(wp), allocatable :: hS(:, :, :), hT(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny, k
      real(wp) :: h_b
      checks: block
         call grid%init(14, 12, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         call alloc3(h, rho, hu, hv, hS, hT, nx, ny, NZ)
         call alloc2(wind_x, wind_y, heat_flux, emp, u_star, buoy_flux, bl_depth, nx, ny)

         h = H_LAYER
         hv = 0.0_wp
         hS = H_LAYER*35.0_wp
         hT = H_LAYER*15.0_wp
         do k = 1, NZ
            if (k > NZ - NZ_ML) then
               hu(:, :, k) = H_LAYER*U_ML
               rho(:, :, k) = 1027.0_wp
            else
               hu(:, :, k) = 0.0_wp
               rho(:, :, k) = 1030.0_wp - 0.5_wp*real(k - 1, wp)
            end if
         end do
         wind_x = TAU
         wind_y = 0.0_wp
         heat_flux = 0.0_wp
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_depth = 0.0_wp

         call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

         h_b = bl_depth(iprobe(grid), jprobe(grid))
         call check(error, h_b > H_LAYER, "stratified BL collapsed to a sub-layer")
         if (allocated(error)) exit checks
         call check(error, h_b < real(NZ, wp)*H_LAYER, "stratified BL reached the bed")
      end block checks
      call dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                       u_star, buoy_flux, bl_depth)
   end subroutine test_strat_bl_depth

   subroutine test_unstrat_fills_column(error)
      !! Uniform density + wind: bulk-Ri never crosses ⇒ BL = total H.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 4
      real(wp), parameter :: H_LAYER = 20.0_wp, DT = 1.0_wp, TAU = 0.05_wp, TOL = 1.0e-9_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :), hu(:, :, :), hv(:, :, :)
      real(wp), allocatable :: hS(:, :, :), hT(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny
      call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      call alloc3(h, rho, hu, hv, hS, hT, nx, ny, NZ)
      call alloc2(wind_x, wind_y, heat_flux, emp, u_star, buoy_flux, bl_depth, nx, ny)

      h = H_LAYER
      rho = 1030.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      hS = H_LAYER*35.0_wp
      hT = H_LAYER*15.0_wp
      wind_x = TAU
      wind_y = 0.0_wp
      heat_flux = 0.0_wp
      emp = 0.0_wp
      u_star = 0.0_wp
      buoy_flux = 0.0_wp
      bl_depth = 0.0_wp

      call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                   wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

      call check(error, abs(bl_depth(iprobe(grid), jprobe(grid)) - real(NZ, wp)*H_LAYER) < TOL, &
                 "unstratified BL did not fill the column")
      call dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                       u_star, buoy_flux, bl_depth)
   end subroutine test_unstrat_fills_column

   subroutine test_overlay_lifts_mixing(error)
      !! Wind on a weakly-stratified, sheared column: the KPP overlay
      !! must mix momentum more strongly than PP81 alone — the spread
      !! of layer velocities (max-min u) drops further under KPP.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 30.0_wp, TAU = 0.2_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :)
      real(wp), allocatable :: hu_a(:, :, :), hv_a(:, :, :), hS_a(:, :, :), hT_a(:, :, :)
      real(wp), allocatable :: hu_b(:, :, :), hv_b(:, :, :), hS_b(:, :, :), hT_b(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny, k, ip, jp
      real(wp) :: spread_a, spread_b
      real(wp) :: u_a(NZ), u_b(NZ)
      checks: block
         call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (h(nx, ny, NZ), rho(nx, ny, NZ))
         allocate (hu_a(nx, ny, NZ), hv_a(nx, ny, NZ), hS_a(nx, ny, NZ), hT_a(nx, ny, NZ))
         allocate (hu_b(nx, ny, NZ), hv_b(nx, ny, NZ), hS_b(nx, ny, NZ), hT_b(nx, ny, NZ))
         allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
         allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_depth(nx, ny))

         h = H_LAYER
         do k = 1, NZ
            rho(:, :, k) = 1028.0_wp - 0.02_wp*real(k - 1, wp)   ! weak stratification
            hu_a(:, :, k) = H_LAYER*0.05_wp*real(k - 1, wp)      ! sheared
            hv_a(:, :, k) = 0.0_wp
            hS_a(:, :, k) = H_LAYER*35.0_wp
            hT_a(:, :, k) = H_LAYER*15.0_wp
         end do
         hu_b = hu_a
         hv_b = hv_a
         hS_b = hS_a
         hT_b = hT_a
         heat_flux = 0.0_wp
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_depth = 0.0_wp

         ! PP81 reference (no wind reaches the interior closure anyway).
         wind_x = 0.0_wp
         wind_y = 0.0_wp
         call run_pp81(grid, NZ, DT, h, hu_a, hv_a, rho, hS_a, hT_a)

         ! KPP with wind.
         wind_x = TAU
         wind_y = 0.0_wp
         call run_kpp(grid, NZ, DT, h, hu_b, hv_b, rho, hS_b, hT_b, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

         ip = iprobe(grid)
         jp = jprobe(grid)
         do k = 1, NZ
            u_a(k) = hu_a(ip, jp, k)/h(ip, jp, k)
            u_b(k) = hu_b(ip, jp, k)/h(ip, jp, k)
         end do
         spread_a = maxval(u_a) - minval(u_a)
         spread_b = maxval(u_b) - minval(u_b)

         call check(error, spread_b < spread_a, &
                    "KPP overlay did not mix momentum more than PP81")
         if (allocated(error)) exit checks
         call check(error, bl_depth(ip, jp) > H_LAYER, "BL depth not established under wind")
      end block checks
      if (allocated(h)) deallocate (h, rho, hu_a, hv_a, hS_a, hT_a, hu_b, hv_b, hS_b, hT_b, &
                                    wind_x, wind_y, heat_flux, emp, u_star, buoy_flux, bl_depth)
   end subroutine test_overlay_lifts_mixing

   subroutine test_convective_no_wind(error)
      !! Surface cooling with NO wind on a neutral (unstratified)
      !! column: the bulk-Ri sweep never crosses ⇒ the BL fills the
      !! column, and the convective velocity scale w_* = (-B_0·h_b)^(1/3)
      !! alone (u_* = 0) must activate the overlay and mix the sheared
      !! momentum more than the quiescent PP81 baseline.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 30.0_wp, Q_COOL = -300.0_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :)
      real(wp), allocatable :: hu_a(:, :, :), hv_a(:, :, :), hS_a(:, :, :), hT_a(:, :, :)
      real(wp), allocatable :: hu_b(:, :, :), hv_b(:, :, :), hS_b(:, :, :), hT_b(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny, k, ip, jp
      real(wp) :: spread_a, spread_b
      real(wp) :: u_a(NZ), u_b(NZ)
      checks: block
         call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         allocate (h(nx, ny, NZ), rho(nx, ny, NZ))
         allocate (hu_a(nx, ny, NZ), hv_a(nx, ny, NZ), hS_a(nx, ny, NZ), hT_a(nx, ny, NZ))
         allocate (hu_b(nx, ny, NZ), hv_b(nx, ny, NZ), hS_b(nx, ny, NZ), hT_b(nx, ny, NZ))
         allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
         allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_depth(nx, ny))

         h = H_LAYER
         do k = 1, NZ
            rho(:, :, k) = 1025.0_wp                            ! neutral: BL fills column
            hu_a(:, :, k) = H_LAYER*0.05_wp*real(k - 1, wp)     ! sheared
            hv_a(:, :, k) = 0.0_wp
            hS_a(:, :, k) = H_LAYER*35.0_wp
            hT_a(:, :, k) = H_LAYER*15.0_wp
         end do
         hu_b = hu_a
         hv_b = hv_a
         hS_b = hS_a
         hT_b = hT_a
         wind_x = 0.0_wp
         wind_y = 0.0_wp
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_depth = 0.0_wp

         heat_flux = 0.0_wp
         call run_pp81(grid, NZ, DT, h, hu_a, hv_a, rho, hS_a, hT_a)

         heat_flux = Q_COOL
         call run_kpp(grid, NZ, DT, h, hu_b, hv_b, rho, hS_b, hT_b, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

         ip = iprobe(grid)
         jp = jprobe(grid)
         do k = 1, NZ
            u_a(k) = hu_a(ip, jp, k)/h(ip, jp, k)
            u_b(k) = hu_b(ip, jp, k)/h(ip, jp, k)
         end do
         spread_a = maxval(u_a) - minval(u_a)
         spread_b = maxval(u_b) - minval(u_b)

         call check(error, buoy_flux(ip, jp) < 0.0_wp, "cooling did not produce B_0 < 0")
         if (allocated(error)) exit checks
         call check(error, spread_b < spread_a, &
                    "convective overlay did not mix momentum without wind")
      end block checks
      if (allocated(h)) deallocate (h, rho, hu_a, hv_a, hS_a, hT_a, hu_b, hv_b, hS_b, hT_b, &
                                    wind_x, wind_y, heat_flux, emp, u_star, buoy_flux, bl_depth)
   end subroutine test_convective_no_wind

   subroutine test_vt2_deepens_bl(error)
      !! V_t² (LMD94 eq 23) adds unresolved shear to the bulk-Ri
      !! denominator ⇒ Ri_bulk crosses ri_crit deeper.  Same wind +
      !! cooling, c_vt2 = 0 vs 1.8: the V_t²-on BL must be ≥ as deep,
      !! and strictly deeper here.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 5.0_wp, DT = 1.0_wp, TAU = 0.05_wp, Q_COOL = -200.0_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :), hu(:, :, :), hv(:, :, :)
      real(wp), allocatable :: hS(:, :, :), hT(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_off(:, :), bl_on(:, :)
      integer :: nx, ny, k, ip, jp
      checks: block
         call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         call alloc3(h, rho, hu, hv, hS, hT, nx, ny, NZ)
         allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
         allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_off(nx, ny), bl_on(nx, ny))

         h = H_LAYER
         hu = 0.0_wp
         hv = 0.0_wp
         hS = H_LAYER*35.0_wp
         hT = H_LAYER*15.0_wp
         do k = 1, NZ
            rho(:, :, k) = 1030.0_wp - 0.2_wp*real(k - 1, wp)
         end do
         wind_x = TAU
         wind_y = 0.0_wp
         heat_flux = Q_COOL
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_off = 0.0_wp
         bl_on = 0.0_wp

         call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                      wind_x, wind_y, heat_flux, emp, 0.0_wp, u_star, buoy_flux, bl_off)
         ! reset prognostics (run_kpp mutated them); BL sweep only reads h/hu/hv/rho/hS
         hu = 0.0_wp
         hv = 0.0_wp
         hS = H_LAYER*35.0_wp
         hT = H_LAYER*15.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_on)

         ip = iprobe(grid)
         jp = jprobe(grid)
         call check(error, bl_on(ip, jp) >= bl_off(ip, jp), "V_t²-on BL shallower than V_t²-off")
         if (allocated(error)) exit checks
         call check(error, bl_on(ip, jp) > bl_off(ip, jp) + 1.0e-6_wp, &
                    "V_t² did not visibly deepen the BL")
      end block checks
      call dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                       u_star, buoy_flux, bl_off)
      if (allocated(bl_on)) deallocate (bl_on)
   end subroutine test_vt2_deepens_bl

   subroutine test_vt2_no_forcing_invariant(error)
      !! No wind, no flux ⇒ u_* = 0, w_* = 0 ⇒ w_s_col = 0 ⇒ V_t² = 0
      !! whatever c_vt2 is.  BL depth must be independent of c_vt2.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 5.0_wp, DT = 1.0_wp, TOL = 1.0e-12_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :), hu(:, :, :), hv(:, :, :)
      real(wp), allocatable :: hS(:, :, :), hT(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_a(:, :), bl_b(:, :)
      integer :: nx, ny, k
      call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      call alloc3(h, rho, hu, hv, hS, hT, nx, ny, NZ)
      allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
      allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_a(nx, ny), bl_b(nx, ny))

      h = H_LAYER
      hu = 0.0_wp
      hv = 0.0_wp
      hS = H_LAYER*35.0_wp
      hT = H_LAYER*15.0_wp
      do k = 1, NZ
         rho(:, :, k) = 1030.0_wp - 0.2_wp*real(k - 1, wp)
      end do
      wind_x = 0.0_wp
      wind_y = 0.0_wp
      heat_flux = 0.0_wp
      emp = 0.0_wp
      u_star = 0.0_wp
      buoy_flux = 0.0_wp
      bl_a = 0.0_wp
      bl_b = 0.0_wp

      call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                   wind_x, wind_y, heat_flux, emp, 0.0_wp, u_star, buoy_flux, bl_a)
      hu = 0.0_wp
      hv = 0.0_wp
      hS = H_LAYER*35.0_wp
      hT = H_LAYER*15.0_wp
      u_star = 0.0_wp
      buoy_flux = 0.0_wp
      call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                   wind_x, wind_y, heat_flux, emp, 10.0_wp, u_star, buoy_flux, bl_b)

      call check(error, maxval(abs(bl_a - bl_b)) < TOL, &
                 "V_t² changed BL depth with no forcing (w_s_col guard missing)")
      call dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                       u_star, buoy_flux, bl_a)
      if (allocated(bl_b)) deallocate (bl_b)
   end subroutine test_vt2_no_forcing_invariant

   subroutine test_nonlocal_redistributes(error)
      !! Destabilizing surface flux (cooling) on a column with UNIFORM
      !! temperature: the gradient-driven diffusive flux is zero, so any
      !! change in hT is purely the KPP non-local counter-gradient term.
      !! It must (a) make the profile non-uniform and (b) conserve the
      !! column-integrated heat (gamma=0 at both boundaries).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 60.0_wp, Q_COOL = -300.0_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :), hu(:, :, :), hv(:, :, :)
      real(wp), allocatable :: hS(:, :, :), hT(:, :, :), hT0(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny, ip, jp
      real(wp) :: col0, col1, spread
      checks: block
         call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         call alloc3(h, rho, hu, hv, hS, hT, nx, ny, NZ)
         allocate (hT0(nx, ny, NZ))
         allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
         allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_depth(nx, ny))

         ! Uniform T and S ⇒ uniform rho ⇒ no buoyancy gradient ⇒ no
         ! diffusive tracer flux.  Quiescent (no shear).
         h = H_LAYER
         hu = 0.0_wp
         hv = 0.0_wp
         rho = 1025.0_wp
         hS = H_LAYER*35.0_wp
         hT = H_LAYER*15.0_wp
         hT0 = hT
         wind_x = 0.0_wp
         wind_y = 0.0_wp
         heat_flux = Q_COOL
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_depth = 0.0_wp

         call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

         ip = iprobe(grid)
         jp = jprobe(grid)
         spread = maxval(hT(ip, jp, :)) - minval(hT(ip, jp, :))
         col0 = sum(hT0(ip, jp, :))
         col1 = sum(hT(ip, jp, :))

         call check(error, spread > 1.0e-6_wp, &
                    "non-local transport left a uniform column unchanged")
         if (allocated(error)) exit checks
         call check(error, abs(col1 - col0) < 1.0e-9_wp*max(abs(col0), 1.0_wp), &
                    "non-local transport did not conserve column heat")
      end block checks
      call dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                       u_star, buoy_flux, bl_depth)
      if (allocated(hT0)) deallocate (hT0)
   end subroutine test_nonlocal_redistributes

   subroutine test_stabilizing_no_nonlocal(error)
      !! Stabilizing surface flux (heating, B_0 > 0) on a uniform-T
      !! column: not destabilizing ⇒ no non-local transport ⇒ a uniform
      !! profile stays uniform (the isolated kernel does not apply the
      !! surface boundary flux itself).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      integer, parameter :: NZ = 6
      real(wp), parameter :: H_LAYER = 10.0_wp, DT = 60.0_wp, Q_WARM = 300.0_wp, TOL = 1.0e-12_wp
      real(wp), allocatable :: h(:, :, :), rho(:, :, :), hu(:, :, :), hv(:, :, :)
      real(wp), allocatable :: hS(:, :, :), hT(:, :, :), hT0(:, :, :)
      real(wp), allocatable :: wind_x(:, :), wind_y(:, :), heat_flux(:, :), emp(:, :)
      real(wp), allocatable :: u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer :: nx, ny, ip, jp
      checks: block
         call grid%init(12, 10, NGHOST, 1.0_wp, 1.0_wp)
         nx = grid%nx_total
         ny = grid%ny_total
         call alloc3(h, rho, hu, hv, hS, hT, nx, ny, NZ)
         allocate (hT0(nx, ny, NZ))
         allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
         allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_depth(nx, ny))

         h = H_LAYER
         hu = 0.0_wp
         hv = 0.0_wp
         rho = 1025.0_wp
         hS = H_LAYER*35.0_wp
         hT = H_LAYER*15.0_wp
         hT0 = hT
         wind_x = 0.0_wp
         wind_y = 0.0_wp
         heat_flux = Q_WARM
         emp = 0.0_wp
         u_star = 0.0_wp
         buoy_flux = 0.0_wp
         bl_depth = 0.0_wp

         call run_kpp(grid, NZ, DT, h, hu, hv, rho, hS, hT, &
                      wind_x, wind_y, heat_flux, emp, C_VT2, u_star, buoy_flux, bl_depth)

         ip = iprobe(grid)
         jp = jprobe(grid)
         call check(error, buoy_flux(ip, jp) > 0.0_wp, "heating did not produce B_0 > 0")
         if (allocated(error)) exit checks
         call check(error, maxval(abs(hT(ip, jp, :) - hT0(ip, jp, :))) < TOL, &
                    "stabilizing flux triggered non-local transport")
      end block checks
      call dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                       u_star, buoy_flux, bl_depth)
      if (allocated(hT0)) deallocate (hT0)
   end subroutine test_stabilizing_no_nonlocal

   ! -----------------------------------------------------------------
   ! Allocation helpers (trim repetition across cases)
   ! -----------------------------------------------------------------

   subroutine alloc3(h, rho, hu, hv, hS, hT, nx, ny, nz)
      real(wp), allocatable, intent(out) :: h(:, :, :), rho(:, :, :), hu(:, :, :)
      real(wp), allocatable, intent(out) :: hv(:, :, :), hS(:, :, :), hT(:, :, :)
      integer, intent(in) :: nx, ny, nz
      allocate (h(nx, ny, nz), rho(nx, ny, nz), hu(nx, ny, nz))
      allocate (hv(nx, ny, nz), hS(nx, ny, nz), hT(nx, ny, nz))
   end subroutine alloc3

   subroutine alloc2(wind_x, wind_y, heat_flux, emp, u_star, buoy_flux, bl_depth, nx, ny)
      real(wp), allocatable, intent(out) :: wind_x(:, :), wind_y(:, :), heat_flux(:, :)
      real(wp), allocatable, intent(out) :: emp(:, :), u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      integer, intent(in) :: nx, ny
      allocate (wind_x(nx, ny), wind_y(nx, ny), heat_flux(nx, ny), emp(nx, ny))
      allocate (u_star(nx, ny), buoy_flux(nx, ny), bl_depth(nx, ny))
   end subroutine alloc2

   subroutine dealloc_all(h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                          u_star, buoy_flux, bl_depth)
      real(wp), allocatable, intent(inout) :: h(:, :, :), rho(:, :, :), hu(:, :, :)
      real(wp), allocatable, intent(inout) :: hv(:, :, :), hS(:, :, :), hT(:, :, :)
      real(wp), allocatable, intent(inout) :: wind_x(:, :), wind_y(:, :), heat_flux(:, :)
      real(wp), allocatable, intent(inout) :: emp(:, :), u_star(:, :), buoy_flux(:, :), bl_depth(:, :)
      if (allocated(h)) deallocate (h, rho, hu, hv, hS, hT, wind_x, wind_y, heat_flux, emp, &
                                    u_star, buoy_flux, bl_depth)
   end subroutine dealloc_all

end module test_coastal_kpp
