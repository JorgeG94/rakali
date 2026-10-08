!! Unit tests for the WENO tracer face-reconstruction ladder.
!!
!! Covers: parse/nghost guards, the rung-degradation ladder, constant-field
!! sanity, polynomial-exactness anchors per rung (would catch coefficient
!! transcription bugs), mirror symmetry (the mandatory u<0 argument-flip
!! guard), kernel-level flux-form conservation, solver-level CWC gates
!! (uniform tracer invariance + closed-basin mass conservation through the
!! full pipeline at nghost=3), square-wave extrema bounds per rung, and a
!! wall-adjacent WENO9 smoke test of the boundary order-reduction.
module test_tracer_recon
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_recon_weno, only: parse_tracer_recon, &
                             tracer_recon_required_nghost, &
                             tracer_recon_support_status, &
                             recon_rung_for_face, &
                             ml_advect_tracer_weno, &
                             TRACER_RECON_PPM, &
                             TRACER_RECON_WENO5, &
                             TRACER_RECON_WENO7, &
                             TRACER_RECON_WENO9, &
                             TRACER_RECON_OK, &
                             TRACER_RECON_REJECT_PQM, &
                             TRACER_RECON_REJECT_UNKNOWN, &
                             TRACER_RECON_REJECT_OCEAN, &
                             TRACER_RECON_REJECT_NO_ML
   use rki_grid, only: hgrid_t
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_tracers, only: ml_init_tracer
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_tracer_recon_tests

contains

   subroutine collect_tracer_recon_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("parse_valid_strings", test_parse_valid), &
                  new_unittest("parse_unknown_returns_minus1", test_parse_unknown), &
                  new_unittest("parse_pqm_returns_minus2", test_parse_pqm), &
                  new_unittest("nghost_ppm", test_nghost_ppm), &
                  new_unittest("nghost_weno5", test_nghost_weno5), &
                  new_unittest("nghost_weno7", test_nghost_weno7), &
                  new_unittest("nghost_weno9", test_nghost_weno9), &
                  new_unittest("support_pqm_rejected_any_regime", test_support_pqm), &
                  new_unittest("support_unknown_rejected", test_support_unknown), &
                  new_unittest("support_weno_ocean_rejected", test_support_ocean), &
                  new_unittest("support_weno_non_multilayer_rejected", test_support_no_ml), &
                  new_unittest("support_weno_coastal_ml_ok", test_support_coastal_ml_ok), &
                  new_unittest("support_ppm_always_ok", test_support_ppm_ok), &
                  new_unittest("rung_boundary_donor", test_rung_boundary_donor), &
                  new_unittest("rung_interior_full_weno9", test_rung_interior_weno9), &
                  new_unittest("rung_degrades_weno5_limited_down", test_rung_degrades), &
                  new_unittest("rung_ladder_table", test_rung_ladder_table), &
                  new_unittest("weno5_constant_field_face", test_weno5_constant_face), &
                  new_unittest("weno9_constant_field_face", test_weno9_constant_face), &
                  new_unittest("advect_weno5_constant_unchanged", test_advect_weno5_constant), &
                  new_unittest("advect_weno9_constant_unchanged", test_advect_weno9_constant), &
                  new_unittest("exact_linear_weno5", test_exact_linear_weno5), &
                  new_unittest("exact_quadratic_weno5", test_exact_quadratic_weno5), &
                  new_unittest("exact_cubic_weno7", test_exact_cubic_weno7), &
                  new_unittest("exact_quartic_weno9", test_exact_quartic_weno9), &
                  new_unittest("mirror_symmetry_weno5", test_mirror_weno5), &
                  new_unittest("mirror_symmetry_weno7", test_mirror_weno7), &
                  new_unittest("mirror_symmetry_weno9", test_mirror_weno9), &
                  new_unittest("conservation_kernel_weno5", test_conservation_weno5), &
                  new_unittest("conservation_kernel_weno9", test_conservation_weno9), &
                  new_unittest("solver_uniform_weno5_nghost3", test_solver_uniform_weno5), &
                  new_unittest("solver_mass_conserved_weno5_nghost3", test_solver_mass_weno5), &
                  new_unittest("solver_uniform_weno7_nghost4", test_solver_uniform_weno7), &
                  new_unittest("solver_mass_conserved_weno7_nghost4", test_solver_mass_weno7), &
                  new_unittest("solver_uniform_weno9_nghost5", test_solver_uniform_weno9), &
                  new_unittest("solver_mass_conserved_weno9_nghost5", test_solver_mass_weno9), &
                  new_unittest("square_wave_extrema_weno5", test_square_wave_weno5), &
                  new_unittest("square_wave_extrema_weno7", test_square_wave_weno7), &
                  new_unittest("square_wave_extrema_weno9", test_square_wave_weno9), &
                  new_unittest("wall_adjacent_weno9_smoke", test_wall_smoke_weno9) &
                  ]
   end subroutine collect_tracer_recon_tests

   ! --------------------------------------------------------------------------
   !  Shared helpers for the physics tests
   ! --------------------------------------------------------------------------

   pure function poly_cell_avg(c, xl, xr) result(avg)
      !! Exact average of P(x) = c(0) + c(1) x + c(2) x^2 + c(3) x^3 + c(4) x^4
      !! over [xl, xr] via the antiderivative.
      real(wp), intent(in) :: c(0:4), xl, xr
      real(wp) :: avg, qr, ql
      qr = c(0)*xr + c(1)*xr**2/2.0_wp + c(2)*xr**3/3.0_wp &
           + c(3)*xr**4/4.0_wp + c(4)*xr**5/5.0_wp
      ql = c(0)*xl + c(1)*xl**2/2.0_wp + c(2)*xl**3/3.0_wp &
           + c(3)*xl**4/4.0_wp + c(4)*xl**5/5.0_wp
      avg = (qr - ql)/(xr - xl)
   end function poly_cell_avg

   subroutine run_exactness_case(error, c, recon, label)
      !! One kernel step on exact polynomial cell averages, uniform u>0,
      !! nz=1, h=1.  In the full-rung interior band the flux-form update
      !! must reproduce the exact advected cell averages to round-off:
      !! the reconstruction reproduces any polynomial inside every
      !! candidate's space, so the swept-average faces are exact and the
      !! update telescopes to the exact translation.  This anchor catches
      !! coefficient transcription errors (it caught the prototype's PQM
      !! slope-stencil bug and would catch any WENO matrix typo).
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in) :: c(0:4)
      integer, intent(in) :: recon
      character(len=*), intent(in) :: label

      integer, parameter :: nghost = 5, ni = 32, nj = 1, nz = 1
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: u0 = 0.4_wp, dt = 1.0_wp, dxv = 1.0_wp
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: x_i, exact, qmax, err
      integer :: i, j

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = dxv; grid%dy = dxv

      h_layer = 1.0_wp
      hu_layer = u0
      hv_layer = 0.0_wp
      mass_flux_x = u0
      mass_flux_y = 0.0_wp

      ! Exact cell averages everywhere (ghosts included).
      qmax = 0.0_wp
      do j = 1, ny
         do i = 1, nx
            x_i = real(i - nx/2, wp)
            hTr_layer(i, j, 1) = poly_cell_avg(c, x_i - 0.5_wp, x_i + 0.5_wp)
            qmax = max(qmax, abs(hTr_layer(i, j, 1)))
         end do
      end do

      call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                 mass_flux_x, mass_flux_y, &
                                 hTr_layer, weno_snap, grid, nz, dt, recon)

      ! Full-rung band: both faces of cell i carry the requested rung for
      ! i in [nghost+5, nx-nghost-3]; use extra margin.
      err = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 6, nx - nghost - 6
            x_i = real(i - nx/2, wp) - u0*dt
            exact = poly_cell_avg(c, x_i - 0.5_wp, x_i + 0.5_wp)
            err = max(err, abs(hTr_layer(i, j, 1) - exact))
         end do
      end do
      call check(error, err < 5.0e-12_wp*qmax, &
                 label//": polynomial advected exactly (rel 5e-12), err/qmax fail")
   end subroutine run_exactness_case

   subroutine run_mirror_case(error, recon, label)
      !! Mirror-symmetry guard for the u<0 argument-flip: advecting a
      !! profile rightward must equal advecting the mirrored profile
      !! leftward, mirrored back, to round-off.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: recon
      character(len=*), intent(in) :: label

      integer, parameter :: nghost = 5, ni = 24, nj = 1, nz = 1
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: u0 = 0.35_wp, dt = 1.0_wp
      real(wp), parameter :: pi = acos(-1.0_wp)
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_a(nx, ny, nz), hTr_b(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: conc, err
      integer :: i, j

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = 1.0_wp; grid%dy = 1.0_wp

      h_layer = 1.0_wp
      hv_layer = 0.0_wp
      mass_flux_y = 0.0_wp

      ! Smooth asymmetric profile (x-only).
      do j = 1, ny
         do i = 1, nx
            conc = 1.0_wp + 0.5_wp*sin(2.0_wp*pi*real(i, wp)/10.0_wp) &
                   + 0.3_wp*cos(2.0_wp*pi*real(i, wp)/7.0_wp)
            hTr_a(i, j, 1) = conc
         end do
      end do
      ! Mirrored profile: b(i) = a(nx+1-i) over the FULL array.
      do j = 1, ny
         do i = 1, nx
            hTr_b(i, j, 1) = hTr_a(nx + 1 - i, j, 1)
         end do
      end do

      ! Run A: rightward.
      hu_layer = u0
      mass_flux_x = u0
      call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                 mass_flux_x, mass_flux_y, &
                                 hTr_a, weno_snap, grid, nz, dt, recon)
      ! Run B: leftward on the mirrored profile.
      hu_layer = -u0
      mass_flux_x = -u0
      call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                 mass_flux_x, mass_flux_y, &
                                 hTr_b, weno_snap, grid, nz, dt, recon)

      err = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            err = max(err, abs(hTr_b(nx + 1 - i, j, 1) - hTr_a(i, j, 1)))
         end do
      end do
      call check(error, err < 1.0e-14_wp, &
                 label//": mirror symmetry (advect-left == flipped advect-right)")
   end subroutine run_mirror_case

   subroutine run_conservation_case(error, recon, label)
      !! Kernel-level flux-form conservation: frozen nonuniform h, smooth
      !! two-layer flow that vanishes near the interior edges (so no flux
      !! crosses into the ghost band), CWC-consistent mass_flux.  Total
      !! interior tracer content per layer must be conserved to round-off
      !! over 50 steps.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: recon
      character(len=*), intent(in) :: label

      integer, parameter :: nghost = 5, ni = 24, nj = 12, nz = 2
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: dt = 0.4_wp
      real(wp), parameter :: pi = acos(-1.0_wp)
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: wcell(nx, ny), vcell(nx, ny)
      real(wp) :: mass0(nz), mass1(nz), taper_x, taper_y
      integer :: i, j, k, step, di, dj

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = 1.0_wp; grid%dy = 1.0_wp

      ! Nonuniform thicknesses, well above the vanishing tolerance.
      do j = 1, ny
         do i = 1, nx
            h_layer(i, j, 1) = 0.6_wp + 0.2_wp*sin(2.0_wp*pi*real(i, wp)/9.0_wp)
            h_layer(i, j, 2) = 0.5_wp + 0.15_wp*cos(2.0_wp*pi*real(j, wp)/11.0_wp)
         end do
      end do

      ! Cell-based transport fields, tapered to zero within 4 cells of the
      ! interior edge so no flux crosses the interior/ghost boundary.
      do j = 1, ny
         do i = 1, nx
            di = min(i - (nghost + 1), nx - nghost - i)
            dj = min(j - (nghost + 1), ny - nghost - j)
            taper_x = 0.0_wp
            taper_y = 0.0_wp
            if (di >= 4 .and. dj >= 4) then
               taper_x = 1.0_wp
               taper_y = 1.0_wp
            end if
            wcell(i, j) = taper_x*0.3_wp*sin(2.0_wp*pi*real(i + j, wp)/13.0_wp)
            vcell(i, j) = taper_y*0.25_wp*cos(2.0_wp*pi*real(i - j, wp)/8.0_wp)
         end do
      end do
      do k = 1, nz
         hu_layer(:, :, k) = 0.5_wp*wcell
         hv_layer(:, :, k) = 0.5_wp*vcell
      end do
      ! CWC-consistent barotropic flux: face value = sum_k face-average of
      ! hu_layer, so the correction term vanishes identically.
      mass_flux_x = 0.0_wp
      mass_flux_y = 0.0_wp
      do j = 1, ny
         do i = 1, nx - 1
            mass_flux_x(i, j) = 0.5_wp*(wcell(i, j) + wcell(i + 1, j))
         end do
      end do
      do j = 1, ny - 1
         do i = 1, nx
            mass_flux_y(i, j) = 0.5_wp*(vcell(i, j) + vcell(i, j + 1))
         end do
      end do

      ! Tracer: distinct smooth bumps per layer.
      do j = 1, ny
         do i = 1, nx
            hTr_layer(i, j, 1) = h_layer(i, j, 1) &
                                 *(30.0_wp + 5.0_wp*exp(-real((i - nx/2)**2 + (j - ny/2)**2, wp)/18.0_wp))
            hTr_layer(i, j, 2) = h_layer(i, j, 2) &
                                 *(10.0_wp + 3.0_wp*exp(-real((i - nx/3)**2 + (j - ny/3)**2, wp)/12.0_wp))
         end do
      end do

      mass0 = 0.0_wp
      do k = 1, nz
         do j = nghost + 1, ny - nghost
            do i = nghost + 1, nx - nghost
               mass0(k) = mass0(k) + hTr_layer(i, j, k)
            end do
         end do
      end do

      do step = 1, 50
         call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                    mass_flux_x, mass_flux_y, &
                                    hTr_layer, weno_snap, grid, nz, dt, recon)
      end do

      mass1 = 0.0_wp
      do k = 1, nz
         do j = nghost + 1, ny - nghost
            do i = nghost + 1, nx - nghost
               mass1(k) = mass1(k) + hTr_layer(i, j, k)
            end do
         end do
      end do

      call check(error, abs(mass1(1) - mass0(1)) < 1.0e-12_wp*abs(mass0(1)), &
                 label//": layer-1 tracer mass conserved to 1e-12 relative")
      if (allocated(error)) return
      call check(error, abs(mass1(2) - mass0(2)) < 1.0e-12_wp*abs(mass0(2)), &
                 label//": layer-2 tracer mass conserved to 1e-12 relative")
   end subroutine run_conservation_case

   subroutine run_square_wave_case(error, recon, bound, label)
      !! Square-wave transit: 50 steps at sigma=0.4 in the full-rung band.
      !! New extrema beyond the IC bounds [0,1] must stay under the
      !! per-rung bound.  The bounds are anchored to the numpy prototype
      !! run in THIS EXACT configuration (11-cell wave, 64 interior cells,
      !! sigma=0.4, 50 steps, zero-gradient BC): measured final overshoot
      !! weno5 9.63e-4, weno7 3.70e-5, weno9 4.48e-4 — the Fortran kernel
      !! reproduced all printed digits, which is the transcription proof.
      !! (The narrow wave overshoots more than the report's 77-cell wave;
      !! weno9 > weno7 here is a real property of the simplified weno9
      !! smoothness indicators at fronts, present in the prototype too.)
      !! All values sit two orders below PPM's accepted 1.8e-2 baseline.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: recon
      real(wp), intent(in) :: bound
      character(len=*), intent(in) :: label

      integer, parameter :: nghost = 5, ni = 64, nj = 1, nz = 1
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: u0 = 0.4_wp, dt = 1.0_wp
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: overshoot, undershoot
      integer :: i, j, g, step

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = 1.0_wp; grid%dy = 1.0_wp

      h_layer = 1.0_wp
      hu_layer = u0
      hv_layer = 0.0_wp
      mass_flux_x = u0
      mass_flux_y = 0.0_wp

      ! Square wave on interior cells 10..20 (absolute nghost+10..nghost+20);
      ! it travels 50*0.4 = 20 cells and stays >= 10 cells inside.
      hTr_layer = 0.0_wp
      do j = 1, ny
         do i = nghost + 10, nghost + 20
            hTr_layer(i, j, 1) = 1.0_wp
         end do
      end do

      do step = 1, 50
         call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                    mass_flux_x, mass_flux_y, &
                                    hTr_layer, weno_snap, grid, nz, dt, recon)
         ! Zero-gradient ghost refresh in x (what the tracer BCs do).
         do j = 1, ny
            do g = 1, nghost
               hTr_layer(g, j, 1) = hTr_layer(nghost + 1, j, 1)
               hTr_layer(nx - nghost + g, j, 1) = hTr_layer(nx - nghost, j, 1)
            end do
         end do
      end do

      overshoot = 0.0_wp
      undershoot = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            overshoot = max(overshoot, hTr_layer(i, j, 1) - 1.0_wp)
            undershoot = max(undershoot, -hTr_layer(i, j, 1))
         end do
      end do
      call check(error, overshoot <= bound .and. undershoot <= bound, &
                 label//": square-wave over/undershoot within per-rung bound")
   end subroutine run_square_wave_case

   ! --------------------------------------------------------------------------
   !  parse_tracer_recon
   ! --------------------------------------------------------------------------

   subroutine test_parse_valid(error)
      !! All recognised strings map to the correct constant.
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, parse_tracer_recon("ppm") == TRACER_RECON_PPM, &
                    "ppm -> TRACER_RECON_PPM")
         if (allocated(error)) exit checks
         call check(error, parse_tracer_recon("weno5") == TRACER_RECON_WENO5, &
                    "weno5 -> TRACER_RECON_WENO5")
         if (allocated(error)) exit checks
         call check(error, parse_tracer_recon("weno7") == TRACER_RECON_WENO7, &
                    "weno7 -> TRACER_RECON_WENO7")
         if (allocated(error)) exit checks
         call check(error, parse_tracer_recon("weno9") == TRACER_RECON_WENO9, &
                    "weno9 -> TRACER_RECON_WENO9")
      end block checks
   end subroutine test_parse_valid

   subroutine test_parse_unknown(error)
      !! Unknown string returns -1.
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_tracer_recon("eno3") == -1, &
                 "unknown string -> -1")
   end subroutine test_parse_unknown

   subroutine test_parse_pqm(error)
      !! "pqm" returns -2 (deferred scheme).
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_tracer_recon("pqm") == -2, &
                 "pqm -> -2 (deferred)")
   end subroutine test_parse_pqm

   ! --------------------------------------------------------------------------
   !  tracer_recon_required_nghost
   ! --------------------------------------------------------------------------

   subroutine test_nghost_ppm(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_required_nghost(TRACER_RECON_PPM) == 2, &
                 "PPM requires nghost = 2")
   end subroutine test_nghost_ppm

   subroutine test_nghost_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_required_nghost(TRACER_RECON_WENO5) == 3, &
                 "WENO5 requires nghost = 3")
   end subroutine test_nghost_weno5

   subroutine test_nghost_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_required_nghost(TRACER_RECON_WENO7) == 4, &
                 "WENO7 requires nghost = 4")
   end subroutine test_nghost_weno7

   subroutine test_nghost_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_required_nghost(TRACER_RECON_WENO9) == 5, &
                 "WENO9 requires nghost = 5")
   end subroutine test_nghost_weno9

   ! --------------------------------------------------------------------------
   !  tracer_recon_support_status  (fail-loud config-guard predicate)
   ! --------------------------------------------------------------------------

   subroutine test_support_pqm(error)
      !! "pqm" is deferred on EVERY config — reject regardless of regime or
      !! multilayer, so the parsed -2 must map to REJECT_PQM whether ocean
      !! or coastal, multilayer or not.
      type(error_type), allocatable, intent(out) :: error
      integer :: pqm
      checks: block
         pqm = parse_tracer_recon("pqm")
         call check(error, tracer_recon_support_status(pqm, .true., .false.) &
                    == TRACER_RECON_REJECT_PQM, "pqm on ocean -> REJECT_PQM")
         if (allocated(error)) exit checks
         call check(error, tracer_recon_support_status(pqm, .false., .true.) &
                    == TRACER_RECON_REJECT_PQM, "pqm on coastal ML -> REJECT_PQM")
         if (allocated(error)) exit checks
         call check(error, tracer_recon_support_status(pqm, .false., .false.) &
                    == TRACER_RECON_REJECT_PQM, "pqm on coastal non-ML -> REJECT_PQM")
      end block checks
   end subroutine test_support_pqm

   subroutine test_support_unknown(error)
      !! An unrecognised string maps to REJECT_UNKNOWN on any config.
      type(error_type), allocatable, intent(out) :: error
      integer :: junk
      junk = parse_tracer_recon("eno3")
      call check(error, tracer_recon_support_status(junk, .false., .true.) &
                 == TRACER_RECON_REJECT_UNKNOWN, "garbage -> REJECT_UNKNOWN")
   end subroutine test_support_unknown

   subroutine test_support_ocean(error)
      !! A weno* rung on sim_type='ocean' is rejected fail-loud (coastal-only
      !! in wave 1) even though ocean has use_multilayer=.false.
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_support_status(TRACER_RECON_WENO5, .true., .false.) &
                 == TRACER_RECON_REJECT_OCEAN, "weno5 on ocean -> REJECT_OCEAN")
   end subroutine test_support_ocean

   subroutine test_support_no_ml(error)
      !! A weno* rung on a coastal run without multilayer tracers (2D
      !! barotropic / NH) is rejected fail-loud.
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_support_status(TRACER_RECON_WENO7, .false., .false.) &
                 == TRACER_RECON_REJECT_NO_ML, "weno7 on coastal non-ML -> REJECT_NO_ML")
   end subroutine test_support_no_ml

   subroutine test_support_coastal_ml_ok(error)
      !! The valid path: a weno* rung on coastal multilayer is accepted (the
      !! nghost / use_ppm_tracer quality guards then apply downstream).
      type(error_type), allocatable, intent(out) :: error
      call check(error, tracer_recon_support_status(TRACER_RECON_WENO9, .false., .true.) &
                 == TRACER_RECON_OK, "weno9 on coastal ML -> OK")
   end subroutine test_support_coastal_ml_ok

   subroutine test_support_ppm_ok(error)
      !! The default ppm is always OK — every config, bit-identical path.
      type(error_type), allocatable, intent(out) :: error
      checks: block
         call check(error, tracer_recon_support_status(TRACER_RECON_PPM, .true., .false.) &
                    == TRACER_RECON_OK, "ppm on ocean -> OK")
         if (allocated(error)) exit checks
         call check(error, tracer_recon_support_status(TRACER_RECON_PPM, .false., .false.) &
                    == TRACER_RECON_OK, "ppm on coastal non-ML -> OK")
      end block checks
   end subroutine test_support_ppm_ok

   ! --------------------------------------------------------------------------
   !  recon_rung_for_face
   ! --------------------------------------------------------------------------

   subroutine test_rung_boundary_donor(error)
      !! A face with only 1 cell available on the upwind side must fall back
      !! to donor-cell (rung 0) regardless of the requested rung_max.
      type(error_type), allocatable, intent(out) :: error
      checks: block
         integer :: rung
         rung = recon_rung_for_face(1, 10, 4)
         call check(error, rung == 0, &
                    "avail_up=1 -> donor-cell even with rung_max=4")
         if (allocated(error)) exit checks
         rung = recon_rung_for_face(1, 1, 2)
         call check(error, rung == 0, &
                    "avail_up=1, avail_down=1 -> donor-cell")
      end block checks
   end subroutine test_rung_boundary_donor

   subroutine test_rung_interior_weno9(error)
      !! Plenty of cells available on both sides: full WENO9 rung.
      type(error_type), allocatable, intent(out) :: error
      integer :: rung
      rung = recon_rung_for_face(10, 10, 4)
      call check(error, rung == 4, &
                 "avail_up=10, avail_down=10, rung_max=4 -> WENO9 (rung 4)")
   end subroutine test_rung_interior_weno9

   subroutine test_rung_degrades(error)
      !! WENO5 needs avail_down >= 2; with avail_down=1 it must degrade to PLM.
      type(error_type), allocatable, intent(out) :: error
      integer :: rung
      rung = recon_rung_for_face(3, 1, 4)
      ! avail_up=3 >= 3 and avail_down=1 < 2 -> WENO5 fails; PLM needs avail_down>=1 -> PLM
      call check(error, rung == 1, &
                 "avail_up=3, avail_down=1, rung_max=4 -> PLM (rung 1)")
   end subroutine test_rung_degrades

   ! --------------------------------------------------------------------------
   !  Direct face-helper checks (accessed through module-level public face
   !  helpers; we test indirectly via ml_advect_tracer_weno with special ICs)
   ! --------------------------------------------------------------------------

   subroutine test_weno5_constant_face(error)
      !! A constant field (all cell averages = c) must give face value = c
      !! for every sigma.  This is a basic sanity check: all candidates
      !! reconstruct c, all betas = 0, omega = d (ideal weights), face = c.
      !!
      !! We verify by running ml_advect_tracer_weno for 1 step on a tiny
      !! uniform-tracer field and checking the result hasn't changed.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nghost = 3, ni = 10, nj = 1, nz = 1
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: tr_val = 3.7_wp
      real(wp), parameter :: h0 = 2.0_wp, u0 = 0.3_wp, dt = 0.5_wp
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: err
      integer :: i, j

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = 100.0_wp; grid%dy = 100.0_wp

      h_layer = h0
      hu_layer = u0*h0
      hv_layer = 0.0_wp
      mass_flux_x = u0*h0
      mass_flux_y = 0.0_wp
      hTr_layer = tr_val*h0  ! hTr = h*Tr

      call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                 mass_flux_x, mass_flux_y, &
                                 hTr_layer, weno_snap, grid, nz, dt, TRACER_RECON_WENO5)

      ! Concentration in interior should still be tr_val
      err = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            err = max(err, abs(hTr_layer(i, j, 1)/h_layer(i, j, 1) - tr_val))
         end do
      end do
      call check(error, err < 1.0e-12_wp, &
                 "WENO5: constant field must remain constant after 1 step")
   end subroutine test_weno5_constant_face

   subroutine test_weno9_constant_face(error)
      !! Same constant-field test with WENO9 (requires nghost=5).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nghost = 5, ni = 10, nj = 1, nz = 1
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: tr_val = 1.23_wp
      real(wp), parameter :: h0 = 1.5_wp, u0 = 0.2_wp, dt = 0.3_wp
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: err
      integer :: i, j

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = 100.0_wp; grid%dy = 100.0_wp

      h_layer = h0
      hu_layer = u0*h0
      hv_layer = 0.0_wp
      mass_flux_x = u0*h0
      mass_flux_y = 0.0_wp
      hTr_layer = tr_val*h0

      call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                 mass_flux_x, mass_flux_y, &
                                 hTr_layer, weno_snap, grid, nz, dt, TRACER_RECON_WENO9)

      err = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            err = max(err, abs(hTr_layer(i, j, 1)/h_layer(i, j, 1) - tr_val))
         end do
      end do
      call check(error, err < 1.0e-12_wp, &
                 "WENO9: constant field must remain constant after 1 step")
   end subroutine test_weno9_constant_face

   ! --------------------------------------------------------------------------
   !  Integration tests: ml_advect_tracer_weno
   ! --------------------------------------------------------------------------

   subroutine test_advect_weno5_constant(error)
      !! 10 steps with uniform tracer and uniform velocity: tracer stays constant.
      !! Uses a larger grid (nghost=3, 20 interior cells) with constant h/hu.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nghost = 3, ni = 20, nj = 3, nz = 2
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: tr_val = 5.0_wp
      real(wp), parameter :: h0 = 3.0_wp, u0 = 0.4_wp, v0 = 0.1_wp
      real(wp), parameter :: dx = 50.0_wp, dy = 50.0_wp
      real(wp), parameter :: dt = 0.5_wp*dx/(abs(u0) + abs(v0) + 1.0e-8_wp)
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: err
      integer :: i, j, step

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = dx; grid%dy = dy

      h_layer = h0
      hu_layer = u0*h0
      hv_layer = v0*h0
      mass_flux_x = u0*h0
      mass_flux_y = v0*h0
      hTr_layer = tr_val*h0

      do step = 1, 10
         call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                    mass_flux_x, mass_flux_y, &
                                    hTr_layer, weno_snap, grid, nz, dt, TRACER_RECON_WENO5)
      end do

      err = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            err = max(err, abs(hTr_layer(i, j, 1)/h_layer(i, j, 1) - tr_val))
            err = max(err, abs(hTr_layer(i, j, 2)/h_layer(i, j, 2) - tr_val))
         end do
      end do
      call check(error, err < 1.0e-11_wp, &
                 "WENO5: constant field stays constant over 10 steps (tol 1e-11)")
   end subroutine test_advect_weno5_constant

   subroutine test_advect_weno9_constant(error)
      !! Same multi-step test with WENO9 (requires nghost=5).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nghost = 5, ni = 20, nj = 3, nz = 2
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: tr_val = 7.777_wp
      real(wp), parameter :: h0 = 1.0_wp, u0 = 0.3_wp, v0 = 0.0_wp
      real(wp), parameter :: dx = 100.0_wp, dy = 100.0_wp
      real(wp), parameter :: dt = 0.5_wp*dx/(abs(u0) + 1.0e-8_wp)
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: err
      integer :: i, j, step

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = dx; grid%dy = dy

      h_layer = h0
      hu_layer = u0*h0
      hv_layer = v0*h0
      mass_flux_x = u0*h0
      mass_flux_y = v0*h0
      hTr_layer = tr_val*h0

      do step = 1, 10
         call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                    mass_flux_x, mass_flux_y, &
                                    hTr_layer, weno_snap, grid, nz, dt, TRACER_RECON_WENO9)
      end do

      err = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            err = max(err, abs(hTr_layer(i, j, 1)/h_layer(i, j, 1) - tr_val))
            err = max(err, abs(hTr_layer(i, j, 2)/h_layer(i, j, 2) - tr_val))
         end do
      end do
      call check(error, err < 1.0e-11_wp, &
                 "WENO9: constant field stays constant over 10 steps (tol 1e-11)")
   end subroutine test_advect_weno9_constant

   ! --------------------------------------------------------------------------
   !  Rung ladder table (exact expected sequence approaching a wall)
   ! --------------------------------------------------------------------------

   subroutine test_rung_ladder_table(error)
      !! The full documented ladder: donor(0) -> plm(1) -> weno5(2) ->
      !! weno7(3) -> weno9(4) as availability grows, the rung_max cap,
      !! the downwind-limited case, and the no-data sentinel.
      type(error_type), allocatable, intent(out) :: error
      call check(error, recon_rung_for_face(0, 99, 4) == -1, "avail_up=0 -> -1")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(1, 99, 4) == 0, "avail 1/99 -> donor")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(2, 99, 4) == 1, "avail 2/99 -> plm")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(3, 2, 4) == 2, "avail 3/2 -> weno5")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(4, 3, 4) == 3, "avail 4/3 -> weno7")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(5, 4, 4) == 4, "avail 5/4 -> weno9")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(5, 3, 4) == 3, "down-limited -> weno7")
      if (allocated(error)) return
      call check(error, recon_rung_for_face(99, 99, 2) == 2, "rung_max caps at weno5")
   end subroutine test_rung_ladder_table

   ! --------------------------------------------------------------------------
   !  Polynomial-exactness anchors (per rung)
   ! --------------------------------------------------------------------------

   subroutine test_exact_linear_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: c(0:4)
      c = [7.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      call run_exactness_case(error, c, TRACER_RECON_WENO5, "weno5/linear")
   end subroutine test_exact_linear_weno5

   subroutine test_exact_quadratic_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: c(0:4)
      c = [5.0_wp, 0.2_wp, 0.05_wp, 0.0_wp, 0.0_wp]
      call run_exactness_case(error, c, TRACER_RECON_WENO5, "weno5/quadratic")
   end subroutine test_exact_quadratic_weno5

   subroutine test_exact_cubic_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: c(0:4)
      c = [5.0_wp, 0.2_wp, 0.05_wp, 0.01_wp, 0.0_wp]
      call run_exactness_case(error, c, TRACER_RECON_WENO7, "weno7/cubic")
   end subroutine test_exact_cubic_weno7

   subroutine test_exact_quartic_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: c(0:4)
      c = [5.0_wp, 0.2_wp, 0.05_wp, 0.01_wp, 0.002_wp]
      call run_exactness_case(error, c, TRACER_RECON_WENO9, "weno9/quartic")
   end subroutine test_exact_quartic_weno9

   ! --------------------------------------------------------------------------
   !  Mirror symmetry (per rung) — the mandatory u<0 argument-flip guard
   ! --------------------------------------------------------------------------

   subroutine test_mirror_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mirror_case(error, TRACER_RECON_WENO5, "weno5")
   end subroutine test_mirror_weno5

   subroutine test_mirror_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mirror_case(error, TRACER_RECON_WENO7, "weno7")
   end subroutine test_mirror_weno7

   subroutine test_mirror_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_mirror_case(error, TRACER_RECON_WENO9, "weno9")
   end subroutine test_mirror_weno9

   ! --------------------------------------------------------------------------
   !  Kernel-level conservation
   ! --------------------------------------------------------------------------

   subroutine test_conservation_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_conservation_case(error, TRACER_RECON_WENO5, "weno5")
   end subroutine test_conservation_weno5

   subroutine test_conservation_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_conservation_case(error, TRACER_RECON_WENO9, "weno9")
   end subroutine test_conservation_weno9

   ! --------------------------------------------------------------------------
   !  Square-wave extrema bounds (per rung)
   ! --------------------------------------------------------------------------

   subroutine test_square_wave_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_square_wave_case(error, TRACER_RECON_WENO5, 2.0e-3_wp, "weno5")
   end subroutine test_square_wave_weno5

   subroutine test_square_wave_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_square_wave_case(error, TRACER_RECON_WENO7, 2.0e-4_wp, "weno7")
   end subroutine test_square_wave_weno7

   subroutine test_square_wave_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_square_wave_case(error, TRACER_RECON_WENO9, 1.0e-3_wp, "weno9")
   end subroutine test_square_wave_weno9

   ! --------------------------------------------------------------------------
   !  Solver-level CWC gates: weno5 through the full pipeline at nghost=3
   !  (twin of test_2d_tracer's uniformity + closed-basin conservation)
   ! --------------------------------------------------------------------------

   subroutine setup_solver_state_recon(state, recon_code, ng)
      !! Closed basin, flat bed, stage gradient driving a sloshing flow;
      !! nz=1 multilayer pipeline with the requested tracer_recon rung and
      !! ghost width (weno5→3, weno7→4, weno9→5).
      type(state_t), intent(inout) :: state
      integer, intent(in) :: recon_code, ng
      integer, parameter :: nxp = 16, nyp = 16
      real(wp), parameter :: dxv = 1.0_wp
      integer :: i, j

      state%use_multilayer = .true.
      state%multilayer%nz_ml = 1
      state%multilayer%initial_salinity = 35.0_wp
      state%multilayer%S_ref = 35.0_wp
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp

      call state%init(nxp, nyp, ng, dxv, dxv, 0.45_wp)

      state%multilayer%use_ppm_tracer = .true.
      state%multilayer%tracer_recon_code = recon_code

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = 1.0_wp + 0.05_wp*real(i, wp)*dxv
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, &
                                          state%barotropic%hv, state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, state%grid, 1)
      call ml_init_tracer(state%multilayer%h_layer, &
                          state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                          state%multilayer%initial_salinity, state%grid, 1)
      call ml_init_tracer(state%multilayer%h_layer, &
                          state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                          15.0_wp, state%grid, 1)
   end subroutine setup_solver_state_recon

   subroutine run_solver_uniform_case(error, recon_code, ng, label)
      !! Uniform tracer stays uniform under the full solver pipeline with
      !! the requested reconstruction (CWC gate through the real barotropic
      !! flux).
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: recon_code, ng
      character(len=*), intent(in) :: label
      type(state_t) :: state
      real(wp) :: s_conc, h, err
      integer :: i, j, n

      call setup_solver_state_recon(state, recon_code, ng)

      call solver_enter_data(state)
      do n = 1, 50
         call solver_step(state)
      end do
      call solver_exit_data(state)

      err = 0.0_wp
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            h = state%multilayer%h_layer(i, j, 1)
            if (h < 1.0e-10_wp) cycle
            s_conc = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1)/h
            err = max(err, abs(s_conc - state%multilayer%initial_salinity))
         end do
      end do
      call state%destroy()
      call check(error, err < 1.0e-10_wp, &
                 label//": uniform S stays uniform (CWC gate)")
   end subroutine run_solver_uniform_case

   subroutine run_solver_mass_case(error, recon_code, ng, label)
      !! Closed-basin total tracer mass conserved through the full solver
      !! pipeline with the requested reconstruction (step-S IC drives real
      !! fronts).
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: recon_code, ng
      character(len=*), intent(in) :: label
      type(state_t) :: state
      real(wp) :: mass0, mass1, drift
      integer :: i, j, n

      call setup_solver_state_recon(state, recon_code, ng)

      ! Step IC: S = 30 west half, S = 35 east half (inside the clamp range).
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            if (i - ng <= 8) then
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1) = &
                  state%multilayer%h_layer(i, j, 1)*30.0_wp
            else
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1) = &
                  state%multilayer%h_layer(i, j, 1)*35.0_wp
            end if
         end do
      end do

      mass0 = 0.0_wp
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            mass0 = mass0 + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1)
         end do
      end do

      call solver_enter_data(state)
      do n = 1, 50
         call solver_step(state)
      end do
      call solver_exit_data(state)

      mass1 = 0.0_wp
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            mass1 = mass1 + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, 1)
         end do
      end do
      drift = abs(mass1 - mass0)/abs(mass0)
      call state%destroy()
      call check(error, drift < 1.0e-12_wp, &
                 label//": closed-basin tracer mass conserved")
   end subroutine run_solver_mass_case

   subroutine test_solver_uniform_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_solver_uniform_case(error, TRACER_RECON_WENO5, 3, "solver/weno5 nghost=3")
   end subroutine test_solver_uniform_weno5

   subroutine test_solver_mass_weno5(error)
      type(error_type), allocatable, intent(out) :: error
      call run_solver_mass_case(error, TRACER_RECON_WENO5, 3, "solver/weno5 nghost=3")
   end subroutine test_solver_mass_weno5

   subroutine test_solver_uniform_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_solver_uniform_case(error, TRACER_RECON_WENO7, 4, "solver/weno7 nghost=4")
   end subroutine test_solver_uniform_weno7

   subroutine test_solver_mass_weno7(error)
      type(error_type), allocatable, intent(out) :: error
      call run_solver_mass_case(error, TRACER_RECON_WENO7, 4, "solver/weno7 nghost=4")
   end subroutine test_solver_mass_weno7

   subroutine test_solver_uniform_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_solver_uniform_case(error, TRACER_RECON_WENO9, 5, "solver/weno9 nghost=5")
   end subroutine test_solver_uniform_weno9

   subroutine test_solver_mass_weno9(error)
      type(error_type), allocatable, intent(out) :: error
      call run_solver_mass_case(error, TRACER_RECON_WENO9, 5, "solver/weno9 nghost=5")
   end subroutine test_solver_mass_weno9

   ! --------------------------------------------------------------------------
   !  Wall-adjacent WENO9 smoke (boundary order-reduction in anger)
   ! --------------------------------------------------------------------------

   subroutine test_wall_smoke_weno9(error)
      !! A front parked against the west interior edge under WENO9: the
      !! rung ladder must engage (donor/plm/weno5/weno7 near the wall),
      !! stay finite, and produce no new extrema beyond the weno5-rung
      !! bound (the widest-overshoot rung the ladder can select).
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nghost = 5, ni = 32, nj = 1, nz = 1
      integer, parameter :: nx = ni + 2*nghost, ny = nj + 2*nghost
      real(wp), parameter :: u0 = 0.3_wp, dt = 1.0_wp
      real(wp), parameter :: bound = 5.0e-4_wp
      type(hgrid_t) :: grid
      real(wp) :: h_layer(nx, ny, nz), hu_layer(nx, ny, nz), hv_layer(nx, ny, nz)
      real(wp) :: mass_flux_x(nx, ny), mass_flux_y(nx, ny)
      real(wp) :: hTr_layer(nx, ny, nz)
      real(wp) :: weno_snap(nx, ny, nz)
      real(wp) :: worst
      integer :: i, j, g, step, dir
      logical :: all_finite

      grid%nx_total = nx; grid%ny_total = ny; grid%nghost = nghost
      grid%dx = 1.0_wp; grid%dy = 1.0_wp

      h_layer = 1.0_wp
      hv_layer = 0.0_wp
      mass_flux_y = 0.0_wp

      do dir = 1, -1, -2
         ! Front at the west interior edge: 1 for the first 4 interior
         ! cells, 0 beyond.
         hTr_layer = 0.0_wp
         do j = 1, ny
            do i = 1, nghost + 4
               hTr_layer(i, j, 1) = 1.0_wp
            end do
         end do
         hu_layer = real(dir, wp)*u0
         mass_flux_x = real(dir, wp)*u0

         do step = 1, 20
            call ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                       mass_flux_x, mass_flux_y, &
                                       hTr_layer, weno_snap, grid, nz, dt, TRACER_RECON_WENO9)
            do j = 1, ny
               do g = 1, nghost
                  hTr_layer(g, j, 1) = hTr_layer(nghost + 1, j, 1)
                  hTr_layer(nx - nghost + g, j, 1) = hTr_layer(nx - nghost, j, 1)
               end do
            end do
         end do

         all_finite = .true.
         worst = 0.0_wp
         do j = nghost + 1, ny - nghost
            do i = nghost + 1, nx - nghost
               if (hTr_layer(i, j, 1) /= hTr_layer(i, j, 1)) all_finite = .false.
               worst = max(worst, hTr_layer(i, j, 1) - 1.0_wp, -hTr_layer(i, j, 1))
            end do
         end do
         call check(error, all_finite, "weno9 wall smoke: all finite")
         if (allocated(error)) return
         call check(error, worst <= bound, &
                    "weno9 wall smoke: no new extrema beyond the ladder bound")
         if (allocated(error)) return
      end do
   end subroutine test_wall_smoke_weno9

end module test_tracer_recon
