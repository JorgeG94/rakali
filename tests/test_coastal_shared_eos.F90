!! Unit tests for the COASTAL opt-in shared nonlinear equation of state
!! (`&tracer_nml eos_form = "wright" | "roquet_spv"` routing the coastal
!! multilayer density through the shared `rki_eos` dispatch).
!!
!! Cases:
!!   * Default-off knob — `ml_parse_eos_form` keeps the legacy linear
!!     kernel for "linear" (and the config default IS "linear"), so
!!     every existing coastal namelist is bit-identical by construction.
!!   * Wright reference seawater — the coastal shared path reproduces
!!     the same physical oracles the ocean's own Wright tests lock:
!!     rho(T=10 degC, S=35 PSU, p=0) within 0.5 of 1027 kg/m^3 and
!!     rho(0, 0, 0) within 0.5 of 999.7 kg/m^3 (Wright 1997, JAOT
!!     14:735-740; also MOM6's reference surface density).
!!   * Nonlinear-vs-linear discriminator — the legacy linear kernel has
!!     an EXACTLY T-independent thermal expansion (equal density steps
!!     for equal T steps); Wright's expansion coefficient grows with T,
!!     so its steps differ measurably.  Also checks the two closures
!!     disagree outright at a warm-water state.
!!   * Roquet wiring — the coastal shared path equals the scalar
!!     `eos_density_point` evaluation cell-for-cell.  The point function
!!     itself is locked against the published-prototype oracle by
!!     `test_ocean_eos_roquet`, so this anchors the coastal plumbing to
!!     that oracle without duplicating the 75-term polynomial here.
!!   * Structured/unstructured agreement — the (nx,ny,nz) shim and the
!!     (nz,ncells) kernel produce the same density for the same (T,S,h).
!!   * Vanishing-layer fallback — h <= H_VANISHED returns eos%rho0
!!     verbatim on both backends (garbage hS/hT must not leak).
module test_coastal_shared_eos
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_config, only: config_t
   use rki_ml_eos, only: ml_equation_of_state, ml_equation_of_state_shared, &
                         ml_parse_eos_form
   use rki_ml_eos_unstr, only: ml_equation_of_state_shared_unstr
   use rki_eos, only: eos_t, eos_density_point, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      EOS_VARIANT_ROQUET_SPV
   implicit none
   private

   public :: collect_coastal_shared_eos_tests

   ! Small ghosted grid shared by the structured tests
   integer, parameter :: NX_PHYS = 4
   integer, parameter :: NY_PHYS = 3
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3

   real(wp), parameter :: RHO_0 = 1025.0_wp

contains

   subroutine collect_coastal_shared_eos_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("eos_form_default_is_legacy_linear", test_default_off), &
                  new_unittest("shared_wright_reference_seawater", test_wright_reference), &
                  new_unittest("shared_wright_differs_from_linear", test_wright_vs_linear), &
                  new_unittest("shared_roquet_matches_point_oracle", test_roquet_point_oracle), &
                  new_unittest("shared_unstr_matches_structured", test_unstr_matches_structured), &
                  new_unittest("shared_vanishing_layer_fallback", test_vanishing_layer) &
                  ]
   end subroutine collect_coastal_shared_eos_tests

   subroutine make_grid(grid)
      type(hgrid_t), intent(out) :: grid
      call grid%init(NX_PHYS, NY_PHYS, NGHOST, 100.0_wp, 100.0_wp)
   end subroutine make_grid

   pure function make_eos(variant) result(eos)
      !! Build a shared-EOS handle the way the coastal wiring does
      !! (rki_state / rki_driver): variant from the knob, rho0 mirroring
      !! the coastal reference density, p_ref left at 0 (surface).
      integer, intent(in) :: variant
      type(eos_t) :: eos
      eos%variant = variant
      eos%rho0 = RHO_0
      eos%is_init = .true.
   end function make_eos

   subroutine fill_uniform(grid, h_layer, hS_layer, hT_layer, h_val, t_val, s_val)
      !! Uniform (h, T, S) fill of the thickness-weighted tracer arrays.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(out) :: h_layer(grid%nx_total, grid%ny_total, NZ)
      real(wp), intent(out) :: hS_layer(grid%nx_total, grid%ny_total, NZ)
      real(wp), intent(out) :: hT_layer(grid%nx_total, grid%ny_total, NZ)
      real(wp), intent(in) :: h_val, t_val, s_val
      h_layer = h_val
      hS_layer = h_val*s_val
      hT_layer = h_val*t_val
   end subroutine fill_uniform

   subroutine test_default_off(error)
      !! The knob's default keeps the legacy linear kernel: config_t
      !! defaults to eos_form = "linear", and ml_parse_eos_form maps
      !! "linear" to use_shared = .false. — the shared path is opt-in
      !! only, so existing namelists are bit-identical by construction.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      logical :: use_shared
      integer :: variant

      call check(error, trim(cfg%eos_form) == "linear", &
                 "config default eos_form must be 'linear' (legacy kernel)")
      if (allocated(error)) return

      call ml_parse_eos_form(cfg%eos_form, use_shared, variant)
      call check(error,.not. use_shared, &
                 "eos_form='linear' must keep the legacy coastal kernel")
      if (allocated(error)) return

      call ml_parse_eos_form("wright", use_shared, variant)
      call check(error, use_shared .and. variant == EOS_VARIANT_WRIGHT_97, &
                 "eos_form='wright' must select the shared Wright EOS")
      if (allocated(error)) return

      call ml_parse_eos_form("roquet_spv", use_shared, variant)
      call check(error, use_shared .and. variant == EOS_VARIANT_ROQUET_SPV, &
                 "eos_form='roquet_spv' must select the shared Roquet EOS")
   end subroutine test_default_off

   subroutine test_wright_reference(error)
      !! Coastal shared path reproduces the Wright (1997) physical
      !! oracles used by the ocean's own EOS tests: standard seawater
      !! rho(10 degC, 35 PSU, 0) ~ 1027 kg/m^3 and fresh water
      !! rho(0, 0, 0) ~ 999.7 kg/m^3.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(eos_t) :: eos
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      integer :: nx, ny

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))
      eos = make_eos(EOS_VARIANT_WRIGHT_97)

      call fill_uniform(grid, h_layer, hS_layer, hT_layer, 2.5_wp, 10.0_wp, 35.0_wp)
      rho_layer = -1.0_wp
      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_layer, &
                                       eos, grid, NZ)
      call check(error, all(abs(rho_layer - 1027.0_wp) < 0.5_wp), &
                 "Wright rho(10 degC, 35 PSU, 0 Pa) off 1027 kg/m^3 by > 0.5")
      if (allocated(error)) return

      call fill_uniform(grid, h_layer, hS_layer, hT_layer, 2.5_wp, 0.0_wp, 0.0_wp)
      rho_layer = -1.0_wp
      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_layer, &
                                       eos, grid, NZ)
      call check(error, all(abs(rho_layer - 999.7_wp) < 0.5_wp), &
                 "Wright rho(0 degC, 0 PSU, 0 Pa) off 999.7 kg/m^3 by > 0.5")
   end subroutine test_wright_reference

   subroutine test_wright_vs_linear(error)
      !! The opt-in path must be a genuinely DIFFERENT closure, not a
      !! relabelled linear one.  Two probes:
      !!  (1) equal T steps at fixed S give exactly equal density steps
      !!      under the legacy linear kernel, but measurably unequal
      !!      steps under Wright (thermal expansion grows with T);
      !!  (2) at a warm-water state the two closures disagree outright.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(eos_t) :: eos
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_w(:, :, :), rho_l(:, :, :)
      real(wp) :: h_val, d1_w, d2_w, d1_l, d2_l
      integer :: nx, ny, k, ii, jj
      real(wp), parameter :: S_REF = 35.0_wp, T_REF = 10.0_wp
      real(wp), parameter :: ALPHA_T = 0.17_wp, BETA_S = 0.78_wp

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ii = NGHOST + 1
      jj = NGHOST + 1
      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_w(nx, ny, NZ), rho_l(nx, ny, NZ))
      eos = make_eos(EOS_VARIANT_WRIGHT_97)

      ! T = 5, 10, 15 degC in layers k = 1, 2, 3 at S = 35.
      h_val = 2.0_wp
      h_layer = h_val
      hS_layer = h_val*S_REF
      do k = 1, NZ
         hT_layer(:, :, k) = h_val*(real(5*k, wp))
      end do

      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_w, &
                                       eos, grid, NZ)
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_l, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, grid, NZ)

      ! (1) Step ratio: linear steps are exactly equal; Wright's differ.
      d1_l = rho_l(ii, jj, 1) - rho_l(ii, jj, 2)
      d2_l = rho_l(ii, jj, 2) - rho_l(ii, jj, 3)
      call check(error, abs(d1_l - d2_l) < 1.0e-10_wp, &
                 "legacy linear EOS must give exactly equal density steps")
      if (allocated(error)) return

      d1_w = rho_w(ii, jj, 1) - rho_w(ii, jj, 2)
      d2_w = rho_w(ii, jj, 2) - rho_w(ii, jj, 3)
      call check(error, abs(d1_w - d2_w) > 0.05_wp, &
                 "Wright density steps must differ (nonlinear thermal expansion)")
      if (allocated(error)) return

      ! (2) Outright disagreement at (T=25, S=35): linear about
      ! (T_ref=10, S_ref=35, rho_0=1025) gives 1025 - 0.17*15 = 1022.45;
      ! Wright gives the real seawater ~1023.3.
      call fill_uniform(grid, h_layer, hS_layer, hT_layer, h_val, 25.0_wp, S_REF)
      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_w, &
                                       eos, grid, NZ)
      call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_l, &
                                RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, grid, NZ)
      call check(error, abs(rho_w(ii, jj, 1) - rho_l(ii, jj, 1)) > 0.3_wp, &
                 "Wright and legacy linear must disagree at a warm-water state")
   end subroutine test_wright_vs_linear

   subroutine test_roquet_point_oracle(error)
      !! The coastal shared path must equal the scalar
      !! `eos_density_point` evaluation cell-for-cell (the point routine
      !! is locked against the published-prototype Roquet oracle by
      !! test_ocean_eos_roquet), and land in the physical seawater band.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(eos_t) :: eos
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      real(wp) :: t_k, s_k, rho_pt
      integer :: nx, ny, k, ii, jj
      real(wp), parameter :: T_OF_K(3) = [4.0_wp, 10.0_wp, 25.0_wp]
      real(wp), parameter :: S_OF_K(3) = [34.0_wp, 35.0_wp, 36.0_wp]

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ii = NGHOST + 2
      jj = NGHOST + 1
      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))
      eos = make_eos(EOS_VARIANT_ROQUET_SPV)

      h_layer = 3.0_wp
      do k = 1, NZ
         hS_layer(:, :, k) = 3.0_wp*S_OF_K(k)
         hT_layer(:, :, k) = 3.0_wp*T_OF_K(k)
      end do

      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_layer, &
                                       eos, grid, NZ)

      do k = 1, NZ
         t_k = T_OF_K(k)
         s_k = S_OF_K(k)
         rho_pt = eos_density_point(eos, t_k, s_k, 0.0_wp)
         call check(error, abs(rho_layer(ii, jj, k) - rho_pt) < 1.0e-9_wp, &
                    "coastal Roquet path must match eos_density_point")
         if (allocated(error)) return
         call check(error, rho_pt > 1020.0_wp .and. rho_pt < 1035.0_wp, &
                    "Roquet density must land in the physical seawater band")
         if (allocated(error)) return
      end do
   end subroutine test_roquet_point_oracle

   subroutine test_unstr_matches_structured(error)
      !! The unstructured (nz, ncells) shared kernel and the structured
      !! (nx, ny, nz) shared shim must agree for the same (T, S, h) —
      !! they reach the same Wright algebra through different code paths
      !! (inlined 3D impl vs the eos_density_point acc-seq routine), so
      !! agreement is to round-off, not asserted bitwise.
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(eos_t) :: eos
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      real(wp), allocatable :: h_u(:, :), hS_u(:, :), hT_u(:, :), rho_u(:, :)
      integer, parameter :: NCELLS = 5
      real(wp), parameter :: T_OF_K(3) = [6.0_wp, 12.0_wp, 22.0_wp]
      real(wp), parameter :: S_OF_K(3) = [33.0_wp, 35.0_wp, 37.0_wp]
      integer :: nx, ny, k, cell, ii, jj

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      ii = NGHOST + 1
      jj = NGHOST + 2
      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))
      allocate (h_u(NZ, NCELLS), hS_u(NZ, NCELLS), hT_u(NZ, NCELLS), rho_u(NZ, NCELLS))
      eos = make_eos(EOS_VARIANT_WRIGHT_97)

      h_layer = 1.5_wp
      do k = 1, NZ
         hS_layer(:, :, k) = 1.5_wp*S_OF_K(k)
         hT_layer(:, :, k) = 1.5_wp*T_OF_K(k)
         h_u(k, :) = 1.5_wp
         hS_u(k, :) = 1.5_wp*S_OF_K(k)
         hT_u(k, :) = 1.5_wp*T_OF_K(k)
      end do

      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_layer, &
                                       eos, grid, NZ)
      call ml_equation_of_state_shared_unstr(h_u, hS_u, hT_u, rho_u, eos, NZ, NCELLS)

      do cell = 1, NCELLS
         do k = 1, NZ
            call check(error, abs(rho_u(k, cell) - rho_layer(ii, jj, k)) < 1.0e-8_wp, &
                       "unstructured shared EOS must match the structured path")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_unstr_matches_structured

   subroutine test_vanishing_layer(error)
      !! h <= H_VANISHED falls back to eos%rho0 verbatim on both
      !! backends — garbage hS/hT in a vanished layer must not leak
      !! into rho (the shared path's D4 vanished-layer guard).
      type(error_type), allocatable, intent(out) :: error
      type(hgrid_t) :: grid
      type(eos_t) :: eos
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :), &
                               rho_layer(:, :, :)
      real(wp), allocatable :: h_u(:, :), hS_u(:, :), hT_u(:, :), rho_u(:, :)
      integer, parameter :: NCELLS = 4
      integer :: nx, ny

      call make_grid(grid)
      nx = grid%nx_total
      ny = grid%ny_total
      allocate (h_layer(nx, ny, NZ), hS_layer(nx, ny, NZ), hT_layer(nx, ny, NZ), &
                rho_layer(nx, ny, NZ))
      allocate (h_u(NZ, NCELLS), hS_u(NZ, NCELLS), hT_u(NZ, NCELLS), rho_u(NZ, NCELLS))
      eos = make_eos(EOS_VARIANT_WRIGHT_97)

      h_layer = 0.0_wp
      hS_layer = 1.0e6_wp
      hT_layer = -1.0e6_wp
      rho_layer = -1.0_wp
      call ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, rho_layer, &
                                       eos, grid, NZ)
      call check(error, all(abs(rho_layer - RHO_0) < 1.0e-12_wp), &
                 "structured shared EOS: vanished layers must fall back to rho0")
      if (allocated(error)) return

      h_u = 0.0_wp
      hS_u = 1.0e6_wp
      hT_u = -1.0e6_wp
      rho_u = -1.0_wp
      call ml_equation_of_state_shared_unstr(h_u, hS_u, hT_u, rho_u, eos, NZ, NCELLS)
      call check(error, all(abs(rho_u - RHO_0) < 1.0e-12_wp), &
                 "unstructured shared EOS: vanished layers must fall back to rho0")
   end subroutine test_vanishing_layer

end module test_coastal_shared_eos
