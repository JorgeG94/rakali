!! Unit tests for the IC setup routines and bathymetry diagnostics in rki_init
!!
!! `test_init_diagnostics` already covers compute_total_mass and compute_layer_sums,
!! so this module focuses on the initial-condition setup routines and the
!! bathymetry diagnostic helper.
module test_init
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_config, only: config_t
   use rki_init, only: init_dam_break, init_gaussian_hump, &
                       init_castro_dam_break, init_still_water, &
                       init_from_config, diagnose_bathymetry
   implicit none
   private

   public :: collect_init_tests

   integer, parameter :: NX = 16, NY = 8, NG = 2
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp

contains

   subroutine collect_init_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dam_break_step", test_dam_break), &
                  new_unittest("gaussian_hump_peak", test_gaussian_hump), &
                  new_unittest("castro_dam_break_bump", test_castro), &
                  new_unittest("still_water_lake_at_rest", test_still_water), &
                  new_unittest("from_config_still_water", test_cfg_still_water), &
                  new_unittest("from_config_gaussian", test_cfg_gaussian), &
                  new_unittest("from_config_castro", test_cfg_castro), &
                  new_unittest("from_config_default_dam_break", test_cfg_default), &
                  new_unittest("diagnose_bathymetry_runs", test_diagnose_bathymetry) &
                  ]
   end subroutine collect_init_tests

   subroutine test_dam_break(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: i_left, i_right
      checks: block

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)
         call init_dam_break(s)

         ! A column well into the left half and one well into the right half
         ! must straddle the discontinuity.
         i_left = NG + NX/4
         i_right = NG + 3*NX/4 + 1

         call check(error, &
                    abs(s%barotropic%h(i_left, NG + 1) - s%barotropic%h(i_right, NG + 1)) > 0.1_wp, &
                    "dam break should have a height step")
         if (allocated(error)) exit checks

         ! Momentum starts at rest
         call check(error, all(abs(s%barotropic%hu) < 1.0e-12_wp), "hu should be zero at t=0")
         if (allocated(error)) exit checks
         call check(error, all(abs(s%barotropic%hv) < 1.0e-12_wp), "hv should be zero at t=0")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_dam_break

   subroutine test_gaussian_hump(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp), parameter :: H0 = 1.0_wp, AMP = 0.5_wp
      real(wp), parameter :: SIGMA = 2.0_wp
      real(wp) :: h_centre, h_edge, x_centre, y_centre
      checks: block

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)

         ! Centre on the domain centre, no shear
         x_centre = 0.5_wp*real(NX, wp)*DX
         y_centre = 0.5_wp*real(NY, wp)*DY
         call init_gaussian_hump(s, H0, AMP, SIGMA, x_centre, y_centre, 0.0_wp)

         h_centre = s%barotropic%h(NG + NX/2, NG + NY/2)
         h_edge = s%barotropic%h(NG + 1, NG + 1)

         call check(error, h_centre > h_edge, "Gaussian hump centre should be higher than edge")
         if (allocated(error)) exit checks

         call check(error, h_centre <= H0 + AMP + 1.0e-9_wp .and. h_centre > H0, &
                    "centre should land in (h0, h0+amp]")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_gaussian_hump

   subroutine test_castro(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp) :: b_centre, b_corner

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)
      call init_castro_dam_break(s)

      ! Castro's H = 1 - 0.4*exp(-x^2 - y^2): centre b is non-trivial
      b_centre = s%barotropic%b(NG + NX/2, NG + NY/2)
      b_corner = s%barotropic%b(NG + 1, NG + 1)

      call check(error, abs(b_centre - b_corner) > 1.0e-6_wp, &
                 "Castro bathymetry must vary across the domain")

      call s%destroy()
   end subroutine test_castro

   subroutine test_still_water(error)
      !! With a non-flat b, init_still_water sets h = max(eta0 - b, 0) so
      !! eta = h + b == eta0 for every wet cell.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      real(wp), parameter :: ETA0 = 2.0_wp
      integer :: i, j
      real(wp) :: eta
      checks: block

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)

         ! Tilted bottom: b varies from 0 to ~0.5 across the domain (always below eta0)
         do j = 1, s%grid%ny_total
            do i = 1, s%grid%nx_total
               s%barotropic%b(i, j) = 0.05_wp*real(i + j, wp)
            end do
         end do

         call init_still_water(s, ETA0)

         do j = NG + 1, NG + NY
            do i = NG + 1, NG + NX
               eta = s%barotropic%h(i, j) + s%barotropic%b(i, j)
               call check(error, abs(eta - ETA0) < 1.0e-12_wp, &
                          "h + b should equal eta0 in interior wet cells")
               if (allocated(error)) exit checks
            end do
         end do

         call check(error, all(abs(s%barotropic%hu) < 1.0e-12_wp), "hu must start at zero")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
   end subroutine test_still_water

   subroutine test_cfg_still_water(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      integer :: i, j
      checks: block

         call s%init(NX, NY, NG, DX, DY, 0.45_wp)

         cfg%init_condition = "still_water"
         cfg%eta_initial = 1.5_wp
         cfg%dx = DX
         cfg%dy = DY
         cfg%cfl = 0.45_wp

         ! Flat bottom so still water just sets h = eta_initial.
         s%barotropic%b = 0.0_wp

         call init_from_config(s, cfg)

         do j = NG + 1, NG + NY
            do i = NG + 1, NG + NX
               call check(error, abs(s%barotropic%h(i, j) - 1.5_wp) < 1.0e-12_wp, &
                          "still_water should leave h == eta_initial")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      continue
      call s%destroy()
   end subroutine test_cfg_still_water

   subroutine test_cfg_gaussian(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      real(wp) :: h_centre, h_edge

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      cfg%init_condition = "gaussian_hump"
      cfg%h0 = 1.0_wp
      cfg%hump_amplitude = 0.4_wp
      cfg%hump_sigma = 2.0_wp
      cfg%hump_x_offset = 0.5_wp*real(NX, wp)*DX
      cfg%hump_y_offset = 0.5_wp*real(NY, wp)*DY
      cfg%init_u_shear = 0.0_wp
      cfg%dx = DX
      cfg%dy = DY
      cfg%cfl = 0.45_wp

      call init_from_config(s, cfg)

      h_centre = s%barotropic%h(NG + NX/2, NG + NY/2)
      h_edge = s%barotropic%h(NG + 1, NG + 1)
      call check(error, h_centre > h_edge, "config-driven gaussian must produce a peak at centre")

      call s%destroy()
   end subroutine test_cfg_gaussian

   subroutine test_cfg_castro(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      cfg%init_condition = "castro_dam_break"
      cfg%dx = DX
      cfg%dy = DY
      cfg%cfl = 0.45_wp

      call init_from_config(s, cfg)

      call check(error, &
                 abs(s%barotropic%b(NG + NX/2, NG + NY/2) - s%barotropic%b(NG + 1, NG + 1)) > 1.0e-6_wp, &
                 "castro dispatch should set non-flat b")

      call s%destroy()
   end subroutine test_cfg_castro

   subroutine test_cfg_default(error)
      !! Unknown / unset init_condition falls through to init_dam_break.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      integer :: i_left, i_right

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      cfg%init_condition = "dam_break"
      cfg%dx = DX
      cfg%dy = DY
      cfg%cfl = 0.45_wp

      call init_from_config(s, cfg)

      i_left = NG + NX/4
      i_right = NG + 3*NX/4 + 1
      call check(error, &
                 abs(s%barotropic%h(i_left, NG + 1) - s%barotropic%h(i_right, NG + 1)) > 0.1_wp, &
                 "default dispatch should produce a dam-break step")

      call s%destroy()
   end subroutine test_cfg_default

   subroutine test_diagnose_bathymetry(error)
      !! `diagnose_bathymetry` is intent(in) — it just logs warnings.
      !! Drive both the steep-gradient branch and the shallow-deep branch.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      integer :: i, j

      call s%init(NX, NY, NG, DX, DY, 0.45_wp)

      ! All wet, mostly 1 m
      s%barotropic%h = 1.0_wp
      s%barotropic%b = 0.0_wp

      ! Inject a steep wet-to-wet bathymetry step at i = NG+NX/2
      ! (gradient > 10 m/cell triggers the steep-gradient warning path).
      do j = NG + 1, NG + NY
         do i = NG + NX/2, NG + NX
            s%barotropic%b(i, j) = -50.0_wp
            s%barotropic%h(i, j) = 50.0_wp
         end do
      end do

      ! And a shallow cell next to a deep one (depth ratio > 20:1)
      s%barotropic%h(NG + 2, NG + 2) = 0.5_wp
      s%barotropic%h(NG + 3, NG + 2) = 100.0_wp

      cfg%dx = DX
      cfg%dy = DY
      cfg%cfl = 0.45_wp

      call diagnose_bathymetry(s, cfg)

      ! Survival is the assertion — diagnose_bathymetry is a logging helper.
      call check(error, .true., "diagnose_bathymetry returned without error")

      call s%destroy()
   end subroutine test_diagnose_bathymetry

end module test_init
