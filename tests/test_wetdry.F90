!! Unit tests for wet/dry handling across kernels
module test_wetdry
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, DESING_EPS, &
                            THIN_LAYER_THRESHOLD
   use rki_kernel_update, only: compute_update
   use rki_kernel_flux, only: compute_flux_hll
   use rki_grid, only: hgrid_t
   use rki_kernel_timestep, only: compute_timestep
   implicit none
   private

   public :: collect_wetdry_tests

contains

   subroutine collect_wetdry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("update_desingularizes_thin_momentum", test_update_thin_layer), &
                  new_unittest("update_zeroes_dry_cell", test_update_dry_cell), &
                  new_unittest("update_preserves_wet_cell", test_update_preserves_wet), &
                  new_unittest("flux_near_dry_front", test_flux_near_dry_front), &
                  new_unittest("timestep_wet_dry_mix", test_timestep_wet_dry_mix), &
                  new_unittest("step_does_not_shrink_dt", test_step_no_dt_shrink) &
                  ]
   end subroutine collect_wetdry_tests

   subroutine test_update_thin_layer(error)
      !! Thin-but-wet cells (DRY_TOLERANCE <= h < eps) keep a Kurganov &
      !! Petrova (2007) desingularized momentum
      !!     hu := h * sqrt(2) h hu / sqrt(h^4 + eps^4)
      !! instead of the old hard zero — thin water must stay mobile (the
      !! hard zero froze wet/dry fronts, see test_dam_break_dry).
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: eps, h_thin, hu_expect, hv_expect, scale
      type(hgrid_t) :: grid

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      flux_h = 0.0_wp
      flux_hu = 0.0_wp
      flux_hv = 0.0_wp

      ! eps = DESING_EPS.  Momenta chosen so the desingularized speed
      ! (~1 m/s) stays below the FROUDE_BACKSTOP bound (~2.7 m/s at this
      ! depth), isolating the KP07 formula.
      eps = DESING_EPS
      h_thin = 0.5_wp*DESING_EPS
      h(4, 4) = h_thin
      hu(4, 4) = 2.0e-4_wp
      hv(4, 4) = 1.0e-4_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copy( h, hu, hv) copyin( flux_h, flux_hu, flux_hv)
      call compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, grid, 0.1_wp)
      !$acc end data

      scale = sqrt(2.0_wp)*h_thin**2/sqrt(h_thin**4 + eps**4)
      hu_expect = scale*2.0e-4_wp
      hv_expect = scale*1.0e-4_wp
      call check(error, abs(hu(4, 4) - hu_expect) < 1.0e-15_wp, &
                 "hu must equal the KP07 desingularized value")
      if (allocated(error)) return
      call check(error, abs(hv(4, 4) - hv_expect) < 1.0e-15_wp, &
                 "hv must equal the KP07 desingularized value")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4)) < 2.0e-4_wp, &
                 "desingularization must reduce thin-cell momentum")
      if (allocated(error)) return
      call check(error, h(4, 4) > 0.0_wp, &
                 "h should be preserved (above DRY_TOLERANCE)")
   end subroutine test_update_thin_layer

   subroutine test_update_dry_cell(error)
      !! Small positive h is preserved for mass conservation;
      !! momentum is zeroed since h < DRY_TOLERANCE.
      !! Only negative h (undershoots) get clamped to zero.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: h_original
      type(hgrid_t) :: grid

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      flux_h = 0.0_wp
      flux_hu = 0.0_wp
      flux_hv = 0.0_wp

      h(4, 4) = 0.5_wp*DRY_TOLERANCE
      h_original = h(4, 4)
      hu(4, 4) = 100.0_wp
      hv(4, 4) = 50.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copy( h, hu, hv) copyin( flux_h, flux_hu, flux_hv)
      call compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, grid, 0.1_wp)
      !$acc end data

      call check(error, abs(h(4, 4) - h_original) < 1.0e-15_wp, &
                 "Small positive h preserved for mass conservation")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4)) < 1.0e-15_wp, &
                 "hu should be zeroed in thin-layer cell")
   end subroutine test_update_dry_cell

   subroutine test_update_preserves_wet(error)
      !! Wet cells above thin-layer threshold should keep their momentum
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 8, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      type(hgrid_t) :: grid

      h = 1.0_wp
      hu = 0.5_wp
      hv = 0.3_wp
      flux_h = 0.0_wp
      flux_hu = 0.0_wp
      flux_hv = 0.0_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copy( h, hu, hv) copyin( flux_h, flux_hu, flux_hv)
      call compute_update(h, hu, hv, flux_h, flux_hu, flux_hv, grid, 0.1_wp)
      !$acc end data

      call check(error, abs(hu(4, 4) - 0.5_wp) < 1.0e-14_wp, &
                 "Wet cell momentum should be unchanged with zero flux")
   end subroutine test_update_preserves_wet

   subroutine test_flux_near_dry_front(error)
      !! Flux at a wet/dry interface should not produce huge values
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 12, NY = 8, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY), b(NX, NY)
      real(wp) :: flux_h(NX, NY), flux_hu(NX, NY), flux_hv(NX, NY)
      real(wp) :: mass_flux_x(NX, NY), mass_flux_y(NX, NY)
      type(hgrid_t) :: grid
      integer :: i, j
      real(wp) :: max_flux

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      b = 0.0_wp
      mass_flux_x = 0.0_wp
      mass_flux_y = 0.0_wp

      ! Wet region on left, dry on right (dam break into dry bed)
      do j = 1, NY
         do i = 1, NX
            if (i <= NX/2) then
               h(i, j) = 1.0_wp
            else
               h(i, j) = 0.0_wp
            end if
         end do
      end do

      call grid%init(NX - 2*NG, NY - 2*NG, NG, 1.0_wp, 1.0_wp)

      !$acc data copyin( h, hu, hv, b) &
      !$acc&     copyout( flux_h, flux_hu, flux_hv, mass_flux_x, mass_flux_y)
      call compute_flux_hll(h, hu, hv, b, flux_h, flux_hu, flux_hv, &
                            mass_flux_x, mass_flux_y, grid)
      !$acc end data

      ! Check that fluxes are bounded
      max_flux = 0.0_wp
      do j = NG + 1, NY - NG
         do i = NG + 1, NX - NG
            max_flux = max(max_flux, abs(flux_h(i, j)))
            max_flux = max(max_flux, abs(flux_hu(i, j)))
            max_flux = max(max_flux, abs(flux_hv(i, j)))
         end do
      end do

      call check(error, max_flux < 100.0_wp, &
                 "Fluxes at wet/dry front should be bounded")
   end subroutine test_flux_near_dry_front

   subroutine test_timestep_wet_dry_mix(error)
      !! dt should be reasonable with a mix of wet, thin, and dry cells
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX = 14, NY = 14, NG = 2
      real(wp) :: h(NX, NY), hu(NX, NY), hv(NX, NY)
      real(wp) :: dt
      real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp
      real(wp) :: dt_expected_min
      type(hgrid_t) :: grid

      h = 0.0_wp
      hu = 0.0_wp
      hv = 0.0_wp

      ! Wet cells with moderate flow
      h(4:8, 4:8) = 1.0_wp
      hu(4:8, 4:8) = 2.0_wp

      ! Thin-layer cells around the wet region
      h(3, 4:8) = 0.5_wp*THIN_LAYER_THRESHOLD
      h(9, 4:8) = 0.5_wp*THIN_LAYER_THRESHOLD
      hu(3, 4:8) = 0.01_wp
      hu(9, 4:8) = 0.01_wp

      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)

      !$acc data copyin( h, hu, hv)
      call compute_timestep(h, hu, hv, grid, CFL, dt)
      !$acc end data

      ! dt should be set by the wet cells: CFL / (u + sqrt(g*h)) * dx
      dt_expected_min = CFL*DX/(2.0_wp + sqrt(GRAVITY*1.0_wp))

      call check(error, dt > 0.01_wp*dt_expected_min, &
                 "dt should not be driven to zero by thin-layer cells")
      if (allocated(error)) return
      call check(error, dt < 10.0_wp*dt_expected_min, &
                 "dt should be in reasonable range")
   end subroutine test_timestep_wet_dry_mix

   subroutine test_step_no_dt_shrink(error)
      !! Running a few steps on a wet/dry problem should not shrink dt to zero
      use rki_state, only: state_t
      use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j, n
      real(wp) :: dt_min

      ! Small grid with a wet/dry front
      call state%init(20, 20, 2, 1.0_wp, 1.0_wp, 0.4_wp)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%b(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            if (i <= state%grid%nghost + 10) then
               state%barotropic%h(i, j) = 1.0_wp
               state%barotropic%hu(i, j) = 0.0_wp
            else
               state%barotropic%h(i, j) = 0.0_wp
               state%barotropic%hu(i, j) = 0.0_wp
            end if
         end do
      end do

      call solver_enter_data(state)

      dt_min = 1.0e10_wp
      do n = 1, 20
         call solver_step(state)
         dt_min = min(dt_min, state%dt)
      end do

      call solver_exit_data(state)

      call check(error, dt_min > 1.0e-4_wp, &
                 "dt should not shrink to near-zero over 20 steps")

      call state%destroy()
   end subroutine test_step_no_dt_shrink

end module test_wetdry
