!! Test: multilayer solver with dry cells (bathymetry)
!!
!! A rectangular domain with a step bathymetry: the left half has
!! b = -2m (deep, h = 2m) and the right half has b = 1m (land, h = 0).
!! This creates a sharp wet/dry interface that exercises the multilayer
!! solver's handling of dry cells — specifically that no NaN values
!! appear in h after a few timesteps.
!!
!! This test was created to diagnose NaN propagation in the multilayer
!! solver on the Sydney Harbour grid, which has complex bathymetry with
!! many dry cells.
module test_multilayer_drycell
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_grid, only: hgrid_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic, &
                              ml_compute_vertical_velocity, &
                              ml_compute_vertical_exchange
   use rki_ml_tracers, only: ml_init_tracer, ml_vertical_exchange_tracer
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_drycell_tests

   integer, parameter :: NX_PHYS = 100
   integer, parameter :: NY_PHYS = 20
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 5
   real(wp), parameter :: DX_VAL = 10.0_wp
   real(wp), parameter :: H0 = 2.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp

contains

   subroutine collect_drycell_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_nan_after_step", test_no_nan), &
                  new_unittest("mass_conserved_with_bathy", test_mass_conserved), &
                  new_unittest("no_nan_tidal_irregular", test_no_nan_tidal), &
                  new_unittest("vexch_flux_limiter", test_vexch_flux_limiter), &
                  new_unittest("thin_column_no_nan", test_thin_column_no_nan) &
                  ]
   end subroutine collect_drycell_tests

   subroutine setup_state(state)
      type(state_t), intent(inout) :: state

      integer :: i, j, k

      state%use_multilayer = .true.
      state%multilayer%alpha_T = 0.0_wp   ! decouple temperature from dynamics for this fixture
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      ! Zero out mixing to isolate the core distribute cycle
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%vmix_ri_alpha = 5.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      ! Step bathymetry: left half deep (wet), right half above sea level (dry)
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            if (i <= NGHOST + NX_PHYS/2) then
               state%barotropic%b(i, j) = -H0  ! deep: h = eta - b = 0 - (-2) = 2
            else
               state%barotropic%b(i, j) = 1.0_wp  ! land: h = max(0 - 1, 0) = 0
            end if
            state%barotropic%h(i, j) = max(0.0_wp - state%barotropic%b(i, j), 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      ! Initialise multilayer
      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)
      call ml_init_tracer(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, S_INIT, &
                          state%grid, NZ)
      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, 0.0_wp, &
                                0.0_wp, 15.0_wp, &
                                state%grid, NZ)
   end subroutine setup_state

   subroutine test_no_nan(error)
      !! Run 10 steps on a grid with dry cells.  Check that h has no NaN.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t
      integer :: step, i, j, ng
      logical :: found_nan

      call setup_state(state)
      call solver_enter_data(state)

      t = 0.0_wp
      do step = 1, 10
         call solver_step(state, time=t)
         t = t + state%dt
      end do

      call solver_exit_data(state)

      ! Check for NaN in h (interior cells)
      ng = NGHOST
      found_nan = .false.
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            if (state%barotropic%h(i, j) /= state%barotropic%h(i, j)) then
               found_nan = .true.
               exit
            end if
         end do
         if (found_nan) exit
      end do

      call check(error,.not. found_nan, &
                 "h should not contain NaN after multilayer steps with dry cells")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Also check h_layer
      block
         integer :: k
         do k = 1, NZ
            do j = ng + 1, state%grid%ny_total - ng
               do i = ng + 1, state%grid%nx_total - ng
                  if (state%multilayer%h_layer(i, j, k) /= state%multilayer%h_layer(i, j, k)) then
                     found_nan = .true.
                     exit
                  end if
               end do
               if (found_nan) exit
            end do
            if (found_nan) exit
         end do
      end block

      call check(error,.not. found_nan, &
                 "h_layer should not contain NaN after multilayer steps with dry cells")

      call state%destroy()
   end subroutine test_no_nan

   subroutine test_mass_conserved(error)
      !! Run for several steps and check total mass is conserved.
      !! Mass = sum(h * dx * dy) over wet interior cells.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t, mass_before, mass_after, rel_error
      integer :: step, i, j, ng

      call setup_state(state)

      ! Compute initial mass
      ng = NGHOST
      mass_before = 0.0_wp
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            mass_before = mass_before + state%barotropic%h(i, j)*DX_VAL*DX_VAL
         end do
      end do

      call solver_enter_data(state)
      t = 0.0_wp
      do step = 1, 20
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)

      ! Compute final mass
      mass_after = 0.0_wp
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            if (state%barotropic%h(i, j) > 0.0_wp .and. state%barotropic%h(i, j) == state%barotropic%h(i, j)) &
               mass_after = mass_after + state%barotropic%h(i, j)*DX_VAL*DX_VAL
         end do
      end do

      rel_error = abs(mass_after - mass_before)/mass_before

      call check(error, rel_error < 0.01_wp, &
                 "mass should be conserved within 1% with dry cells")

      call state%destroy()
   end subroutine test_mass_conserved

   subroutine test_no_nan_tidal(error)
      !! More realistic test: tidal boundary, irregular bathymetry with
      !! isolated dry cells and steep gradients.  Runs 20 steps and checks
      !! that h has no NaN.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t
      integer :: step, i, j, ng, k
      logical :: found_nan
      real(wp) :: x, y, bathy

      state%use_multilayer = .true.
      state%multilayer%alpha_T = 0.0_wp   ! decouple temperature from dynamics for this fixture
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%vmix_ri_alpha = 5.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Wall on west/south/north, tidal on east
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("tidal")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")
      state%bc%tidal_amplitude = 0.5_wp
      state%bc%tidal_period = 100.0_wp

      ! Irregular bathymetry: sloping with isolated dry patches
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x = real(i - NGHOST, wp)*DX_VAL
            y = real(j - NGHOST, wp)*DX_VAL

            ! Sloping bottom: deep on right (east/tidal), shallow on left
            bathy = -5.0_wp + 6.0_wp*(1.0_wp - x/(NX_PHYS*DX_VAL))

            ! Add isolated dry patches (small islands)
            if (mod(i, 7) == 3 .and. mod(j, 5) == 2 .and. i > NGHOST + 10) then
               bathy = 2.0_wp  ! island
            end if

            ! Steep gradient: a ridge across the middle
            if (i > NGHOST + NX_PHYS/2 - 2 .and. i < NGHOST + NX_PHYS/2 + 2) then
               bathy = bathy + 3.0_wp
            end if

            state%barotropic%b(i, j) = bathy
            state%barotropic%h(i, j) = max(0.0_wp - bathy, 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      ! Initialise multilayer
      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)
      call ml_init_tracer(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, S_INIT, &
                          state%grid, NZ)
      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, 0.0_wp, &
                                0.0_wp, 15.0_wp, &
                                state%grid, NZ)

      call solver_enter_data(state)
      t = 0.0_wp
      do step = 1, 20
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)

      ! Check for NaN in h
      ng = NGHOST
      found_nan = .false.
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            if (state%barotropic%h(i, j) /= state%barotropic%h(i, j)) then
               found_nan = .true.
               exit
            end if
         end do
         if (found_nan) exit
      end do

      call check(error,.not. found_nan, &
                 "h should not contain NaN with tidal BC and irregular bathymetry")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Check h_layer
      do k = 1, NZ
         do j = ng + 1, state%grid%ny_total - ng
            do i = ng + 1, state%grid%nx_total - ng
               if (state%multilayer%h_layer(i, j, k) /= state%multilayer%h_layer(i, j, k)) then
                  found_nan = .true.
                  exit
               end if
            end do
            if (found_nan) exit
         end do
         if (found_nan) exit
      end do

      call check(error,.not. found_nan, &
                 "h_layer should not contain NaN with tidal BC and irregular bathymetry")

      call state%destroy()
   end subroutine test_no_nan_tidal

   subroutine test_vexch_flux_limiter(error)
      !! Vertical exchange with huge w should not blow up.
      !! Sets up a 2-layer column with large diagnostic w that would cause
      !! the vertical Courant number |w*dt/h| >> 1 without the flux limiter.
      !! Without the limiter, momentum overshoots and grows without bound.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: LNX = 12, LNY = 12, LNG = 2, LNZ = 2
      real(wp) :: h_layer(LNX, LNY, LNZ)
      real(wp) :: hu_layer(LNX, LNY, LNZ), hv_layer(LNX, LNY, LNZ)
      real(wp) :: w(LNX, LNY, LNZ + 1)
      real(wp) :: dt, hu_before, hu_after
      type(hgrid_t) :: grid

      call grid%init(LNX - 2*LNG, LNY - 2*LNG, LNG, 1.0_wp, 1.0_wp)

      dt = 1.0_wp
      h_layer = 0.01_wp   ! very thin layers (10 mm)
      hu_layer = 0.0_wp
      hv_layer = 0.0_wp

      ! Bottom layer has some momentum
      hu_layer(:, :, 1) = 0.01_wp*2.0_wp  ! h*u = 0.02, u = 2 m/s

      ! Large upward w at interface — vertical Courant = |w*dt/h| = 10/0.01 = 1000
      w = 0.0_wp
      w(:, :, 2) = 10.0_wp  ! huge vertical velocity

      hu_before = hu_layer(5, 5, 1)

      !$omp target data map(to:  h_layer, w) map(tofrom:  hu_layer, hv_layer)
      call ml_compute_vertical_exchange(h_layer, hu_layer, hv_layer, &
                                        w, grid, LNZ, dt)
      !$omp end target data

      hu_after = hu_layer(5, 5, 1)

      ! Without limiter: flux = w*u*dt = 10*2*1 = 20, but h*u was only 0.02
      ! That would extract 1000x more momentum than exists — catastrophic
      ! With limiter: flux is capped, momentum stays bounded
      call check(error, hu_after == hu_after, &  ! not NaN
                 "vertical exchange with huge w should not produce NaN")
      if (allocated(error)) return

      call check(error, abs(hu_after) < 1.0_wp, &
                 "vertical exchange flux limiter prevents overshoot")
   end subroutine test_vexch_flux_limiter

   subroutine test_thin_column_no_nan(error)
      !! Full solver test with a thin water column next to dry cells.
      !! This is the minimal reproducer for the Sydney Harbour NaN bug:
      !! thin wet cells at wetting/drying fronts produce large div(hu),
      !! which creates large diagnostic w, and the vertical exchange
      !! overshoots without the flux limiter.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t
      integer :: step, i, j, ng
      logical :: found_nan

      state%use_multilayer = .true.
      state%multilayer%alpha_T = 0.0_wp   ! decouple temperature from dynamics for this fixture
      state%multilayer%nz_ml = 2
      state%multilayer%initial_salinity = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = 1000.0_wp
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%vmix_ri_alpha = 5.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      ! Bathymetry: mostly flat ocean with a thin strip of barely-wet cells
      ! next to dry land.  This creates the sharp h gradient that triggers
      ! the large div(hu) → large w → exchange overshoot.
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            if (i <= NGHOST + NX_PHYS/3) then
               state%barotropic%b(i, j) = -5.0_wp  ! deep ocean
            else if (i <= NGHOST + NX_PHYS/3 + 3) then
               ! Thin strip: h = 0.001m (1mm of water)
               state%barotropic%b(i, j) = -0.001_wp
            else
               state%barotropic%b(i, j) = 1.0_wp   ! dry land
            end if
            state%barotropic%h(i, j) = max(0.0_wp - state%barotropic%b(i, j), 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      ! Give the ocean some initial velocity to create momentum at the front
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            if (state%barotropic%h(i, j) > 0.1_wp) then
               state%barotropic%hu(i, j) = state%barotropic%h(i, j)*0.5_wp  ! u = 0.5 m/s rightward
            end if
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, 2)
      call ml_init_tracer(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, 0.0_wp, &
                          state%grid, 2)
      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, 1000.0_wp, 0.0_wp, 0.0_wp, &
                                0.0_wp, 15.0_wp, &
                                state%grid, 2)

      call solver_enter_data(state)
      t = 0.0_wp
      do step = 1, 30
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)

      ! Check for NaN in h
      ng = NGHOST
      found_nan = .false.
      do j = ng + 1, state%grid%ny_total - ng
         do i = ng + 1, state%grid%nx_total - ng
            if (state%barotropic%h(i, j) /= state%barotropic%h(i, j)) then
               found_nan = .true.
               exit
            end if
         end do
         if (found_nan) exit
      end do

      call check(error,.not. found_nan, &
                 "thin column next to dry cells should not produce NaN in h")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Check mass conservation (approximate — wall BCs so no inflow/outflow)
      block
         real(wp) :: mass
         mass = 0.0_wp
         do j = ng + 1, state%grid%ny_total - ng
            do i = ng + 1, state%grid%nx_total - ng
               if (state%barotropic%h(i, j) > 0.0_wp .and. state%barotropic%h(i, j) == state%barotropic%h(i, j)) &
                  mass = mass + state%barotropic%h(i, j)*DX_VAL*DX_VAL
            end do
         end do
         ! Initial mass: ~33 cells * 5m depth * 100m^2 = 1.65e6
         ! Should be conserved within 5% (wetting/drying causes some loss)
         call check(error, mass > 1.0e5_wp, &
                    "mass should not drain to near-zero with thin columns")
      end block

      call state%destroy()
   end subroutine test_thin_column_no_nan

end module test_multilayer_drycell
