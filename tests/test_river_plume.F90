!! Analytical test: buoyant river plume spreading
!!
!! A rectangular channel is initialised with freshwater (S=0) on the left
!! 20% and saltwater (S=35 PSU) on the right 80%.  The fresh water is
!! buoyant and spreads rightward as a surface-trapped plume over the denser
!! saltwater.  Dense saltwater undercuts leftward along the bottom.
!!
!! Unlike the lock exchange test (which validates horizontal gravity current
!! speed), this test focuses on vertical structure:
!!   1. Surface trapping: the plume stays in the upper layers
!!   2. Pycnocline mixing: PP81 vertical mixing creates intermediate salinity
!!      at the interface between fresh and salt layers
!!   3. Salt conservation: total hS is preserved through transport + mixing
!!
!! This validates the full coupled machinery:
!!   salinity → EOS → BPG → differential layer flow → advection → PP mixing
module test_river_plume
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_river_plume_tests

   ! Domain parameters
   real(wp), parameter :: H0 = 1.0_wp        !! Undisturbed depth (m)
   real(wp), parameter :: DX_VAL = 0.2_wp    !! Grid spacing (m)
   integer, parameter :: NX_PHYS = 100       !! 20 m channel
   integer, parameter :: NY_PHYS = 4         !! Narrow in y
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 5              !! Five sigma layers

   ! Density / salinity parameters
   real(wp), parameter :: S_OCEAN = 35.0_wp  !! Ambient ocean salinity (PSU)
   real(wp), parameter :: S_RIVER = 0.0_wp   !! Freshwater salinity (PSU)
   real(wp), parameter :: BETA_S = 0.78_wp   !! Haline contraction (kg/m^3/PSU)
   real(wp), parameter :: RHO_0 = 1000.0_wp  !! Reference density

   ! Freshwater region: left 20% of physical domain
   integer, parameter :: I_INTERFACE = NX_PHYS/5

contains

   subroutine collect_river_plume_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("plume_surface_trapped", test_plume_surface_trapped), &
                  new_unittest("pycnocline_mixing", test_pycnocline_mixing), &
                  new_unittest("total_salt_conserved", test_total_salt_conserved) &
                  ]
   end subroutine collect_river_plume_tests

   subroutine init_river_plume(state)
      !! Fresh water on left 20%, salt water on right 80%, at rest
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_iface

      ! Flat bottom, uniform depth, at rest
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            state%barotropic%b(i, j) = 0.0_wp
         end do
      end do

      ! Initialise layers from barotropic (uniform sigma)
      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Set salinity: fresh on left, salty on right
      i_iface = NGHOST + I_INTERFACE
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (i <= i_iface) then
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_RIVER
               else
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_OCEAN
               end if
            end do
         end do
      end do

      ! Compute initial density from EOS (temperature disabled here: alpha_T=0)
      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, 0.0_wp, &
                                0.0_wp, 15.0_wp, &
                                state%grid, NZ)

   end subroutine init_river_plume

   subroutine setup_state(state)
      !! Create and initialise the river plume state
      type(state_t), intent(inout) :: state

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_OCEAN
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0
      state%multilayer%alpha_T = 0.0_wp   ! decouple temperature from dynamics for this fixture

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Wall BCs on all sides (closed channel)
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      call init_river_plume(state)
   end subroutine setup_state

   subroutine run_plume(state, t_final)
      !! Run the solver up to t_final
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: t_final

      real(wp) :: t

      call solver_enter_data(state)

      t = 0.0_wp
      do while (t < t_final)
         call solver_step(state, time=t)
         t = t + state%dt
      end do

      call solver_exit_data(state)
   end subroutine run_plume

   function compute_total_salt(state) result(total)
      !! Integrate hS over all layers and interior cells
      type(state_t), intent(in) :: state
      real(wp) :: total

      integer :: i, j, k

      total = 0.0_wp
      do k = 1, NZ
         do j = NGHOST + 1, state%grid%ny_total - NGHOST
            do i = NGHOST + 1, state%grid%nx_total - NGHOST
               total = total + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k)
            end do
         end do
      end do
      total = total*DX_VAL*DX_VAL
   end function compute_total_salt

   subroutine test_plume_surface_trapped(error)
      !! After the buoyant plume spreads, the top layer in the salty region
      !! should have reduced salinity (fresh water arrived) while the bottom
      !! layer should still be mostly salty (plume is surface-trapped).
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: g_prime, u_theory, t_final, L
      real(wp) :: S_top, S_bot, h_top, h_bot
      integer :: i_probe, j_mid

      call setup_state(state)

      ! Reduced gravity and internal wave speed
      g_prime = GRAVITY*BETA_S*(S_OCEAN - S_RIVER)/RHO_0
      u_theory = 0.5_wp*sqrt(g_prime*H0)

      ! Run long enough for the plume to spread well past the interface
      ! Front travels ~u_theory * t_final; probe is at 30% (interface at 20%)
      L = NX_PHYS*DX_VAL
      t_final = 0.30_wp*L/u_theory

      call run_plume(state, t_final)

      ! Probe point: 30% into the domain (10% right of the 20% interface)
      i_probe = NGHOST + (NX_PHYS*3)/10
      j_mid = NGHOST + NY_PHYS/2

      ! Top layer salinity (k=NZ is top)
      h_top = state%multilayer%h_layer(i_probe, j_mid, NZ)
      if (h_top > 1.0e-6_wp) then
         S_top = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i_probe, j_mid, NZ)/h_top
      else
         S_top = S_OCEAN
      end if

      ! Bottom layer salinity (k=1 is bottom)
      h_bot = state%multilayer%h_layer(i_probe, j_mid, 1)
      if (h_bot > 1.0e-6_wp) then
         S_bot = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i_probe, j_mid, 1)/h_bot
      else
         S_bot = S_OCEAN
      end if

      ! Top layer should have reduced salinity (fresh plume arrived)
      call check(error, S_top < 30.0_wp, &
                 "Top layer should have reduced salinity from plume")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Bottom layer should still be mostly salty (plume is surface-trapped)
      call check(error, S_bot > 30.0_wp, &
                 "Bottom layer should remain salty (surface-trapped plume)")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call state%destroy()
   end subroutine test_plume_surface_trapped

   subroutine test_pycnocline_mixing(error)
      !! At the plume front, PP mixing should create intermediate salinity
      !! in the mid-depth layers.  Without mixing, the salinity jump would
      !! be perfectly sharp (only top or bottom layer has anomalous S).
      !! With PP mixing active, at least one mid-layer should have
      !! 2 < S < 33 (well between fresh and ambient).
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: g_prime, u_theory, t_final, L
      real(wp) :: S_k, h_k
      integer :: i_probe, j_mid, k
      logical :: found_intermediate

      call setup_state(state)

      g_prime = GRAVITY*BETA_S*(S_OCEAN - S_RIVER)/RHO_0
      u_theory = 0.5_wp*sqrt(g_prime*H0)
      L = NX_PHYS*DX_VAL
      t_final = 0.30_wp*L/u_theory

      call run_plume(state, t_final)

      ! Probe just right of the initial interface (25% into domain, interface at 20%)
      ! where the plume has been active longest and mixing has had most time
      i_probe = NGHOST + (NX_PHYS*25)/100
      j_mid = NGHOST + NY_PHYS/2

      ! Search mid-depth layers for intermediate salinity
      found_intermediate = .false.
      do k = 2, NZ - 1  ! exclude top and bottom layer
         h_k = state%multilayer%h_layer(i_probe, j_mid, k)
         if (h_k > 1.0e-6_wp) then
            S_k = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i_probe, j_mid, k)/h_k
            if (S_k > 2.0_wp .and. S_k < 33.0_wp) then
               found_intermediate = .true.
               exit
            end if
         end if
      end do

      call check(error, found_intermediate, &
                 "Mid-depth layers should have intermediate salinity from PP mixing")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call state%destroy()
   end subroutine test_pycnocline_mixing

   subroutine test_total_salt_conserved(error)
      !! Total salt (integral of hS over domain) should be conserved
      !! through the transport and mixing operations.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: g_prime, u_theory, t_final, L
      real(wp) :: salt_before, salt_after, rel_error

      call setup_state(state)

      salt_before = compute_total_salt(state)

      g_prime = GRAVITY*BETA_S*(S_OCEAN - S_RIVER)/RHO_0
      u_theory = 0.5_wp*sqrt(g_prime*H0)
      L = NX_PHYS*DX_VAL
      t_final = 0.30_wp*L/u_theory

      call run_plume(state, t_final)

      salt_after = compute_total_salt(state)

      ! Relative error in total salt
      rel_error = abs(salt_after - salt_before)/salt_before

      call check(error, rel_error < 0.01_wp, &
                 "Total salt should be conserved within 1%")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call state%destroy()
   end subroutine test_total_salt_conserved

end module test_river_plume
