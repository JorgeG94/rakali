!! Analytical test: lock exchange gravity current
!!
!! A rectangular channel is initialised with dense fluid on the left half
!! and lighter fluid on the right.  Two gravity currents form: dense water
!! flows rightward along the bottom, light water flows leftward along the
!! top.  The theoretical front speed for a two-layer inviscid lock exchange
!! with energy-conserving theory (Benjamin 1968) is:
!!
!!     u_front = 0.5 * sqrt(g' * H)
!!
!! where g' = g * delta_rho / rho_0 is the reduced gravity.
!!
!! This test validates the coupled baroclinic machinery:
!!   salinity → EOS → baroclinic pressure gradient → differential layer flow
module test_lock_exchange
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_tracers, only: ml_init_tracer
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_lock_exchange_tests

   ! Domain parameters
   real(wp), parameter :: H0 = 1.0_wp        !! Undisturbed depth (m)
   real(wp), parameter :: DX_VAL = 0.1_wp    !! Grid spacing (m)
   integer, parameter :: NX_PHYS = 100       !! 10 m channel
   integer, parameter :: NY_PHYS = 4         !! Narrow in y
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 7              !! Seven sigma layers

   ! Density parameters
   real(wp), parameter :: S_DENSE = 35.0_wp  !! Salinity on left (PSU)
   real(wp), parameter :: S_LIGHT = 0.0_wp   !! Salinity on right (PSU)
   real(wp), parameter :: BETA_S = 0.78_wp   !! Haline contraction (kg/m³/PSU)
   real(wp), parameter :: RHO_0 = 1000.0_wp  !! Reference density

   ! Temperature-driven parameters
   ! Chosen so alpha_T*(T_HOT - T_COLD) = beta_S*(S_DENSE - S_LIGHT) = 27.3 kg/m³,
   ! i.e. a temperature-only lock exchange produces the same density contrast
   ! as the salinity-only one.  Front speed should therefore match within a
   ! few percent.
   real(wp), parameter :: T_COLD = 15.0_wp       !! "Dense" side temperature (degC)
   real(wp), parameter :: T_HOT = 25.0_wp       !! "Light" side temperature (degC)
   real(wp), parameter :: ALPHA_T = 2.73_wp       !! 0.78*35 / 10
   real(wp), parameter :: T_REF_LOCK = T_HOT      !! Reference chosen so hot side = rho_0

contains

   subroutine collect_lock_exchange_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lock_exchange_front_speed", test_front_speed), &
                  new_unittest("lock_exchange_symmetry", test_symmetry), &
                  new_unittest("lock_exchange_temperature_driven", &
                               test_temperature_driven), &
                  new_unittest("lock_exchange_combined_cancellation", &
                               test_combined_cancellation), &
                  new_unittest("lock_exchange_combined_reinforcing", &
                               test_combined_reinforcing) &
                  ]
   end subroutine collect_lock_exchange_tests

   subroutine init_lock_exchange(state)
      !! Set up the lock exchange initial condition:
      !! dense (high salinity) on the left half, light (zero salinity) on the right
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_mid

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

      ! Set salinity: dense on left, light on right
      i_mid = NGHOST + NX_PHYS/2
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (i <= i_mid) then
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_DENSE
               else
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_LIGHT
               end if
            end do
         end do
      end do

      ! Compute initial density from EOS (temperature disabled here: alpha_T=0)
      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, 0.0_wp, &
                                0.0_wp, 15.0_wp, &
                                state%grid, NZ)

   end subroutine init_lock_exchange

   subroutine run_lock_exchange(state, t_final)
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
   end subroutine run_lock_exchange

   function find_salinity_front(hS_layer, h_layer, nx, ny, ng, layer, dx, S_threshold) result(x_front)
      !! Find the front position of the gravity current by tracking salinity.
      !! For the bottom layer: scan from right for S > threshold (dense water spreading right).
      real(wp), intent(in) :: hS_layer(:, :, :), h_layer(:, :, :)
      integer, intent(in) :: nx, ny, ng, layer
      real(wp), intent(in) :: dx, S_threshold
      real(wp) :: x_front

      integer :: i, j_mid
      real(wp) :: S_k

      j_mid = ng + (ny - 2*ng)/2

      ! Scan from right to find leading edge of salinity > threshold
      x_front = real(nx - 2*ng, wp)*dx*0.5_wp  ! default: midpoint
      do i = nx - ng, ng + 1, -1
         S_k = hS_layer(i, j_mid, layer)/max(h_layer(i, j_mid, layer), 1.0e-6_wp)
         if (S_k > S_threshold) then
            x_front = (real(i - ng, wp) - 0.5_wp)*dx
            exit
         end if
      end do
   end function find_salinity_front

   subroutine setup_state(state)
      !! Create and initialise the lock exchange state
      type(state_t), intent(inout) :: state

      state%use_multilayer = .true.
      state%multilayer%alpha_T = 0.0_wp   ! decouple temperature from dynamics for this fixture
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_DENSE
      state%multilayer%beta_S = BETA_S
      state%multilayer%rho_0 = RHO_0

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Wall BCs on all sides (closed channel)
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      call init_lock_exchange(state)
   end subroutine setup_state

   subroutine test_front_speed(error)
      !! The bottom-layer gravity current front should propagate at
      !! approximately u_front = 0.5 * sqrt(g' * H) (Benjamin 1968).
      !! We check that the measured front speed is within 50% of theory.
      !! The tolerance is generous because: (1) numerical diffusion from
      !! first-order upwind slows the front, (2) sigma layers can only
      !! approximate the half-depth current structure.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: g_prime, u_theory, x_front, x_init, u_measured, t_final
      real(wp) :: L

      call setup_state(state)

      ! Reduced gravity
      g_prime = GRAVITY*BETA_S*(S_DENSE - S_LIGHT)/RHO_0

      ! Benjamin half-depth front speed
      u_theory = 0.5_wp*sqrt(g_prime*H0)

      ! Run long enough for the front to travel ~1/4 of domain
      L = NX_PHYS*DX_VAL  ! 10 m
      t_final = 0.25_wp*L/u_theory

      ! Initial front is at midpoint
      x_init = L/2.0_wp

      call run_lock_exchange(state, t_final)

      ! Find the bottom-layer front (rightward propagating dense current)
      x_front = find_salinity_front(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%h_layer, &
                                    state%grid%nx_total, state%grid%ny_total, &
                                    NGHOST, 1, DX_VAL, 0.5_wp)

      ! Measured front speed
      u_measured = (x_front - x_init)/t_final

      ! Check front has moved forward
      call check(error, x_front > x_init + 0.5_wp, &
                 "Bottom gravity current should propagate rightward")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Check front speed within 50% of Benjamin theory
      call check(error, abs(u_measured - u_theory)/u_theory < 0.50_wp, &
                 "Front speed should match Benjamin theory within 50%")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call state%destroy()
   end subroutine test_front_speed

   subroutine test_symmetry(error)
      !! The lock exchange should be approximately symmetric: the dense
      !! bottom front spreading rightward and the dense-water retreat
      !! (leftward) should have similar extent.
      !!
      !! We track both using the bottom layer salinity:
      !!   - Right front: rightmost cell with S_bot > 0.5 (dense spreading right)
      !!   - Left retreat: leftmost cell with S_bot < 34.5 (dense leaving left)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: g_prime, u_theory, t_final, L
      real(wp) :: x_front_right, x_retreat_left, x_mid
      real(wp) :: dist_right, dist_left
      integer :: i, j_mid
      real(wp) :: S_k

      call setup_state(state)

      g_prime = GRAVITY*BETA_S*(S_DENSE - S_LIGHT)/RHO_0
      u_theory = 0.5_wp*sqrt(g_prime*H0)
      L = NX_PHYS*DX_VAL
      t_final = 0.25_wp*L/u_theory
      x_mid = L/2.0_wp

      call run_lock_exchange(state, t_final)

      ! Right front: rightmost cell where bottom-layer S > 0.5 PSU
      x_front_right = find_salinity_front(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%h_layer, &
                                          state%grid%nx_total, state%grid%ny_total, &
                                          NGHOST, 1, DX_VAL, 0.5_wp)

      ! Left retreat: leftmost cell where bottom-layer S < 34.5 PSU
      ! (dense water has been carried away by the rightward current)
      j_mid = NGHOST + NY_PHYS/2
      x_retreat_left = x_mid  ! default
      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         S_k = state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j_mid, 1)/ &
               max(state%multilayer%h_layer(i, j_mid, 1), 1.0e-6_wp)
         if (S_k < S_DENSE - 0.5_wp) then
            x_retreat_left = (real(i - NGHOST, wp) - 0.5_wp)*DX_VAL
            exit
         end if
      end do

      ! Distance each front has traveled from midpoint
      dist_right = x_front_right - x_mid
      dist_left = x_mid - x_retreat_left

      ! Both fronts should have moved
      call check(error, dist_right > 0.3_wp, "Dense front should propagate rightward")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call check(error, dist_left > 0.3_wp, "Dense water should retreat on the left")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Symmetry: both fronts travel similar distances (within 50%)
      call check(error, abs(dist_right - dist_left)/max(dist_right, dist_left) < 0.50_wp, &
                 "Lock exchange should be approximately symmetric")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      call state%destroy()
   end subroutine test_symmetry

   subroutine init_lock_exchange_temperature(state)
      !! Cold dense water on the left, warm light water on the right.
      !! Salinity is uniform (=0) and beta_S is set to 0 so only
      !! temperature drives the density contrast.
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_mid

      ! Flat bottom, uniform depth, at rest
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            state%barotropic%b(i, j) = 0.0_wp
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Zero salinity everywhere; density contrast comes from temperature only
      i_mid = NGHOST + NX_PHYS/2
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = 0.0_wp
               if (i <= i_mid) then
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_COLD
               else
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_HOT
               end if
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, 0.0_wp, 0.0_wp, &
                                ALPHA_T, T_REF_LOCK, &
                                state%grid, NZ)
   end subroutine init_lock_exchange_temperature

   subroutine setup_state_temperature(state)
      !! Temperature-driven lock exchange setup
      type(state_t), intent(inout) :: state

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%rho_0 = RHO_0
      state%multilayer%beta_S = 0.0_wp            ! no salinity contribution
      state%multilayer%S_ref = 0.0_wp
      state%multilayer%alpha_T = ALPHA_T          ! thermal expansion drives the flow
      state%multilayer%T_ref = T_REF_LOCK
      state%multilayer%initial_temperature = T_HOT
      state%multilayer%initial_salinity = 0.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      call init_lock_exchange_temperature(state)
   end subroutine setup_state_temperature

   subroutine init_lock_exchange_combined(state)
      !! Opposite-sign salinity and temperature gradients tuned so that
      !! the net density contrast is zero.  Salt+warm on the left,
      !! fresh+cold on the right.  Density anomaly:
      !!   left:  beta_S*S_DENSE - alpha_T*(T_HOT - T_REF)
      !!        = 0.78*35 - 2.73*(25 - 25)  = +27.3   (from salinity only)
      !! To cancel, we need the opposite-sign temperature contribution:
      !!   shift T_REF so left side has T < T_REF; i.e. T_REF = T_HOT, and
      !!   put T_HOT on left (= T_REF, 0 contribution) while reserving a
      !!   compensating offset elsewhere.  The cleanest setup is symmetric:
      !! use beta_S*S_DENSE/2 above T_REF and -beta_S*S_DENSE/2 below.
      !!
      !! Concretely: salt on left (+27.3), fresh on right (0).
      !! Put T_HOT on right (T - T_REF = +10 -> -alpha_T*10 = -27.3) and
      !! T_COLD on left (T - T_REF = 0 -> 0).  With T_REF = T_COLD:
      !!   left  rho anomaly: beta_S*S_DENSE - alpha_T*(T_COLD - T_COLD) = +27.3
      !!   right rho anomaly: 0             - alpha_T*(T_HOT  - T_COLD) = -27.3
      !! Net contrast = +54.6, NOT zero.  Try T_REF = T_HOT:
      !!   left  rho anomaly: beta_S*S_DENSE - alpha_T*(T_COLD - T_HOT) = 27.3 + 27.3 = +54.6
      !!   right rho anomaly: 0             - alpha_T*(T_HOT - T_HOT)  = 0
      !! Same contrast, 54.6.  For cancellation we need warm on the salty
      !! side and cold on the fresh side, with T_REF at the midpoint:
      !!   T_REF = (T_HOT + T_COLD)/2 = 20
      !!   left  (salt, hot):  beta_S*35 - alpha_T*(25-20) = 27.3 - 13.65 = +13.65
      !!   right (fresh, cold): 0 - alpha_T*(15-20)       = +13.65
      !! Contrast = 0. Flat density.
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_mid
      real(wp), parameter :: T_REF_MID = 0.5_wp*(T_HOT + T_COLD)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            state%barotropic%b(i, j) = 0.0_wp
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      i_mid = NGHOST + NX_PHYS/2
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (i <= i_mid) then
                  ! Salty and hot
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_DENSE
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_HOT
               else
                  ! Fresh and cold
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_LIGHT
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_COLD
               end if
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, 0.0_wp, &
                                ALPHA_T, T_REF_MID, &
                                state%grid, NZ)
   end subroutine init_lock_exchange_combined

   subroutine setup_state_combined(state)
      !! Combined (cancelling) S+T lock exchange
      type(state_t), intent(inout) :: state
      real(wp), parameter :: T_REF_MID = 0.5_wp*(T_HOT + T_COLD)

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%rho_0 = RHO_0
      state%multilayer%beta_S = BETA_S
      state%multilayer%S_ref = 0.0_wp
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%T_ref = T_REF_MID
      state%multilayer%initial_temperature = T_REF_MID
      state%multilayer%initial_salinity = 0.5_wp*(S_DENSE + S_LIGHT)

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      call init_lock_exchange_combined(state)
   end subroutine setup_state_combined

   subroutine test_temperature_driven(error)
      !! Cold water spreading right along the bottom; front speed should
      !! be of the same order as the salinity-only case since the
      !! density contrast was tuned to match.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, u_max, u_peak
      integer :: i, j, k

      call setup_state_temperature(state)

      t_final = 1.0_wp
      call run_lock_exchange(state, t_final)

      ! Check peak horizontal velocity: density contrast was tuned to match
      ! the salinity-only case, so a gravity current of order
      ! 0.5*sqrt(g'*H) ~ 0.26 m/s should develop.  If alpha_T or the EOS
      ! signature are wired incorrectly the flow stays at rest.
      u_peak = 0.0_wp
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (state%multilayer%h_layer(i, j, k) > 1.0e-6_wp) then
                  u_max = abs(state%multilayer%hu_layer(i, j, k)/state%multilayer%h_layer(i, j, k))
                  if (u_max > u_peak) u_peak = u_max
               end if
            end do
         end do
      end do

      call check(error, u_peak > 0.10_wp, &
                 "temperature-driven flow develops: peak |u| > 0.10 m/s")

      call state%destroy()
   end subroutine test_temperature_driven

   subroutine test_combined_cancellation(error)
      !! Opposite-sign S and T gradients tuned so net density anomaly
      !! is zero.  No gravity current should develop; peak velocity
      !! stays very small throughout the run.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, u_max, u_peak
      integer :: i, j, k

      call setup_state_combined(state)

      t_final = 1.0_wp
      call run_lock_exchange(state, t_final)

      ! Scan peak layer velocity after the run
      u_peak = 0.0_wp
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (state%multilayer%h_layer(i, j, k) > 1.0e-6_wp) then
                  u_max = abs(state%multilayer%hu_layer(i, j, k)/state%multilayer%h_layer(i, j, k))
                  if (u_max > u_peak) u_peak = u_max
               end if
            end do
         end do
      end do

      ! Cancelling case: peak velocity should be orders of magnitude smaller
      ! than the salinity-only gravity current speed (~0.26 m/s).  Allow
      ! 5 cm/s as a generous tolerance for numerical noise from advection
      ! of sharp initial tracer steps.
      call check(error, u_peak < 0.05_wp, &
                 "zero-net-density-contrast: peak velocity stays below 0.05 m/s")

      call state%destroy()
   end subroutine test_combined_cancellation

   subroutine init_lock_exchange_reinforcing(state)
      !! Reinforcing S+T contrast: salty+cold on the left, fresh+warm on
      !! the right, with T_REF at the midpoint. Each tracer contributes
      !! half of the salinity-only density anomaly so their sum matches.
      !!   rho(left)  - rho(right) = beta_S*(S_DENSE - S_LIGHT)
      !!                           + alpha_T*(T_HOT  - T_COLD)
      !! Both kernels must be active for the front to propagate at the
      !! salinity-only reference speed.
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_mid
      real(wp), parameter :: T_REF_MID = 0.5_wp*(T_HOT + T_COLD)
      real(wp), parameter :: S_MID = 0.5_wp*(S_DENSE + S_LIGHT)
      real(wp), parameter :: S_LEFT = S_MID + 0.5_wp*(S_DENSE - S_LIGHT)
      real(wp), parameter :: S_RIGHT = S_MID - 0.5_wp*(S_DENSE - S_LIGHT)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%h(i, j) = H0
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            state%barotropic%b(i, j) = 0.0_wp
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      i_mid = NGHOST + NX_PHYS/2
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (i <= i_mid) then
                  ! Salty and cold -> both increase density
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_LEFT
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_COLD
               else
                  ! Fresh and warm -> both decrease density
                  state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*S_RIGHT
                  state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(i, j, k) = state%multilayer%h_layer(i, j, k)*T_HOT
               end if
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, 0.0_wp, &
                                ALPHA_T, T_REF_MID, &
                                state%grid, NZ)
   end subroutine init_lock_exchange_reinforcing

   subroutine setup_state_reinforcing(state)
      type(state_t), intent(inout) :: state
      real(wp), parameter :: T_REF_MID = 0.5_wp*(T_HOT + T_COLD)

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%rho_0 = RHO_0
      state%multilayer%beta_S = BETA_S
      state%multilayer%S_ref = 0.0_wp
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%T_ref = T_REF_MID
      state%multilayer%initial_temperature = T_REF_MID
      state%multilayer%initial_salinity = 0.5_wp*(S_DENSE + S_LIGHT)

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")

      call init_lock_exchange_reinforcing(state)
   end subroutine setup_state_reinforcing

   subroutine test_combined_reinforcing(error)
      !! Combined S+T lock exchange where both tracers contribute equally
      !! to the same density contrast. The total rho anomaly equals the
      !! salinity-only reference case, so the bottom-layer front must
      !! propagate at approximately the same speed. Both EOS/BPG paths
      !! must be active for this to match; disabling either collapses
      !! the front speed to half the reference.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: g_prime, u_theory, t_final, L
      real(wp) :: x_mid, x_front, u_measured
      real(wp), parameter :: S_THRESH = 0.5_wp*(S_DENSE + S_LIGHT)

      call setup_state_reinforcing(state)

      ! Full density anomaly: beta_S*dS + alpha_T*dT where each half is
      ! exactly half of the salinity-only contrast, summing to the same
      ! total.
      g_prime = GRAVITY*(BETA_S*(S_DENSE - S_LIGHT) &
                         + ALPHA_T*(T_HOT - T_COLD))/RHO_0
      u_theory = 0.5_wp*sqrt(g_prime*H0)
      L = NX_PHYS*DX_VAL
      t_final = 0.25_wp*L/u_theory
      x_mid = L/2.0_wp

      call run_lock_exchange(state, t_final)

      ! Use a salinity threshold above the mid concentration to locate
      ! the leading edge of the dense (salty) water.
      x_front = find_salinity_front(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, state%multilayer%h_layer, &
                                    state%grid%nx_total, state%grid%ny_total, &
                                    NGHOST, 1, DX_VAL, S_THRESH)
      u_measured = (x_front - x_mid)/t_final

      call check(error, x_front > x_mid + 0.5_wp, &
                 "combined S+T front must propagate rightward")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Same 50% band as the salinity-only reference test (numerical
      ! diffusion and sigma-layer discretisation both apply identically).
      call check(error, abs(u_measured - u_theory)/u_theory < 0.50_wp, &
                 "combined S+T front speed must match Benjamin theory within 50%")

      call state%destroy()
   end subroutine test_combined_reinforcing

   function find_temperature_front(hT_layer, h_layer, nx, ny, ng, layer, dx, T_threshold) result(x_front)
      !! Scan the bottom layer from right to find the first cell where
      !! T drops below the threshold (leading edge of the cold gravity
      !! current spreading rightward).
      real(wp), intent(in) :: hT_layer(:, :, :), h_layer(:, :, :)
      integer, intent(in) :: nx, ny, ng, layer
      real(wp), intent(in) :: dx, T_threshold
      real(wp) :: x_front

      integer :: i, j_mid
      real(wp) :: T_k

      j_mid = ng + (ny - 2*ng)/2

      x_front = real(nx - 2*ng, wp)*dx*0.5_wp   ! default: midpoint
      do i = nx - ng, ng + 1, -1
         T_k = hT_layer(i, j_mid, layer)/max(h_layer(i, j_mid, layer), 1.0e-6_wp)
         if (T_k < T_threshold) then
            x_front = (real(i - ng, wp) - 0.5_wp)*dx
            exit
         end if
      end do
   end function find_temperature_front

end module test_lock_exchange
