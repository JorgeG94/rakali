!! Analytical test: per-layer inertial oscillation, parallelism check
!!
!! Closed basin, uniform depth, uniform density (no horizontal pressure
!! gradient, no buoyancy), spatially-uniform velocity field with
!! BAROCLINIC vertical shear (different velocity magnitude per layer,
!! all initially aligned with +x).
!!
!! Under per-layer Coriolis every layer's velocity vector rotates by
!! the same angle each step (the rotation depends only on `f` and `dt`,
!! not on the velocity itself), so all layer velocity vectors remain
!! PARALLEL throughout the run.
!!
!! Under the old depth-averaged Coriolis only the depth-averaged
!! momentum rotates; the per-layer baroclinic perturbations stay frozen
!! in their initial direction.  After even a small rotation the layer
!! velocity vectors diverge in direction (each layer ends up as
!! perturbation_x + rotated_mean), so the per-layer cross product
!! becomes non-zero.
!!
!! The discriminator: after running for a non-trivial fraction of the
!! inertial period, the cross product of any two layer velocity vectors
!! must be ≈ 0 (parallel).  This assertion is INSENSITIVE to the exact
!! rotation magnitude (so RK2 phase errors don't affect it) — it
!! purely checks the per-layer rotation property that the new kernel
!! provides and the old one does not.
module test_inertial_oscillation
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   implicit none
   private

   public :: collect_inertial_oscillation_tests

   ! Domain — large closed basin so wall reflections (gravity waves
   ! generated when the wall mirrors the spatially-uniform velocity
   ! into the ghost cells) don't reach the centre probe within the
   ! test horizon.  At c = sqrt(g*H0) ≈ 9.9 m/s and T_i/4 ≈ 157 s,
   ! waves can travel ~1.6 km — so a 6 km basin (probe 3 km from walls)
   ! gives ~2× margin.
   integer, parameter :: NX_PHYS = 30
   integer, parameter :: NY_PHYS = 30
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 4
   real(wp), parameter :: H0 = 10.0_wp        !! 10 m water column
   real(wp), parameter :: DX_VAL = 200.0_wp   !! 200 m cells → 6 km domain

   ! Coriolis parameter: pick large f to keep T_i short so the test
   ! runs quickly.  Physical mid-latitude is ~1e-4; we use 1e-2 →
   ! T_i = 2*pi/f ≈ 628 s, T_i/4 ≈ 157 s.  Per-step rotation angle
   ! f*dt at the barotropic CFL is small (~0.1 rad) so many small
   ! rotations compose to the full pi/2 turn.
   real(wp), parameter :: F_COR = 1.0e-2_wp

   ! Density (uniform — no buoyancy)
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: T_INIT = 15.0_wp

contains

   subroutine collect_inertial_oscillation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("baroclinic_layers_stay_parallel", &
                               test_baroclinic_layers_stay_parallel) &
                  ]
   end subroutine collect_inertial_oscillation_tests

   subroutine setup_state(state)
      !! Build a closed-basin state with uniform density and a
      !! spatially-uniform baroclinic velocity profile (different
      !! velocity in each layer, all aligned with +x at t=0).
      type(state_t), intent(inout) :: state

      integer :: i, j, k

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp     ! decouple T from dynamics
      state%multilayer%beta_S = 0.0_wp      ! decouple S from dynamics (no BPG)
      state%multilayer%rho_0 = RHO_0
      ! Disable PP81 vertical viscosity: with baroclinic IC + uniform
      ! density it would mix the per-layer velocities (Ri = 0), breaking
      ! the per-layer parallelism assertion.  This test is specifically
      ! about isolating the per-layer Coriolis rotation from any other
      ! kernel that touches per-layer momentum.
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp   ! belt + braces (default 0)
      state%multilayer%smagorinsky_cs = 0.0_wp

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.45_wp)

      ! Closed basin
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      ! Coriolis on, no friction, no wind
      state%barotropic%coriolis_f = F_COR
      state%barotropic%manning_n = 0.0_wp

      ! Flat bathymetry, uniform depth
      state%barotropic%h = H0
      state%barotropic%b = 0.0_wp
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      ! Build uniform sigma layers
      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Per-layer initial velocity: U_k = real(k, wp) along +x.  Every
      ! layer at every interior cell carries this velocity, so spatial
      ! gradients are zero and the only thing that can change the layer
      ! velocity is the Coriolis kernel.  (Wall BCs reflect at the
      ! boundary — the interior cells stay spatially uniform.)
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               state%multilayer%hu_layer(i, j, k) = state%multilayer%h_layer(i, j, k)*real(k, wp)
               state%multilayer%hv_layer(i, j, k) = 0.0_wp
            end do
         end do
      end do

      ! Tracers: uniform, just so they exist (alpha_T = beta_S = 0 means
      ! density is just rho_0 everywhere, so no BPG).
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = state%multilayer%h_layer*S_INIT
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = state%multilayer%h_layer*T_INIT

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_state

   subroutine run_to(state, t_final)
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
   end subroutine run_to

   subroutine test_baroclinic_layers_stay_parallel(error)
      !! After running for a non-trivial fraction of the inertial period,
      !! all layer velocity vectors must remain mutually PARALLEL.  The
      !! cross product (u_a * v_b - u_b * v_a) is the parallelism
      !! discriminant — it is identically zero when both vectors point
      !! in the same direction (or opposite directions), regardless of
      !! magnitude.
      !!
      !! The initial velocities are (U_k, 0) for k = 1..NZ — all
      !! parallel along +x.  Per-layer Coriolis rotates every layer by
      !! the same f*dt per step, preserving parallelism for all time.
      !! Depth-averaged Coriolis would rotate only the mean and leave
      !! the perturbations in +x, so layer 1 (perturbation negative
      !! along x) and layer NZ (perturbation positive along x) end up
      !! at different angles after even one step, and the cross product
      !! grows linearly with the rotation angle.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      real(wp) :: T_run, h_top, h_bot
      real(wp) :: u_top, v_top, u_bot, v_bot
      real(wp) :: cross, mag_top, mag_bot, mag_product
      integer :: ic, jc

      call setup_state(state)

      ! Run for an eighth of the inertial period — long enough for the
      ! rotation to be substantial (~45 degrees) but well within the
      ! window before wall-reflected gravity waves reach the centre.
      T_run = 0.25_wp*PI/F_COR     ! T_i/8 = pi/(4 f) ≈ 79 s for f=0.01

      ic = NGHOST + NX_PHYS/2
      jc = NGHOST + NY_PHYS/2

      call run_to(state, T_run)

      h_top = state%multilayer%h_layer(ic, jc, NZ)
      h_bot = state%multilayer%h_layer(ic, jc, 1)
      if (h_top <= 0.0_wp .or. h_bot <= 0.0_wp) then
         call check(error, .false., "Probe layer thickness vanished — bad IC")
         call state%destroy()
         return
      end if
      u_top = state%multilayer%hu_layer(ic, jc, NZ)/h_top
      v_top = state%multilayer%hv_layer(ic, jc, NZ)/h_top
      u_bot = state%multilayer%hu_layer(ic, jc, 1)/h_bot
      v_bot = state%multilayer%hv_layer(ic, jc, 1)/h_bot

      ! Cross product, normalised by the magnitude product to give
      ! sin(angle_between).  A value near 0 means the two vectors are
      ! parallel (or anti-parallel); under the per-layer Coriolis they
      ! must point the same direction (sin ≈ 0, ~1 % tolerance).
      cross = u_top*v_bot - u_bot*v_top
      mag_top = sqrt(u_top*u_top + v_top*v_top)
      mag_bot = sqrt(u_bot*u_bot + v_bot*v_bot)
      mag_product = mag_top*mag_bot

      if (mag_product <= 0.0_wp) then
         call check(error, .false., "Layer velocities collapsed to zero — bad IC")
         call state%destroy()
         return
      end if

      call check(error, abs(cross)/mag_product < 1.0e-2_wp, &
                 "Layer 1 and layer NZ velocities stay parallel under per-layer Coriolis")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Sanity: rotation has actually occurred (otherwise parallelism is
      ! trivially preserved by doing nothing).  Initial v was zero; after
      ! T_i/8 the v components must be substantially non-zero with the
      ! correct sign (negative — clockwise rotation of +x velocity).
      call check(error, v_top < -0.1_wp*mag_top, &
                 "Layer NZ has rotated substantially (v_top < 0 from initial 0)")

      call state%destroy()
   end subroutine test_baroclinic_layers_stay_parallel

end module test_inertial_oscillation
