!! Gaussian hump gravity wave propagation test
module test_gaussian_hump
   !! Runs a small-amplitude Gaussian surface perturbation on a flat
   !! bottom with open boundaries. Tests wave dispersion, circular
   !! symmetry preservation, and clean open boundary radiation.
   !!
   !! Spec reference: Section 7.5
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_boundary, only: bc_params_t, BC_OPEN
   implicit none
   private

   public :: collect_gaussian_hump_tests

   integer, parameter :: NX = 50
      !! Physical cells in x
   integer, parameter :: NY = 50
      !! Physical cells in y
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DEPTH = 10.0_wp
      !! Flat bottom depth (m)
   real(wp), parameter :: DOMAIN_SIZE = 500.0_wp
      !! Domain extent in each direction (m)
   real(wp), parameter :: DX = DOMAIN_SIZE/real(NX, wp)
   real(wp), parameter :: DY = DOMAIN_SIZE/real(NY, wp)
   real(wp), parameter :: CFL_NUM = 0.45_wp
   real(wp), parameter :: AMP = 0.01_wp
      !! Gaussian amplitude (m) — small for linear regime
   real(wp), parameter :: SIGMA = 50.0_wp
      !! Gaussian width (m)
   real(wp), parameter :: WAVE_SPEED = 9.9045_wp
      !! sqrt(g * H) ~ sqrt(9.81 * 10)

contains

   subroutine collect_gaussian_hump_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gaussian_circular_symmetry", test_symmetry), &
                  new_unittest("gaussian_mass_conserved", test_mass), &
                  new_unittest("gaussian_no_negative_depth", test_no_negative), &
                  new_unittest("gaussian_boundary_exit", test_boundary_exit) &
                  ]
   end subroutine collect_gaussian_hump_tests

   subroutine setup_gaussian(state)
      !! Initialise state with Gaussian hump at domain centre, open BCs
      type(state_t), intent(inout) :: state

      integer :: i, j
      real(wp) :: x_pos, y_pos, x_center, y_center, r_sq

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      ! Open boundaries on all sides
      state%bc%west = BC_OPEN
      state%bc%east = BC_OPEN
      state%bc%south = BC_OPEN
      state%bc%north = BC_OPEN

      x_center = DOMAIN_SIZE/2.0_wp
      y_center = DOMAIN_SIZE/2.0_wp

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x_pos = (real(i - NGHOST, wp) - 0.5_wp)*DX
            y_pos = (real(j - NGHOST, wp) - 0.5_wp)*DY

            state%barotropic%b(i, j) = -DEPTH
            r_sq = (x_pos - x_center)**2 + (y_pos - y_center)**2
            state%barotropic%h(i, j) = DEPTH + AMP*exp(-r_sq/(2.0_wp*SIGMA*SIGMA))
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

   end subroutine setup_gaussian

   function compute_mass(s) result(mass)
      type(state_t), intent(in) :: s
      real(wp) :: mass
      integer :: i, j

      mass = 0.0_wp
      do j = s%grid%nghost + 1, s%grid%ny_total - s%grid%nghost
         do i = s%grid%nghost + 1, s%grid%nx_total - s%grid%nghost
            mass = mass + s%barotropic%h(i, j)*s%grid%dx*s%grid%dy
         end do
      end do
   end function compute_mass

   subroutine test_symmetry(error)
      !! After wave propagation, check that circular symmetry is preserved.
      !! Compare eta perturbation along x-axis vs y-axis at same radius.
      use, intrinsic :: iso_fortran_env, only: error_unit
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t, t_end
      integer :: i_mid, j_mid, k, n_check, n_steps
      real(wp) :: eta_x, eta_y, max_asym, max_pert
      real(wp) :: eta_center

      call setup_gaussian(state)

      write (error_unit, '(a,f8.2,a,f8.2)') '  dx=', DX, '  dy=', DY
      write (error_unit, '(a,i6,a,i6)') '  nx_total=', state%grid%nx_total, &
         '  ny_total=', state%grid%ny_total
      write (error_unit, '(a,es12.4)') '  initial h at center=', &
         state%barotropic%h(NGHOST + NX/2 + 1, NGHOST + NY/2 + 1)
      write (error_unit, '(a,es12.4)') '  initial b at center=', &
         state%barotropic%b(NGHOST + NX/2 + 1, NGHOST + NY/2 + 1)

      ! Run until wave front is about halfway to boundary
      ! t ~ (domain/4) / c = 125 / 9.9 ~ 12.5 s
      t_end = 12.0_wp

      call solver_enter_data(state)
      t = 0.0_wp
      n_steps = 0
      do while (t < t_end)
         call solver_step(state, dt_max=(t_end - t))
         t = t + state%dt
         n_steps = n_steps + 1
         if (n_steps <= 3 .or. mod(n_steps, 50) == 0) then
            write (error_unit, '(a,i6,a,es12.4,a,es12.4)') &
               '  step=', n_steps, '  dt=', state%dt, '  t=', t
         end if
      end do
      call solver_exit_data(state)

      write (error_unit, '(a,i6,a,es12.4)') '  total steps=', n_steps, &
         '  final t=', t

      i_mid = NGHOST + NX/2 + 1
      j_mid = NGHOST + NY/2 + 1

      eta_center = state%barotropic%h(i_mid, j_mid) + state%barotropic%b(i_mid, j_mid)
      write (error_unit, '(a,es12.4)') '  eta at center after run=', eta_center

      ! Check symmetry along +x vs +y at several radii
      max_asym = 0.0_wp
      max_pert = 0.0_wp
      n_check = min(NX/4, NY/4)
      do k = 1, n_check
         eta_x = state%barotropic%h(i_mid + k, j_mid) + state%barotropic%b(i_mid + k, j_mid)
         eta_y = state%barotropic%h(i_mid, j_mid + k) + state%barotropic%b(i_mid, j_mid + k)
         max_asym = max(max_asym, abs(eta_x - eta_y))
         max_pert = max(max_pert, abs(eta_x), abs(eta_y))
         if (k <= 5 .or. k == n_check) then
            write (error_unit, '(a,i4,a,es12.4,a,es12.4,a,es12.4)') &
               '  k=', k, '  eta_x=', eta_x, '  eta_y=', eta_y, &
               '  diff=', abs(eta_x - eta_y)
         end if
      end do

      write (error_unit, '(a,es12.4,a,es12.4)') &
         '  max_asym=', max_asym, '  max_pert=', max_pert

      ! Symmetry should hold to near machine precision (dx = dy, same grid)
      call check(error, max_asym < 1.0e-12_wp, &
                 "Gaussian hump circular symmetry violated")

      call state%destroy()
   end subroutine test_symmetry

   subroutine test_mass(error)
      !! Mass should be conserved as wave propagates and exits.
      !! With open BCs, mass will decrease as waves leave. We check
      !! conservation during the early phase (wave still interior).
      use, intrinsic :: iso_fortran_env, only: error_unit
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t, t_end, mass_initial, mass_after
      real(wp) :: mass_error

      call setup_gaussian(state)
      mass_initial = compute_mass(state)

      ! Run for short time while wave is fully interior
      ! Wave needs ~25s to reach boundary (250m / 9.9 m/s)
      t_end = 10.0_wp

      call solver_enter_data(state)
      t = 0.0_wp
      do while (t < t_end)
         call solver_step(state, dt_max=(t_end - t))
         t = t + state%dt
      end do
      call solver_exit_data(state)

      mass_after = compute_mass(state)
      mass_error = abs(mass_after - mass_initial)/mass_initial

      write (error_unit, '(a,es16.8,a,es16.8,a,es12.4)') &
         '  mass_init=', mass_initial, '  mass_after=', mass_after, &
         '  rel_err=', mass_error

      ! The wave front hasn't reached the boundary at t = 10 s (250 m
      ! to go at ~10 m/s), but the Gaussian IC has a non-zero tail at
      ! the boundary: amp * exp(-r^2 / 2sigma^2) = 0.01 * exp(-12.5) ~
      ! 4e-8 at the wall. Flather BCs interact with this tail and
      ! leak mass at ~1e-7 per step in relative terms. Tighter than
      ! the original 1e-4 by two orders of magnitude, but not as
      ! tight as a truly closed basin — this IC isn't one.
      call check(error, mass_error < 1.0e-6_wp, &
                 "Gaussian hump: mass loss while wave is interior must be <1e-6")

      call state%destroy()
   end subroutine test_mass

   subroutine test_no_negative(error)
      !! No negative depths should occur during propagation
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t, t_end, min_h
      integer :: i, j

      call setup_gaussian(state)

      t_end = 12.0_wp

      call solver_enter_data(state)
      t = 0.0_wp
      do while (t < t_end)
         call solver_step(state, dt_max=(t_end - t))
         t = t + state%dt
      end do
      call solver_exit_data(state)

      min_h = huge(1.0_wp)
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            min_h = min(min_h, state%barotropic%h(i, j))
         end do
      end do

      call check(error, min_h >= -1.0e-10_wp, &
                 "Negative depth detected in Gaussian hump test")

      call state%destroy()
   end subroutine test_no_negative

   subroutine test_boundary_exit(error)
      !! After wave exits domain, perturbation should be negligible.
      !! Run long enough for wave to fully exit, then check max |eta - 0|.
      use, intrinsic :: iso_fortran_env, only: error_unit
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t, t_end, max_pert, eta_pert
      integer :: i, j, n_steps

      call setup_gaussian(state)

      ! Wave needs ~25s to reach boundary (250m / 9.9 m/s); run ~50s for full exit
      t_end = 50.0_wp

      call solver_enter_data(state)
      t = 0.0_wp
      n_steps = 0
      do while (t < t_end)
         call solver_step(state, dt_max=(t_end - t))
         t = t + state%dt
         n_steps = n_steps + 1
         if (n_steps <= 5 .or. mod(n_steps, 20) == 0) then
            write (error_unit, '(a,i6,a,es12.4,a,es12.4)') &
               '  exit: step=', n_steps, '  dt=', state%dt, '  t=', t
         end if
         ! Bail out early if timestep collapses
         if (state%dt < 1.0e-6_wp) then
            write (error_unit, '(a,es12.4,a,i6)') &
               '  EXIT: dt collapsed to ', state%dt, ' at step ', n_steps
            ! Dump boundary region diagnostics
            call solver_exit_data(state)
            write (error_unit, '(a)') '  Boundary region h, hu, hv (west edge):'
            do j = NGHOST + 1, min(NGHOST + 3, state%grid%ny_total - NGHOST)
               do i = 1, NGHOST + 2
                  write (error_unit, '(a,i3,a,i3,a,es12.4,a,es12.4,a,es12.4)') &
                     '    i=', i, ' j=', j, ' h=', state%barotropic%h(i, j), ' hu=', state%barotropic%hu(i, j), &
                     ' hv=', state%barotropic%hv(i, j)
               end do
            end do
            write (error_unit, '(a)') '  Center region:'
            do i = NX/2, NX/2 + 2
               write (error_unit, '(a,i3,a,es12.4,a,es12.4)') &
                  '    i=', i, ' h=', state%barotropic%h(NGHOST + i, NGHOST + NY/2 + 1), &
                  ' hu=', state%barotropic%hu(NGHOST + i, NGHOST + NY/2 + 1)
            end do
            ! Find max |hu| and where
            block
               real(wp) :: max_hu, max_hv
               integer :: mi, mj
               max_hu = 0.0_wp
               max_hv = 0.0_wp
               mi = 0; mj = 0
               do j = 1, state%grid%ny_total
                  do i = 1, state%grid%nx_total
                     if (abs(state%barotropic%hu(i, j)) > max_hu) then
                        max_hu = abs(state%barotropic%hu(i, j))
                        mi = i; mj = j
                     end if
                     max_hv = max(max_hv, abs(state%barotropic%hv(i, j)))
                  end do
               end do
               write (error_unit, '(a,es12.4,a,i4,a,i4)') &
                  '  max |hu|=', max_hu, ' at i=', mi, ' j=', mj
               write (error_unit, '(a,es12.4)') '  max |hv|=', max_hv
               if (mi > 0) then
                  write (error_unit, '(a,es12.4,a,es12.4)') &
                     '  h at max hu=', state%barotropic%h(mi, mj), &
                     ' b at max hu=', state%barotropic%b(mi, mj)
               end if
            end block
            call check(error, .false., &
                       "Timestep collapsed — open BC instability")
            call state%destroy()
            return
         end if
      end do
      call solver_exit_data(state)

      write (error_unit, '(a,i6,a,es12.4)') '  exit: total steps=', n_steps, &
         '  final t=', t

      ! Check that surface perturbation is small everywhere
      max_pert = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            eta_pert = abs(state%barotropic%h(i, j) + state%barotropic%b(i, j))
            max_pert = max(max_pert, eta_pert)
         end do
      end do

      write (error_unit, '(a,es12.4,a,es12.4)') &
         '  max_pert=', max_pert, '  threshold=', 0.5_wp*AMP

      ! After full exit, perturbation should be small fraction of initial amp
      ! Flather BCs are not perfectly non-reflective; allow 50% of initial amp
      call check(error, max_pert < 0.5_wp*AMP, &
                 "Significant perturbation remains after wave should have exited")

      call state%destroy()
   end subroutine test_boundary_exit

end module test_gaussian_hump
