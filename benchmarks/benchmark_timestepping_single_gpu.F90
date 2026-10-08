!! Single-GPU time-integrator benchmark (performance + convergence)
program benchmark_timestepping_single_gpu
   !! Companion to bench_single_gpu, focused on the selectable explicit
   !! time integrator (`time_scheme` = euler / rk2 / rk3 / rk43) on the 2D
   !! barotropic coastal path.  Runs identically on CPU (serial / multicore)
   !! and GPU (OpenACC) — the backend is selected at build time, exactly
   !! like bench_single_gpu.
   !!
   !! Two sections, both CSV to stdout (human-readable to stderr):
   !!
   !!   [performance] per-scheme step throughput (Mcells/s) on a dam-break
   !!     grid driven through the production solver_step.  euler is one
   !!     flux/stage, rk43 is four — so raw Mcells/s drops with stage count,
   !!     while the `ssp_coeff` column shows the larger stable dt that buys
   !!     it back (rk43 covers 2x the simulated time per step).
   !!
   !!   [convergence] per-scheme work-precision: a smooth standing wave
   !!     integrated to a fixed time at several fixed dt, L2 error against a
   !!     fine reference vs the measured wall time.  This is the error-vs-
   !!     walltime data for the TUFLOW HPC / SCHISM comparison — the fair
   !!     metric when the schemes (and their stable dt) differ.
   !!
   !! Usage:
   !!   benchmark_timestepping_single_gpu            default sweep
   !!   benchmark_timestepping_single_gpu <N>        single grid N x N (perf)
   !!   benchmark_timestepping_single_gpu <N> <steps>
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_timestepper, only: timestepper_t, make_timestepper, &
                              TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   implicit none

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: CFL = 0.45_wp

   character(len=*), parameter :: BACKEND = &
#if defined(RAKALI_GPU_OFFLOAD)
                                  "openacc"
#else
   "cpu_serial"
#endif

   integer, parameter :: N_SCHEMES = 4
   integer, parameter :: SCHEMES(N_SCHEMES) = &
                         [TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43]
   character(len=5), parameter :: SCHEME_NAMES(N_SCHEMES) = &
                                  [character(len=5) :: "euler", "rk2", "rk3", "rk43"]

   integer :: ig, n_sizes, timed_steps, warmup_steps, isc
   integer, allocatable :: grid_sizes(:)
   character(len=32) :: arg

   if (command_argument_count() >= 1) then
      allocate (grid_sizes(1))
      call get_command_argument(1, arg)
      read (arg, *) grid_sizes(1)
      n_sizes = 1
   else
#if defined(RAKALI_GPU_OFFLOAD)
      grid_sizes = [512, 1024, 2048]
#else
      grid_sizes = [64, 128, 256]
#endif
      n_sizes = size(grid_sizes)
   end if

   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      read (arg, *) timed_steps
   else
#if defined(RAKALI_GPU_OFFLOAD)
      timed_steps = 500
#else
      timed_steps = 100
#endif
   end if

   warmup_steps = max(10, timed_steps/10)

   ! ---- Performance section ----
   write (*, '(a)') "# section=performance"
   write (*, '(a)') "backend,scheme,n_stages,ssp_coeff,nx,ny,warmup,"// &
      "timed_steps,total_s,mcells_s,us_per_cell_step"
   do ig = 1, n_sizes
      do isc = 1, N_SCHEMES
         call run_perf(SCHEMES(isc), SCHEME_NAMES(isc), &
                       grid_sizes(ig), grid_sizes(ig), warmup_steps, timed_steps)
      end do
   end do

   ! ---- Convergence / work-precision section ----
   write (*, '(a)') "# section=convergence"
   write (*, '(a)') "backend,scheme,dt,n_steps,l2_error,walltime_s"
   call run_convergence()

contains

   subroutine run_perf(scheme, name, nx_phys, ny_phys, n_warmup, n_timed)
      !! Time n_timed production solver_step calls for one scheme.
      integer, intent(in) :: scheme
      character(len=*), intent(in) :: name
      integer, intent(in) :: nx_phys, ny_phys, n_warmup, n_timed

      type(state_t) :: state
      type(timestepper_t) :: ts
      type(timer_type) :: t_total
      integer :: n, i, j, nx, ny
      real(wp) :: elapsed, mcells, us_per_cell_step
      character(len=256) :: csv_line

      ts = make_timestepper(scheme)

      call state%init(nx_phys, ny_phys, NGHOST, 1.0_wp, 1.0_wp, CFL)
      state%time_scheme = scheme

      nx = state%grid%nx_total
      ny = state%grid%ny_total

      ! Dam-break initial condition (same as bench_single_gpu).
      do j = 1, ny
         do i = 1, nx
            state%barotropic%b(i, j) = 0.0_wp
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
            if (i <= NGHOST + nx_phys/2) then
               state%barotropic%h(i, j) = 2.0_wp
            else
               state%barotropic%h(i, j) = 1.0_wp
            end if
         end do
      end do

      call solver_enter_data(state)

      do n = 1, n_warmup
         call solver_step(state)
      end do

      call t_total%start()
      do n = 1, n_timed
         call solver_step(state)
      end do
      call t_total%stop()
      elapsed = t_total%get_elapsed_time()

      call solver_exit_data(state)

      mcells = real(nx_phys, wp)*real(ny_phys, wp)*real(n_timed, wp) &
               /elapsed/1.0e6_wp
      us_per_cell_step = elapsed*1.0e6_wp &
                         /(real(nx_phys, wp)*real(ny_phys, wp)*real(n_timed, wp))

      call logger%info("[perf] "//trim(name)//" "//to_string(nx_phys)//"^2 ["// &
                       BACKEND//"]  "//to_string(mcells)//" Mcells/s, "// &
                       to_string(ts%n_stages)//" stages, ssp_coeff="// &
                       to_string(ts%ssp_coeff))

      write (csv_line, '(a,",",a,",",i0,",",f6.3,",",i0,",",i0,",",i0,",",'// &
             'i0,",",es12.5,",",es12.5,",",es12.5)') &
         BACKEND, trim(name), ts%n_stages, ts%ssp_coeff, nx_phys, ny_phys, &
         n_warmup, n_timed, elapsed, mcells, us_per_cell_step
      write (*, '(a)') trim(csv_line)

      call state%destroy()
   end subroutine run_perf

   subroutine run_convergence()
      !! Smooth standing wave -> per-scheme L2 error vs a fine reference,
      !! with measured wall time, over a sweep of fixed timesteps.
      integer, parameter :: NXC = 128, NYC = 4
      real(wp), parameter :: DXC = 100.0_wp, H0 = 10.0_wp, AMP = 0.01_wp
      real(wp), parameter :: T_FINAL = 60.0_wp
      integer, parameter :: N_DT = 4
      ! Coarsest step per scheme then halved N_DT times.
      real(wp), parameter :: DT_BASE(N_SCHEMES) = &
                             [1.0_wp, 2.0_wp, 4.0_wp, 5.0_wp]

      real(wp) :: ref(NXC, NYC), sol(NXC, NYC)
      real(wp) :: dt, err, wall
      integer :: isc, idt, nsteps
      character(len=256) :: csv_line

      ! Fine reference: rk43 with a small step to T_FINAL.
      call wave_run(TIME_RK43, NXC, NYC, DXC, H0, AMP, 0.125_wp, &
                    nint(T_FINAL/0.125_wp), ref, wall)

      do isc = 1, N_SCHEMES
         dt = DT_BASE(isc)
         do idt = 1, N_DT
            nsteps = nint(T_FINAL/dt)
            call wave_run(SCHEMES(isc), NXC, NYC, DXC, H0, AMP, dt, nsteps, &
                          sol, wall)
            err = l2_diff(sol, ref, NXC, NYC)

            call logger%info("[conv] "//trim(SCHEME_NAMES(isc))//" dt="// &
                             to_string(dt)//" err="//to_string(err)// &
                             " wall="//to_string(wall)//"s")

            write (csv_line, '(a,",",a,",",es12.5,",",i0,",",es12.5,",",es12.5)') &
               BACKEND, trim(SCHEME_NAMES(isc)), dt, nsteps, err, wall
            write (*, '(a)') trim(csv_line)

            dt = 0.5_wp*dt
         end do
      end do
   end subroutine run_convergence

   subroutine wave_run(scheme, nxc, nyc, dxc, h0, amp, dt, nsteps, h_out, wall)
      !! Integrate the standing wave with fixed dt; return interior h + wall time.
      integer, intent(in) :: scheme, nxc, nyc, nsteps
      real(wp), intent(in) :: dxc, h0, amp, dt
      real(wp), intent(out) :: h_out(nxc, nyc)
      real(wp), intent(out) :: wall

      type(state_t) :: s
      type(timer_type) :: t
      integer :: i, j, n
      real(wp) :: xc, kwave

      call s%init(nxc, nyc, NGHOST, dxc, dxc, CFL)
      s%time_scheme = scheme
      s%dt_fixed = dt

      kwave = acos(-1.0_wp)/(real(nxc, wp)*dxc)
      s%barotropic%b = 0.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            xc = (real(i - s%grid%nghost, wp) - 0.5_wp)*dxc
            s%barotropic%h(i, j) = h0 + amp*cos(kwave*xc)
         end do
      end do

      call solver_enter_data(s)
      call t%start()
      do n = 1, nsteps
         call solver_step(s)
      end do
      call t%stop()
      wall = t%get_elapsed_time()
      call solver_exit_data(s)

      do j = 1, nyc
         do i = 1, nxc
            h_out(i, j) = s%barotropic%h(i + NGHOST, j + NGHOST)
         end do
      end do

      call s%destroy()
   end subroutine wave_run

   function l2_diff(a, b, nxc, nyc) result(d)
      integer, intent(in) :: nxc, nyc
      real(wp), intent(in) :: a(nxc, nyc), b(nxc, nyc)
      real(wp) :: d
      integer :: i, j

      d = 0.0_wp
      do j = 1, nyc
         do i = 1, nxc
            d = d + (a(i, j) - b(i, j))**2
         end do
      end do
      d = sqrt(d/real(nxc*nyc, wp))
   end function l2_diff

end program benchmark_timestepping_single_gpu
