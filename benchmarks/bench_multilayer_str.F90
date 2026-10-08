!! Full-physics multilayer benchmark, structured grid (single GPU)
program bench_multilayer_str
   !! Structured-grid sister of `bench_multilayer`.  Drives `solver_step`
   !! in multilayer mode on a Cartesian Nx × Ny mesh with a stratified
   !! initial condition.  Exercises the full structured kernel chain:
   !! barotropic HLL flux + update, layer distribute, S/T tracer
   !! advection (PPM), Smagorinsky horizontal viscosity, vertical
   !! exchange, PP81 vertical mixing, EOS, BPG (Jacobian), bottom
   !! friction, MOM6 velocity truncation, and Coriolis.
   !!
   !! Bonus over the unstr sister: the structured `solver_step` is
   !! already wrapped in `profiler_start/stop` regions, so this bench
   !! prints a per-kernel breakdown to stderr after the timed loop.
   !!
   !! Usage:
   !!   bench_multilayer_str                              Default sweep
   !!   bench_multilayer_str <nx> <ny>                    Single grid
   !!   bench_multilayer_str <nx> <ny> <nz>               Single config
   !!   bench_multilayer_str <nx> <ny> <nz> <steps>       Custom step count
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_tracers, only: ml_init_tracer
   use rki_ml_eos, only: ml_equation_of_state
   use rki_profiler, only: profiler_init, profiler_report, profiler_reset, profiler_end
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   implicit none

   character(len=*), parameter :: BACKEND = &
#if defined(RAKALI_GPU_OFFLOAD)
                                  "openacc"
#else
   "cpu_serial"
#endif

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DOMAIN_LX = 4000.0_wp
   real(wp), parameter :: DOMAIN_LY = 4000.0_wp
   real(wp), parameter :: BED_DEPTH = 12.0_wp
   real(wp), parameter :: CFL = 0.45_wp
   real(wp), parameter :: S_BOTTOM = 35.0_wp
   real(wp), parameter :: S_SURFACE = 30.0_wp
   real(wp), parameter :: T_BOTTOM = 15.0_wp
   real(wp), parameter :: T_SURFACE = 20.0_wp

   integer, allocatable :: nx_sweep(:), ny_sweep(:), nz_sweep(:)
   integer :: timed_steps, warmup_steps, ig, inz
   character(len=32) :: arg

   if (command_argument_count() >= 2) then
      allocate (nx_sweep(1), ny_sweep(1))
      call get_command_argument(1, arg); read (arg, *) nx_sweep(1)
      call get_command_argument(2, arg); read (arg, *) ny_sweep(1)
   else
#if defined(RAKALI_GPU_OFFLOAD)
      nx_sweep = [256, 512, 1024]
      ny_sweep = [256, 512, 1024]
#else
      nx_sweep = [64, 128]
      ny_sweep = [64, 128]
#endif
   end if

   if (command_argument_count() >= 3) then
      allocate (nz_sweep(1))
      call get_command_argument(3, arg); read (arg, *) nz_sweep(1)
   else
      nz_sweep = [5, 10, 20]
   end if

   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg); read (arg, *) timed_steps
   else
#if defined(RAKALI_GPU_OFFLOAD)
      timed_steps = 200
#else
      timed_steps = 25
#endif
   end if
   warmup_steps = max(5, timed_steps/10)

   write (*, '(a)') "backend,nx,ny,ncells,nz,warmup,timed_steps,"// &
      "total_s,mcells_s,us_per_step"

   do ig = 1, size(nx_sweep)
      do inz = 1, size(nz_sweep)
         call run_benchmark(nx_sweep(ig), ny_sweep(ig), nz_sweep(inz), &
                            warmup_steps, timed_steps)
      end do
   end do

contains

   subroutine run_benchmark(nx_phys, ny_phys, nz, n_warmup, n_timed)
      integer, intent(in) :: nx_phys, ny_phys, nz, n_warmup, n_timed

      type(state_t)    :: state
      type(timer_type) :: t_total
      integer :: n, ncells
      real(wp) :: dx, dy, elapsed, mcells, us_per_step
      character(len=256) :: csv

      dx = DOMAIN_LX/real(nx_phys, wp)
      dy = DOMAIN_LY/real(ny_phys, wp)

      ! Multilayer must be flagged before state%init so the multilayer
      ! sub-object allocates the per-layer arrays + tracer registry.
      state%use_multilayer = .true.
      state%multilayer%nz_ml = nz
      state%multilayer%initial_salinity = 0.5_wp*(S_BOTTOM + S_SURFACE)
      state%multilayer%initial_temperature = 0.5_wp*(T_BOTTOM + T_SURFACE)
      state%multilayer%S_ref = 0.0_wp
      state%multilayer%T_ref = 15.0_wp
      state%multilayer%beta_S = 0.78_wp
      state%multilayer%alpha_T = 0.17_wp
      state%multilayer%rho_0 = 1025.0_wp
      state%multilayer%hdiff_kappa = 1.0_wp
      state%multilayer%smagorinsky_cs = 0.15_wp
      state%multilayer%use_ppm_tracer = .true.
      state%multilayer%use_bpg_jacobian = .true.

      call state%init(nx_phys, ny_phys, NGHOST, dx, dy, CFL)
      ncells = nx_phys*ny_phys

      call set_bowl_bathymetry(state)
      call init_stratified(state, nz)

      call solver_enter_data(state)

      call profiler_init()
      call profiler_reset()

      call logger%info("--- ML-str "//to_string(nx_phys)//"x"//to_string(ny_phys) &
                       //" nz="//to_string(nz) &
                       //" ["//BACKEND//"] " &
                       //to_string(n_warmup)//"w+"//to_string(n_timed)//"t ---")

      do n = 1, n_warmup
         call solver_step(state)
      end do
      call gpu_sync()

      ! Reset profiler so warmup time isn't counted.
      call profiler_reset()

      call t_total%start()
      do n = 1, n_timed
         call solver_step(state)
      end do
      call gpu_sync()
      call t_total%stop()
      elapsed = t_total%get_elapsed_time()

      mcells = real(ncells, wp)*real(n_timed, wp)/elapsed/1.0e6_wp
      us_per_step = elapsed/real(n_timed, wp)*1.0e6_wp

      call logger%info("  Total time:  "//to_string(elapsed)//" s")
      call logger%info("  Mcells/s:    "//to_string(mcells))
      call logger%info("  us / step:   "//to_string(us_per_step))

      call profiler_report("Per-kernel breakdown")
      call profiler_end()

      write (csv, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",'// &
             'es12.5,",",es12.5,",",es12.5)') &
         BACKEND, nx_phys, ny_phys, ncells, nz, n_warmup, n_timed, &
         elapsed, mcells, us_per_step
      write (*, '(a)') trim(csv)

      call solver_exit_data(state)
      call state%destroy()
   end subroutine run_benchmark

   subroutine set_bowl_bathymetry(state)
      !! Same gentle radial bowl as the unstr sister: 12 m at the
      !! centre, 6 m at the edges.  Keeps the domain fully wet so the
      !! protect kernel doesn't dominate.
      type(state_t), intent(inout) :: state

      integer :: i, j, ng
      real(wp) :: cx, cy, r2, r2_max, frac, depth

      ng = state%grid%nghost
      r2_max = (0.5_wp*DOMAIN_LX)**2 + (0.5_wp*DOMAIN_LY)**2

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            cx = (real(i - ng, wp) - 0.5_wp)*state%grid%dx - 0.5_wp*DOMAIN_LX
            cy = (real(j - ng, wp) - 0.5_wp)*state%grid%dy - 0.5_wp*DOMAIN_LY
            r2 = cx*cx + cy*cy
            frac = min(r2/r2_max, 1.0_wp)
            depth = BED_DEPTH - 0.5_wp*BED_DEPTH*frac
            ! Convention: b is the bed elevation (negative below sea level)
            state%barotropic%b(i, j) = -depth
            state%barotropic%h(i, j) = depth        ! still water at z = 0
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do
   end subroutine set_bowl_bathymetry

   subroutine init_stratified(state, nz)
      !! Seed layer arrays from the barotropic IC, then fill the tracer
      !! registry with constant-S / constant-T columns (same convention
      !! as the river-plume test: scalar `initial_*` is the column
      !! average; vertical structure is not seeded explicitly here).
      !! `solver_step` will compute the EOS each stage so initial
      !! `rho_layer` only matters for the very first BPG.
      type(state_t), intent(inout) :: state
      integer, intent(in) :: nz

      call ml_init_layers_from_barotropic(state%barotropic%h, state%barotropic%hu, &
                                          state%barotropic%hv, state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, nz)

      call ml_init_tracer(state%multilayer%h_layer, &
                          state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                          state%multilayer%initial_salinity, state%grid, nz)
      call ml_init_tracer(state%multilayer%h_layer, &
                          state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                          state%multilayer%initial_temperature, state%grid, nz)

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                state%multilayer%rho_0, state%multilayer%beta_S, state%multilayer%S_ref, &
                                state%multilayer%alpha_T, state%multilayer%T_ref, &
                                state%grid, nz)
   end subroutine init_stratified

   subroutine gpu_sync()
   end subroutine gpu_sync

end program bench_multilayer_str
