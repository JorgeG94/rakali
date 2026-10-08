!! Single-GPU performance benchmark
program bench_single_gpu
   !! Measures cell-updates/second for the shallow water solver
   !!
   !! Usage:
   !!   bench_single_gpu              Run default sweep for the backend
   !!   bench_single_gpu <N>          Run single grid size N x N
   !!   bench_single_gpu <N> <steps>  Run N x N with custom timed step count
   !!
   !! Human-readable output goes to stderr (via logger).
   !! CSV output goes to stdout.
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_boundary, only: apply_boundary_conditions
   use rki_kernel_flux, only: compute_flux_hllc
   use rki_kernel_update, only: compute_update, compute_rk2_average
   use rki_kernel_timestep, only: compute_timestep
   use rki_kernel_friction, only: apply_friction
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   implicit none

   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = 1.0_wp
   real(wp), parameter :: DY = 1.0_wp
   real(wp), parameter :: CFL = 0.45_wp

   character(len=*), parameter :: BACKEND = &
#if defined(RAKALI_GPU_OFFLOAD)
                                  "openacc"
#else
   "cpu_serial"
#endif

   integer :: ig, n_sizes, timed_steps, warmup_steps
   integer, allocatable :: grid_sizes(:)
   character(len=32) :: arg

   ! Parse command-line arguments
   if (command_argument_count() >= 1) then
      allocate (grid_sizes(1))
      call get_command_argument(1, arg)
      read (arg, *) grid_sizes(1)
      n_sizes = 1
   else
#if defined(RAKALI_GPU_OFFLOAD)
      grid_sizes = [512, 1024, 2048, 4096]
#else
      grid_sizes = [64, 128, 256, 512]
#endif
      n_sizes = size(grid_sizes)
   end if

   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      read (arg, *) timed_steps
   else
#if defined(RAKALI_GPU_OFFLOAD)
      timed_steps = 1000
#else
      timed_steps = 100
#endif
   end if

   warmup_steps = max(10, timed_steps/10)

   ! CSV header
   write (*, '(a)') "backend,nx,ny,warmup,timed_steps,"// &
      "total_s,mcells_s,"// &
      "flux_s,update_s,timestep_s,friction_s,bc_s,save_s,rk2_s,"// &
      "flux_pct,update_pct,timestep_pct,friction_pct,bc_pct,save_pct,rk2_pct"

   do ig = 1, n_sizes
      call run_benchmark(grid_sizes(ig), grid_sizes(ig), warmup_steps, timed_steps)
   end do

contains

   subroutine gpu_sync()
      !! Synchronize GPU before reading timer (no-op, synchronous)
   end subroutine gpu_sync

   subroutine run_benchmark(nx_phys, ny_phys, n_warmup, n_timed)
      integer, intent(in) :: nx_phys, ny_phys, n_warmup, n_timed

      type(state_t) :: state
      type(timer_type) :: t_total, t_flux, t_update, t_timestep
      type(timer_type) :: t_friction, t_bc, t_save, t_rk2
      integer :: n, i, j, nx, ny, ng
      real(wp) :: elapsed, mcells
      real(wp) :: s_flux, s_update, s_timestep, s_friction, s_bc, s_save, s_rk2
      character(len=256) :: csv_line

      call state%init(nx_phys, ny_phys, NGHOST, DX, DY, CFL)

      nx = state%grid%nx_total
      ny = state%grid%ny_total
      ng = state%grid%nghost

      ! Dam break initial condition
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

      call logger%info("--- "//to_string(nx_phys)//" x "//to_string(ny_phys) &
                       //" ["//BACKEND//"] " &
                       //to_string(n_warmup)//"w+"//to_string(n_timed)//"t ---")

      ! Warmup (use solver_step, no timing needed)
      do n = 1, n_warmup
         call solver_step(state)
      end do
      call gpu_sync()

      ! Reset per-kernel accumulators
      s_flux = 0.0_wp
      s_update = 0.0_wp
      s_timestep = 0.0_wp
      s_friction = 0.0_wp
      s_bc = 0.0_wp
      s_save = 0.0_wp
      s_rk2 = 0.0_wp

      ! Timed run with per-kernel instrumentation
      call t_total%start()

      do n = 1, n_timed

         ! ---- BCs + timestep from U^n ----
         call t_bc%start()
         call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                        state%grid, state%bc)
         call gpu_sync()
         call t_bc%stop()
         s_bc = s_bc + t_bc%get_elapsed_time()

         call t_timestep%start()
         call compute_timestep(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                               state%grid, state%cfl, state%dt)
         call gpu_sync()
         call t_timestep%stop()
         s_timestep = s_timestep + t_timestep%get_elapsed_time()

         ! ---- Save U^n (OpenACC only; CUDA fuses save into update) ----
         call t_save%start()
         call save_state_bench(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                               state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, nx, ny)
         call gpu_sync()
         call t_save%stop()
         s_save = s_save + t_save%get_elapsed_time()

         ! ---- Stage 1: U* = U^n + dt * L(U^n) ----
         call t_flux%start()
         call compute_flux_hllc(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                                state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                                state%grid)
         call gpu_sync()
         call t_flux%stop()
         s_flux = s_flux + t_flux%get_elapsed_time()

         call t_update%start()
         call compute_update(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                             state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                             state%grid, state%dt)
         call gpu_sync()
         call t_update%stop()
         s_update = s_update + t_update%get_elapsed_time()

         call t_friction%start()
         call apply_friction(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                             state%grid, state%dt, &
                             state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                             state%forcing%atm_pressure, state%forcing%precipitation)
         call gpu_sync()
         call t_friction%stop()
         s_friction = s_friction + t_friction%get_elapsed_time()

         ! ---- Stage 2: U** = U* + dt * L(U*) ----
         call t_bc%start()
         call apply_boundary_conditions(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                        state%grid, state%bc)
         call gpu_sync()
         call t_bc%stop()
         s_bc = s_bc + t_bc%get_elapsed_time()

         call t_flux%start()
         call compute_flux_hllc(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, state%barotropic%b, &
                                state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                                state%barotropic%mass_flux_x, state%barotropic%mass_flux_y, &
                                state%grid)
         call gpu_sync()
         call t_flux%stop()
         s_flux = s_flux + t_flux%get_elapsed_time()

         call t_rk2%start()
         call compute_update(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                             state%barotropic%flux_h, state%barotropic%flux_hu, state%barotropic%flux_hv, &
                             state%grid, state%dt)

         call apply_friction(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                             state%grid, state%dt, &
                             state%barotropic%manning_n, state%forcing%wind_x, state%forcing%wind_y, &
                             state%forcing%atm_pressure, state%forcing%precipitation)

         call compute_rk2_average(state%barotropic%h, state%barotropic%hu, state%barotropic%hv, &
                                  state%barotropic%h0, state%barotropic%hu0, state%barotropic%hv0, &
                                  state%grid)
         call gpu_sync()
         call t_rk2%stop()
         s_rk2 = s_rk2 + t_rk2%get_elapsed_time()

      end do

      call gpu_sync()
      call t_total%stop()
      elapsed = t_total%get_elapsed_time()

      call solver_exit_data(state)

      ! Compute metrics
      mcells = real(nx_phys, wp)*real(ny_phys, wp)*real(n_timed, wp) &
               /elapsed/1.0e6_wp

      ! Human-readable summary (to stderr via logger)
      call logger%info("  Total time:  "//to_string(elapsed)//" s")
      call logger%info("  Mcells/s:    "//to_string(mcells))
      call logger%info("  Kernel breakdown (% of total):")
      call logger%info("    flux:      "//to_string(100.0_wp*s_flux/elapsed)//"%" &
                       //"  ("//to_string(s_flux*1.0e3_wp/real(2*n_timed, wp))//" ms/call, HLLC)")
      call logger%info("    update:    "//to_string(100.0_wp*s_update/elapsed)//"%")
      call logger%info("    timestep:  "//to_string(100.0_wp*s_timestep/elapsed)//"%")
      call logger%info("    friction:  "//to_string(100.0_wp*s_friction/elapsed)//"%")
      call logger%info("    bc:        "//to_string(100.0_wp*s_bc/elapsed)//"%")
      call logger%info("    save:      "//to_string(100.0_wp*s_save/elapsed)//"%")
      call logger%info("    rk2_avg:   "//to_string(100.0_wp*s_rk2/elapsed)//"%")

      ! CSV line (to stdout)
      write (csv_line, '(a,",",i0,",",i0,",",i0,",",i0,",",'// &
             'es12.5,",",es12.5,",",'// &
             '7(es12.5,","),'// &
             '6(f6.2,","),f6.2)') &
         BACKEND, nx_phys, ny_phys, n_warmup, n_timed, &
         elapsed, mcells, &
         s_flux, s_update, s_timestep, s_friction, s_bc, s_save, s_rk2, &
         100.0_wp*s_flux/elapsed, &
         100.0_wp*s_update/elapsed, &
         100.0_wp*s_timestep/elapsed, &
         100.0_wp*s_friction/elapsed, &
         100.0_wp*s_bc/elapsed, &
         100.0_wp*s_save/elapsed, &
         100.0_wp*s_rk2/elapsed
      write (*, '(a)') trim(csv_line)

      call state%destroy()

   end subroutine run_benchmark

   subroutine save_state_bench(h, hu, hv, h0, hu0, hv0, nx, ny)
      !! Copy current state into saved arrays (on device)
      integer, intent(in) :: nx, ny
      real(wp), intent(in) :: h(nx, ny), hu(nx, ny), hv(nx, ny)
      real(wp), intent(out) :: h0(nx, ny), hu0(nx, ny), hv0(nx, ny)

      integer :: i, j

      do concurrent(j=1:ny, i=1:nx)
         h0(i, j) = h(i, j)
         hu0(i, j) = hu(i, j)
         hv0(i, j) = hv(i, j)
      end do

   end subroutine save_state_bench

end program bench_single_gpu
