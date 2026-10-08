!! Unstructured central-upwind solver benchmark (single GPU)
program bench_unstr
   !! Measures cell-updates/second for the unstructured KNP solver on
   !! a programmatically built regular right-triangle mesh.  Each
   !! Nx × Ny rectangle is split into two CCW triangles, yielding
   !! 2*Nx*Ny owned cells with the same valence everywhere -- a
   !! deliberately simple stress test for the cell-centric flux
   !! pattern, with no irregular cache misses.
   !!
   !! What the numbers reflect:
   !!   * Per-step the solver visits each cell's three edges twice
   !!     (interior edges are computed from both sides), runs the
   !!     central-upwind flux, and writes back to its own arrays only
   !!     -- no atomics.  That double-visit cost is real and shows up
   !!     in the bandwidth column.
   !!   * On a regular grid every triangle has identical geometry, so
   !!     the cell loop is fully coalesced if the renumberer placed
   !!     adjacent cells contiguously.  Compare against an unrenumbered
   !!     mesh by editing `apply_renumber` below.
   !!
   !! Usage:
   !!   bench_unstr                        Default sweep
   !!   bench_unstr <nx> <ny>              Single grid
   !!   bench_unstr <nx> <ny> <steps>      Custom step count
   !!
   !! Human-readable lines on stderr (logger).  CSV on stdout.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_enter_data, mesh_unstr_exit_data, &
                             mesh_unstr_destroy
   use rki_mesh_gen_unstr, only: gen_regular_triangulation
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_destroy, state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_solver_unstr, only: solver_unstr_step
   use rki_kernel_timestep_unstr, only: compute_timestep_unstr
   use rki_kernel_extrapolate_unstr, only: extrapolate_second_order_edge_unstr
   use rki_kernel_flux_unstr, only: compute_fluxes_central_unstr
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   implicit none

   real(wp), parameter :: DOMAIN_LX = 1000.0_wp
   real(wp), parameter :: DOMAIN_LY = 1000.0_wp
   real(wp), parameter :: BED_DEPTH = 10.0_wp
   real(wp), parameter :: CFL = 0.45_wp

   character(len=*), parameter :: BACKEND = &
#if defined(RAKALI_GPU_OFFLOAD)
                                  "openacc"
#else
   "cpu_serial"
#endif

   integer, allocatable :: nx_sweep(:), ny_sweep(:)
   integer :: timed_steps, warmup_steps, ig
   character(len=32) :: arg

   if (command_argument_count() >= 2) then
      allocate (nx_sweep(1), ny_sweep(1))
      call get_command_argument(1, arg)
      read (arg, *) nx_sweep(1)
      call get_command_argument(2, arg)
      read (arg, *) ny_sweep(1)
   else
#if defined(RAKALI_GPU_OFFLOAD)
      nx_sweep = [128, 256, 512, 1024]
      ny_sweep = [128, 256, 512, 1024]
#else
      nx_sweep = [32, 64, 128]
      ny_sweep = [32, 64, 128]
#endif
   end if

   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read (arg, *) timed_steps
   else
#if defined(RAKALI_GPU_OFFLOAD)
      timed_steps = 500
#else
      timed_steps = 50
#endif
   end if
   warmup_steps = max(10, timed_steps/10)

   write (*, '(a)') "backend,nx_rect,ny_rect,ncells,warmup,timed_steps,"// &
      "total_s,mcells_s,us_per_step"

   do ig = 1, size(nx_sweep)
      call run_benchmark(nx_sweep(ig), ny_sweep(ig), warmup_steps, timed_steps)
   end do

contains

   subroutine run_benchmark(nx_rect, ny_rect, n_warmup, n_timed)
      integer, intent(in) :: nx_rect, ny_rect, n_warmup, n_timed

      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(timer_type) :: t_total
      real(wp), allocatable :: vert_x(:), vert_y(:)
      integer, allocatable :: tri_v(:, :)
      integer :: ncells, n
      real(wp) :: elapsed, mcells, us_per_step
      character(len=256) :: csv

      call gen_regular_triangulation(nx_rect, ny_rect, DOMAIN_LX, DOMAIN_LY, vert_x, vert_y, tri_v)
      call mesh_unstr_init_from_triangles(mesh, vert_x, vert_y, tri_v)
      ncells = mesh%ncells_owned

      ! State allocations (barotropic only -- the multilayer machinery
      ! has its own benches via the Python river-plume etc.)
      call state_unstr_init(state, mesh)

      ! Still-water IC: bed = -10 m, stage = 0 -> uniform 10 m depth.
      ! Tiny perturbation in one corner so the wave speeds are non-zero
      ! and the timestep is meaningful.
      state%barotropic%bed_cv = -BED_DEPTH
      state%barotropic%stage_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%barotropic%height_cv = BED_DEPTH
      state%cfl = CFL
      state%barotropic%stage_cv(1) = 0.05_wp ! seed perturbation so dt is finite

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      ! Bootstrap: populate max_speed for the first timestep
      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)
      call compute_timestep_unstr(state, mesh)

      call logger%info("--- "//to_string(nx_rect)//"x"//to_string(ny_rect)// &
                       " rect ("//to_string(ncells)//" tris) [" &
                       //BACKEND//"] " &
                       //to_string(n_warmup)//"w+"//to_string(n_timed)//"t ---")
      call logger%info("  Initial dt = "//to_string(state%dt)//" s")

      ! Warmup
      do n = 1, n_warmup
         call solver_unstr_step(state, mesh, state%dt)
      end do
      call gpu_sync()

      call t_total%start()
      do n = 1, n_timed
         call solver_unstr_step(state, mesh, state%dt)
      end do
      call gpu_sync()
      call t_total%stop()
      elapsed = t_total%get_elapsed_time()

      mcells = real(ncells, wp)*real(n_timed, wp)/elapsed/1.0e6_wp
      us_per_step = elapsed/real(n_timed, wp)*1.0e6_wp

      call logger%info("  Total time:  "//to_string(elapsed)//" s")
      call logger%info("  Mcells/s:    "//to_string(mcells))
      call logger%info("  us / step:   "//to_string(us_per_step))

      write (csv, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",'// &
             'es12.5,",",es12.5,",",es12.5)') &
         BACKEND, nx_rect, ny_rect, ncells, n_warmup, n_timed, &
         elapsed, mcells, us_per_step
      write (*, '(a)') trim(csv)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      deallocate (vert_x, vert_y, tri_v)
   end subroutine run_benchmark

   subroutine gpu_sync()
      !! No-op (acc directives synchronous under stdpar).
   end subroutine gpu_sync

end program bench_unstr
