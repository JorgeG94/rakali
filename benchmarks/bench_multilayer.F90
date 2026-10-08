!! Full-physics multilayer benchmark, unstructured grid (single GPU)
program bench_multilayer
   !! Drives `solver_unstr_step` in multilayer mode on a regular
   !! triangulated mesh with a stratified initial condition.  Exercises
   !! the full kernel chain that's active in a real estuary run:
   !!   extrapolate -> central-upwind flux -> fused update
   !!     -> tracer advection (MUSCL) -> distribute barotropic
   !!     -> Smagorinsky horizontal viscosity
   !!     -> vertical exchange (momentum + tracers)
   !!     -> PP81 vertical viscosity + tracer diffusion
   !!     -> RK2 average -> clamp -> EOS -> BPG (Jacobian)
   !!     -> recompute barotropic momentum
   !!
   !! Runs each configured (ncells, nz) combination, reports
   !! cell-updates/s and us/step.  No NetCDF, no MPI.
   !!
   !! For per-kernel breakdowns, the unstructured solver currently has
   !! no `profiler_start/stop` calls inside its step (the structured
   !! sister bench `bench_multilayer_str` gets one for free).  Add the
   !! profiler regions to `solver_unstr_step` if you need finer detail
   !! out of this bench.
   !!
   !! Usage:
   !!   bench_multilayer                          Default sweep
   !!   bench_multilayer <nx_rect> <ny_rect>      Single grid, default nz
   !!   bench_multilayer <nx_rect> <ny_rect> <nz> Single config
   !!   bench_multilayer <nx_rect> <ny_rect> <nz> <steps>  Custom step count
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_enter_data, mesh_unstr_exit_data, &
                             mesh_unstr_destroy
   use rki_mesh_gen_unstr, only: gen_regular_triangulation
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_destroy, state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_solver_unstr, only: solver_unstr_step, solver_unstr_init_multilayer
   use rki_kernel_timestep_unstr, only: compute_timestep_unstr
   use rki_kernel_extrapolate_unstr, only: extrapolate_second_order_edge_unstr
   use rki_kernel_flux_unstr, only: compute_fluxes_central_unstr
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

   ! ---- Args ----
   if (command_argument_count() >= 2) then
      allocate (nx_sweep(1), ny_sweep(1))
      call get_command_argument(1, arg); read (arg, *) nx_sweep(1)
      call get_command_argument(2, arg); read (arg, *) ny_sweep(1)
   else
#if defined(RAKALI_GPU_OFFLOAD)
      nx_sweep = [128, 256, 512]
      ny_sweep = [128, 256, 512]
#else
      nx_sweep = [32, 64]
      ny_sweep = [32, 64]
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

   write (*, '(a)') "backend,nx_rect,ny_rect,ncells,nz,warmup,timed_steps,"// &
      "total_s,mcells_s,us_per_step"

   do ig = 1, size(nx_sweep)
      do inz = 1, size(nz_sweep)
         call run_benchmark(nx_sweep(ig), ny_sweep(ig), nz_sweep(inz), &
                            warmup_steps, timed_steps)
      end do
   end do

contains

   subroutine run_benchmark(nx_rect, ny_rect, nz, n_warmup, n_timed)
      integer, intent(in) :: nx_rect, ny_rect, nz, n_warmup, n_timed

      type(mesh_unstr_t)  :: mesh
      type(unstr_state_t) :: state
      type(timer_type)    :: t_total
      real(wp), allocatable :: vert_x(:), vert_y(:)
      integer, allocatable :: tri_v(:, :)
      integer :: ncells, n
      real(wp) :: elapsed, mcells, us_per_step
      character(len=256) :: csv

      call gen_regular_triangulation(nx_rect, ny_rect, DOMAIN_LX, DOMAIN_LY, vert_x, vert_y, tri_v)
      call mesh_unstr_init_from_triangles(mesh, vert_x, vert_y, tri_v)
      ncells = mesh%ncells_owned

      ! Configure ML knobs BEFORE state_unstr_init so they take effect
      ! during multilayer allocation / tracer-registry setup.
      call state_unstr_init(state, mesh, nz_ml=nz)
      state%cfl = CFL
      state%multilayer%use_bpg_jacobian = .true.
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

      ! Bathymetry: smooth radial bowl, 12 m at the centre, 6 m at the
      ! edges -- enough variation to stress BPG and BC handling without
      ! introducing wet/dry cells.
      call set_bowl_bathymetry(state, mesh)

      ! Stage = 0 (still water surface)
      state%barotropic%stage_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%barotropic%height_cv = state%barotropic%stage_cv - state%barotropic%bed_cv

      ! Push to device and seed the multilayer fields (h_layer, hu_layer,
      ! hv_layer + tracers).  `solver_unstr_init_multilayer` reads
      ! `initial_salinity` / `initial_temperature` and fills hS/hT.
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_init_multilayer(state, mesh)

      ! Bootstrap: max_speed must be populated before the CFL kernel.
      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)
      call compute_timestep_unstr(state, mesh)

      call logger%info("--- ML "//to_string(nx_rect)//"x"//to_string(ny_rect) &
                       //" rect ("//to_string(ncells)//" tris)" &
                       //" nz="//to_string(nz) &
                       //" ["//BACKEND//"] " &
                       //to_string(n_warmup)//"w+"//to_string(n_timed)//"t ---")
      call logger%info("  Initial dt = "//to_string(state%dt)//" s")

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

      write (csv, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",'// &
             'es12.5,",",es12.5,",",es12.5)') &
         BACKEND, nx_rect, ny_rect, ncells, nz, n_warmup, n_timed, &
         elapsed, mcells, us_per_step
      write (*, '(a)') trim(csv)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
      deallocate (vert_x, vert_y, tri_v)
   end subroutine run_benchmark

   subroutine set_bowl_bathymetry(state, mesh)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in)    :: mesh

      real(wp) :: cx, cy, r2, r2_max, frac, depth
      integer :: k

      r2_max = (0.5_wp*DOMAIN_LX)**2 + (0.5_wp*DOMAIN_LY)**2
      do k = 1, mesh%ncells_owned
         cx = mesh%centroid_x(k) - 0.5_wp*DOMAIN_LX
         cy = mesh%centroid_y(k) - 0.5_wp*DOMAIN_LY
         r2 = cx*cx + cy*cy
         frac = min(r2/r2_max, 1.0_wp)
         depth = BED_DEPTH - 0.5_wp*BED_DEPTH*frac  ! 12 m centre -> 6 m edges
         state%barotropic%bed_cv(k) = -depth
      end do
   end subroutine set_bowl_bathymetry

   subroutine gpu_sync()
   end subroutine gpu_sync

end program bench_multilayer
