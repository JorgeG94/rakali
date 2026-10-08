!! Driver: full run lifecycle for the Rakali solver.
module rki_driver
   !! Owns the end-to-end run path that used to live inline in
   !! `app/main.F90`: setup (state allocation / BC /
   !! forcing / nesting / output / GPU enter_data), the time-stepping
   !! loop (forcing update, solver_step, status print, snapshot output,
   !! restart, gauges), and finalize (exit_data, final snapshot,
   !! diagnostics, profiler report, resource cleanup).
   !!
   !! Exposed as `driver_run(cfg)` so:
   !!   * `app/main.F90` stays thin (argv parsing, dispatch).
   !!   * Unit tests can drive a full compute lifecycle from an in-code
   !!     `config_t`, catching integration-level bugs (OpenACC partial
   !!     presence, init ordering, missing attaches) that per-kernel
   !!     unit tests miss.
   use rki_constants, only: wp, &
                            KEPS_STAB_CONSTANT, KEPS_STAB_GALPERIN, KEPS_STAB_CANUTO
   use rki_vcoord, only: parse_vcoord_type, parse_remap_method
   use rki_timestepper, only: parse_time_scheme, TIME_RK2, TIME_SEMI_IMPLICIT
   use rki_config, only: config_t
   use rki_state, only: state_t, register_default_tracers, &
                        load_bathymetry, fill_bathymetry_ghosts
   use rki_grid, only: hgrid_t
   use rki_boundary_types, only: bc_init_from_config, bc_init_chapman
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_mesh_gen_unstr, only: gen_regular_triangulation, gen_delaunay_triangulation
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_enter_data, mesh_unstr_exit_data, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_enter_data, state_unstr_exit_data, &
                              state_unstr_destroy, state_unstr_enable_si
   use rki_si_unstr_state, only: parse_si_precond_unstr
   use rki_si_common, only: parse_si_bsub_fill
   use rki_solver_unstr, only: solver_unstr_advance, solver_unstr_init_multilayer
   use rki_console_stats, only: console_stats_t
   use rki_coastal_console_stats, only: coastal_console_stats_report
   ! (load_bathymetry/fill_bathymetry_ghosts now re-exported from rki_state)
   use rki_forcing, only: forcing_t, forcing_init_from_config, forcing_update, forcing_destroy
   use rki_nesting, only: nesting_t, nesting_init_from_config, nesting_update, nesting_destroy
   use rki_output, only: output_state_t, output_init, output_write_snapshot, &
                         output_close, restart_read, restart_write, &
                         output_rank_filename, &
                         write_ascii_output
   use rki_gauges, only: gauge_output_t, gauges_init, gauges_write, gauges_destroy
   use rki_init, only: init_from_config, compute_total_mass, &
                       init_multilayer_from_barotropic, restore_multilayer_after_restart, &
                       compute_layer_sums, diagnose_bathymetry
   use pic_logger, only: logger => global_logger, &
                         debug_level, verbose_level, info_level, &
                         performance_level, warning_level, error_level
   use pic_strings, only: to_string
   use pic_io, only: to_char_count
   use pic_types, only: int64
   use pic_timer, only: timer_type
   use rki_profiler, only: profiler_init, profiler_start, profiler_stop, &
                           profiler_report, profiler_end
   use rki_banner, only: print_banner
   use rki_mem_report, only: mem_host_rss_bytes, mem_log_state_budget, &
                             mem_log_device_actuals, mem_log_device_growth, &
                             mem_set_counted_budget, mem_log_computed_budget, &
                             mem_log_computed_line
   use pic_knowledge, only: get_knowledge
   implicit none
   private

   public :: driver_run
   public :: configure_log_level

contains

   subroutine driver_run(cfg)
      !! Compute-rank entry point.  Dispatches to the per-backend helper
      !! for the coastal HLL/HLLC + KNP dynamical core.  Any non-coastal
      !! `sim_type` is rejected fail-loud by `validate_config`; guard
      !! again here so a programmatic caller that skips validation can
      !! never silently fall through to the wrong dynamical core.
      type(config_t), intent(inout) :: cfg

      if (trim(cfg%sim_type) /= "coastal") then
         call logger%error("driver_run: sim_type='"//trim(cfg%sim_type)// &
                           "' is not supported by this build (coastal only)")
         error stop "driver_run: unsupported sim_type"
      end if
      if (trim(cfg%grid_type) == "unstructured") then
         call driver_run_coastal_unstr(cfg)
      else
         call driver_run_coastal(cfg)
      end if
   end subroutine driver_run

   subroutine driver_run_coastal_unstr(cfg)
      !! End-to-end unstructured coastal run on a native regular-triangulation
      !! mesh (Python-free).  Routed from `driver_run` when
      !! `grid_type=="unstructured"`; reuses `&grid_nml` (nx,ny,dx,dy) for the
      !! mesh and `nz_layers` for the column.  Flat bed, lake at rest,
      !! horizontally-uniform IC.
      type(config_t), intent(inout) :: cfg

      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      type(timer_type) :: wall_timer
      real(wp), allocatable :: vert_x(:), vert_y(:)
      integer, allocatable :: tri_v(:, :)
      integer :: ncells, nz
      real(wp) :: lx, ly, depth, elapsed
      logical :: si_ok

      nz = max(cfg%nz_layers, 1)
      lx = real(cfg%nx, wp)*cfg%dx
      ly = real(cfg%ny, wp)*cfg%dy
      depth = cfg%h0   ! flat bed at the configured still-water depth (h0)

      ! The semi-implicit TPFA operator needs an acute-Delaunay mesh: the
      ! plain right-triangle split has coincident circumcentres per quad
      ! (edge_delta = 0) and the SI quality gate rejects it.  The explicit
      ! path keeps the regular split (bit-identity for existing runs).
      if (parse_time_scheme(cfg%time_scheme) == TIME_SEMI_IMPLICIT) then
         call gen_delaunay_triangulation(cfg%nx, cfg%ny, lx, ly, vert_x, vert_y, tri_v)
      else
         call gen_regular_triangulation(cfg%nx, cfg%ny, lx, ly, vert_x, vert_y, tri_v)
      end if

      call mesh_unstr_init_from_triangles(mesh, vert_x, vert_y, tri_v)
      ncells = mesh%ncells_owned

      call logger%info("Coastal unstructured run: "//to_string(cfg%nx)//" x "// &
                       to_string(cfg%ny)//" regular mesh ("// &
                       to_string(2*cfg%nx*cfg%ny)//" triangles), nz = "// &
                       to_string(nz))

      ! Vertical coordinate is set before init -- it drives the multilayer
      ! vcoord setup inside state_unstr_init.
      state%multilayer%vcoord%coord_type = parse_vcoord_type(cfg%vcoord_type)
      ! Closure knobs that GATE allocations inside state_unstr_init must be
      ! set BEFORE it (the SI k-eps Picard save slots allocate only when
      ! vmix_use_keps .and. keps_si_picard > 1 — a post-init sync would
      ! leave them unallocated and the Picard loop would hit unallocated
      ! arrays).  Mirrors the structured init_from_config ordering.
      state%multilayer%vmix_use_keps = cfg%vmix_use_keps
      state%multilayer%vmix_use_kpp = cfg%vmix_use_kpp
      state%multilayer%keps_si_picard = cfg%keps_si_picard
      call state_unstr_init(state, mesh, nz_ml=nz)

      ! Physics from the namelist: the structured coastal knobs mapped onto the
      ! unstructured state (EOS + tracers + vmix + lateral mixing + forcing).
      state%cfl = cfg%cfl
      state%time_scheme = parse_time_scheme(cfg%time_scheme)

      ! The selectable explicit integrator (euler/rk3/rk43) is wired into
      ! the 2D barotropic path only; multilayer keeps SSP-RK2 (per-layer
      ! BPG / mixing stages sit between the Euler update and the combine).
      ! Semi-implicit is exempt: the SI step has its own 3D column solve.
      if (state%time_scheme /= TIME_RK2 .and. state%time_scheme /= TIME_SEMI_IMPLICIT &
          .and. state%use_multilayer) then
         error stop "time_scheme: euler/rk3/rk43 require the 2D barotropic "// &
            "unstructured path (nz_layers <= 1); multilayer runs use rk2."
      end if

      ! Semi-implicit (Casulli) mode: run the fail-loud TPFA mesh-quality
      ! gate and allocate the si slot now — after state_unstr_init (nz_ml
      ! known) and before state_unstr_enter_data (the slot maps with the
      ! state).  The knob envelope (dt_fixed > 0, wall/tidal/clamped edges,
      ! single rank, precond) was already validated at configure time.
      if (state%time_scheme == TIME_SEMI_IMPLICIT) then
         call state_unstr_enable_si(state, mesh, &
                                    parse_si_precond_unstr(cfg%si_precond), &
                                    cfg%si_rtol, cfg%dt_fixed, si_ok, &
                                    wetdry=cfg%si_wetdry, n_sub=cfg%si_wetdry_nsub, &
                                    bsub_fill=parse_si_bsub_fill(cfg%si_wetdry_bsub_fill))
         if (.not. si_ok) then
            error stop "driver: mesh failed the semi-implicit TPFA quality "// &
               "gate (non-Delaunay / obtuse cell pairs); refine the mesh or "// &
               "use the explicit path."
         end if
      end if
      state%multilayer%rho_0 = cfg%rho_0
      state%multilayer%alpha_T = cfg%alpha_T
      state%multilayer%beta_S = cfg%beta_S
      ! Coastal EOS selector (&tracer_nml eos_form) — mirrors the
      ! structured wiring in rki_state: default "linear" keeps the
      ! legacy kernel; "wright" / "roquet_spv" opt in to the shared
      ! nonlinear EOS (handle is a flat POD, no device map needed).
      block
         use rki_ml_eos, only: ml_parse_eos_form
         call ml_parse_eos_form(cfg%eos_form, &
                                state%multilayer%eos_use_shared, &
                                state%multilayer%eos%variant)
      end block
      state%multilayer%eos%rho0 = cfg%rho_0
      state%multilayer%eos%S_ref = cfg%S_ref
      state%multilayer%eos%T_ref = cfg%T_ref
      state%multilayer%eos%alpha_T = cfg%alpha_T
      state%multilayer%eos%beta_S = cfg%beta_S
      state%multilayer%eos%is_init = .true.
      state%multilayer%initial_salinity = cfg%initial_salinity
      state%multilayer%S_init_surface = cfg%S_init_surface
      state%multilayer%S_init_bottom = cfg%S_init_bottom
      state%multilayer%initial_temperature = cfg%initial_temperature
      state%multilayer%T_init_surface = cfg%T_init_surface
      state%multilayer%T_init_bottom = cfg%T_init_bottom
      state%multilayer%kappa_S_bg = cfg%kappa_S_bg
      state%multilayer%kappa_T_bg = cfg%kappa_T_bg
      state%multilayer%use_ppm_tracer = cfg%use_ppm_tracer
      ! (vmix_use_keps / vmix_use_kpp / keps_si_picard synced BEFORE
      ! state_unstr_init above — they gate init allocations.)
      state%multilayer%keps_stability = keps_stability_from_string(cfg%keps_stability)
      state%multilayer%smagorinsky_cs = cfg%smagorinsky_cs
      state%multilayer%evap_mass_sink = cfg%evap_mass_sink
      state%multilayer%evap_latent_cooling = cfg%evap_latent_cooling
      state%multilayer%hdiff_kappa = cfg%hdiff_kappa
      state%barotropic%manning_n = cfg%manning_n
      state%barotropic%wind_tau_x = cfg%wind_stress_x
      state%barotropic%wind_tau_y = cfg%wind_stress_y
      state%barotropic%coriolis_f = cfg%coriolis_f

      ! Flat bed, lake at rest.
      state%barotropic%bed_cv = -depth
      state%barotropic%stage_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      state%barotropic%height_cv = state%barotropic%stage_cv - state%barotropic%bed_cv

      call mesh_unstr_enter_data(mesh)
      call profiler_start("unstr_enter_data")
      call state_unstr_enter_data(state)
      call profiler_stop("unstr_enter_data")
      call solver_unstr_init_multilayer(state, mesh)

      call logger%info("Starting time evolution: target "//to_string(cfg%t_end)//" s")
      call wall_timer%start()
      call solver_unstr_advance(state, mesh, cfg%t_end)
      call wall_timer%stop()
      elapsed = real(wall_timer%get_elapsed_time(), wp)

      call profiler_start("unstr_exit_data")
      call state_unstr_exit_data(state)
      call profiler_stop("unstr_exit_data")
      call mesh_unstr_exit_data(mesh)

      call logger%info("Unstructured run complete: "//to_string(state%step_count)// &
                       " steps in "//to_string(elapsed)//" s")
      call report_throughput(real(ncells, wp), state%step_count, elapsed)

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine driver_run_coastal_unstr

   pure function keps_stability_from_string(name) result(code)
      !! Map the `keps_stability` namelist string onto the KEPS_STAB_* enum
      !! (mirrors the structured rki_state mapping).
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(name))
      case ("galperin")
         code = KEPS_STAB_GALPERIN
      case ("canuto")
         code = KEPS_STAB_CANUTO
      case default
         code = KEPS_STAB_CONSTANT
      end select
   end function keps_stability_from_string

   subroutine report_throughput(n_cells, n_steps, elapsed, total_mcells)
      !! Log horizontal throughput (cells x steps / wall time).  Throughput is
      !! horizontal cells x steps, NOT x nz (the codebase convention).
      !! Optionally returns the figure so a caller can thread it into a later
      !! summary.
      real(wp), intent(in) :: n_cells
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: elapsed
      real(wp), intent(out), optional :: total_mcells

      real(wp) :: total

      total = 0.0_wp
      if (elapsed > 0.0_wp) then
         total = n_cells*real(n_steps, wp)/elapsed/1.0e6_wp
      end if
      if (present(total_mcells)) total_mcells = total
      call logger%info("  throughput: "//to_string(total)//" Mcells/s")
   end subroutine report_throughput

   subroutine driver_run_coastal(cfg)
      !! Run the full lifecycle for the structured coastal regime.
      !!
      !! Allocates and initialises state, configures BCs / forcing / nesting / output / gauges, places
      !! data on GPU, runs the time-stepping loop with periodic
      !! diagnostics and output, then finalises (D2H, final snapshot,
      !! diagnostics, profiler, resource cleanup).
      type(config_t), intent(in) :: cfg

      type(state_t) :: state
      type(forcing_t) :: forc
      type(nesting_t) :: nest
      type(output_state_t) :: ostate
      type(timer_type) :: wall_timer, setup_timer
      type(gauge_output_t) :: gauges
      type(console_stats_t) :: console_stats
         !! Persistent MOM6-style console reporter (latches the t=0
         !! reference for the drift columns).

      character(len=512) :: rank_filename
      real(wp) :: t_current, total_mass_initial, total_mass_final
      real(wp) :: next_output_time, next_restart_time, next_status_time
      real(wp) :: next_gauge_time
      real(wp) :: wall_elapsed, mcells_per_sec
      integer :: n_steps, n_snapshots
      integer(int64) :: rss0, rss1
      logical :: use_netcdf_output, use_restart_output
      logical :: use_async_output

      call setup_timer%start()

      ! Host RSS before state allocation — the growth across init is an
      ! ADVISORY host-side figure + a conservative upper bound for the
      ! pre-flight OOM warning; the true device footprint is measured
      ! across enter_data (see rki_mem_report).
      rss0 = mem_host_rss_bytes()

      call init_local_state(cfg, state, t_current, total_mass_initial)

      ! Set boundary condition parameters and allocate Chapman arrays
      call init_bc_from_config(state, cfg)
      call init_bc_chapman(state)

      ! Initialise atmospheric forcing from NetCDF (if configured)
      call forcing_init_from_config(forc, cfg, state)
      state%forcing%active = forc%active

      ! Initialise one-way nesting from parent domain (if configured)
      call nesting_init_from_config(nest, cfg, state)

      ! Initialise NetCDF output
      use_netcdf_output = cfg%output_to_file .and. &
                          trim(cfg%output_format) == "netcdf"
      if (use_netcdf_output) then
         rank_filename = output_rank_filename(trim(cfg%output_dir), "output", 0)
         call output_init(ostate, trim(rank_filename), state, cfg)
      end if
      if (use_netcdf_output) then
         next_output_time = t_current + cfg%output_interval
      end if

      ! Initialise gauge point output (high-frequency time series)
      call gauges_init(gauges, cfg, state, 0)
      if (gauges%active) then
         if (cfg%gauge_output_interval > 0.0_wp) then
            next_gauge_time = t_current + cfg%gauge_output_interval
         else
            next_gauge_time = t_current
         end if
      end if

      use_async_output = .false.

      ! Restart output scheduling
      use_restart_output = (cfg%restart_interval > 0.0_wp)
      if (use_restart_output) then
         next_restart_time = t_current + cfg%restart_interval
      end if

      ! Status output scheduling
      if (cfg%status_interval > 0.0_wp) then
         next_status_time = t_current + cfg%status_interval
      end if

      ! Initialise profiler (NVTX annotations + wall-clock timers)
      call profiler_init()

      ! Memory budget: host growth since rss0 is the upper-bound estimate
      ! of what enter_data will map.  Log it (+ device free/total compare)
      ! BEFORE the mapping so an over-budget run reports a readable line
      ! instead of dying opaquely inside the CUDA runtime.
      rss1 = mem_host_rss_bytes()
      if (rss0 >= 0_int64 .and. rss1 >= 0_int64) then
         call mem_log_state_budget("coastal state", max(rss1 - rss0, 0_int64))
      end if

      ! Latch the counted state-array footprint for the enter_data
      ! reconciliation (mem_log_device_actuals reads it); the human-readable
      ! budget prints in the run-header banner (log_run_header).
      call mem_set_counted_budget(state%bytes())

      ! Place data on GPU
      call profiler_start("coastal_enter_data")
      call solver_enter_data(state)
      call profiler_stop("coastal_enter_data")

      call mem_log_device_actuals("coastal state mapped")

      call setup_timer%stop()
      call logger%info("Setup time: "//to_string(setup_timer%get_elapsed_time())//" s")

      ! Warmup step to trigger lazy GPU initialisation (JIT, context setup)
      if (t_current < cfg%t_end) then
         call solver_step(state, dt_max=(cfg%t_end - t_current), time=t_current)
         t_current = t_current + state%dt
      end if

      ! Banner before time evolution
      call log_time_evolution_banner(cfg, use_netcdf_output)

      ! Salt / Temp columns + drift lines only make sense with the
      ! multilayer tracers; a 2D-barotropic run reports Mass + En only.
      console_stats%report_thermodynamics = state%use_multilayer
      ! Coastal stability is gravity-wave-limited, so a high *advective*
      ! CFL is normal — keep the NaN guard, drop the CFL panic.
      console_stats%panic_on_cfl = .false.

      ! Main time-stepping loop
      n_steps = 1
      n_snapshots = 0
      call wall_timer%start()
      call profiler_start("time_loop")

      do while (t_current < cfg%t_end)
         ! Update atmospheric forcing arrays for current time
         if (forc%active) then
            call forcing_update(forc, state, t_current)
            !$omp target update to(state%forcing%wind_x, state%forcing%wind_y, &
            !$omp&              state%forcing%atm_pressure, &
            !$omp&              state%forcing%precipitation, &
            !$omp&              state%forcing%heat_flux, &
            !$omp&              state%forcing%evap_minus_precip)
         end if

         ! Update nested boundary ghost cells from parent domain
         if (nest%active) then
            !$omp target update from(state%barotropic%h, state%barotropic%hu, state%barotropic%hv)
            call nesting_update(nest, state, t_current)
            !$omp target update to(state%barotropic%h, state%barotropic%hu, state%barotropic%hv)
         end if

         call profiler_start("solver_step")
         call solver_step(state, dt_max=(cfg%t_end - t_current), time=t_current)
         call profiler_stop("solver_step")
         t_current = t_current + state%dt
         n_steps = n_steps + 1

         ! Periodic status output (table format)
         if (cfg%status_interval > 0.0_wp) then
            if (t_current >= next_status_time) then
               call emit_periodic_status(state, cfg, n_steps, t_current, &
                                         real(wall_timer%get_elapsed_time(), wp), &
                                         console_stats)
               ! Quiet drift check: one-shot warning if device memory grew
               ! > 15% past the post-enter_data baseline (lazy workspaces /
               ! leak) — no info lines at status cadence.
               call mem_log_device_growth("coastal status", quiet=.true.)
               next_status_time = next_status_time + cfg%status_interval
            end if
         else if (mod(n_steps, 100) == 0) then
            call print_status_line(n_steps, t_current, state%dt, &
                                   real(wall_timer%get_elapsed_time(), wp), cfg%t_end)
         end if

         ! NetCDF snapshot output
         if (use_netcdf_output .and. t_current >= next_output_time) then
            call flush_snapshot(state, ostate, cfg, t_current, next_output_time, &
                                n_snapshots)
         end if

         ! Restart file output
         if (use_restart_output .and. t_current >= next_restart_time) then
            call flush_restart(state, cfg, t_current, next_restart_time)
         end if

         ! Gauge point output (high-frequency time series)
         if (gauges%active) then
            if (cfg%gauge_output_interval <= 0.0_wp .or. t_current >= next_gauge_time) then
               call write_gauges(gauges, state, cfg, t_current, next_gauge_time)
            end if
         end if
      end do
      call profiler_stop("time_loop")
      wall_elapsed = wall_timer%get_elapsed_time()
      call logger%info("")
      call logger%info("Solver complete: "//to_string(n_steps)//" steps in "// &
                       to_string(wall_elapsed)//" s")
      call report_throughput(real(cfg%nx, wp)*real(cfg%ny, wp), n_steps, &
                             wall_elapsed, mcells_per_sec)
      ! End-of-run device-memory truth: state is still device-resident
      ! here (exit_data below) — the growth since enter_data attributes
      ! the lazily allocated kernel workspaces even on short runs.
      call mem_log_device_growth("coastal end of run")

      ! --- Finalisation (output flush, D2H, gather, diagnostics) ---
      call setup_timer%start()

      ! Compute final mass on device (GPU reduction, transfers only 1 scalar)
      total_mass_final = compute_total_mass(state)

      ! Retrieve data from GPU
      call profiler_start("coastal_exit_data")
      call solver_exit_data(state)
      call profiler_stop("coastal_exit_data")

      ! Write the final snapshot
      call finalize_output(state, ostate, t_current, use_netcdf_output)

      call setup_timer%stop()
      call logger%info("Finalisation time: "//to_string(setup_timer%get_elapsed_time())//" s")

      ! Final diagnostics
      call wall_timer%stop()
      call log_final_diagnostics(state, cfg, n_steps, wall_elapsed, &
                                 real(wall_timer%get_elapsed_time(), wp), &
                                 mcells_per_sec, total_mass_final, &
                                 total_mass_initial, use_netcdf_output)

      ! Print profiling summary
      call profiler_report("Compute", root_region="time_loop")
      call profiler_end()

      call get_knowledge()

      call nesting_destroy(nest)
      call forcing_destroy(forc)
      call gauges_destroy(gauges)
      call state%destroy()
   end subroutine driver_run_coastal

   subroutine log_run_header(cfg, total_mass_initial, state)
      !! Startup summary: banner + grid / spacing / layers / duration /
      !! mass / fixed-dt note / counted memory
      !! footprint.  Pure logger output.
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: total_mass_initial
      type(state_t), intent(in) :: state
         !! State the footprint is counted from.

      character(len=32) :: mass_str

      call print_banner(cfg%sim_type)
      call logger%info("")
      call logger%info("  Grid:     "//to_string(cfg%nx)//" x "//to_string(cfg%ny)// &
                       " ("//to_char_count(int(cfg%nx, int64)*int(cfg%ny, int64))//")")
      call logger%info("  Spacing:  dx = "//to_string(cfg%dx)//" m, dy = "// &
                       to_string(cfg%dy)//" m")
      if (cfg%use_multilayer) then
         call logger%info("  Layers:   "//to_string(cfg%nz_layers)//" (multilayer hydrostatic)")
      else if (cfg%use_nonhydrostatic) then
         call logger%info("  Layers:   "//to_string(cfg%nz_layers)//" (non-hydrostatic)")
      else
         call logger%info("  Layers:   1 (2D barotropic)")
      end if
      call logger%info("  Duration: "//to_string(cfg%t_end)//" s")
      write (mass_str, "(ES14.6)") total_mass_initial
      call logger%info("  Mass:     "//trim(adjustl(mass_str)))
      if (cfg%dt_fixed > 0.0_wp) then
         call logger%warning("  Fixed dt: "//to_string(cfg%dt_fixed)//" s (CFL disabled)")
      end if
      ! Counted state-array footprint (exact allocatable sum; gated-off
      ! sub-slots count 0, so the breakdown reflects conditional allocation).
      block
         integer(int64) :: b_total, b_baro, b_forcing
         b_total = state%bytes()
         b_baro = state%barotropic%bytes()
         b_forcing = state%forcing%bytes()
         call mem_log_computed_budget("coastal state", b_total)
         call mem_log_computed_line("barotropic", b_baro)
         call mem_log_computed_line("surface forcing", b_forcing)
         if (state%use_multilayer) then
            call mem_log_computed_line("multilayer + tracers", state%multilayer%bytes())
         end if
         if (state%use_nonhydrostatic) then
            call mem_log_computed_line("non-hydrostatic", state%nh%bytes())
         end if
      end block
      call logger%info("")
   end subroutine log_run_header

   subroutine log_time_evolution_banner(cfg, use_netcdf_output)
      !! "Starting time evolution" banner + status-table header.  On
      !! non-GPU builds, detects the multicore thread count from ACC_NUM_CORES
      !! (then OMP_NUM_THREADS).
      type(config_t), intent(in) :: cfg
      logical, intent(in) :: use_netcdf_output

      call logger%info("")
#ifdef RAKALI_GPU_OFFLOAD
      call logger%info("Starting time evolution: target "//to_string(cfg%t_end)// &
                       " s, GPU")
#else
      block
         integer :: n_threads
         character(len=64) :: env_val
         character(len=128) :: ios_msg
         integer :: env_len, env_stat, ios
         n_threads = 0
         ! NVHPC -stdpar=multicore / -acc=multicore reads ACC_NUM_CORES first,
         ! then falls back to OMP_NUM_THREADS. Gcc -ftree-parallelize-loops
         ! bakes the count in at compile time (RAKALI_NPROC) and ignores both.
         call get_environment_variable("ACC_NUM_CORES", env_val, env_len, env_stat)
         if (env_stat == 0 .and. env_len > 0) then
            read (env_val, *, iostat=ios, iomsg=ios_msg) n_threads
            if (ios /= 0) n_threads = 0
         end if
         if (n_threads <= 0) then
            call get_environment_variable("OMP_NUM_THREADS", env_val, env_len, env_stat)
            if (env_stat == 0 .and. env_len > 0) then
               read (env_val, *, iostat=ios, iomsg=ios_msg) n_threads
               if (ios /= 0) n_threads = 0
            end if
         end if
         if (n_threads > 1) then
            call logger%info("Starting time evolution: target "//to_string(cfg%t_end)// &
                             " s, "//to_string(n_threads)//" CPU threads")
         else
            call logger%info("Starting time evolution: target "//to_string(cfg%t_end)// &
                             " s, CPU")
         end if
      end block
#endif
      if (use_netcdf_output) then
         call logger%info("  Output every "//to_string(cfg%output_interval)//" s")
      end if
      call logger%info("")
      ! Slim table; the conserved quantities are
      ! reported below each row in the shared [stats] block.
      call logger%info("     Step         t (s)         dt (s)     Wall (s)  Remaining (s)")
      call logger%info("  -------  ------------  ------------  -----------  -------------")
   end subroutine log_time_evolution_banner

   subroutine emit_periodic_status(state, cfg, n_steps, t_current, &
                                   wall_elapsed, console_stats)
      !! Periodic status: the slim Step / t / dt / Wall / Remaining
      !! row, then the shared MOM6-style
      !! `[stats]` block (Mass / Salt / Heat drift + energy), computed +
      !! all-reduced by `coastal_console_stats_report`.  `state` is
      !! `intent(inout)` because the device reductions require it.
      type(state_t), intent(inout) :: state
      type(config_t), intent(in) :: cfg
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: t_current
      real(wp), intent(in) :: wall_elapsed
      type(console_stats_t), intent(inout) :: console_stats

      call print_status_line(n_steps, t_current, state%dt, &
                             wall_elapsed, cfg%t_end)
      call coastal_console_stats_report(console_stats, state, state%dt, &
                                        t_current, n_steps)
   end subroutine emit_periodic_status

   subroutine init_local_state(cfg, state, t_current, total_mass_initial)
      !! Build the ready-to-run state: allocate, load bathymetry, read a
      !! restart or set the initial conditions, log the run header and the
      !! bathymetry diagnostics, and fill the ghost bathymetry.
      type(config_t), intent(in) :: cfg
      type(state_t), intent(inout) :: state
      real(wp), intent(out) :: t_current
      real(wp), intent(out) :: total_mass_initial

      logical :: loaded_ml

      t_current = 0.0_wp
      total_mass_initial = 0.0_wp

      call state%init_from_config(cfg)

      ! Load bathymetry from file or use flat bottom
      if (len_trim(cfg%bathymetry_file) > 0) then
         call load_bathymetry(trim(cfg%bathymetry_file), state)
      end if

      ! Load restart file or set up initial conditions
      if (len_trim(cfg%restart_file) > 0) then
         call restart_read(trim(cfg%restart_file), state, t_current, &
                           loaded_multilayer=loaded_ml)
         call logger%info("Warm start from restart at t = "//to_string(t_current)//" s")
         if (state%use_multilayer) then
            if (loaded_ml) then
               call restore_multilayer_after_restart(state)
            else
               ! Barotropic-only restart under a multilayer config: lay the
               ! layers out from the restored barotropic state.
               call init_multilayer_from_barotropic(state, cfg)
            end if
         end if
      else
         call init_from_config(state, cfg)
      end if

      total_mass_initial = compute_total_mass(state)
      call diagnose_bathymetry(state, cfg)
      call log_run_header(cfg, total_mass_initial, state)

      ! Extrapolate domain-edge ghost bathymetry
      call fill_bathymetry_ghosts(state)
   end subroutine init_local_state
   subroutine log_final_diagnostics(state, cfg, n_steps, wall_elapsed, wall_total, &
                                    mcells_per_sec, total_mass_final, &
                                    total_mass_initial, use_netcdf_output)
      !! Rank-0 end-of-run summary: step / time / throughput / mass-change log
      !! lines, the ASCII fallback dump when NetCDF output is off, and the
      !! merge hint.  Caller stops the wall timer and passes its elapsed time.
      type(state_t), intent(in) :: state
      type(config_t), intent(in) :: cfg
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: wall_elapsed
      real(wp), intent(in) :: wall_total
      real(wp), intent(in) :: mcells_per_sec
      real(wp), intent(in) :: total_mass_final
      real(wp), intent(in) :: total_mass_initial
      logical, intent(in) :: use_netcdf_output

      call logger%info("Total steps: "//to_string(n_steps))
      call logger%info("Solver time: "//to_string(wall_elapsed)//" s")
      call logger%info("Total wall time: "//to_string(wall_total)//" s")
      call logger%info("Mcells/s: "//to_string(mcells_per_sec))
      call logger%info("Final mass: "//to_string(total_mass_final))
      call logger%info("Mass change: "// &
                       to_string(abs(total_mass_final - total_mass_initial)))

      ! ASCII fallback output
      if (.not. use_netcdf_output) then
         call write_ascii_output(state, cfg%output_dir)
      end if

   end subroutine log_final_diagnostics

   subroutine flush_snapshot(state, ostate, cfg, t_current, next_output_time, &
                             n_snapshots)
      !! Flush one NetCDF snapshot mid-run: D->H, write locally, H->D, bump
      !! the snapshot counter + cadence cursor.  Caller
      !! gates on `use_netcdf_output .and. t_current >= next_output_time`.
      type(state_t), intent(inout) :: state
      type(output_state_t), intent(inout) :: ostate
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: t_current
      real(wp), intent(inout) :: next_output_time
      integer, intent(inout) :: n_snapshots

      call profiler_start("sync_output")
      call solver_exit_data(state)
      call output_write_snapshot(ostate, state, t_current)
      call solver_enter_data(state)
      call profiler_stop("sync_output")
      n_snapshots = n_snapshots + 1
      call logger%info("  >> Snapshot "//to_string(n_snapshots)// &
                       " at t = "//to_string(t_current)//" s")
      next_output_time = next_output_time + cfg%output_interval
   end subroutine flush_snapshot

   subroutine finalize_output(state, ostate, t_current, use_netcdf_output)
      !! End-of-run output: write the final snapshot and close the file.
      !! `state` has already been pulled D->H by the caller.
      type(state_t), intent(inout) :: state
      type(output_state_t), intent(inout) :: ostate
      real(wp), intent(in) :: t_current
      logical, intent(in) :: use_netcdf_output

      if (use_netcdf_output) then
         call output_write_snapshot(ostate, state, t_current)
         call output_close(ostate)
      end if
   end subroutine finalize_output

   subroutine flush_restart(state, cfg, t_current, next_restart_time)
      !! Write a restart file: D->H, restart_write, H->D, then advance the
      !! restart cadence cursor.  Caller gates on
      !! `use_restart_output .and. t_current >= next_restart_time`.
      type(state_t), intent(inout) :: state
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: t_current
      real(wp), intent(inout) :: next_restart_time

      character(len=512) :: rank_filename

      call solver_exit_data(state)
      rank_filename = output_rank_filename(trim(cfg%output_dir), "restart", 0)
      call restart_write(trim(rank_filename), state, t_current)
      call solver_enter_data(state)
      next_restart_time = next_restart_time + cfg%restart_interval
   end subroutine flush_restart
   subroutine write_gauges(gauges, state, cfg, t_current, next_gauge_time)
      !! Pull the gauge cells D->H and append one gauge time-series record,
      !! advancing the cadence cursor.  Caller gates on gauges%active and the
      !! gauge cadence.
      type(gauge_output_t), intent(inout) :: gauges
      type(state_t), intent(inout) :: state
      type(config_t), intent(in) :: cfg
      real(wp), intent(in) :: t_current
      real(wp), intent(inout) :: next_gauge_time

      integer :: g, i0, j0

      do g = 1, gauges%n_local
         i0 = gauges%ix0(g)
         j0 = gauges%iy0(g)
         !$omp target update from(state%barotropic%h(i0:i0 + 1, j0:j0 + 1))
         !$omp target update from(state%barotropic%hu(i0:i0 + 1, j0:j0 + 1))
         !$omp target update from(state%barotropic%hv(i0:i0 + 1, j0:j0 + 1))
      end do
      call gauges_write(gauges, state, t_current)
      if (cfg%gauge_output_interval > 0.0_wp) then
         next_gauge_time = next_gauge_time + cfg%gauge_output_interval
      end if
   end subroutine write_gauges

   ! ================================================================
   ! Helper subroutines (private to the driver module)
   ! ================================================================

   subroutine print_status_line(step, t, dt, wall, t_end, mass, mean_S, mean_T)
      !! Print a table-formatted status line with estimated remaining time.
      !! When mean_S / mean_T are present, two extra columns are written
      !! between Mass and Wall so the multilayer header lines up.
      integer, intent(in) :: step
      real(wp), intent(in) :: t, dt, wall, t_end
      real(wp), intent(in), optional :: mass
      real(wp), intent(in), optional :: mean_S, mean_T

      real(wp) :: remaining
      character(len=256) :: line

      if (t > 0.0_wp) then
         remaining = wall*(t_end - t)/t
      else
         remaining = 0.0_wp
      end if

      if (present(mass) .and. present(mean_S) .and. present(mean_T)) then
         write (line, "(I9,F14.2,F14.4,ES13.3,F13.4,F13.4,F13.2,F15.2)") &
            step, t, dt, mass, mean_S, mean_T, wall, remaining
      else if (present(mass)) then
         write (line, "(I9,F14.2,F14.4,ES13.3,F13.2,F15.2)") &
            step, t, dt, mass, wall, remaining
      else
         write (line, "(I9,F14.2,F14.4,13X,F13.2,F15.2)") &
            step, t, dt, wall, remaining
      end if
      call logger%info("  "//trim(line))
   end subroutine print_status_line

   subroutine configure_log_level(level_str)
      !! Set logger verbosity from string
      character(len=*), intent(in) :: level_str

      select case (trim(level_str))
      case ("debug")
         call logger%configure(debug_level)
      case ("verbose")
         call logger%configure(verbose_level)
      case ("info")
         call logger%configure(info_level)
      case ("performance")
         call logger%configure(performance_level)
      case ("warning")
         call logger%configure(warning_level)
      case ("error")
         call logger%configure(error_level)
      case default
         call logger%configure(info_level)
      end select
   end subroutine configure_log_level

   subroutine init_bc_from_config(state, cfg)
      !! Initialise boundary condition parameters from config
      type(state_t), intent(inout) :: state
      type(config_t), intent(in) :: cfg

      call bc_init_from_config(state%bc, cfg)
   end subroutine init_bc_from_config

   subroutine init_bc_chapman(state)
      !! Allocate and initialise Chapman BC arrays from state
      type(state_t), intent(inout) :: state

      call bc_init_chapman(state%bc, state%barotropic%h, state%barotropic%b, &
                           state%grid%nx_total, state%grid%ny_total, state%grid%nghost)
   end subroutine init_bc_chapman

end module rki_driver
