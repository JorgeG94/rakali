!! End-to-end test that drives the real coastal compute-rank lifecycle.
!!
!! Unlike `test_driver_smoke` — which replicates the
!! init -> enter_data -> step -> diagnostic -> exit_data *sequence* by
!! calling the `solver_*` routines directly — this test invokes the actual
!! `driver_run` entry point with a coastal `config_t`.  That exercises the
!! orchestration in `driver_run_coastal` as a unit: domain decomposition,
!! global init + scatter, BC / forcing / nesting / output / gauge setup,
!! GPU `enter_data`, the time-stepping loop with snapshot scheduling, and
!! the finalisation (final snapshot, `exit_data`, diagnostics).  It is the
!! regression net for refactoring that routine — a control-flow or
!! init-ordering mistake surfaces here even when per-kernel tests stay green.
!!
!! Setup: a tiny closed dam break (wall BCs, flat bed) on a 2D barotropic
!! grid, writing NetCDF snapshots to the working directory.  The assertions
!! re-open that file and check:
!!   * the unlimited `time` dimension has >= 2 records — proves the loop ran,
!!     the snapshot-scheduling branch fired, and the final snapshot wrote;
!!   * total water mass `sum(h)*dx*dy` is invariant between the first and
!!     last recorded snapshot — proves the conservative solver actually
!!     advanced state through the driver and nothing leaked at the walls.
module test_driver_coastal
   use rki_constants, only: wp
   use rki_config, only: config_t
   use rki_driver, only: driver_run, configure_log_level
   use rki_output, only: output_rank_filename
   use rki_io_netcdf, only: nc_open_read, nc_close, nc_get_dim_len, &
                            nc_get_varid, nc_check
   use netcdf, only: nf90_get_var
   use testdrive, only: error_type, check, new_unittest, unittest_type
   implicit none
   private

   public :: collect_driver_coastal_tests

   integer, parameter :: NX = 16, NY = 16
   real(wp), parameter :: DX = 100.0_wp
   real(wp), parameter :: T_END = 60.0_wp
   real(wp), parameter :: DT_OUT = 20.0_wp

   ! driver_run needs the MPI comm-env initialised (the real app does
   ! this in app/main); testdrive runs all cases in one process, so we
   ! init once via this saved flag.

contains

   subroutine collect_driver_coastal_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("driver_run_coastal_e2e", test_coastal_e2e), &
                  new_unittest("driver_run_coastal_restart", test_coastal_restart), &
                  new_unittest("driver_run_coastal_gauges", test_coastal_gauges) &
                  ]
   end subroutine collect_driver_coastal_tests

   subroutine make_cfg(cfg)
      !! Minimal coastal config: small 2D barotropic dam break with NetCDF
      !! snapshots to the cwd.  Everything else (wall BCs, no forcing /
      !! nesting / gauges / restart / I/O server, CFL-based dt) falls out of
      !! the `config_t` defaults.
      type(config_t), intent(out) :: cfg
      cfg%sim_type = "coastal"
      cfg%nx = NX
      cfg%ny = NY
      cfg%dx = DX
      cfg%dy = DX
      cfg%nghost = 2
      cfg%t_end = T_END
      cfg%init_condition = "dam_break"
      cfg%output_to_file = .true.
      cfg%output_format = "netcdf"
      cfg%output_dir = "."
      cfg%output_interval = DT_OUT
   end subroutine make_cfg

   subroutine cleanup_file(filename)
      character(len=*), intent(in) :: filename
      integer :: u
      logical :: exists
      inquire (file=filename, exist=exists)
      if (exists) then
         open (newunit=u, file=filename, status='old')
         close (u, status='delete')
      end if
   end subroutine cleanup_file

   subroutine write_gauge_file(filename, name, x, y)
      !! Write a one-station gauge locations file ("name  x  y", # comments).
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: x, y
      integer :: u
      open (newunit=u, file=filename, status='replace', action='write')
      write (u, '(A)') '# gauge locations for test_driver_coastal'
      write (u, '(A,2(2X,F12.4))') trim(name), x, y
      close (u)
   end subroutine write_gauge_file

   subroutine test_coastal_e2e(error)
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=512) :: fn
      integer :: ncid, varid, nt
      real(wp) :: h_first(NX, NY), h_last(NX, NY)
      real(wp) :: mass_first, mass_last, tol

      ! Keep ctest output to the failure path; the driver is chatty at info.
      call configure_log_level("warning")
      call make_cfg(cfg)
      fn = output_rank_filename(".", "output", 0)
      call cleanup_file(trim(fn))

      ! Tolerance adapts to working precision: tight in double, loose enough
      ! in single that file round-trip round-off doesn't trip a false leak.
      tol = max(1.0e-9_wp, 1.0e2_wp*epsilon(1.0_wp))

      checks: block
         call driver_run(cfg)

         call nc_open_read(trim(fn), ncid)

         call nc_get_dim_len(ncid, "time", nt)
         call check(error, nt >= 2, &
                    "driver should write >= 2 snapshots (scheduled + final)")
         if (allocated(error)) then; call nc_close(ncid); exit checks; end if

         call nc_get_varid(ncid, "h", varid)
         call nc_check(nf90_get_var(ncid, varid, h_first, &
                                    start=[1, 1, 1], count=[NX, NY, 1]), &
                       "reading first h snapshot")
         call nc_check(nf90_get_var(ncid, varid, h_last, &
                                    start=[1, 1, nt], count=[NX, NY, 1]), &
                       "reading last h snapshot")
         call nc_close(ncid)

         mass_first = sum(h_first)*DX*DX
         mass_last = sum(h_last)*DX*DX

         call check(error, mass_first > 0.0_wp, &
                    "first recorded snapshot has positive total mass")
         if (allocated(error)) exit checks

         call check(error, abs(mass_last - mass_first)/mass_first < tol, &
                    "closed dam break conserves mass across the driver run")
      end block checks

      call cleanup_file(trim(fn))
   end subroutine test_coastal_e2e

   subroutine test_coastal_restart(error)
      !! Restart round-trip through the driver.  A cold run writes a restart
      !! file mid-loop (exercising the restart-write loop block: D->H,
      !! restart_write_local, H->D); a second run warm-starts from it
      !! (exercising the restart_read path in init_local_state).  Asserts the
      !! restart captured a positive time + mass, the warm run resumed after
      !! that time (not t=0), and mass is preserved across the round-trip.
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=512) :: out_fn, restart_fn
      integer :: ncid, varid, nt
      real(wp) :: t_restart, mass_restart, mass_warm, first_out_time, tol
      real(wp) :: h2d(NX, NY), tbuf(1)

      call configure_log_level("warning")
      tol = max(1.0e-9_wp, 1.0e2_wp*epsilon(1.0_wp))
      out_fn = output_rank_filename(".", "output", 0)
      restart_fn = output_rank_filename(".", "restart", 0)
      call cleanup_file(trim(out_fn))
      call cleanup_file(trim(restart_fn))

      checks: block
         ! --- Phase 1: cold run that writes a restart file mid-loop ---
         call make_cfg(cfg)
         cfg%restart_interval = DT_OUT
         call driver_run(cfg)

         call nc_open_read(trim(restart_fn), ncid)
         call nc_get_varid(ncid, "time", varid)
         call nc_check(nf90_get_var(ncid, varid, t_restart), "reading restart time")
         call nc_get_varid(ncid, "h", varid)
         call nc_check(nf90_get_var(ncid, varid, h2d), "reading restart h")
         call nc_close(ncid)
         mass_restart = sum(h2d)*DX*DX

         call check(error, t_restart > 0.0_wp, &
                    "restart file should record a positive simulation time")
         if (allocated(error)) exit checks
         call check(error, mass_restart > 0.0_wp, &
                    "restart file should hold positive total mass")
         if (allocated(error)) exit checks

         ! --- Phase 2: warm start from that restart ---
         call make_cfg(cfg)
         cfg%restart_file = trim(restart_fn)
         cfg%t_end = t_restart + 3.0_wp*DT_OUT
         call driver_run(cfg)

         call nc_open_read(trim(out_fn), ncid)
         call nc_get_dim_len(ncid, "time", nt)
         call check(error, nt >= 1, "warm-start run should write >= 1 snapshot")
         if (allocated(error)) then; call nc_close(ncid); exit checks; end if
         call nc_get_varid(ncid, "time", varid)
         call nc_check(nf90_get_var(ncid, varid, tbuf, start=[1], count=[1]), &
                       "reading first warm-start output time")
         first_out_time = tbuf(1)
         call nc_get_varid(ncid, "h", varid)
         call nc_check(nf90_get_var(ncid, varid, h2d, start=[1, 1, 1], count=[NX, NY, 1]), &
                       "reading first warm-start h")
         call nc_close(ncid)
         mass_warm = sum(h2d)*DX*DX

         call check(error, first_out_time > t_restart, &
                    "warm start should resume from the restart time, not t=0")
         if (allocated(error)) exit checks
         call check(error, abs(mass_warm - mass_restart)/mass_restart < tol, &
                    "mass preserved across the restart round-trip")
      end block checks

      call cleanup_file(trim(out_fn))
      call cleanup_file(trim(restart_fn))
   end subroutine test_coastal_restart

   subroutine test_coastal_gauges(error)
      !! Run with a gauge locations file so gauges are active, exercising the
      !! driver's per-step gauge-write loop block (incl. the D->H pulls of the
      !! gauge cells).  Asserts the gauge NetCDF is written with >= 1 time
      !! record and a positive interpolated surface elevation (wet cell).
      type(error_type), allocatable, intent(out) :: error
      type(config_t) :: cfg
      character(len=512) :: out_fn
      character(len=*), parameter :: gauge_fn = "test_driver_coastal_gauges.txt"
      character(len=*), parameter :: gauge_nc = "./gauges.nc"
      integer :: ncid, varid, nt
      real(wp) :: eta_last, ebuf(1)

      call configure_log_level("warning")
      out_fn = output_rank_filename(".", "output", 0)
      call cleanup_file(trim(out_fn))
      call cleanup_file(gauge_nc)
      call cleanup_file(gauge_fn)

      ! One gauge at the domain centre (domain is [0, NX*DX] x [0, NY*DX])
      call write_gauge_file(gauge_fn, "centre", 0.5_wp*NX*DX, 0.5_wp*NY*DX)

      checks: block
         call make_cfg(cfg)
         cfg%gauge_file = gauge_fn
         cfg%gauge_output_interval = DT_OUT
         call driver_run(cfg)

         call nc_open_read(gauge_nc, ncid)
         call nc_get_dim_len(ncid, "time", nt)
         call check(error, nt >= 1, "gauge output should have >= 1 time record")
         if (allocated(error)) then; call nc_close(ncid); exit checks; end if
         call nc_get_varid(ncid, "eta", varid)
         call nc_check(nf90_get_var(ncid, varid, ebuf, start=[1, nt], count=[1, 1]), &
                       "reading last gauge eta")
         eta_last = ebuf(1)
         call nc_close(ncid)

         call check(error, eta_last > 0.0_wp, &
                    "interpolated gauge eta should be positive in a wet cell")
      end block checks

      call cleanup_file(trim(out_fn))
      call cleanup_file(gauge_nc)
      call cleanup_file(gauge_fn)
   end subroutine test_coastal_gauges

end module test_driver_coastal
