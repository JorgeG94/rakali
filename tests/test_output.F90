!! Unit tests for the structured solution output and restart writer
module test_output
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_config, only: config_t
   use rki_init, only: init_from_config, restore_multilayer_after_restart
   use rki_output, only: output_state_t, output_init, output_write_snapshot, &
                         output_close, restart_write, restart_read, &
                         output_rank_filename, output_write_snapshot_from_arrays, &
                         write_ascii_output
   use netcdf, only: nf90_open, nf90_close, nf90_nowrite, nf90_noerr, &
                     nf90_inq_varid, nf90_get_var, nf90_inq_dimid, &
                     nf90_inquire_dimension
   implicit none
   private

   public :: collect_output_tests

   integer, parameter :: NX = 8, NY = 6, NG = 2, NZ_ML = 3
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp, CFL = 0.45_wp

contains

   subroutine collect_output_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("barotropic_round_trip", test_barotropic_round_trip), &
                  new_unittest("multilayer_round_trip", test_multilayer_round_trip), &
                  new_unittest("bed_sediment_output", test_bed_sediment_output), &
                  new_unittest("snapshot_from_arrays", test_snapshot_from_arrays), &
                  new_unittest("restart_round_trip_barotropic", test_restart_barotropic), &
                  new_unittest("restart_round_trip_multilayer", test_restart_multilayer), &
                  new_unittest("restart_restores_derived_sigma", test_restart_derived_sigma), &
                  new_unittest("restart_restores_derived_zstar_full", test_restart_derived_zstar), &
                  new_unittest("output_rank_filename_pattern", test_rank_filename), &
                  new_unittest("write_ascii_output_runs", test_write_ascii) &
                  ]
   end subroutine collect_output_tests

   ! ----------------------------------------------------------------------
   ! Helpers
   ! ----------------------------------------------------------------------

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: u
      logical :: exists
      inquire (file=trim(filename), exist=exists)
      if (exists) then
         open (newunit=u, file=trim(filename), status='old')
         close (u, status='delete')
      end if
   end subroutine delete_file

   subroutine fill_barotropic(s, time)
      type(state_t), intent(inout) :: s
      real(wp), intent(in) :: time
      integer :: i, j

      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            s%barotropic%h(i, j) = 2.0_wp + 0.1_wp*real(i, wp) + time
            s%barotropic%hu(i, j) = 0.5_wp*real(j, wp)
            s%barotropic%hv(i, j) = -0.25_wp*real(i, wp)
            s%barotropic%b(i, j) = -10.0_wp
         end do
      end do
   end subroutine fill_barotropic

   subroutine make_multilayer_state(s)
      type(state_t), intent(inout) :: s
      integer :: i, j, k

      s%use_multilayer = .true.
      s%multilayer%nz_ml = NZ_ML
      call s%init(NX, NY, NG, DX, DY, CFL)

      ! sigma layers — the multilayer init defaults to uniform sigma
      ! when state%init runs, so we just populate the prognostic fields.
      do k = 1, NZ_ML
         do j = 1, s%grid%ny_total
            do i = 1, s%grid%nx_total
               s%multilayer%h_layer(i, j, k) = 1.0_wp/real(NZ_ML, wp)*(2.0_wp + 0.1_wp*real(i, wp))
               s%multilayer%hu_layer(i, j, k) = 0.0_wp
               s%multilayer%hv_layer(i, j, k) = 0.0_wp
               s%multilayer%tracers(s%multilayer%idx_salinity)%hTr(i, j, k) = &
                  s%multilayer%h_layer(i, j, k)*35.0_wp
               s%multilayer%tracers(s%multilayer%idx_temperature)%hTr(i, j, k) = &
                  s%multilayer%h_layer(i, j, k)*15.0_wp
            end do
         end do
      end do

      call fill_barotropic(s, 0.0_wp)
   end subroutine make_multilayer_state

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_barotropic_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      type(output_state_t) :: ostate
      character(len=*), parameter :: fname = "test_output_baro.nc"
      integer :: ncid, dimid, varid, status, ntime
      real(wp) :: time_vals(2)
      checks: block

         call s%init(NX, NY, NG, DX, DY, CFL)
         call fill_barotropic(s, 0.0_wp)

         cfg%dx = DX
         cfg%dy = DY
         cfg%nghost = NG
         cfg%cfl = CFL

         call output_init(ostate, fname, s, cfg)
         call output_write_snapshot(ostate, s, 0.0_wp)
         call output_write_snapshot(ostate, s, 0.5_wp)
         call output_close(ostate)

         ! Reopen and check the time dimension grew to 2
         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open output file")
         if (allocated(error)) exit checks

         status = nf90_inq_dimid(ncid, "time", dimid)
         call check(error, status == nf90_noerr, "find time dim")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_inquire_dimension(ncid, dimid, len=ntime)
         call check(error, ntime == 2, "expected 2 time steps")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_inq_varid(ncid, "time", varid)
         status = nf90_get_var(ncid, varid, time_vals)
         call check(error, abs(time_vals(1) - 0.0_wp) < 1.0e-10_wp, "time(1) == 0")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if
         call check(error, abs(time_vals(2) - 0.5_wp) < 1.0e-10_wp, "time(2) == 0.5")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_close(ncid)

      end block checks
      continue
      call s%destroy()
      call delete_file(fname)
   end subroutine test_barotropic_round_trip

   subroutine test_multilayer_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      type(output_state_t) :: ostate
      character(len=*), parameter :: fname = "test_output_ml.nc"
      integer :: ncid, varid, status
      checks: block

         call make_multilayer_state(s)

         cfg%dx = DX
         cfg%dy = DY
         cfg%nghost = NG
         cfg%cfl = CFL
         cfg%use_multilayer = .true.
         cfg%nz_layers = NZ_ML

         call output_init(ostate, fname, s, cfg)
         call output_write_snapshot(ostate, s, 0.0_wp)
         call output_close(ostate)

         ! The multilayer file should expose a layer-thickness variable
         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open ml output")
         if (allocated(error)) exit checks

         status = nf90_inq_varid(ncid, "h_layer", varid)
         call check(error, status == nf90_noerr, "h_layer variable should exist in ML output")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_close(ncid)

      end block checks
      continue
      call s%destroy()
      call delete_file(fname)
   end subroutine test_multilayer_round_trip

   subroutine test_bed_sediment_output(error)
      !! With sediment on, the bed reservoir is written as a 2D (x, y, time)
      !! NetCDF variable so deposition / erosion can be visualised.  Check the
      !! variable exists and round-trips its (uniform) initial bed mass.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      type(output_state_t) :: ostate
      character(len=*), parameter :: fname = "test_output_bed.nc"
      real(wp), parameter :: BED0 = 4.0_wp
      integer :: ncid, varid, status
      real(wp) :: bed_read(NX, NY, 1)
      checks: block

         ! use_sediment before init so bed_sediment allocates (= BED0).
         s%use_multilayer = .true.
         s%multilayer%nz_ml = NZ_ML
         s%multilayer%use_sediment = .true.
         s%multilayer%sediment_bed_init = BED0
         call s%init(NX, NY, NG, DX, DY, CFL)
         call fill_barotropic(s, 0.0_wp)

         cfg%dx = DX
         cfg%dy = DY
         cfg%nghost = NG
         cfg%cfl = CFL
         cfg%use_multilayer = .true.
         cfg%nz_layers = NZ_ML

         call output_init(ostate, fname, s, cfg)
         call output_write_snapshot(ostate, s, 0.0_wp)
         call output_close(ostate)

         status = nf90_open(fname, nf90_nowrite, ncid)
         call check(error, status == nf90_noerr, "open bed output")
         if (allocated(error)) exit checks

         status = nf90_inq_varid(ncid, "bed_sediment", varid)
         call check(error, status == nf90_noerr, &
                    "bed_sediment variable should exist when sediment is on")
         if (allocated(error)) then
            status = nf90_close(ncid); exit checks
         end if

         status = nf90_get_var(ncid, varid, bed_read)
         call check(error, abs(bed_read(3, 2, 1) - BED0) < 1.0e-10_wp, &
                    "bed_sediment should round-trip its initial mass")

         status = nf90_close(ncid)

      end block checks
      continue
      call s%destroy()
      call delete_file(fname)
   end subroutine test_bed_sediment_output

   subroutine test_snapshot_from_arrays(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      type(config_t) :: cfg
      type(output_state_t) :: ostate
      character(len=*), parameter :: fname = "test_output_arrays.nc"
      real(wp) :: h_arr(NX, NY), hu_arr(NX, NY), hv_arr(NX, NY), b_arr(NX, NY)
      integer :: i, j

      call s%init(NX, NY, NG, DX, DY, CFL)
      call fill_barotropic(s, 0.0_wp)

      cfg%dx = DX
      cfg%dy = DY
      cfg%nghost = NG
      cfg%cfl = CFL

      call output_init(ostate, fname, s, cfg)

      ! Build interior-only host arrays (matches the async I/O contract)
      do j = 1, NY
         do i = 1, NX
            h_arr(i, j) = 3.0_wp
            hu_arr(i, j) = 0.0_wp
            hv_arr(i, j) = 0.0_wp
            b_arr(i, j) = -10.0_wp
         end do
      end do

      call output_write_snapshot_from_arrays(ostate, h_arr, hu_arr, hv_arr, b_arr, 0.0_wp)
      call output_close(ostate)

      call check(error, .true., "snapshot_from_arrays should run without error")

      call s%destroy()
      call delete_file(fname)
   end subroutine test_snapshot_from_arrays

   subroutine test_restart_barotropic(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s, s2
      character(len=*), parameter :: fname = "test_restart_baro.nc"
      real(wp) :: t_in, t_out
      logical :: loaded_ml
      integer :: i, j, ierr
      checks: block

         call s%init(NX, NY, NG, DX, DY, CFL)
         call fill_barotropic(s, 0.0_wp)
         t_in = 12.5_wp

         call restart_write(fname, s, t_in)

         call s2%init(NX, NY, NG, DX, DY, CFL)
         call restart_read(fname, s2, t_out, ierr=ierr, loaded_multilayer=loaded_ml)

         call check(error, ierr == 0, "restart_read ierr should be 0")
         if (allocated(error)) exit checks
         call check(error,.not. loaded_ml, "barotropic restart should not flag multilayer")
         if (allocated(error)) exit checks
         call check(error, abs(t_out - t_in) < 1.0e-12_wp, "time should round-trip")
         if (allocated(error)) exit checks

         do j = NG + 1, NG + NY
            do i = NG + 1, NG + NX
               call check(error, abs(s%barotropic%h(i, j) - s2%barotropic%h(i, j)) < 1.0e-12_wp, &
                          "h interior should round-trip")
               if (allocated(error)) exit checks
               call check(error, abs(s%barotropic%hu(i, j) - s2%barotropic%hu(i, j)) < 1.0e-12_wp, &
                          "hu interior should round-trip")
               if (allocated(error)) exit checks
            end do
         end do

      end block checks
      continue
      call s%destroy()
      call s2%destroy()
      call delete_file(fname)
   end subroutine test_restart_barotropic

   subroutine test_restart_multilayer(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s, s2
      character(len=*), parameter :: fname = "test_restart_ml.nc"
      real(wp) :: t_in, t_out
      logical :: loaded_ml
      integer :: i, j, k, ierr
      checks: block

         call make_multilayer_state(s)
         t_in = 7.0_wp

         call restart_write(fname, s, t_in)

         ! Fresh multilayer state for read-back
         s2%use_multilayer = .true.
         s2%multilayer%nz_ml = NZ_ML
         call s2%init(NX, NY, NG, DX, DY, CFL)

         call restart_read(fname, s2, t_out, ierr=ierr, loaded_multilayer=loaded_ml)

         call check(error, ierr == 0, "restart_read ierr should be 0")
         if (allocated(error)) exit checks
         call check(error, loaded_ml, "multilayer flag should be set on restart_read")
         if (allocated(error)) exit checks

         ! Spot-check one layer interior
         i = NG + 3; j = NG + 2; k = 2
         call check(error, abs(s%multilayer%h_layer(i, j, k) - s2%multilayer%h_layer(i, j, k)) < 1.0e-12_wp, &
                    "h_layer round-trip")
         if (allocated(error)) exit checks

         ! Every multilayer prognostic field, interior cells only
         call check(error, interior_maxdiff(s%multilayer%h_layer, s2%multilayer%h_layer) < 1.0e-12_wp, &
                    "h_layer full round-trip")
         if (allocated(error)) exit checks
         call check(error, interior_maxdiff(s%multilayer%hu_layer, s2%multilayer%hu_layer) < 1.0e-12_wp, &
                    "hu_layer full round-trip")
         if (allocated(error)) exit checks
         call check(error, interior_maxdiff(s%multilayer%hv_layer, s2%multilayer%hv_layer) < 1.0e-12_wp, &
                    "hv_layer full round-trip")
         if (allocated(error)) exit checks
         call check(error, interior_maxdiff( &
                    s%multilayer%tracers(s%multilayer%idx_salinity)%hTr, &
                    s2%multilayer%tracers(s2%multilayer%idx_salinity)%hTr) < 1.0e-12_wp, &
                    "hS_layer full round-trip")
         if (allocated(error)) exit checks
         call check(error, interior_maxdiff( &
                    s%multilayer%tracers(s%multilayer%idx_temperature)%hTr, &
                    s2%multilayer%tracers(s2%multilayer%idx_temperature)%hTr) < 1.0e-12_wp, &
                    "hT_layer full round-trip")
         if (allocated(error)) exit checks

      end block checks
      continue
      call s%destroy()
      call s2%destroy()
      call delete_file(fname)
   end subroutine test_restart_multilayer

   function interior_maxdiff(a, b) result(d)
      !! Max |a - b| over the interior (non-ghost) cells of two 3D fields.
      real(wp), intent(in) :: a(:, :, :), b(:, :, :)
      real(wp) :: d
      d = maxval(abs(a(NG + 1:NG + NX, NG + 1:NG + NY, :) &
                     - b(NG + 1:NG + NX, NG + 1:NG + NY, :)))
   end function interior_maxdiff

   subroutine restart_derived_case(vcoord, error)
      !! A restart file restores the prognostic multilayer fields but not the
      !! derived ones.  `restore_multilayer_after_restart` must rebuild the
      !! density, the reference-density profile (and, for z*-full, the
      !! per-column reference levels) so a warm start sees the same state as
      !! the cold start it was written from.  Stratified, sloping bed.
      character(len=*), intent(in) :: vcoord
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: s, s2
      type(config_t) :: cfg
      character(len=*), parameter :: fname = "test_restart_derived.nc"
      real(wp) :: t_out, drho, dref
      integer :: i, j, ierr
      logical :: loaded_ml
      checks: block

         cfg%nx = NX; cfg%ny = NY; cfg%nghost = NG
         cfg%dx = 100.0_wp; cfg%dy = 100.0_wp; cfg%cfl = CFL
         cfg%use_multilayer = .true.
         cfg%nz_layers = NZ_ML
         cfg%vcoord_type = vcoord
         cfg%init_condition = "still_water"
         cfg%eta_initial = 0.0_wp
         cfg%S_init_surface = 30.0_wp
         cfg%S_init_bottom = 35.0_wp
         cfg%T_init_surface = 20.0_wp
         cfg%T_init_bottom = 10.0_wp
         cfg%alpha_T = 0.2_wp
         cfg%beta_S = 0.8_wp

         call s%init_from_config(cfg)
         do j = 1, s%grid%ny_total
            do i = 1, s%grid%nx_total
               s%barotropic%b(i, j) = -(20.0_wp + 2.0_wp*real(i, wp))
            end do
         end do
         call init_from_config(s, cfg)

         call restart_write(fname, s, 11.0_wp)

         call s2%init_from_config(cfg)
         s2%barotropic%b = s%barotropic%b
         call restart_read(fname, s2, t_out, ierr=ierr, loaded_multilayer=loaded_ml)
         call check(error, ierr == 0 .and. loaded_ml, "restart_read should load the multilayer fields")
         if (allocated(error)) exit checks

         ! Before the restore the derived density is unset
         call restore_multilayer_after_restart(s2)

         drho = interior_maxdiff(s%multilayer%rho_layer, s2%multilayer%rho_layer)
         call check(error, drho < 1.0e-9_wp, &
                    "rho_layer after restart restore differs from the cold start")
         if (allocated(error)) exit checks

         dref = maxval(abs(s%multilayer%rho_ref - s2%multilayer%rho_ref))
         call check(error, maxval(abs(s%multilayer%rho_ref)) > 0.0_wp, &
                    "cold-start rho_ref should be populated")
         if (allocated(error)) exit checks
         call check(error, dref < 1.0e-9_wp, &
                    "rho_ref after restart restore differs from the cold start")
         if (allocated(error)) exit checks
         call check(error, abs(s%multilayer%rho_ref_dz - s2%multilayer%rho_ref_dz) < 1.0e-12_wp, &
                    "rho_ref_dz after restart restore differs from the cold start")
         if (allocated(error)) exit checks

         if (trim(vcoord) == "zstar_full") then
            call check(error, maxval(abs(s%multilayer%vcoord_z_ref_col)) > 0.0_wp, &
                       "cold-start z*-full reference levels should be populated")
            if (allocated(error)) exit checks
            call check(error, maxval(abs(s%multilayer%vcoord_z_ref_col &
                                         - s2%multilayer%vcoord_z_ref_col)) < 1.0e-12_wp, &
                       "z*-full reference levels after restart restore differ from the cold start")
            if (allocated(error)) exit checks
         end if

      end block checks
      continue
      call s%destroy()
      call s2%destroy()
      call delete_file(fname)
   end subroutine restart_derived_case

   subroutine test_restart_derived_sigma(error)
      type(error_type), allocatable, intent(out) :: error
      call restart_derived_case("sigma", error)
   end subroutine test_restart_derived_sigma

   subroutine test_restart_derived_zstar(error)
      type(error_type), allocatable, intent(out) :: error
      call restart_derived_case("zstar_full", error)
   end subroutine test_restart_derived_zstar

   subroutine test_rank_filename(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=512) :: fname

      fname = output_rank_filename(".", "out", 5)
      call check(error, index(trim(fname), "out_rank_000005.nc") > 0, &
                 "rank filename should embed 6-digit zero-padded rank")
      if (allocated(error)) return

      fname = output_rank_filename("/tmp/abc", "restart", 42)
      call check(error, index(trim(fname), "/tmp/abc/restart_rank_000042.nc") > 0, &
                 "rank filename should join base_dir, prefix and rank")
   end subroutine test_rank_filename

   subroutine test_write_ascii(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      character(len=*), parameter :: out_dir = "."
      character(len=*), parameter :: out_file = "./final_h.dat"
      logical :: exists

      call s%init(NX, NY, NG, DX, DY, CFL)
      call fill_barotropic(s, 0.0_wp)

      call write_ascii_output(s, out_dir)

      inquire (file=out_file, exist=exists)
      call check(error, exists, "ASCII output file should exist on disk")

      call s%destroy()
      call delete_file(out_file)
   end subroutine test_write_ascii

end module test_output
