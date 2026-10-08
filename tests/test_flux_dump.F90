!! Unit tests for the optional NetCDF flux I/O snapshot (rki_flux_dump)
module test_flux_dump
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_kernel_flux, only: compute_flux_hll, compute_flux_hllc
   use rki_solver, only: solver_enter_data, solver_exit_data
   use rki_flux_dump, only: dump_flux_io
   use netcdf, only: nf90_open, nf90_close, nf90_nowrite, nf90_noerr, &
                     nf90_inq_varid, nf90_get_var, nf90_inq_dimid, &
                     nf90_inquire_dimension, nf90_get_att, nf90_global
   implicit none
   private

   public :: collect_flux_dump_tests

contains

   subroutine collect_flux_dump_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dam_break_hll_snapshot", test_dam_break_hll), &
                  new_unittest("dam_break_hllc_snapshot", test_dam_break_hllc), &
                  new_unittest("lake_at_rest_snapshot", test_lake_at_rest_dump) &
                  ]
   end subroutine collect_flux_dump_tests

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

   subroutine fill_dam_break(s)
      !! Same input as test_flux's test_symmetric_dam_break: h_L=2, h_R=1,
      !! u=v=0, flat bed, dam at i = nx/2.
      type(state_t), intent(inout) :: s
      integer :: i, j, i_mid, nx, ny

      nx = s%grid%nx_total
      ny = s%grid%ny_total
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      s%barotropic%b = 0.0_wp
      i_mid = nx/2
      do j = 1, ny
         do i = 1, nx
            if (i <= i_mid) then
               s%barotropic%h(i, j) = 2.0_wp
            else
               s%barotropic%h(i, j) = 1.0_wp
            end if
         end do
      end do
   end subroutine fill_dam_break

   subroutine compute_flux_on_device(s)
      !! Mimic the per-RK-stage flux dispatch from rki_solver.F90:99/:107
      !! so the test exercises the same data path the dump consumes.
      type(state_t), intent(inout) :: s

      if (s%barotropic%use_hllc) then
         call compute_flux_hllc(s%barotropic%h, s%barotropic%hu, s%barotropic%hv, s%barotropic%b, &
                                s%barotropic%flux_h, s%barotropic%flux_hu, s%barotropic%flux_hv, &
                                s%barotropic%mass_flux_x, s%barotropic%mass_flux_y, s%grid)
      else
         call compute_flux_hll(s%barotropic%h, s%barotropic%hu, s%barotropic%hv, s%barotropic%b, &
                               s%barotropic%flux_h, s%barotropic%flux_hu, s%barotropic%flux_hv, &
                               s%barotropic%mass_flux_x, s%barotropic%mass_flux_y, s%grid)
      end if
   end subroutine compute_flux_on_device

   subroutine check_input_var(error, ncid, varname, expected, ok)
      !! Read a full-grid input variable from `ncid` and confirm it matches.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: varname
      real(wp), intent(in) :: expected(:, :)
      logical, intent(out) :: ok

      integer :: varid, status
      real(wp), allocatable :: buf(:, :)

      ok = .false.
      status = nf90_inq_varid(ncid, varname, varid)
      call check(error, status == nf90_noerr, "find input variable "//varname)
      if (allocated(error)) return

      allocate (buf(size(expected, 1), size(expected, 2)))
      status = nf90_get_var(ncid, varid, buf)
      call check(error, status == nf90_noerr, "read input variable "//varname)
      if (allocated(error)) return

      call check(error, maxval(abs(buf - expected)) < 1.0e-14_wp, &
                 "input "//varname//" must round-trip exactly")
      if (allocated(error)) return
      ok = .true.
   end subroutine check_input_var

   subroutine check_output_var(error, ncid, varname, expected, tol, ok)
      !! Read an interior-only output variable and compare to the value
      !! the kernel produced in-state.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: varname
      real(wp), intent(in) :: expected(:, :)
      real(wp), intent(in) :: tol
      logical, intent(out) :: ok

      integer :: varid, status
      real(wp), allocatable :: buf(:, :)

      ok = .false.
      status = nf90_inq_varid(ncid, varname, varid)
      call check(error, status == nf90_noerr, "find output variable "//varname)
      if (allocated(error)) return

      allocate (buf(size(expected, 1), size(expected, 2)))
      status = nf90_get_var(ncid, varid, buf)
      call check(error, status == nf90_noerr, "read output variable "//varname)
      if (allocated(error)) return

      call check(error, maxval(abs(buf - expected)) < tol, &
                 "output "//varname//" must match the kernel result")
      if (allocated(error)) return
      ok = .true.
   end subroutine check_output_var

   subroutine test_dam_break_hll(error)
      !! Dam-break state, HLL solver, snapshotted between flux and update.
      !! Verify: (a) input arrays round-trip exactly, (b) output arrays
      !! match the in-state kernel result bit-exactly, (c) the analytic
      !! HLL mass flux sqrt(g/2) appears at the dam cells, (d) global
      !! metadata is correct.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX_PHYS = 12, NY_PHYS = 4, NG = 2
      integer, parameter :: NX = NX_PHYS + 2*NG, NY = NY_PHYS + 2*NG
      type(state_t) :: s
      character(len=*), parameter :: fname = "test_flux_dump_hll.nc"
      character(len=8) :: solver_attr
      integer :: ncid, dimid, status, nghost_attr, i_mid, j
      integer :: dim_xt, dim_yt, dim_xi, dim_yi
      real(wp) :: expected_mass_flux
      logical :: ok

      call s%init(NX_PHYS, NY_PHYS, NG, 1.0_wp, 1.0_wp, 0.45_wp)
      s%barotropic%use_hllc = .false.
      call fill_dam_break(s)

      ! Mirror the solver: arrays go on device, kernel populates flux_*
      ! and mass_flux_* there.  dump_flux_io then snapshots the same
      ! resident bytes — no second kernel call.
      call solver_enter_data(s)
      call compute_flux_on_device(s)
      call dump_flux_io(fname, s)
      call solver_exit_data(s)

      status = nf90_open(fname, nf90_nowrite, ncid)
      call check(error, status == nf90_noerr, "open dump file")
      if (allocated(error)) goto 100

      status = nf90_inq_dimid(ncid, "x_total", dimid)
      status = nf90_inquire_dimension(ncid, dimid, len=dim_xt)
      call check(error, dim_xt == NX, "x_total dim should be NX")
      if (allocated(error)) goto 100
      status = nf90_inq_dimid(ncid, "y_total", dimid)
      status = nf90_inquire_dimension(ncid, dimid, len=dim_yt)
      call check(error, dim_yt == NY, "y_total dim should be NY")
      if (allocated(error)) goto 100
      status = nf90_inq_dimid(ncid, "x", dimid)
      status = nf90_inquire_dimension(ncid, dimid, len=dim_xi)
      call check(error, dim_xi == NX_PHYS, "x dim should be nx_phys")
      if (allocated(error)) goto 100
      status = nf90_inq_dimid(ncid, "y", dimid)
      status = nf90_inquire_dimension(ncid, dimid, len=dim_yi)
      call check(error, dim_yi == NY_PHYS, "y dim should be ny_phys")
      if (allocated(error)) goto 100

      call check_input_var(error, ncid, "h", s%barotropic%h, ok); if (.not. ok) goto 100
      call check_input_var(error, ncid, "hu", s%barotropic%hu, ok); if (.not. ok) goto 100
      call check_input_var(error, ncid, "hv", s%barotropic%hv, ok); if (.not. ok) goto 100
      call check_input_var(error, ncid, "b", s%barotropic%b, ok); if (.not. ok) goto 100

      call check_output_var(error, ncid, "flux_h", &
                            s%barotropic%flux_h(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 100
      call check_output_var(error, ncid, "flux_hu", &
                            s%barotropic%flux_hu(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 100
      call check_output_var(error, ncid, "flux_hv", &
                            s%barotropic%flux_hv(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 100
      call check_output_var(error, ncid, "mass_flux_x", &
                            s%barotropic%mass_flux_x(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 100
      call check_output_var(error, ncid, "mass_flux_y", &
                            s%barotropic%mass_flux_y(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 100

      ! Sanity: the in-state flux_h at the dam matches the analytic HLL
      ! value — same property test_symmetric_dam_break asserts in test_flux.
      expected_mass_flux = sqrt(0.5_wp*GRAVITY)
      i_mid = NX/2
      do j = NG + 1, NY - NG
         call check(error, &
                    abs(s%barotropic%flux_h(i_mid, j) + expected_mass_flux) < 1.0e-12_wp, &
                    "flux_h at dam must equal -sqrt(g/2)")
         if (allocated(error)) goto 100
      end do

      status = nf90_get_att(ncid, nf90_global, "nghost", nghost_attr)
      call check(error, status == nf90_noerr .and. nghost_attr == NG, &
                 "nghost global attribute should equal NG")
      if (allocated(error)) goto 100

      solver_attr = ""
      status = nf90_get_att(ncid, nf90_global, "solver", solver_attr)
      call check(error, status == nf90_noerr, "read solver global attribute")
      if (allocated(error)) goto 100
      call check(error, trim(solver_attr) == "HLL", &
                 "solver global attribute should be HLL")

100   continue
      status = nf90_close(ncid)
      call s%destroy()
      call delete_file(fname)
   end subroutine test_dam_break_hll

   subroutine test_dam_break_hllc(error)
      !! Same input but HLLC solver — solver attribute and outputs differ.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX_PHYS = 12, NY_PHYS = 4, NG = 2
      integer, parameter :: NX = NX_PHYS + 2*NG, NY = NY_PHYS + 2*NG
      type(state_t) :: s
      character(len=*), parameter :: fname = "test_flux_dump_hllc.nc"
      character(len=8) :: solver_attr
      integer :: ncid, status
      logical :: ok

      call s%init(NX_PHYS, NY_PHYS, NG, 1.0_wp, 1.0_wp, 0.45_wp)
      s%barotropic%use_hllc = .true.
      call fill_dam_break(s)

      call solver_enter_data(s)
      call compute_flux_on_device(s)
      call dump_flux_io(fname, s)
      call solver_exit_data(s)

      status = nf90_open(fname, nf90_nowrite, ncid)
      call check(error, status == nf90_noerr, "open HLLC dump file")
      if (allocated(error)) goto 200

      solver_attr = ""
      status = nf90_get_att(ncid, nf90_global, "solver", solver_attr)
      call check(error, status == nf90_noerr, "read solver global attribute")
      if (allocated(error)) goto 200
      call check(error, trim(solver_attr) == "HLLC", &
                 "solver attribute should be HLLC when use_hllc=.true.")
      if (allocated(error)) goto 200

      call check_output_var(error, ncid, "flux_h", &
                            s%barotropic%flux_h(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 200
      call check_output_var(error, ncid, "flux_hu", &
                            s%barotropic%flux_hu(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 200
      call check_output_var(error, ncid, "flux_hv", &
                            s%barotropic%flux_hv(NG + 1:NX - NG, NG + 1:NY - NG), &
                            1.0e-14_wp, ok); if (.not. ok) goto 200

200   continue
      status = nf90_close(ncid)
      call s%destroy()
      call delete_file(fname)
   end subroutine test_dam_break_hllc

   subroutine test_lake_at_rest_dump(error)
      !! Lake at rest over variable bathymetry — the well-balanced
      !! property must survive the snapshot: every interior output cell
      !! must read back as near-zero.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NX_PHYS = 10, NY_PHYS = 10, NG = 2
      integer, parameter :: NX = NX_PHYS + 2*NG, NY = NY_PHYS + 2*NG
      type(state_t) :: s
      real(wp), allocatable :: buf(:, :)
      character(len=*), parameter :: fname = "test_flux_dump_lake.nc"
      real(wp) :: eta_0, max_flux
      integer :: ncid, varid, status, i, j

      call s%init(NX_PHYS, NY_PHYS, NG, 1.0_wp, 1.0_wp, 0.45_wp)
      s%barotropic%use_hllc = .false.

      eta_0 = 2.0_wp
      s%barotropic%hu = 0.0_wp
      s%barotropic%hv = 0.0_wp
      do j = 1, NY
         do i = 1, NX
            s%barotropic%b(i, j) = 0.5_wp* &
                                   exp(-0.5_wp*real((i - NX/2)**2 + (j - NY/2)**2, wp))
            s%barotropic%h(i, j) = max(0.0_wp, eta_0 - s%barotropic%b(i, j))
         end do
      end do

      call solver_enter_data(s)
      call compute_flux_on_device(s)
      call dump_flux_io(fname, s)
      call solver_exit_data(s)

      status = nf90_open(fname, nf90_nowrite, ncid)
      call check(error, status == nf90_noerr, "open lake dump file")
      if (allocated(error)) goto 300

      allocate (buf(NX_PHYS, NY_PHYS))
      max_flux = 0.0_wp

      status = nf90_inq_varid(ncid, "flux_h", varid)
      status = nf90_get_var(ncid, varid, buf)
      max_flux = max(max_flux, maxval(abs(buf)))

      status = nf90_inq_varid(ncid, "flux_hu", varid)
      status = nf90_get_var(ncid, varid, buf)
      max_flux = max(max_flux, maxval(abs(buf)))

      status = nf90_inq_varid(ncid, "flux_hv", varid)
      status = nf90_get_var(ncid, varid, buf)
      max_flux = max(max_flux, maxval(abs(buf)))

      call check(error, max_flux < 1.0e-10_wp, &
                 "lake-at-rest snapshot must report near-zero fluxes (well-balanced)")

300   continue
      status = nf90_close(ncid)
      call s%destroy()
      call delete_file(fname)
   end subroutine test_lake_at_rest_dump

end module test_flux_dump
