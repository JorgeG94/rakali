!! Optional NetCDF snapshot of the barotropic flux-kernel per-cell I/O
module rki_flux_dump
   !! Diagnostic snapshot of the barotropic state arrays as they sit
   !! between the flux call and the update call inside an RK stage.
   !!
   !! Reads the existing `state%barotropic%{h, hu, hv, b, flux_*,
   !! mass_flux_*}` arrays that the solver has already populated on
   !! device — does **not** re-run the kernel.  The dump must be
   !! inserted between `solver_compute_barotropic_flux` and
   !! `compute_update` so that `(h, hu, hv, b)` are still the flux
   !! kernel's inputs and `flux_*` / `mass_flux_*` are the matching
   !! outputs.  After `compute_update`, `h/hu/hv` are overwritten.
   !!
   !! Input fields are written over the full ghost-padded array shape
   !! `(nx_total, ny_total)`; outputs are written over the interior
   !! `(nx_phys, ny_phys)`.  Combined with the `nghost` global
   !! attribute, that lets a reader reconstruct the 9-point cross
   !! stencil `(i-2..i+2, j) + (i, j-2..j+2)` that each interior cell
   !! consumes.
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_io_netcdf, only: nc_create_file, nc_close, nc_def_dim, &
                            nc_def_var_2d, nc_put_att, &
                            nc_put_att_global, nc_put_att_global_int, &
                            nc_enddef, nc_put_var_2d, &
                            rki_def_var_1d, rki_put_var_1d
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: dump_flux_io

contains

   subroutine dump_flux_io(filename, state)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Snapshot the flux-kernel inputs and outputs from `state` to a
      !! NetCDF file at `filename`.
      character(len=*), intent(in) :: filename
      type(state_t), intent(in) :: state
         !! Simulation state with the flux arrays already populated by
         !! `solver_compute_barotropic_flux`.

      real(wp), allocatable :: x_total(:), y_total(:)
      integer :: ncid
      integer :: x_tot_dimid, y_tot_dimid, x_int_dimid, y_int_dimid
      integer :: x_tot_varid, y_tot_varid, x_int_varid, y_int_varid
      integer :: h_varid, hu_varid, hv_varid, b_varid
      integer :: fh_varid, fhu_varid, fhv_varid, mx_varid, my_varid
      integer :: i, j, lo_i, hi_i, lo_j, hi_j
      integer :: nx, ny, ng

      nx = state%grid%nx_total
      ny = state%grid%ny_total
      ng = state%grid%nghost
      lo_i = ng + 1
      hi_i = nx - ng
      lo_j = ng + 1
      hi_j = ny - ng

      ! Pull the slices we are about to write back from device.  No
      ! data region is opened — the state arrays remain device-resident
      ! after this call.  `if_present` makes the routine usable on a
      ! host-only state too (e.g. in unit tests where the caller never
      ! enters a data region).
      !$acc update self(state%barotropic%h, state%barotropic%hu, &
      !$acc&            state%barotropic%hv, state%barotropic%b) if_present
      !$acc update self(state%barotropic%flux_h(lo_i:hi_i, lo_j:hi_j), &
      !$acc&            state%barotropic%flux_hu(lo_i:hi_i, lo_j:hi_j), &
      !$acc&            state%barotropic%flux_hv(lo_i:hi_i, lo_j:hi_j), &
      !$acc&            state%barotropic%mass_flux_x(lo_i:hi_i, lo_j:hi_j), &
      !$acc&            state%barotropic%mass_flux_y(lo_i:hi_i, lo_j:hi_j)) if_present

      ! Cell-centre coordinates with the interior origin at x=y=0; ghost
      ! cells land at negative indices.
      allocate (x_total(nx), y_total(ny))
      do i = 1, nx
         x_total(i) = (real(i - ng, wp) - 0.5_wp)*state%grid%dx
      end do
      do j = 1, ny
         y_total(j) = (real(j - ng, wp) - 0.5_wp)*state%grid%dy
      end do

      call logger%info("Writing flux I/O dump: "//trim(filename))
      call nc_create_file(filename, ncid)

      call nc_def_dim(ncid, "x_total", nx, x_tot_dimid)
      call nc_def_dim(ncid, "y_total", ny, y_tot_dimid)
      call nc_def_dim(ncid, "x", state%grid%nx_phys, x_int_dimid)
      call nc_def_dim(ncid, "y", state%grid%ny_phys, y_int_dimid)

      call rki_def_var_1d(ncid, "x_total", x_tot_dimid, x_tot_varid)
      call nc_put_att(ncid, x_tot_varid, "units", "m")
      call nc_put_att(ncid, x_tot_varid, "long_name", &
                      "cell-centre x including ghost layers")

      call rki_def_var_1d(ncid, "y_total", y_tot_dimid, y_tot_varid)
      call nc_put_att(ncid, y_tot_varid, "units", "m")
      call nc_put_att(ncid, y_tot_varid, "long_name", &
                      "cell-centre y including ghost layers")

      call rki_def_var_1d(ncid, "x", x_int_dimid, x_int_varid)
      call nc_put_att(ncid, x_int_varid, "units", "m")
      call nc_put_att(ncid, x_int_varid, "long_name", &
                      "cell-centre x of interior cells")
      call nc_put_att(ncid, x_int_varid, "axis", "X")

      call rki_def_var_1d(ncid, "y", y_int_dimid, y_int_varid)
      call nc_put_att(ncid, y_int_varid, "units", "m")
      call nc_put_att(ncid, y_int_varid, "long_name", &
                      "cell-centre y of interior cells")
      call nc_put_att(ncid, y_int_varid, "axis", "Y")

      ! Inputs: full ghost-padded arrays. Each interior cell (i, j) reads
      ! the 9-point cross stencil (i-2..i+2, j) and (i, j-2..j+2).
      call nc_def_var_2d(ncid, "h", [x_tot_dimid, y_tot_dimid], h_varid)
      call nc_put_att(ncid, h_varid, "units", "m")
      call nc_put_att(ncid, h_varid, "long_name", "input water depth")

      call nc_def_var_2d(ncid, "hu", [x_tot_dimid, y_tot_dimid], hu_varid)
      call nc_put_att(ncid, hu_varid, "units", "m2 s-1")
      call nc_put_att(ncid, hu_varid, "long_name", "input x-momentum (h*u)")

      call nc_def_var_2d(ncid, "hv", [x_tot_dimid, y_tot_dimid], hv_varid)
      call nc_put_att(ncid, hv_varid, "units", "m2 s-1")
      call nc_put_att(ncid, hv_varid, "long_name", "input y-momentum (h*v)")

      call nc_def_var_2d(ncid, "b", [x_tot_dimid, y_tot_dimid], b_varid)
      call nc_put_att(ncid, b_varid, "units", "m")
      call nc_put_att(ncid, b_varid, "long_name", &
                      "input bottom elevation (positive up)")

      ! Outputs: interior only.
      call nc_def_var_2d(ncid, "flux_h", [x_int_dimid, y_int_dimid], fh_varid)
      call nc_put_att(ncid, fh_varid, "units", "m s-1")
      call nc_put_att(ncid, fh_varid, "long_name", &
                      "net flux divergence of h")

      call nc_def_var_2d(ncid, "flux_hu", [x_int_dimid, y_int_dimid], fhu_varid)
      call nc_put_att(ncid, fhu_varid, "units", "m2 s-2")
      call nc_put_att(ncid, fhu_varid, "long_name", &
                      "net flux divergence of hu (with well-balanced correction)")

      call nc_def_var_2d(ncid, "flux_hv", [x_int_dimid, y_int_dimid], fhv_varid)
      call nc_put_att(ncid, fhv_varid, "units", "m2 s-2")
      call nc_put_att(ncid, fhv_varid, "long_name", &
                      "net flux divergence of hv (with well-balanced correction)")

      call nc_def_var_2d(ncid, "mass_flux_x", [x_int_dimid, y_int_dimid], mx_varid)
      call nc_put_att(ncid, mx_varid, "units", "m2 s-1")
      call nc_put_att(ncid, mx_varid, "long_name", &
                      "HLL h-flux at right face (i+1/2)")

      call nc_def_var_2d(ncid, "mass_flux_y", [x_int_dimid, y_int_dimid], my_varid)
      call nc_put_att(ncid, my_varid, "units", "m2 s-1")
      call nc_put_att(ncid, my_varid, "long_name", &
                      "HLL h-flux at top face (j+1/2)")

      call nc_put_att_global(ncid, "title", &
                             "Per-interior-cell I/O snapshot for the barotropic flux kernel")
      call nc_put_att_global(ncid, "source", "Rakali rki_flux_dump")
      if (state%barotropic%use_hllc) then
         call nc_put_att_global(ncid, "solver", "HLLC")
      else
         call nc_put_att_global(ncid, "solver", "HLL")
      end if
      call nc_put_att_global_int(ncid, "nghost", ng)
      call nc_put_att_global_int(ncid, "nx_phys", state%grid%nx_phys)
      call nc_put_att_global_int(ncid, "ny_phys", state%grid%ny_phys)

      call nc_enddef(ncid)

      call rki_put_var_1d(ncid, x_tot_varid, x_total)
      call rki_put_var_1d(ncid, y_tot_varid, y_total)
      call rki_put_var_1d(ncid, x_int_varid, x_total(lo_i:hi_i))
      call rki_put_var_1d(ncid, y_int_varid, y_total(lo_j:hi_j))

      call nc_put_var_2d(ncid, h_varid, state%barotropic%h)
      call nc_put_var_2d(ncid, hu_varid, state%barotropic%hu)
      call nc_put_var_2d(ncid, hv_varid, state%barotropic%hv)
      call nc_put_var_2d(ncid, b_varid, state%barotropic%b)

      call nc_put_var_2d(ncid, fh_varid, &
                         state%barotropic%flux_h(lo_i:hi_i, lo_j:hi_j))
      call nc_put_var_2d(ncid, fhu_varid, &
                         state%barotropic%flux_hu(lo_i:hi_i, lo_j:hi_j))
      call nc_put_var_2d(ncid, fhv_varid, &
                         state%barotropic%flux_hv(lo_i:hi_i, lo_j:hi_j))
      call nc_put_var_2d(ncid, mx_varid, &
                         state%barotropic%mass_flux_x(lo_i:hi_i, lo_j:hi_j))
      call nc_put_var_2d(ncid, my_varid, &
                         state%barotropic%mass_flux_y(lo_i:hi_i, lo_j:hi_j))

      call nc_close(ncid)

   end subroutine dump_flux_io

end module rki_flux_dump
