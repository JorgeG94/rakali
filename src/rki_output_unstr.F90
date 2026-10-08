!! UGRID NetCDF output for the unstructured solver
module rki_output_unstr
   !! Writes solution snapshots to CF UGRID 1.0 NetCDF files.
   !!
   !! Topology (nodes, face-node connectivity, centroids, areas) is
   !! written once at file creation. Time-varying fields (stage, xmom,
   !! ymom) are appended as 2D slices `(face, time)` on each snapshot.
   !! Bathymetry is written once as a static 1D field `(face)`.
   !!
   !! One file with `ncells_owned` faces. The `cell_global_id` variable is the
   !! cell index `1..ncells_owned`.
   !!
   !! Two init paths:
   !!   * `output_unstr_init` — takes `mesh_unstr_t` + `unstr_state_t`.
   !!   * `output_unstr_init_from_arrays` — takes raw arrays.
   !!
   !! Two write paths:
   !!   * `output_unstr_write_snapshot` — reads from `unstr_state_t`.
   !!   * `output_unstr_write_snapshot_from_arrays` — takes raw arrays.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   use rki_io_netcdf, only: nc_check, nc_create_file, nc_close, &
                            nc_def_dim, rki_def_var_1d, nc_def_var_2d, &
                            nc_put_att, nc_put_att_real, nc_put_att_int, &
                            nc_put_att_global, nc_put_att_global_int, &
                            nc_enddef, rki_put_var_1d, rki_put_var_2d_slice, &
                            NC_WP
   use rki_constants, only: DRY_TOLERANCE
   use netcdf, only: nf90_unlimited, nf90_int, nf90_int64, nf90_def_var, &
                     nf90_put_var, nf90_sync, nf90_global, &
                     nf90_def_var_deflate
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: output_unstr_state_t
   public :: output_unstr_init
   public :: output_unstr_init_from_arrays
   public :: output_unstr_write_snapshot
   public :: output_unstr_write_snapshot_from_arrays
   public :: output_unstr_close

   type :: output_unstr_state_t
      !! Tracks the state of an open UGRID output file
      integer :: ncid = -1
         !! NetCDF file ID (-1 = not open)
      integer :: time_dimid = -1
         !! Time dimension ID
      integer :: time_varid = -1
         !! Time variable ID
      integer :: stage_varid = -1
         !! Free-surface elevation variable ID
      integer :: xmom_varid = -1
         !! x-momentum variable ID
      integer :: ymom_varid = -1
         !! y-momentum variable ID
      integer :: time_index = 0
         !! Current time record index
      integer :: ncells = 0
         !! Number of owned cells written per snapshot

      ! ---- multilayer (only populated when state%use_multilayer) ----
      logical :: has_ml = .false.
         !! True when 3D layer / tracer variables were defined.
      integer :: nz_ml = 0
         !! Number of sigma layers (cached for snapshot writes).
      integer :: ml_h_layer_varid = -1
         !! Layer thickness variable ID
      integer :: ml_u_varid = -1
         !! Layer x-velocity variable ID (derived: hu_layer / h_layer)
      integer :: ml_v_varid = -1
         !! Layer y-velocity variable ID
      integer :: ml_rho_varid = -1
         !! Layer density variable ID
      integer :: ml_w_varid = -1
         !! Vertical velocity at interfaces variable ID
      integer, allocatable :: ml_tracer_varids(:)
         !! Per-tracer NetCDF variable IDs, one per
         !! `state%multilayer%tracers(:)` entry.  `-1` for tracers that
         !! aren't being written.  Adding a sediment / passive scalar
         !! requires no I/O changes — this loop picks it up.
   end type output_unstr_state_t

contains

   ! ================================================================
   ! Private core: creates the UGRID file from raw arrays
   ! ================================================================

   subroutine create_ugrid_file(ostate, filename, nverts, ncells, &
                                vert_x, vert_y, tri_v, &
                                centroid_x, centroid_y, area, bed, &
                                cell_global_id, dlev, state)
      !! Shared implementation: create file, define UGRID topology, write
      !! static data. Both public init routines delegate here.  When
      !! `state` is present and `state%use_multilayer` is true, also
      !! defines the 3D layer / tracer variables before exiting define
      !! mode.
      type(output_unstr_state_t), intent(out) :: ostate
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nverts, ncells
      real(wp), intent(in) :: vert_x(nverts), vert_y(nverts)
      integer, intent(in) :: tri_v(3, ncells)
      real(wp), intent(in) :: centroid_x(ncells), centroid_y(ncells)
      real(wp), intent(in) :: area(ncells), bed(ncells)
      integer(8), intent(in), optional :: cell_global_id(ncells)
      integer, intent(in) :: dlev
      type(unstr_state_t), intent(in), optional :: state

      integer :: ncid
      integer :: node_dimid, face_dimid, three_dimid, time_dimid
      integer :: topo_varid
      integer :: node_x_varid, node_y_varid
      integer :: conn_varid
      integer :: face_x_varid, face_y_varid, face_area_varid
      integer :: bed_varid
      integer :: stage_varid, xmom_varid, ymom_varid
      integer :: time_varid
      integer :: cgid_varid

      ! ---- create file ----
      call nc_create_file(trim(filename), ncid)

      ! ---- global attributes ----
      call nc_put_att_global(ncid, "Conventions", "CF-1.11 UGRID-1.0")
      call nc_put_att_global(ncid, "title", "Rakali unstructured solver output")
      call nc_put_att_global(ncid, "source", "Rakali")
      call nc_put_att_global(ncid, "institution", "NCI")
      call nc_put_att_global_int(ncid, "ncells_owned", ncells)
      call nc_put_att_global_int(ncid, "nverts", nverts)

      ! ---- dimensions ----
      call nc_def_dim(ncid, "node", nverts, node_dimid)
      call nc_def_dim(ncid, "face", ncells, face_dimid)
      call nc_def_dim(ncid, "three", 3, three_dimid)
      call nc_def_dim(ncid, "time", nf90_unlimited, time_dimid)

      ! ---- mesh topology dummy variable ----
      call nc_check(nf90_def_var(ncid, "mesh2d", nf90_int, topo_varid), &
                    "defining mesh2d")
      call nc_put_att(ncid, topo_varid, "cf_role", "mesh_topology")
      call nc_put_att_int(ncid, topo_varid, "topology_dimension", 2)
      call nc_put_att(ncid, topo_varid, "node_coordinates", &
                      "mesh_node_x mesh_node_y")
      call nc_put_att(ncid, topo_varid, "face_node_connectivity", &
                      "mesh_face_nodes")
      call nc_put_att(ncid, topo_varid, "face_coordinates", &
                      "mesh_face_x mesh_face_y")

      ! ---- node coordinates ----
      call rki_def_var_1d(ncid, "mesh_node_x", node_dimid, node_x_varid, dlev)
      call nc_put_att(ncid, node_x_varid, "standard_name", &
                      "projection_x_coordinate")
      call nc_put_att(ncid, node_x_varid, "long_name", &
                      "x-coordinate of mesh nodes")
      call nc_put_att(ncid, node_x_varid, "units", "m")

      call rki_def_var_1d(ncid, "mesh_node_y", node_dimid, node_y_varid, dlev)
      call nc_put_att(ncid, node_y_varid, "standard_name", &
                      "projection_y_coordinate")
      call nc_put_att(ncid, node_y_varid, "long_name", &
                      "y-coordinate of mesh nodes")
      call nc_put_att(ncid, node_y_varid, "units", "m")

      ! ---- face-node connectivity (integer, shape three x face) ----
      call nc_check(nf90_def_var(ncid, "mesh_face_nodes", nf90_int, &
                                 [three_dimid, face_dimid], conn_varid), &
                    "defining mesh_face_nodes")
      call nc_put_att(ncid, conn_varid, "cf_role", "face_node_connectivity")
      call nc_put_att_int(ncid, conn_varid, "start_index", 1)
      if (dlev > 0) then
         call nc_check(nf90_def_var_deflate(ncid, conn_varid, &
                                            shuffle=1, deflate=1, &
                                            deflate_level=dlev), &
                       "enabling compression for mesh_face_nodes")
      end if

      ! ---- face centroids ----
      call rki_def_var_1d(ncid, "mesh_face_x", face_dimid, face_x_varid, dlev)
      call nc_put_att(ncid, face_x_varid, "long_name", &
                      "x-coordinate of face centroids")
      call nc_put_att(ncid, face_x_varid, "units", "m")

      call rki_def_var_1d(ncid, "mesh_face_y", face_dimid, face_y_varid, dlev)
      call nc_put_att(ncid, face_y_varid, "long_name", &
                      "y-coordinate of face centroids")
      call nc_put_att(ncid, face_y_varid, "units", "m")

      ! ---- face area ----
      call rki_def_var_1d(ncid, "mesh_face_area", face_dimid, &
                          face_area_varid, dlev)
      call nc_put_att(ncid, face_area_varid, "long_name", "Triangle area")
      call nc_put_att(ncid, face_area_varid, "units", "m2")

      ! ---- bathymetry (static, face) ----
      call rki_def_var_1d(ncid, "bed", face_dimid, bed_varid, dlev)
      call nc_put_att(ncid, bed_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, bed_varid, "location", "face")
      call nc_put_att(ncid, bed_varid, "long_name", "Bed elevation")
      call nc_put_att(ncid, bed_varid, "units", "m")

      ! ---- cell global id ----
      if (present(cell_global_id)) then
         call nc_check(nf90_def_var(ncid, "cell_global_id", nf90_int64, &
                                    [face_dimid], cgid_varid), &
                       "defining cell_global_id")
         call nc_put_att(ncid, cgid_varid, "long_name", &
                         "Cell index")
      end if

      ! ---- time ----
      call rki_def_var_1d(ncid, "time", time_dimid, time_varid)
      call nc_put_att(ncid, time_varid, "units", &
                      "seconds since simulation start")
      call nc_put_att(ncid, time_varid, "calendar", "none")
      call nc_put_att(ncid, time_varid, "axis", "T")

      ! ---- time-varying solution fields (face, time) ----
      call nc_def_var_2d(ncid, "stage", [face_dimid, time_dimid], &
                         stage_varid, dlev)
      call nc_put_att(ncid, stage_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, stage_varid, "location", "face")
      call nc_put_att(ncid, stage_varid, "long_name", &
                      "Free-surface elevation")
      call nc_put_att(ncid, stage_varid, "units", "m")

      call nc_def_var_2d(ncid, "xmom", [face_dimid, time_dimid], &
                         xmom_varid, dlev)
      call nc_put_att(ncid, xmom_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, xmom_varid, "location", "face")
      call nc_put_att(ncid, xmom_varid, "long_name", "x-momentum (h*u)")
      call nc_put_att(ncid, xmom_varid, "units", "m2/s")

      call nc_def_var_2d(ncid, "ymom", [face_dimid, time_dimid], &
                         ymom_varid, dlev)
      call nc_put_att(ncid, ymom_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, ymom_varid, "location", "face")
      call nc_put_att(ncid, ymom_varid, "long_name", "y-momentum (h*v)")
      call nc_put_att(ncid, ymom_varid, "units", "m2/s")

      ! ---- multilayer outputs (only when state%use_multilayer) ----
      ! Defined here while still in define mode so we get one nc_enddef
      ! for the whole file.  Variables go in (face, layer, time) order
      ! to match the unstructured convention `(nz_ml, ncells)` flipped
      ! so layer is the slow axis on disk.
      if (present(state)) then
         if (state%use_multilayer) then
            call define_multilayer_outputs(ostate, ncid, face_dimid, &
                                           time_dimid, dlev, state)
         end if
      end if

      ! ---- end define mode, write static data ----
      call nc_enddef(ncid)

      call rki_put_var_1d(ncid, node_x_varid, vert_x)
      call rki_put_var_1d(ncid, node_y_varid, vert_y)
      call nc_check(nf90_put_var(ncid, conn_varid, tri_v), &
                    "writing mesh_face_nodes")
      call rki_put_var_1d(ncid, face_x_varid, centroid_x)
      call rki_put_var_1d(ncid, face_y_varid, centroid_y)
      call rki_put_var_1d(ncid, face_area_varid, area)
      call rki_put_var_1d(ncid, bed_varid, bed)

      if (present(cell_global_id)) then
         call nc_check(nf90_put_var(ncid, cgid_varid, cell_global_id), &
                       "writing cell_global_id")
      end if

      ! Stash handles for snapshot writes
      ostate%ncid = ncid
      ostate%time_dimid = time_dimid
      ostate%time_varid = time_varid
      ostate%stage_varid = stage_varid
      ostate%xmom_varid = xmom_varid
      ostate%ymom_varid = ymom_varid
      ostate%time_index = 0
      ostate%ncells = ncells

      call logger%info("UGRID output initialised: "//trim(filename)// &
                       " ("//to_string(ncells)//" faces, "// &
                       to_string(nverts)//" nodes)")

   end subroutine create_ugrid_file

   subroutine define_multilayer_outputs(ostate, ncid, face_dimid, time_dimid, &
                                        dlev, state)
      !! Define the 3D layer / tracer variables on a UGRID file that's
      !! still in define mode.  Per-tracer variables come from the
      !! registry — one NetCDF variable per `state%multilayer%tracers(:)`
      !! entry, metadata pulled directly from each `unstr_tracer_t`.
      !! Adding a sediment / passive scalar produces a correctly-
      !! labelled output variable with no I/O changes.
      !!
      !! Dim order `[nz_dimid, face_dimid, time_dimid]` matches the
      !! in-memory `(nz_ml, ncells_total)` layout so each snapshot is a
      !! single contiguous `nf90_put_var` (no transpose buffer needed).
      type(output_unstr_state_t), intent(inout) :: ostate
      integer, intent(in)    :: ncid, face_dimid, time_dimid, dlev
      type(unstr_state_t), intent(in)    :: state

      integer :: nz_dimid, nzi_dimid
      integer :: dimids_3d(3), dimids_3d_w(3)
      integer :: it, varid

      ostate%has_ml = .true.
      ostate%nz_ml = state%multilayer%nz_ml

      call nc_def_dim(ncid, "nz_ml", state%multilayer%nz_ml, nz_dimid)
      call nc_def_dim(ncid, "nz_ml_interface", state%multilayer%nz_ml + 1, nzi_dimid)

      dimids_3d = [nz_dimid, face_dimid, time_dimid]
      dimids_3d_w = [nzi_dimid, face_dimid, time_dimid]

      ! Layer thickness (m)
      call nc_check(nf90_def_var(ncid, "h_layer", NC_WP, dimids_3d, &
                                 ostate%ml_h_layer_varid), &
                    "defining h_layer")
      call nc_put_att(ncid, ostate%ml_h_layer_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, ostate%ml_h_layer_varid, "location", "face")
      call nc_put_att(ncid, ostate%ml_h_layer_varid, "long_name", "Layer thickness")
      call nc_put_att(ncid, ostate%ml_h_layer_varid, "units", "m")
      if (dlev > 0) then
         call nc_check(nf90_def_var_deflate(ncid, ostate%ml_h_layer_varid, &
                                            shuffle=1, deflate=1, deflate_level=dlev), &
                       "compress h_layer")
      end if

      ! Layer x-velocity (derived: hu_layer / h_layer)
      call nc_check(nf90_def_var(ncid, "u_ml", NC_WP, dimids_3d, ostate%ml_u_varid), &
                    "defining u_ml")
      call nc_put_att(ncid, ostate%ml_u_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, ostate%ml_u_varid, "location", "face")
      call nc_put_att(ncid, ostate%ml_u_varid, "long_name", "Layer x-velocity")
      call nc_put_att(ncid, ostate%ml_u_varid, "units", "m s-1")
      if (dlev > 0) then
         call nc_check(nf90_def_var_deflate(ncid, ostate%ml_u_varid, &
                                            shuffle=1, deflate=1, deflate_level=dlev), &
                       "compress u_ml")
      end if

      ! Layer y-velocity
      call nc_check(nf90_def_var(ncid, "v_ml", NC_WP, dimids_3d, ostate%ml_v_varid), &
                    "defining v_ml")
      call nc_put_att(ncid, ostate%ml_v_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, ostate%ml_v_varid, "location", "face")
      call nc_put_att(ncid, ostate%ml_v_varid, "long_name", "Layer y-velocity")
      call nc_put_att(ncid, ostate%ml_v_varid, "units", "m s-1")
      if (dlev > 0) then
         call nc_check(nf90_def_var_deflate(ncid, ostate%ml_v_varid, &
                                            shuffle=1, deflate=1, deflate_level=dlev), &
                       "compress v_ml")
      end if

      ! Layer density
      call nc_check(nf90_def_var(ncid, "rho_layer", NC_WP, dimids_3d, &
                                 ostate%ml_rho_varid), &
                    "defining rho_layer")
      call nc_put_att(ncid, ostate%ml_rho_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, ostate%ml_rho_varid, "location", "face")
      call nc_put_att(ncid, ostate%ml_rho_varid, "standard_name", "sea_water_density")
      call nc_put_att(ncid, ostate%ml_rho_varid, "long_name", "Layer density")
      call nc_put_att(ncid, ostate%ml_rho_varid, "units", "kg m-3")
      if (dlev > 0) then
         call nc_check(nf90_def_var_deflate(ncid, ostate%ml_rho_varid, &
                                            shuffle=1, deflate=1, deflate_level=dlev), &
                       "compress rho_layer")
      end if

      ! Vertical velocity at layer interfaces (nz_ml + 1)
      call nc_check(nf90_def_var(ncid, "w_interface", NC_WP, dimids_3d_w, &
                                 ostate%ml_w_varid), &
                    "defining w_interface")
      call nc_put_att(ncid, ostate%ml_w_varid, "mesh", "mesh2d")
      call nc_put_att(ncid, ostate%ml_w_varid, "location", "face")
      call nc_put_att(ncid, ostate%ml_w_varid, "long_name", &
                      "Vertical velocity at layer interfaces")
      call nc_put_att(ncid, ostate%ml_w_varid, "units", "m s-1")
      if (dlev > 0) then
         call nc_check(nf90_def_var_deflate(ncid, ostate%ml_w_varid, &
                                            shuffle=1, deflate=1, deflate_level=dlev), &
                       "compress w_interface")
      end if

      ! ---- Tracer registry: one variable per registered tracer ----
      if (allocated(state%multilayer%tracers)) then
         allocate (ostate%ml_tracer_varids(size(state%multilayer%tracers)))
         ostate%ml_tracer_varids = -1
         do it = 1, size(state%multilayer%tracers)
            associate (tr => state%multilayer%tracers(it))
               if (len_trim(tr%name) == 0) cycle
               call nc_check(nf90_def_var(ncid, trim(tr%name), NC_WP, &
                                          dimids_3d, varid), &
                             "defining tracer "//trim(tr%name))
               ostate%ml_tracer_varids(it) = varid
               call nc_put_att(ncid, varid, "mesh", "mesh2d")
               call nc_put_att(ncid, varid, "location", "face")
               if (len_trim(tr%units) > 0) then
                  call nc_put_att(ncid, varid, "units", trim(tr%units))
               end if
               if (len_trim(tr%long_name) > 0) then
                  call nc_put_att(ncid, varid, "long_name", trim(tr%long_name))
               end if
               if (len_trim(tr%standard_name) > 0) then
                  call nc_put_att(ncid, varid, "standard_name", trim(tr%standard_name))
               end if
               if (dlev > 0) then
                  call nc_check(nf90_def_var_deflate(ncid, varid, &
                                                     shuffle=1, deflate=1, &
                                                     deflate_level=dlev), &
                                "compress "//trim(tr%name))
               end if
            end associate
         end do
      end if
   end subroutine define_multilayer_outputs

   ! ================================================================
   ! Public init from mesh/state objects (compute-rank path)
   ! ================================================================

   subroutine output_unstr_init(ostate, filename, state, mesh, deflate_level)
      !! Create a CF UGRID 1.0 NetCDF file from mesh + state objects.
      type(output_unstr_state_t), intent(out) :: ostate
      character(len=*), intent(in) :: filename
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in), optional :: deflate_level

      integer :: dlev, n

      dlev = 0
      if (present(deflate_level)) dlev = deflate_level
      n = mesh%ncells_owned

      if (allocated(mesh%cell_global_id)) then
         call create_ugrid_file(ostate, filename, mesh%nverts, n, &
                                mesh%vert_x(1:mesh%nverts), &
                                mesh%vert_y(1:mesh%nverts), &
                                mesh%tri_v(:, 1:n), &
                                mesh%centroid_x(1:n), mesh%centroid_y(1:n), &
                                mesh%area(1:n), state%barotropic%bed_cv(1:n), &
                                cell_global_id=mesh%cell_global_id(1:n), &
                                dlev=dlev, state=state)
      else
         call create_ugrid_file(ostate, filename, mesh%nverts, n, &
                                mesh%vert_x(1:mesh%nverts), &
                                mesh%vert_y(1:mesh%nverts), &
                                mesh%tri_v(:, 1:n), &
                                mesh%centroid_x(1:n), mesh%centroid_y(1:n), &
                                mesh%area(1:n), state%barotropic%bed_cv(1:n), &
                                dlev=dlev, state=state)
      end if

   end subroutine output_unstr_init

   ! ================================================================
   ! Public init from raw arrays (I/O server path)
   ! ================================================================

   subroutine output_unstr_init_from_arrays(ostate, filename, nverts, ncells, &
                                            vert_x, vert_y, tri_v, &
                                            centroid_x, centroid_y, area, bed, &
                                            cell_global_id, deflate_level)
      !! Create a CF UGRID 1.0 NetCDF file from raw arrays.
      !! For callers that hold raw topology arrays instead of a mesh.
      type(output_unstr_state_t), intent(out) :: ostate
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nverts, ncells
      real(wp), intent(in) :: vert_x(nverts), vert_y(nverts)
      integer, intent(in) :: tri_v(3, ncells)
      real(wp), intent(in) :: centroid_x(ncells), centroid_y(ncells)
      real(wp), intent(in) :: area(ncells), bed(ncells)
      integer(8), intent(in), optional :: cell_global_id(ncells)
      integer, intent(in), optional :: deflate_level

      integer :: dlev

      dlev = 0
      if (present(deflate_level)) dlev = deflate_level

      if (present(cell_global_id)) then
         call create_ugrid_file(ostate, filename, nverts, ncells, &
                                vert_x, vert_y, tri_v, &
                                centroid_x, centroid_y, area, bed, &
                                cell_global_id=cell_global_id, dlev=dlev)
      else
         call create_ugrid_file(ostate, filename, nverts, ncells, &
                                vert_x, vert_y, tri_v, &
                                centroid_x, centroid_y, area, bed, &
                                dlev=dlev)
      end if

   end subroutine output_unstr_init_from_arrays

   ! ================================================================
   ! Public snapshot writers
   ! ================================================================

   subroutine output_unstr_write_snapshot(ostate, state, time)
      !! Append a snapshot from state object. Caller must ensure host sync.
      type(output_unstr_state_t), intent(inout) :: ostate
      type(unstr_state_t), intent(in) :: state
      real(wp), intent(in) :: time

      integer :: n
      n = ostate%ncells
      call output_unstr_write_snapshot_from_arrays( &
         ostate, state%barotropic%stage_cv(1:n), state%barotropic%xmom_cv(1:n), &
         state%barotropic%ymom_cv(1:n), time)

      ! Multilayer + tracer fields (only when the file was opened with
      ! state%use_multilayer; otherwise has_ml is false).
      if (ostate%has_ml) then
         call write_multilayer_snapshot(ostate, state)
      end if

   end subroutine output_unstr_write_snapshot

   subroutine write_multilayer_snapshot(ostate, state)
      !! Write the 3D layer fields + per-tracer fields for the current
      !! time index.  Velocities are derived from `hu_layer / h_layer`
      !! with a dry-cell guard; tracers from `tracers(it)%hTr / h_layer`.
      !!
      !! Memory layout `(nz_ml, ncells_total)` matches the file layout
      !! `[nz_dimid, face_dimid, time_dimid]` so each field is one
      !! contiguous `nf90_put_var` — no transpose buffer.
      type(output_unstr_state_t), intent(inout) :: ostate
      type(unstr_state_t), intent(in)    :: state

      integer :: ti, n, nz, kk, cc, it
      real(wp), allocatable :: buf(:, :)
      real(wp) :: h_k

      ti = ostate%time_index
      n = ostate%ncells
      nz = ostate%nz_ml

      ! Layer thickness — write directly from memory.
      call nc_check(nf90_put_var(ostate%ncid, ostate%ml_h_layer_varid, &
                                 state%multilayer%h_layer(1:nz, 1:n), &
                                 start=[1, 1, ti], &
                                 count=[nz, n, 1]), &
                    "writing h_layer slice")

      ! Velocity components — derive u = hu / h with dry-cell guard.
      allocate (buf(nz, n))
      do cc = 1, n
         do kk = 1, nz
            h_k = state%multilayer%h_layer(kk, cc)
            if (h_k > DRY_TOLERANCE) then
               buf(kk, cc) = state%multilayer%hu_layer(kk, cc)/h_k
            else
               buf(kk, cc) = 0.0_wp
            end if
         end do
      end do
      call nc_check(nf90_put_var(ostate%ncid, ostate%ml_u_varid, buf, &
                                 start=[1, 1, ti], count=[nz, n, 1]), &
                    "writing u_ml slice")

      do cc = 1, n
         do kk = 1, nz
            h_k = state%multilayer%h_layer(kk, cc)
            if (h_k > DRY_TOLERANCE) then
               buf(kk, cc) = state%multilayer%hv_layer(kk, cc)/h_k
            else
               buf(kk, cc) = 0.0_wp
            end if
         end do
      end do
      call nc_check(nf90_put_var(ostate%ncid, ostate%ml_v_varid, buf, &
                                 start=[1, 1, ti], count=[nz, n, 1]), &
                    "writing v_ml slice")

      ! Density (already a primitive variable — direct write).
      call nc_check(nf90_put_var(ostate%ncid, ostate%ml_rho_varid, &
                                 state%multilayer%rho_layer(1:nz, 1:n), &
                                 start=[1, 1, ti], count=[nz, n, 1]), &
                    "writing rho_layer slice")

      ! Vertical velocity at interfaces (nz + 1 levels).
      call nc_check(nf90_put_var(ostate%ncid, ostate%ml_w_varid, &
                                 state%multilayer%w_interface(1:nz + 1, 1:n), &
                                 start=[1, 1, ti], count=[nz + 1, n, 1]), &
                    "writing w_interface slice")

      ! ---- Per-tracer variables ----
      ! Divide the conserved hTr by h_layer to recover the primitive
      ! variable in the same units the registry advertises.
      if (allocated(state%multilayer%tracers) .and. allocated(ostate%ml_tracer_varids)) then
         do it = 1, size(state%multilayer%tracers)
            if (ostate%ml_tracer_varids(it) < 0) cycle
            do cc = 1, n
               do kk = 1, nz
                  h_k = state%multilayer%h_layer(kk, cc)
                  if (h_k > DRY_TOLERANCE) then
                     buf(kk, cc) = state%multilayer%tracers(it)%hTr(kk, cc)/h_k
                  else
                     buf(kk, cc) = 0.0_wp
                  end if
               end do
            end do
            call nc_check(nf90_put_var(ostate%ncid, ostate%ml_tracer_varids(it), &
                                       buf, start=[1, 1, ti], count=[nz, n, 1]), &
                          "writing tracer slice")
         end do
      end if

      deallocate (buf)
   end subroutine write_multilayer_snapshot

   subroutine output_unstr_write_snapshot_from_arrays(ostate, stage, xmom, &
                                                      ymom, time)
      !! Append a snapshot from raw arrays. Used by both compute-rank
      !! and I/O server paths.
      type(output_unstr_state_t), intent(inout) :: ostate
      real(wp), intent(in) :: stage(:), xmom(:), ymom(:)
      real(wp), intent(in) :: time

      integer :: ti

      if (ostate%ncid < 0) return

      ostate%time_index = ostate%time_index + 1
      ti = ostate%time_index

      call nc_check(nf90_put_var(ostate%ncid, ostate%time_varid, time, &
                                 start=[ti]), "writing time")

      call rki_put_var_2d_slice(ostate%ncid, ostate%stage_varid, stage, ti)
      call rki_put_var_2d_slice(ostate%ncid, ostate%xmom_varid, xmom, ti)
      call rki_put_var_2d_slice(ostate%ncid, ostate%ymom_varid, ymom, ti)

      call nc_check(nf90_sync(ostate%ncid), "syncing output file")

      call logger%debug("UGRID snapshot "//to_string(ti)// &
                        " at t = "//to_string(time)//" s")

   end subroutine output_unstr_write_snapshot_from_arrays

   subroutine output_unstr_close(ostate)
      !! Close the output file.
      type(output_unstr_state_t), intent(inout) :: ostate

      if (ostate%ncid < 0) return

      call nc_close(ostate%ncid)
      if (allocated(ostate%ml_tracer_varids)) deallocate (ostate%ml_tracer_varids)
      call logger%info("UGRID output closed ("// &
                       to_string(ostate%time_index)//" snapshots written)")
      ostate%ncid = -1

   end subroutine output_unstr_close

end module rki_output_unstr
