!! Integration test: VCOORD_ZSTAR_FULL on a triangular mesh with a
!! hard 10 m → 200 m bathymetry step.  F3 mirror of the structured
!! `test_zstar_full_shelf_break`.
!!
!! What we validate:
!!   1. Solver runs through several SSP-RK2 steps without NaN/blowup
!!      under VCOORD_ZSTAR_FULL on triangles.
!!   2. On BOTH sides of the shelf, the surface ROMS layer (k=NZ)
!!      lands close to `zstar_h_surf_target` — full z*'s headline
!!      property (fixed physical surface-layer thickness regardless
!!      of total column depth).  z*-lite cannot deliver this.
!!   3. Per-column mass conservation: sum(h_layer(:, cell)) ≈ h(cell)
!!      at every interior cell.
module test_zstar_full_shelf_break_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, &
                             mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, &
                              state_unstr_destroy, &
                              state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_solver_unstr, only: solver_unstr_step, &
                               solver_unstr_init_multilayer, &
                               solver_unstr_update_host
   use rki_vcoord, only: STRETCH_UNIFORM
   implicit none
   private

   public :: collect_zstar_full_shelf_break_unstr_tests

   integer, parameter :: NZ = 10
   real(wp), parameter :: H_SHALLOW = 10.0_wp
   real(wp), parameter :: H_DEEP = 200.0_wp
   real(wp), parameter :: H_SURF_TARGET = 2.0_wp
   real(wp), parameter :: H_MIN = 1.0e-4_wp
   integer, parameter :: N_SURF = 3
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: T_INIT = 15.0_wp

   ! Seiche params for the passive-tracer registry test (flat bed + SSH bump).
   real(wp), parameter :: H_FLAT = 20.0_wp
   real(wp), parameter :: ETA_AMP = 0.5_wp

contains

   subroutine collect_zstar_full_shelf_break_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("unstr_solver_runs_no_nan", &
                               test_unstr_solver_runs_no_nan), &
                  new_unittest("unstr_deep_surface_at_target", &
                               test_unstr_deep_surface_at_target), &
                  new_unittest("unstr_shallow_surface_at_target", &
                               test_unstr_shallow_surface_at_target), &
                  new_unittest("unstr_per_column_mass_conservation", &
                               test_unstr_per_column_mass), &
                  new_unittest("unstr_passive_tracer_tracks_salinity", &
                               test_unstr_passive_tracer_tracks_salinity) &
                  ]
   end subroutine collect_zstar_full_shelf_break_unstr_tests

   subroutine build_shelf_mesh(mesh)
      !! 4x2 strip of unit squares (x in [0, 4], y in [0, 1]), each
      !! split into 4 triangles fanning from a centre vertex.  Cells
      !! whose centroid x < 2 sit on the shallow side; cells with
      !! centroid x >= 2 on the deep side.  16 cells total.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(14), vy(14)
      integer :: tv(3, 16)
      integer :: i, sq

      ! Corners of 4 squares + 4 centre vertices = 8 corners (5x2 grid) + 4 centres = 14? Let me just do 4 vertical columns of squares.
      ! Layout:
      !   y=1:  v6  v7  v8  v9  v10
      !   y=0:  v1  v2  v3  v4  v5
      !  centres: v11 v12 v13 v14 (at y=0.5)
      vx = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, &
            0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, &
            0.5_wp, 1.5_wp, 2.5_wp, 3.5_wp]
      vy = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, &
            0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp]

      ! Each square split into 4 fans from its centre vertex.
      ! Square 1: corners v1, v2, v7, v6, centre v11.
      ! Square 2: corners v2, v3, v8, v7, centre v12.
      ! Square 3: corners v3, v4, v9, v8, centre v13.
      ! Square 4: corners v4, v5, v10, v9, centre v14.
      do sq = 1, 4
         i = (sq - 1)*4
         ! sq's centre vertex
         block
            integer :: c, a, b, c2, d
            c = 10 + sq
            a = sq         ! bot-left
            b = sq + 1     ! bot-right
            c2 = sq + 6    ! top-right (5 + sq + 1)
            d = sq + 5     ! top-left  (5 + sq)
            tv(:, i + 1) = [a, b, c]      ! bottom
            tv(:, i + 2) = [b, c2, c]     ! right
            tv(:, i + 3) = [c2, d, c]     ! top
            tv(:, i + 4) = [d, a, c]      ! left
         end block
      end do
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_shelf_mesh

   subroutine setup_shelf_state(state, mesh)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      integer :: cell

      call build_shelf_mesh(mesh)

      ! Set ZSTAR_FULL parameters BEFORE state_unstr_init so the
      ! multilayer init allocates the per-cell z_ref_col array.
      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      call state_unstr_init(state, mesh, nz_ml=NZ)

      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      ! Bathymetry: hard step at x=2.  bed_cv = bed elevation
      ! (negative for submerged).
      do cell = 1, state%ncells_total
         if (mesh%centroid_x(cell) < 2.0_wp) then
            state%barotropic%bed_cv(cell) = -H_SHALLOW
         else
            state%barotropic%bed_cv(cell) = -H_DEEP
         end if
         state%barotropic%stage_cv(cell) = 0.0_wp
         state%barotropic%xmom_cv(cell) = 0.0_wp
         state%barotropic%ymom_cv(cell) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_init_multilayer(state, mesh)
   end subroutine setup_shelf_state

   logical function any_nan_in_h_layer(state)
      type(unstr_state_t), intent(in) :: state
      integer :: cell, k
      any_nan_in_h_layer = .false.
      do cell = 1, state%ncells_total
         do k = 1, NZ
            if (state%multilayer%h_layer(k, cell) &
                /= state%multilayer%h_layer(k, cell)) then
               any_nan_in_h_layer = .true.
               return
            end if
         end do
      end do
   end function any_nan_in_h_layer

   subroutine run_to(state, mesh, t_final, dt)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: t_final, dt
      real(wp) :: t

      t = 0.0_wp
      do while (t < t_final)
         call solver_unstr_step(state, mesh, dt)
         t = t + dt
      end do
   end subroutine run_to

   subroutine test_unstr_solver_runs_no_nan(error)
      type(error_type), allocatable, intent(out) :: error
      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp), parameter :: T_END = 0.1_wp
      real(wp), parameter :: DT = 0.005_wp

      call setup_shelf_state(state, mesh)
      call run_to(state, mesh, T_END, DT)

      call check(error,.not. any_nan_in_h_layer(state), &
                 "ZSTAR_FULL unstr shelf-break: no NaN in h_layer")

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_solver_runs_no_nan

   subroutine test_unstr_deep_surface_at_target(error)
      type(error_type), allocatable, intent(out) :: error
      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: cell, cell_deep
      real(wp) :: dz_surf
      checks: block

         call setup_shelf_state(state, mesh)

         ! Find a deep cell
         cell_deep = 0
         do cell = 1, state%ncells_owned
            if (mesh%centroid_x(cell) > 2.5_wp) then
               cell_deep = cell
               exit
            end if
         end do
         call check(error, cell_deep > 0, "deep cell located")
         if (cell_deep == 0) exit checks

         dz_surf = state%multilayer%h_layer(NZ, cell_deep)
         call check(error, dz_surf > 0.5_wp*H_SURF_TARGET .and. &
                    dz_surf < 2.0_wp*H_SURF_TARGET, &
                    "deep-side surface ROMS layer thickness near H_SURF_TARGET")

      end block checks
      continue
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_deep_surface_at_target

   subroutine test_unstr_shallow_surface_at_target(error)
      type(error_type), allocatable, intent(out) :: error
      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: cell, cell_shallow
      real(wp) :: dz_surf
      checks: block

         call setup_shelf_state(state, mesh)

         cell_shallow = 0
         do cell = 1, state%ncells_owned
            if (mesh%centroid_x(cell) < 1.5_wp) then
               cell_shallow = cell
               exit
            end if
         end do
         call check(error, cell_shallow > 0, "shallow cell located")
         if (cell_shallow == 0) exit checks

         dz_surf = state%multilayer%h_layer(NZ, cell_shallow)
         call check(error, dz_surf > 0.5_wp*H_SURF_TARGET .and. &
                    dz_surf < 2.0_wp*H_SURF_TARGET, &
                    "shallow-side surface ROMS layer thickness near H_SURF_TARGET")

      end block checks
      continue
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_shallow_surface_at_target

   subroutine test_unstr_per_column_mass(error)
      type(error_type), allocatable, intent(out) :: error
      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp), parameter :: T_END = 0.05_wp
      real(wp), parameter :: DT = 0.005_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      integer :: cell, k
      real(wp) :: sum_h, h_cell, diff

      call setup_shelf_state(state, mesh)
      call run_to(state, mesh, T_END, DT)

      do cell = 1, state%ncells_owned
         sum_h = 0.0_wp
         do k = 1, NZ
            sum_h = sum_h + state%multilayer%h_layer(k, cell)
         end do
         h_cell = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
         diff = abs(sum_h - h_cell)
         call check(error, diff < TOL*max(h_cell, 1.0_wp), &
                    "unstr per-column sum(h_layer) = h after solver step")
         if (allocated(error)) exit
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_per_column_mass

   subroutine setup_seiche_state(state, mesh)
      !! Flat-bed closed basin under VCOORD_ZSTAR_FULL with a Gaussian SSH
      !! bump that drives a barotropic seiche, so advection + the per-cell
      !! ml_distribute_with_remap_full_unstr + ml_remap_tracer_unstr
      !! pipeline all do real work.  beta_S = alpha_T = 0 (salinity
      !! effectively passive), no mixing, widened tracer bounds so the
      !! clamp never fires.  Registers a passive dye seeded identically to
      !! salinity so it must track S bit-for-bit through generic transport.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      integer :: cell
      real(wp) :: xc, yc

      call build_shelf_mesh(mesh)

      ! ZSTAR_FULL params BEFORE state_unstr_init (allocates z_ref_col).
      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      call state_unstr_init(state, mesh, nz_ml=NZ)

      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp
      state%multilayer%S_min = -1000.0_wp
      state%multilayer%S_max = 1000.0_wp
      state%multilayer%T_min = -1000.0_wp
      state%multilayer%T_max = 1000.0_wp

      ! Append a passive dye with salinity's IC + transport config BEFORE
      ! enter_data so tracers_enter_data maps its hTr; init_multilayer
      ! then seeds it from tr_init = S_INIT.
      call state%multilayer%register_passive_tracer(mesh, "dye", S_INIT, &
                                                    -1000.0_wp, 1000.0_wp, &
                                                    0.0_wp, 0.0_wp)

      ! Flat bed + Gaussian SSH bump at the basin centre (x=2, y=0.5).
      do cell = 1, state%ncells_total
         state%barotropic%bed_cv(cell) = -H_FLAT
         xc = mesh%centroid_x(cell)
         yc = mesh%centroid_y(cell)
         state%barotropic%stage_cv(cell) = ETA_AMP &
                                           *exp(-((xc - 2.0_wp)**2 + (yc - 0.5_wp)**2)/0.5_wp)
         state%barotropic%xmom_cv(cell) = 0.0_wp
         state%barotropic%ymom_cv(cell) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_init_multilayer(state, mesh)
   end subroutine setup_seiche_state

   subroutine test_unstr_passive_tracer_tracks_salinity(error)
      !! Tracer-registry generalization (unstructured): register a 3rd
      !! PASSIVE tracer with salinity's IC + transport config and run a
      !! seiche.  With beta_S = 0 salinity is effectively passive, so the
      !! dye must track it bit-for-bit through the generic transport
      !! pipeline (advection + ml_remap_tracer_unstr + RK2 + clamp).  Any
      !! kernel that silently skips the 3rd tracer (e.g. an un-generalized
      !! remap) makes the two diverge — strictly stronger than a bare
      !! sum-conservation check, which a skipped tracer would falsely pass.
      !! F3 mirror of the structured
      !! test_zstar_full_tracer_conservation::passive_tracer_tracks_salinity.
      type(error_type), allocatable, intent(out) :: error
      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp), parameter :: T_END = 0.5_wp
      real(wp), parameter :: DT = 0.005_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      integer :: idx_dye, idx_s
      real(wp) :: max_diff

      call setup_seiche_state(state, mesh)
      call run_to(state, mesh, T_END, DT)
      call solver_unstr_update_host(state)

      idx_dye = size(state%multilayer%tracers)
      idx_s = state%multilayer%idx_salinity
      max_diff = maxval(abs(state%multilayer%tracers(idx_dye)%hTr &
                            - state%multilayer%tracers(idx_s)%hTr))

      call check(error, max_diff < TOL, &
                 "unstr passive tracer must track salinity through the "// &
                 "generic transport pipeline (advection + remap + RK2 + clamp)")

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_passive_tracer_tracks_salinity

end module test_zstar_full_shelf_break_unstr
