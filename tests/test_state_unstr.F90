!! Unit tests for the unstructured state type
module test_state_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, &
                              state_unstr_destroy, &
                              LOW_FROUDE_OFF
   use rki_boundary_unstr, only: BC_WALL, MAX_BC_GROUPS
   implicit none
   private

   public :: collect_state_unstr_tests

contains

   subroutine collect_state_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_allocates_arrays", test_init_allocates), &
                  new_unittest("init_sets_sizes", test_init_sets_sizes), &
                  new_unittest("init_zeroes_arrays", test_init_zeroes), &
                  new_unittest("edge_arrays_have_shape_3xN", test_edge_arrays_shape), &
                  new_unittest("destroy_deallocates", test_destroy_deallocates), &
                  new_unittest("double_init_safe", test_double_init), &
                  new_unittest("default_physics_knobs", test_default_knobs), &
                  new_unittest("default_no_multilayer", test_default_no_multilayer), &
                  new_unittest("multilayer_init_allocates", test_multilayer_init_allocates), &
                  new_unittest("multilayer_destroy_resets", test_multilayer_destroy_resets), &
                  new_unittest("multilayer_reinit_as_2d", test_multilayer_reinit_as_2d) &
                  ]
   end subroutine collect_state_unstr_tests

   ! ----------------------------------------------------------------------
   ! Helper -- the same 4-triangle square mesh used in test_mesh_unstr
   ! ----------------------------------------------------------------------

   subroutine build_test_mesh(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(5), vy(5)
      integer :: tv(3, 4)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.5_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.5_wp]
      tv(:, 1) = [1, 2, 5]
      tv(:, 2) = [2, 3, 5]
      tv(:, 3) = [3, 4, 5]
      tv(:, 4) = [4, 1, 5]

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_test_mesh

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_init_allocates(error)
      !! Every allocatable must be allocated after init.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)

      call check(error, allocated(state%barotropic%stage_cv), "stage_cv allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%xmom_cv), "xmom_cv allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%ymom_cv), "ymom_cv allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%bed_cv), "bed_cv allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%height_cv), "height_cv allocated")
      if (allocated(error)) return

      call check(error, allocated(state%barotropic%stage_ev), "stage_ev allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%xmom_ev), "xmom_ev allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%ymom_ev), "ymom_ev allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%bed_ev), "bed_ev allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%height_ev), "height_ev allocated")
      if (allocated(error)) return

      call check(error, allocated(state%barotropic%stage_eu), "stage_eu allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%xmom_eu), "xmom_eu allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%ymom_eu), "ymom_eu allocated")
      if (allocated(error)) return

      call check(error, allocated(state%barotropic%xmom_siu), "xmom_siu allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%ymom_siu), "ymom_siu allocated")
      if (allocated(error)) return

      call check(error, allocated(state%barotropic%stage_backup), "stage_backup allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%xmom_backup), "xmom_backup allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%ymom_backup), "ymom_backup allocated")
      if (allocated(error)) return

      call check(error, allocated(state%barotropic%x_centroid_work), "x_centroid_work allocated")
      if (allocated(error)) return
      call check(error, allocated(state%barotropic%y_centroid_work), "y_centroid_work allocated")
      if (allocated(error)) return

      call check(error, allocated(state%barotropic%max_speed), "max_speed allocated")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_init_allocates

   subroutine test_init_sets_sizes(error)
      !! Cell counts on the state must mirror the mesh.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)

      call check(error, state%ncells_owned == mesh%ncells_owned, &
                 "ncells_owned should mirror mesh")
      if (allocated(error)) return
      call check(error, state%ncells_total == mesh%ncells_total, &
                 "ncells_total should mirror mesh")
      if (allocated(error)) return

      call check(error, size(state%barotropic%stage_cv) == mesh%ncells_total, &
                 "stage_cv length should be ncells_total")
      if (allocated(error)) return
      call check(error, size(state%barotropic%bed_cv) == mesh%ncells_total, &
                 "bed_cv length should be ncells_total")
      if (allocated(error)) return
      call check(error, size(state%barotropic%max_speed) == mesh%ncells_total, &
                 "max_speed length should be ncells_total")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_init_sets_sizes

   subroutine test_init_zeroes(error)
      !! Every array must be exactly zero after init.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)

      call check(error, maxval(abs(state%barotropic%stage_cv)) < 1.0e-15_wp, "stage_cv zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%xmom_cv)) < 1.0e-15_wp, "xmom_cv zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%ymom_cv)) < 1.0e-15_wp, "ymom_cv zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%bed_cv)) < 1.0e-15_wp, "bed_cv zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%height_cv)) < 1.0e-15_wp, "height_cv zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%stage_ev)) < 1.0e-15_wp, "stage_ev zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%stage_eu)) < 1.0e-15_wp, "stage_eu zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%xmom_eu)) < 1.0e-15_wp, "xmom_eu zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%xmom_siu)) < 1.0e-15_wp, "xmom_siu zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%stage_backup)) < 1.0e-15_wp, "backup zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%max_speed)) < 1.0e-15_wp, "max_speed zero")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_init_zeroes

   subroutine test_edge_arrays_shape(error)
      !! Per-edge-of-cell arrays must have shape (3, ncells_total).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)

      call check(error, size(state%barotropic%stage_ev, 1) == 3, "stage_ev dim 1 = 3")
      if (allocated(error)) return
      call check(error, size(state%barotropic%stage_ev, 2) == mesh%ncells_total, &
                 "stage_ev dim 2 = ncells_total")
      if (allocated(error)) return
      call check(error, size(state%barotropic%xmom_ev, 1) == 3, "xmom_ev dim 1 = 3")
      if (allocated(error)) return
      call check(error, size(state%barotropic%height_ev, 1) == 3, "height_ev dim 1 = 3")
      if (allocated(error)) return
      call check(error, size(state%barotropic%bed_ev, 2) == mesh%ncells_total, &
                 "bed_ev dim 2 = ncells_total")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_edge_arrays_shape

   subroutine test_destroy_deallocates(error)
      !! Every allocatable must be released after destroy and counters reset.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)
      call state_unstr_destroy(state)

      call check(error,.not. allocated(state%barotropic%stage_cv), "stage_cv deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%barotropic%bed_cv), "bed_cv deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%barotropic%stage_ev), "stage_ev deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%barotropic%stage_eu), "stage_eu deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%barotropic%xmom_siu), "xmom_siu deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%barotropic%stage_backup), "backup deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%barotropic%max_speed), "max_speed deallocated")
      if (allocated(error)) return
      call check(error, state%ncells_owned == 0, "ncells_owned reset")
      if (allocated(error)) return
      call check(error, state%ncells_total == 0, "ncells_total reset")

      call mesh_unstr_destroy(mesh)
   end subroutine test_destroy_deallocates

   subroutine test_double_init(error)
      !! Calling init twice on the same state must not leak: it should
      !! destroy the old allocation, then resize to the new mesh.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh1, mesh2
      type(unstr_state_t) :: state

      ! First mesh: the 4-triangle square (4 cells).
      call build_test_mesh(mesh1)
      call state_unstr_init(state, mesh1)

      call check(error, state%ncells_owned == 4, "first init: 4 cells")
      if (allocated(error)) return

      ! Second mesh: a single right triangle (1 cell).
      block
         real(wp) :: vx(3), vy(3)
         integer :: tv(3, 1)
         vx = [0.0_wp, 1.0_wp, 0.0_wp]
         vy = [0.0_wp, 0.0_wp, 1.0_wp]
         tv(:, 1) = [1, 2, 3]
         call mesh_unstr_init_from_triangles(mesh2, vx, vy, tv)
      end block

      call state_unstr_init(state, mesh2)

      call check(error, state%ncells_owned == 1, "second init: 1 cell")
      if (allocated(error)) return
      call check(error, size(state%barotropic%stage_cv) == 1, "stage_cv resized to 1")
      if (allocated(error)) return
      call check(error, size(state%barotropic%stage_ev, 2) == 1, "stage_ev resized to (3,1)")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh1)
      call mesh_unstr_destroy(mesh2)
   end subroutine test_double_init

   subroutine test_default_knobs(error)
      !! Default physics knobs should match the ANUGA defaults.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)

      call check(error, &
                 abs(state%minimum_allowed_height - 1.0e-3_wp) < 1.0e-15_wp, &
                 "default minimum_allowed_height should be 1e-3")
      if (allocated(error)) return
      call check(error, abs(state%cfl - 0.5_wp) < 1.0e-15_wp, &
                 "default CFL should be 0.5 (ANUGA convention)")
      if (allocated(error)) return
      call check(error, state%step_count == 0, "step_count starts at 0")
      if (allocated(error)) return
      call check(error, abs(state%time) < 1.0e-15_wp, "time starts at 0")
      if (allocated(error)) return
      call check(error, abs(state%mass_error) < 1.0e-15_wp, "mass_error starts at 0")
      if (allocated(error)) return

      ! ANUGA-default slope-limiter betas
      call check(error, abs(state%beta_w - 1.0_wp) < 1.0e-15_wp, &
                 "default beta_w should be 1.0")
      if (allocated(error)) return
      call check(error, abs(state%beta_w_dry - 0.2_wp) < 1.0e-15_wp, &
                 "default beta_w_dry should be 0.2")
      if (allocated(error)) return
      call check(error, abs(state%beta_uh - 0.8_wp) < 1.0e-15_wp, &
                 "default beta_uh should be 0.8")
      if (allocated(error)) return
      call check(error, abs(state%beta_uh_dry - 0.2_wp) < 1.0e-15_wp, &
                 "default beta_uh_dry should be 0.2")
      if (allocated(error)) return
      call check(error, abs(state%beta_vh - 0.8_wp) < 1.0e-15_wp, &
                 "default beta_vh should be 0.8")
      if (allocated(error)) return
      call check(error, abs(state%beta_vh_dry - 0.2_wp) < 1.0e-15_wp, &
                 "default beta_vh_dry should be 0.2")
      if (allocated(error)) return
      call check(error, state%extrapolate_velocity_second_order, &
                 "default extrapolate_velocity_second_order should be .true.")
      if (allocated(error)) return
      call check(error, state%low_froude_mode == LOW_FROUDE_OFF, &
                 "default low_froude_mode should be LOW_FROUDE_OFF")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%manning_n) < 1.0e-15_wp, &
                 "default manning_n should be 0")
      if (allocated(error)) return
      ! Every boundary group defaults to BC_WALL so a fresh state with
      ! no explicit BC setup behaves like an all-walls domain.
      call check(error, all(state%bc_type(1:MAX_BC_GROUPS) == BC_WALL), &
                 "default bc_type(:) should be all BC_WALL")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_default_knobs

   ! ----------------------------------------------------------------------
   ! Multilayer scaffolding (Phase 0 of the 3D unstructured plan).
   ! ----------------------------------------------------------------------

   subroutine test_default_no_multilayer(error)
      !! Calling init without nz_ml leaves the state in pure 2D mode
      !! and does not allocate any layer arrays.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh)

      call check(error,.not. state%use_multilayer, &
                 "use_multilayer should default to false")
      if (allocated(error)) return
      call check(error, state%multilayer%nz_ml == 1, "nz_ml should default to 1")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%h_layer), "h_layer unallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%tracers), "tracer registry unallocated")
      if (allocated(error)) return
      call check(error, state%multilayer%idx_salinity == 0, "idx_salinity unset")
      if (allocated(error)) return
      call check(error, state%multilayer%idx_temperature == 0, "idx_temperature unset")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%dsig_ml), "dsig_ml unallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%w_interface), "w_interface unallocated")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_default_no_multilayer

   subroutine test_multilayer_init_allocates(error)
      !! Passing nz_ml enables multilayer mode: every layer array is
      !! allocated with layer as the fast index and zero-initialised;
      !! dsig_ml sums to 1.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer, parameter :: NZ = 5

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      call check(error, state%use_multilayer, "use_multilayer should be true")
      if (allocated(error)) return
      call check(error, state%multilayer%nz_ml == NZ, "nz_ml mirrors optional arg")
      if (allocated(error)) return

      ! Shape checks -- layer is dim 1 (fast index).
      call check(error, size(state%multilayer%h_layer, 1) == NZ, "h_layer dim 1 = nz_ml")
      if (allocated(error)) return
      call check(error, size(state%multilayer%h_layer, 2) == mesh%ncells_total, &
                 "h_layer dim 2 = ncells_total")
      if (allocated(error)) return
      call check(error, size(state%multilayer%hu_layer, 1) == NZ, "hu_layer dim 1 = nz_ml")
      if (allocated(error)) return
      call check(error, size(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, 1) == NZ, "hS_layer dim 1 = nz_ml")
      if (allocated(error)) return
      call check(error, size(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, 1) == NZ, "hT_layer dim 1 = nz_ml")
      if (allocated(error)) return
      call check(error, size(state%multilayer%rho_layer, 1) == NZ, "rho_layer dim 1 = nz_ml")
      if (allocated(error)) return
      call check(error, size(state%multilayer%w_interface, 1) == NZ + 1, &
                 "w_interface dim 1 = nz_ml + 1")
      if (allocated(error)) return
      call check(error, size(state%multilayer%hu_layer0, 1) == NZ, "hu_layer0 backup shape")
      if (allocated(error)) return

      ! dsig sums to 1 (uniform fractions for now)
      call check(error, abs(sum(state%multilayer%dsig_ml) - 1.0_wp) < 1.0e-14_wp, &
                 "dsig_ml must sum to 1")
      if (allocated(error)) return

      ! Zero-init
      call check(error, maxval(abs(state%multilayer%h_layer)) < 1.0e-15_wp, "h_layer zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr)) < 1.0e-15_wp, "hS_layer zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%multilayer%tracers(state%multilayer%idx_temperature)%hTr)) < 1.0e-15_wp, "hT_layer zero")
      if (allocated(error)) return
      call check(error, maxval(abs(state%multilayer%w_interface)) < 1.0e-15_wp, "w_interface zero")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_multilayer_init_allocates

   subroutine test_multilayer_destroy_resets(error)
      !! After destroy, layer arrays must be released and the
      !! use_multilayer flag reset.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=3)
      call state_unstr_destroy(state)

      call check(error,.not. state%use_multilayer, "use_multilayer reset")
      if (allocated(error)) return
      call check(error, state%multilayer%nz_ml == 1, "nz_ml reset to 1")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%h_layer), "h_layer released")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%tracers), "tracer registry released")
      if (allocated(error)) return
      call check(error, state%multilayer%idx_salinity == 0, "idx_salinity reset")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%dsig_ml), "dsig_ml released")

      call mesh_unstr_destroy(mesh)
   end subroutine test_multilayer_destroy_resets

   subroutine test_multilayer_reinit_as_2d(error)
      !! Init with nz_ml, then init again without: the second call must
      !! drop the layer arrays and restore 2D-only mode.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_test_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=4)
      call check(error, state%use_multilayer, "first init: multilayer on")
      if (allocated(error)) return

      call state_unstr_init(state, mesh)
      call check(error,.not. state%use_multilayer, &
                 "second init without nz_ml: multilayer off")
      if (allocated(error)) return
      call check(error,.not. allocated(state%multilayer%h_layer), &
                 "h_layer dropped on re-init")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_multilayer_reinit_as_2d

end module test_state_unstr
