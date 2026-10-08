!! Behavioural tests for the unstructured extrapolation kernel
module test_extrapolate_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, &
                              state_unstr_destroy, &
                              state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_mesh_unstr, only: mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_kernel_extrapolate_unstr, only: extrapolate_second_order_edge_unstr
   implicit none
   private

   public :: collect_extrapolate_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_extrapolate_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("isolated_triangle_constant_edges", test_isolated_triangle_constant), &
                  new_unittest("two_tri_constant_field", test_two_tri_constant), &
                  new_unittest("two_tri_lake_at_rest_flat", test_two_tri_lake_flat), &
                  new_unittest("two_tri_lake_at_rest_sloped_bed", test_two_tri_lake_sloped), &
                  new_unittest("four_tri_lake_at_rest_sloped_bed", test_four_tri_lake_sloped), &
                  new_unittest("well_balanced_identity_holds", test_wb_identity), &
                  new_unittest("dry_cell_centroid_zeroed", test_dry_centroid_zeroed), &
                  new_unittest("velocity_form_centroid_recovered", test_velocity_form_recovery), &
                  new_unittest("constant_momentum_passes_through", test_constant_momentum), &
                  new_unittest("flag_off_does_not_touch_centroid", test_flag_off_centroid_untouched) &
                  ]
   end subroutine collect_extrapolate_tests

   ! ----------------------------------------------------------------------
   ! Mesh builders (mirror the test_mesh_unstr fixtures)
   ! ----------------------------------------------------------------------

   subroutine build_isolated_triangle(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(3), vy(3)
      integer :: tv(3, 1)

      vx = [0.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_isolated_triangle

   subroutine build_two_triangle_square(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      tv(:, 2) = [1, 3, 4]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_two_triangle_square

   subroutine build_four_triangle_square(mesh)
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
   end subroutine build_four_triangle_square

   subroutine set_lake_at_rest(state, eta_const, bed_fn_x, bed_fn_y, bed_const)
      !! Initialise stage to a constant `eta_const` everywhere and bed
      !! to a linear function `bed_const + bed_fn_x*x + bed_fn_y*y`
      !! evaluated at the cell centroids. Momenta zero.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: eta_const, bed_fn_x, bed_fn_y, bed_const
      integer :: k

      ! Caller must have called state_unstr_init already; we lean on the
      ! fact that the state arrays are sized to ncells_total. We don't
      ! have a mesh handle here, so we just write into all elements.
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      do k = 1, state%ncells_total
         state%barotropic%stage_cv(k) = eta_const
         state%barotropic%bed_cv(k) = bed_const + bed_fn_x*0.0_wp + bed_fn_y*0.0_wp
      end do
      ! Bed needs the actual centroid coordinates -- the caller can pass
      ! a mesh and refill bed_cv after this if it wants a non-flat bed.
   end subroutine set_lake_at_rest

   subroutine run_kernel(state, mesh)
      !! Helper that wraps the extrapolation kernel with the device
      !! enter/exit data dance. On CPU builds the enter/exit calls are
      !! no-ops, so this is identical to calling the kernel directly.
      !! On GPU builds it ensures the state and mesh derived-type
      !! components are mapped to device memory before the kernel runs
      !! and the host-visible results are copied back afterwards.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      call extrapolate_second_order_edge_unstr(state, mesh)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_kernel

   subroutine set_bed_linear(state, mesh, b0, bx, by)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: b0, bx, by
      integer :: k

      do k = 1, mesh%ncells_total
         state%barotropic%bed_cv(k) = b0 + bx*mesh%centroid_x(k) + by*mesh%centroid_y(k)
      end do
   end subroutine set_bed_linear

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_isolated_triangle_constant(error)
      !! n_boundaries == 3 -> every edge value should equal the
      !! centroid (the constant-edge-values branch).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: i

      call build_isolated_triangle(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv(1) = 0.0_wp
      state%barotropic%stage_cv(1) = 2.0_wp     ! depth = 2
      state%barotropic%xmom_cv(1) = 0.5_wp
      state%barotropic%ymom_cv(1) = -0.3_wp

      call run_kernel(state, mesh)

      do i = 1, 3
         call check(error, abs(state%barotropic%stage_ev(i, 1) - 2.0_wp) < TOL, &
                    "stage edge equals centroid")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%height_ev(i, 1) - 2.0_wp) < TOL, &
                    "height edge equals centroid")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_ev(i, 1) - 0.5_wp) < TOL, &
                    "xmom edge equals centroid (constant case)")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_ev(i, 1) + 0.3_wp) < TOL, &
                    "ymom edge equals centroid (constant case)")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%bed_ev(i, 1) - 0.0_wp) < TOL, &
                    "bed edge from stage_ev - height_ev = 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_isolated_triangle_constant

   subroutine test_two_tri_constant(error)
      !! Constant stage and constant momentum on the 2-triangle square
      !! (n_boundaries == 2 case). Every edge should match the centroid.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, i

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.5_wp
      state%barotropic%xmom_cv = 0.7_wp
      state%barotropic%ymom_cv = -0.2_wp

      call run_kernel(state, mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            call check(error, abs(state%barotropic%stage_ev(i, k) - 1.5_wp) < TOL, &
                       "constant stage edge")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%xmom_ev(i, k) - 0.7_wp) < TOL, &
                       "constant xmom edge")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%ymom_ev(i, k) + 0.2_wp) < TOL, &
                       "constant ymom edge")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%height_ev(i, k) - 1.5_wp) < TOL, &
                       "constant height edge")
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_two_tri_constant

   subroutine test_two_tri_lake_flat(error)
      !! Lake at rest with a flat bed -- the simplest possible case.
      !! Every reconstructed stage_ev must equal stage_cv to machine
      !! precision and every momentum edge value must be exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, i

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 3.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call run_kernel(state, mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            call check(error, abs(state%barotropic%stage_ev(i, k) - 3.0_wp) < TOL, &
                       "lake at rest: stage_ev should equal centroid")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%xmom_ev(i, k)) < TOL, &
                       "lake at rest: xmom_ev should be zero")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%ymom_ev(i, k)) < TOL, &
                       "lake at rest: ymom_ev should be zero")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_two_tri_lake_flat

   subroutine test_two_tri_lake_sloped(error)
      !! Lake at rest with a sloped (linear) bed: stage = constant
      !! everywhere, bed = b0 + bx*x + by*y. The water depth varies
      !! linearly with position but the free-surface elevation (stage)
      !! is exactly constant. Stage_ev must therefore equal stage_cv,
      !! and momenta must remain zero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, i

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 5.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 1.0_wp, 0.5_wp, -0.3_wp)

      call run_kernel(state, mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            call check(error, abs(state%barotropic%stage_ev(i, k) - 5.0_wp) < TOL, &
                       "lake at rest (sloped bed): stage_ev = stage_cv")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%xmom_ev(i, k)) < TOL, &
                       "lake at rest (sloped bed): xmom_ev = 0")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%ymom_ev(i, k)) < TOL, &
                       "lake at rest (sloped bed): ymom_ev = 0")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_two_tri_lake_sloped

   subroutine test_four_tri_lake_sloped(error)
      !! Same well-balancing test on the 4-triangle square. This puts
      !! every cell through the n_boundaries == 1 (full 3-neighbour)
      !! branch instead of the 2-neighbour branch.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, i

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 4.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 0.5_wp, 0.7_wp, 0.2_wp)

      call run_kernel(state, mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            call check(error, abs(state%barotropic%stage_ev(i, k) - 4.0_wp) < TOL, &
                       "4-tri lake at rest: stage_ev = stage_cv")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%xmom_ev(i, k)) < TOL, &
                       "4-tri lake at rest: xmom_ev = 0")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%ymom_ev(i, k)) < TOL, &
                       "4-tri lake at rest: ymom_ev = 0")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_four_tri_lake_sloped

   subroutine test_wb_identity(error)
      !! Across all cells and edges, the well-balanced identity
      !! `stage_ev = bed_ev + height_ev` must hold to machine precision.
      !! This is what the kernel enforces by computing
      !! `bed_ev = stage_ev - height_ev` rather than reconstructing bed
      !! independently.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, i
      real(wp) :: residual

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Some non-trivial state -- doesn't have to be physical, just
      ! something with non-zero gradients in stage and bed.
      state%barotropic%stage_cv = [4.0_wp, 4.5_wp, 4.2_wp, 4.1_wp]
      call set_bed_linear(state, mesh, 0.0_wp, 0.5_wp, -0.5_wp)
      state%barotropic%xmom_cv = [0.1_wp, -0.2_wp, 0.3_wp, 0.0_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.1_wp, -0.1_wp, 0.2_wp]

      call run_kernel(state, mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            residual = state%barotropic%stage_ev(i, k) - state%barotropic%bed_ev(i, k) - state%barotropic%height_ev(i, k)
            call check(error, abs(residual) < TOL, &
                       "stage_ev should equal bed_ev + height_ev exactly")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wb_identity

   subroutine test_dry_centroid_zeroed(error)
      !! A cell whose depth is below `minimum_allowed_height` should
      !! have its centroid momenta zeroed by the kernel (after the
      !! velocity-form restore pass). The work arrays were stashed as
      !! zero on the way in, so the restore pass writes zero back.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Cell 1 dry, cell 2 wet.
      state%barotropic%bed_cv = [10.0_wp, 0.0_wp]
      state%barotropic%stage_cv = [10.0_wp, 5.0_wp]    ! depth: 0, 5
      state%barotropic%xmom_cv = [3.0_wp, 0.4_wp]
      state%barotropic%ymom_cv = [-2.0_wp, 0.1_wp]

      call run_kernel(state, mesh)

      call check(error, abs(state%barotropic%xmom_cv(1)) < TOL, &
                 "dry cell xmom_cv should be zeroed")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1)) < TOL, &
                 "dry cell ymom_cv should be zeroed")
      if (allocated(error)) return
      ! Wet cell: original momentum should be restored from the work array.
      call check(error, abs(state%barotropic%xmom_cv(2) - 0.4_wp) < TOL, &
                 "wet cell xmom_cv should be restored")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.1_wp) < TOL, &
                 "wet cell ymom_cv should be restored")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dry_centroid_zeroed

   subroutine test_velocity_form_recovery(error)
      !! When velocity-form extrapolation is on (default), all wet
      !! centroid momenta must be byte-identical before and after the
      !! kernel runs. The work-array round-trip is the contract.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: xmom_in(4), ymom_in(4)
      integer :: k

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 3.0_wp
      state%barotropic%xmom_cv = [0.5_wp, -0.3_wp, 0.7_wp, 0.1_wp]
      state%barotropic%ymom_cv = [0.2_wp, 0.4_wp, -0.5_wp, -0.1_wp]

      xmom_in = state%barotropic%xmom_cv(1:4)
      ymom_in = state%barotropic%ymom_cv(1:4)

      call run_kernel(state, mesh)

      do k = 1, 4
         call check(error, abs(state%barotropic%xmom_cv(k) - xmom_in(k)) < TOL, &
                    "centroid xmom should round-trip")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_cv(k) - ymom_in(k)) < TOL, &
                    "centroid ymom should round-trip")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_velocity_form_recovery

   subroutine test_constant_momentum(error)
      !! With constant stage AND constant momentum, the limiter should
      !! see zero gradients everywhere and pass through unchanged. Edge
      !! momenta must equal centroid momenta on every cell.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k, i

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 1.5_wp
      state%barotropic%ymom_cv = -0.8_wp

      call run_kernel(state, mesh)

      do k = 1, mesh%ncells_owned
         do i = 1, 3
            call check(error, abs(state%barotropic%stage_ev(i, k) - 2.0_wp) < TOL, &
                       "constant momentum: stage_ev = 2")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%xmom_ev(i, k) - 1.5_wp) < TOL, &
                       "constant momentum: xmom_ev = 1.5")
            if (allocated(error)) return
            call check(error, abs(state%barotropic%ymom_ev(i, k) + 0.8_wp) < TOL, &
                       "constant momentum: ymom_ev = -0.8")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_constant_momentum

   subroutine test_flag_off_centroid_untouched(error)
      !! With velocity-form extrapolation disabled, the kernel should
      !! still produce sensible edge values and the centroid momenta
      !! must NOT be touched (no division-by-height or restore pass).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: xmom_in(2)

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%extrapolate_velocity_second_order = .false.
      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 4.0_wp
      state%barotropic%xmom_cv = [0.6_wp, -0.4_wp]
      state%barotropic%ymom_cv = 0.0_wp
      xmom_in = state%barotropic%xmom_cv(1:2)

      call run_kernel(state, mesh)

      call check(error, abs(state%barotropic%xmom_cv(1) - xmom_in(1)) < TOL, &
                 "flag off: xmom_cv(1) untouched")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - xmom_in(2)) < TOL, &
                 "flag off: xmom_cv(2) untouched")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_flag_off_centroid_untouched

end module test_extrapolate_unstr
