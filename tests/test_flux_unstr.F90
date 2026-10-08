!! Behavioural tests for the unstructured central-upwind flux kernel
module test_flux_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
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
   use rki_kernel_extrapolate_unstr, only: extrapolate_second_order_edge_unstr
   use rki_kernel_flux_unstr, only: compute_fluxes_central_unstr
   implicit none
   private

   public :: collect_flux_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_flux_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("isolated_triangle_lake_at_rest", test_isolated_lake_at_rest), &
                  new_unittest("two_tri_lake_at_rest_flat", test_two_tri_lake_flat), &
                  new_unittest("two_tri_lake_at_rest_sloped", test_two_tri_lake_sloped), &
                  new_unittest("four_tri_lake_at_rest_sloped", test_four_tri_lake_sloped), &
                  new_unittest("isolated_wall_zero_mass_flux", test_isolated_wall_zero_mass_flux), &
                  new_unittest("dry_state_no_nan", test_dry_state_no_nan), &
                  new_unittest("max_speed_recorded_lake_at_rest", test_max_speed_lake_at_rest), &
                  new_unittest("eu_arrays_overwritten_each_call", test_eu_overwritten) &
                  ]
   end subroutine collect_flux_tests

   ! ----------------------------------------------------------------------
   ! Mesh + state helpers
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

   subroutine set_bed_linear(state, mesh, b0, bx, by)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: b0, bx, by
      integer :: k

      do k = 1, mesh%ncells_total
         state%barotropic%bed_cv(k) = b0 + bx*mesh%centroid_x(k) + by*mesh%centroid_y(k)
      end do
   end subroutine set_bed_linear

   subroutine run_flux_pipeline(state, mesh)
      !! Run the extrapolation kernel and then the flux kernel,
      !! wrapped with the device data dance. The two kernels chained
      !! together is exactly what the SSP-RK2 driver will do later --
      !! this is the first place we exercise that pattern on GPU.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      call extrapolate_second_order_edge_unstr(state, mesh)
      call compute_fluxes_central_unstr(state, mesh)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_flux_pipeline

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_isolated_lake_at_rest(error)
      !! Single triangle, all walls, lake at rest. Wall mirroring with
      !! zero momentum is a no-op, and equal stages on both sides give
      !! zero flux. Every accumulator must be exactly zero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_isolated_triangle(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call run_flux_pipeline(state, mesh)

      call check(error, abs(state%barotropic%stage_eu(1)) < TOL, &
                 "isolated lake-at-rest: stage_eu = 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_eu(1)) < TOL, &
                 "isolated lake-at-rest: xmom_eu = 0")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_eu(1)) < TOL, &
                 "isolated lake-at-rest: ymom_eu = 0")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_isolated_lake_at_rest

   subroutine test_two_tri_lake_flat(error)
      !! Two triangles sharing the diagonal. Flat bed, constant stage,
      !! zero momentum. The interior edge sees identical states on
      !! both sides, the wall edges see mirrored zero-momentum states,
      !! so every flux contribution is zero.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.5_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call run_flux_pipeline(state, mesh)

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_eu(k)) < TOL, &
                    "2-tri flat lake: stage_eu = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_eu(k)) < TOL, &
                    "2-tri flat lake: xmom_eu = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_eu(k)) < TOL, &
                    "2-tri flat lake: ymom_eu = 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_two_tri_lake_flat

   subroutine test_two_tri_lake_sloped(error)
      !! THE WELL-BALANCED TEST. Sloped bed, constant stage. The
      !! depth varies but the free surface is exactly horizontal.
      !! All explicit-update accumulators must be zero to machine
      !! precision -- this exercises the bed-slope source telescoping
      !! against the central-upwind pressure flux.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 5.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 1.0_wp, 0.5_wp, -0.3_wp)

      call run_flux_pipeline(state, mesh)

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_eu(k)) < TOL, &
                    "2-tri sloped lake: stage_eu = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_eu(k)) < TOL, &
                    "2-tri sloped lake: xmom_eu = 0 (well-balanced)")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_eu(k)) < TOL, &
                    "2-tri sloped lake: ymom_eu = 0 (well-balanced)")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_two_tri_lake_sloped

   subroutine test_four_tri_lake_sloped(error)
      !! Same well-balancing test on the 4-triangle square. Each cell
      !! has 1 boundary edge, so the extrapolation goes through the
      !! 3-neighbour branch with hfactor blending. This is a tougher
      !! test for well-balancing than the 2-tri case.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k

      call build_four_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%stage_cv = 4.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      call set_bed_linear(state, mesh, 0.5_wp, 0.7_wp, 0.2_wp)

      call run_flux_pipeline(state, mesh)

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_eu(k)) < TOL, &
                    "4-tri sloped lake: stage_eu = 0")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_eu(k)) < TOL, &
                    "4-tri sloped lake: xmom_eu = 0 (well-balanced)")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_eu(k)) < TOL, &
                    "4-tri sloped lake: ymom_eu = 0 (well-balanced)")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_four_tri_lake_sloped

   subroutine test_isolated_wall_zero_mass_flux(error)
      !! Isolated triangle with non-zero momentum. Wall mirroring
      !! reflects the normal velocity, so the central-upwind mass
      !! flux through every wall edge is exactly zero (s_max + s_min
      !! cancels). The momentum accumulators may be nonzero (the
      !! walls exert a pressure force) but mass conservation must
      !! hold.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_isolated_triangle(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 1.0_wp
      state%barotropic%xmom_cv = 0.3_wp
      state%barotropic%ymom_cv = -0.2_wp

      call run_flux_pipeline(state, mesh)

      call check(error, abs(state%barotropic%stage_eu(1)) < 1.0e-10_wp, &
                 "wall reflection: mass flux through walls = 0")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_isolated_wall_zero_mass_flux

   subroutine test_dry_state_no_nan(error)
      !! Entire mesh dry: depth at the dry threshold, zero momentum.
      !! Every flux call hits the dry/dry branch which returns zeros
      !! cleanly -- there must be no NaN or Inf in the explicit-update
      !! arrays at the end.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k
      logical :: any_nan

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Bed = 0, stage = 0 -> depth = 0 everywhere.
      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call run_flux_pipeline(state, mesh)

      any_nan = .false.
      do k = 1, mesh%ncells_owned
         if (state%barotropic%stage_eu(k) /= state%barotropic%stage_eu(k)) any_nan = .true.
         if (state%barotropic%xmom_eu(k) /= state%barotropic%xmom_eu(k)) any_nan = .true.
         if (state%barotropic%ymom_eu(k) /= state%barotropic%ymom_eu(k)) any_nan = .true.
      end do
      call check(error,.not. any_nan, "dry state must not produce NaN")
      if (allocated(error)) return

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_eu(k)) < TOL, &
                    "dry state: stage_eu = 0")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dry_state_no_nan

   subroutine test_max_speed_lake_at_rest(error)
      !! For lake at rest with depth h, the central-upwind wave speeds
      !! are s_max = +c, s_min = -c with c = sqrt(g*h). The kernel
      !! tracks max(s_max, -s_min) = c per edge and per cell. For a
      !! flat constant lake every cell sees the same c on every edge,
      !! so state%barotropic%max_speed(k) must equal sqrt(g*h) on every owned cell.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: expected_c
      integer :: k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 3.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call run_flux_pipeline(state, mesh)

      expected_c = sqrt(GRAVITY*3.0_wp)
      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%max_speed(k) - expected_c) < 1.0e-10_wp, &
                    "lake at rest: max_speed = sqrt(g*h)")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_max_speed_lake_at_rest

   subroutine test_eu_overwritten(error)
      !! The kernel must overwrite (stage_eu, xmom_eu, ymom_eu) on each
      !! call, not accumulate into whatever was there before. Verifying
      !! this directly by pre-filling the host arrays does not work
      !! under our GPU data-residency model: the *_eu arrays are mapped
      !! `alloc:` on enter_data and `delete:` on exit_data, so a host
      !! pre-fill never reaches the device and a device write never
      !! reaches the host. The contract is real, but we have to check
      !! it without crossing the host/device boundary.
      !!
      !! Trick: run the kernel twice on the same state and compare the
      !! two results. If the kernel were accumulating, the second run
      !! would double every accumulator. If it's overwriting (the
      !! correct contract), the second run produces an identical
      !! result. We use a non-trivial state so the accumulators are
      !! non-zero -- a doubling would be unmistakable.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp), allocatable :: stage_eu_first(:), xmom_eu_first(:), ymom_eu_first(:)
      integer :: k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Non-trivial state: sloped bed, finite momentum -- the flux
      ! accumulators will be nonzero, so any "accumulate instead of
      ! overwrite" bug would show up as exactly 2x on the second run.
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.5_wp
      state%barotropic%ymom_cv = -0.3_wp
      call set_bed_linear(state, mesh, 0.2_wp, 0.4_wp, 0.1_wp)

      ! First run -- snapshot the result.
      call run_flux_pipeline(state, mesh)

      allocate (stage_eu_first(mesh%ncells_owned))
      allocate (xmom_eu_first(mesh%ncells_owned))
      allocate (ymom_eu_first(mesh%ncells_owned))
      stage_eu_first = state%barotropic%stage_eu(1:mesh%ncells_owned)
      xmom_eu_first = state%barotropic%xmom_eu(1:mesh%ncells_owned)
      ymom_eu_first = state%barotropic%ymom_eu(1:mesh%ncells_owned)

      ! Reset the centroid state to identical inputs (the previous run
      ! left height_cv populated and may have shifted things), so the
      ! second flux pass starts from the same conserved variables.
      state%barotropic%stage_cv = 2.0_wp
      state%barotropic%xmom_cv = 0.5_wp
      state%barotropic%ymom_cv = -0.3_wp
      call set_bed_linear(state, mesh, 0.2_wp, 0.4_wp, 0.1_wp)

      ! Second run -- must produce the same accumulators, not double them.
      call run_flux_pipeline(state, mesh)

      do k = 1, mesh%ncells_owned
         call check(error, abs(state%barotropic%stage_eu(k) - stage_eu_first(k)) < 1.0e-12_wp, &
                    "second run: stage_eu must equal first (not double)")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%xmom_eu(k) - xmom_eu_first(k)) < 1.0e-12_wp, &
                    "second run: xmom_eu must equal first (not double)")
         if (allocated(error)) return
         call check(error, abs(state%barotropic%ymom_eu(k) - ymom_eu_first(k)) < 1.0e-12_wp, &
                    "second run: ymom_eu must equal first (not double)")
         if (allocated(error)) return
      end do

      deallocate (stage_eu_first, xmom_eu_first, ymom_eu_first)
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_eu_overwritten

end module test_flux_unstr
