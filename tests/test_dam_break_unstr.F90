!! Behavioural test: Ritter dam break on a triangulated channel
module test_dam_break_unstr
   !! The first end-to-end PHYSICS test for the unstructured solver.
   !!
   !! Lake-at-rest only proves the solver does *nothing* correctly --
   !! necessary but insufficient. The dam break has a known analytic
   !! solution (Ritter / Toro Riemann solver) and is the canonical
   !! shallow-water shock-capturing test. Getting the rarefaction
   !! and the shock in roughly the right place at roughly the right
   !! amplitude is a real signal that the central-upwind flux,
   !! second-order extrapolation, well-balanced bed-slope source,
   !! SSP-RK2 driver, protect kernel, and BC dispatch all work
   !! together correctly.
   !!
   !! Setup: a thin triangulated channel along x in [0, 10]. The
   !! dam sits at x = 5. Left side has depth 2 m, right side has
   !! depth 1 m, both initially at rest. We advance with the SSP-RK2
   !! driver for ~0.1 s and compare the centroid stage along the
   !! channel against the exact Ritter solution evaluated at the
   !! same time. All boundaries are walls (the default), but the
   !! disturbance does not reach them in 0.1 s so they don't matter.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
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
   use rki_solver_unstr, only: solver_unstr_advance
   use ritter_dambreak, only: exact_riemann_dambreak, &
                              exact_riemann_sample
   implicit none
   private

   public :: collect_dam_break_tests

   ! Channel geometry
   real(wp), parameter :: LX = 10.0_wp
   real(wp), parameter :: LY = 0.5_wp
   integer, parameter :: NX = 100
   integer, parameter :: NY = 2
   real(wp), parameter :: DAM_X = 5.0_wp
   real(wp), parameter :: H_LEFT = 2.0_wp
   real(wp), parameter :: H_RIGHT = 1.0_wp
   real(wp), parameter :: T_END = 0.1_wp

contains

   subroutine collect_dam_break_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ritter_dam_break_l1_error_below_threshold", test_dam_break_l1), &
                  new_unittest("ritter_dam_break_ymomentum_bounded", test_dam_break_ymom_bounded), &
                  new_unittest("ritter_dam_break_no_negative_depth", test_dam_break_positive), &
                  new_unittest("ritter_dam_break_undisturbed_endpoints", test_dam_break_endpoints), &
                  new_unittest("ritter_dam_break_convergence_with_refinement", test_dam_break_convergence) &
                  ]
   end subroutine collect_dam_break_tests

   ! ----------------------------------------------------------------------
   ! Mesh builder: rectangular channel split into triangles
   ! ----------------------------------------------------------------------

   subroutine build_channel_mesh(mesh, nx_in, ny_in)
      !! Triangulate the rectangle [0, LX] x [0, LY] with `nx_in x ny_in`
      !! quad cells, each split into two CCW right triangles. Total
      !! `2 * nx_in * ny_in` triangles. Defaults to the module-level
      !! NX and NY when no resolution is specified, so existing tests
      !! get the same 100 x 2 channel they did before.
      !!
      !! For quad (i, j) with i in 0..nx-1, j in 0..ny-1, the four
      !! corners are
      !!   v00 = (i*dx,     j*dy)
      !!   v10 = ((i+1)*dx, j*dy)
      !!   v11 = ((i+1)*dx, (j+1)*dy)
      !!   v01 = (i*dx,     (j+1)*dy)
      !! and the two triangles are
      !!   T1 = (v00, v10, v11)   ! lower-right, CCW
      !!   T2 = (v00, v11, v01)   ! upper-left,  CCW
      type(mesh_unstr_t), intent(out) :: mesh
      integer, intent(in), optional :: nx_in
      integer, intent(in), optional :: ny_in

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: nx_use, ny_use
      integer :: i, j, nv, ntri
      integer :: v00, v10, v11, v01, t
      real(wp) :: dx, dy

      nx_use = NX
      ny_use = NY
      if (present(nx_in)) nx_use = nx_in
      if (present(ny_in)) ny_use = ny_in

      dx = LX/real(nx_use, wp)
      dy = LY/real(ny_use, wp)

      ! Vertex grid: (nx+1) * (ny+1) nodes, row-major (j outer, i inner).
      nv = (nx_use + 1)*(ny_use + 1)
      allocate (vx(nv), vy(nv))
      do j = 0, ny_use
         do i = 0, nx_use
            nv = j*(nx_use + 1) + i + 1
            vx(nv) = real(i, wp)*dx
            vy(nv) = real(j, wp)*dy
         end do
      end do

      ! Triangles: 2 per quad, ordered T1 then T2 within each quad,
      ! quads in (j, i) order.
      ntri = 2*nx_use*ny_use
      allocate (tv(3, ntri))
      do j = 0, ny_use - 1
         do i = 0, nx_use - 1
            v00 = j*(nx_use + 1) + i + 1
            v10 = j*(nx_use + 1) + (i + 1) + 1
            v11 = (j + 1)*(nx_use + 1) + (i + 1) + 1
            v01 = (j + 1)*(nx_use + 1) + i + 1

            t = 2*(j*nx_use + i) + 1
            tv(:, t) = [v00, v10, v11]
            tv(:, t + 1) = [v00, v11, v01]
         end do
      end do

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)

      deallocate (vx, vy, tv)
   end subroutine build_channel_mesh

   subroutine set_dam_break_ic(state, mesh)
      !! Initial condition: stage = h_left for cells with centroid_x
      !! <= dam_x, stage = h_right otherwise. Bed flat at zero,
      !! momenta zero everywhere.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer :: k

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      do k = 1, mesh%ncells_total
         if (mesh%centroid_x(k) <= DAM_X) then
            state%barotropic%stage_cv(k) = H_LEFT
         else
            state%barotropic%stage_cv(k) = H_RIGHT
         end if
      end do
   end subroutine set_dam_break_ic

   subroutine run_dam_break(state, mesh)
      !! Run the SSP-RK2 driver from t=0 to T_END inside the device
      !! data dance.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_advance(state, mesh, T_END)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_dam_break

   function relative_l1_error_vs_ritter(state, mesh) result(rel_l1)
      !! Cell-area-weighted L1 error between the solver's stage field
      !! and the exact Ritter solution evaluated at each cell centroid
      !! at the current time. Normalised by the initial jump magnitude
      !! so the result is dimensionless and easy to compare across
      !! resolutions.
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp) :: rel_l1

      real(wp) :: h_star, u_star
      real(wp) :: h_exact, u_exact
      real(wp) :: l1_sum, area_sum
      integer :: k

      call exact_riemann_dambreak(H_LEFT, H_RIGHT, h_star, u_star)

      l1_sum = 0.0_wp
      area_sum = 0.0_wp
      do k = 1, mesh%ncells_owned
         call exact_riemann_sample(H_LEFT, H_RIGHT, h_star, u_star, &
                                   DAM_X, state%time, mesh%centroid_x(k), &
                                   h_exact, u_exact)
         l1_sum = l1_sum + abs(state%barotropic%stage_cv(k) - h_exact)*mesh%area(k)
         area_sum = area_sum + mesh%area(k)
      end do

      rel_l1 = (l1_sum/area_sum)/(H_LEFT - H_RIGHT)
   end function relative_l1_error_vs_ritter

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_dam_break_l1(error)
      !! L1 error against the exact Ritter solution. With ~100 cells
      !! along the channel and second-order KNP, the scheme typically
      !! lands around 3% L1. The earlier 5% bound accepted a doubling
      !! of numerical dissipation — 3.5% keeps a ~15% headroom over
      !! the measured noise floor.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: rel_l1

      call build_channel_mesh(mesh)
      call state_unstr_init(state, mesh)
      call set_dam_break_ic(state, mesh)

      call run_dam_break(state, mesh)

      rel_l1 = relative_l1_error_vs_ritter(state, mesh)

      call check(error, rel_l1 < 0.035_wp, &
                 "L1 error / jump should be below 3.5%")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dam_break_l1

   subroutine test_dam_break_ymom_bounded(error)
      !! The dam break is purely 1D in x: the *physical* y-momentum
      !! is exactly zero everywhere. Any non-zero `ymom_cv` in the
      !! solution is a pure discretisation artifact -- ideally also
      !! zero, but in practice small because the triangulation is not
      !! y-symmetric (same-direction diagonals in every quad) and the
      !! wall mirroring on the bottom wall is not the same operator
      !! as the interior flux on the row-to-row edge.
      !!
      !! This is the right diagnostic for the dam break. A clean
      !! relative bound:
      !!
      !!   max |ymom_cv| <  RELATIVE_TOL * max |xmom_cv|
      !!
      !! is dimensionless, scales with the actual flow magnitude, and
      !! catches gross y-asymmetry (e.g. a sign error in n2) which
      !! would produce |ymom| comparable to |xmom|, while tolerating
      !! the small natural asymmetry from the triangulation. We use
      !! a loose 10% tolerance because the per-cell wall asymmetry
      !! near the rarefaction is non-trivial; a tighter version would
      !! need a smoother triangulation (e.g. alternating diagonals).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: max_xmom, max_ymom
      integer :: k

      call build_channel_mesh(mesh)
      call state_unstr_init(state, mesh)
      call set_dam_break_ic(state, mesh)

      call run_dam_break(state, mesh)

      max_xmom = 0.0_wp
      max_ymom = 0.0_wp
      do k = 1, mesh%ncells_owned
         max_xmom = max(max_xmom, abs(state%barotropic%xmom_cv(k)))
         max_ymom = max(max_ymom, abs(state%barotropic%ymom_cv(k)))
      end do

      ! Sanity: the dam break must produce non-trivial x-flow.
      call check(error, max_xmom > 0.1_wp, &
                 "dam break should produce non-trivial x-momentum")
      if (allocated(error)) return

      ! The actual y-momentum check.
      call check(error, max_ymom < 0.10_wp*max_xmom, &
                 "max |ymom| should be below 10% of max |xmom|")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dam_break_ymom_bounded

   subroutine test_dam_break_positive(error)
      !! No cell may have stage below the bed. Both endpoints stay
      !! near their initial depths and the rarefaction is bounded
      !! between H_RIGHT and H_LEFT, so any depth below H_RIGHT (or
      !! above H_LEFT) is a sign of overshoot.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k
      real(wp) :: depth, depth_min, depth_max

      call build_channel_mesh(mesh)
      call state_unstr_init(state, mesh)
      call set_dam_break_ic(state, mesh)

      call run_dam_break(state, mesh)

      depth_min = huge(0.0_wp)
      depth_max = -huge(0.0_wp)
      do k = 1, mesh%ncells_owned
         depth = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)
         depth_min = min(depth_min, depth)
         depth_max = max(depth_max, depth)
      end do

      call check(error, depth_min > 0.0_wp, &
                 "no negative depths anywhere")
      if (allocated(error)) return
      ! Allow up to 1% overshoot above H_LEFT (slope-limited monotone
      ! schemes are TVD only at first order; second-order can have
      ! small bounded oscillations).
      call check(error, depth_max < 1.01_wp*H_LEFT, &
                 "depth shouldn't overshoot H_LEFT by more than 1%")
      if (allocated(error)) return
      call check(error, depth_min > 0.99_wp*H_RIGHT, &
                 "depth shouldn't undershoot H_RIGHT by more than 1%")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dam_break_positive

   subroutine test_dam_break_endpoints(error)
      !! At t = 0.1 s the disturbance has not reached the channel
      !! ends (left edge of fan moves at -sqrt(g*H_LEFT) ~ -4.4 m/s,
      !! shock moves at ~3 m/s). After 0.1 s the disturbance reaches
      !! roughly x in [4.5, 5.4]. Cells at x < 4 or x > 6 should
      !! still have their initial stage to high accuracy.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: k

      call build_channel_mesh(mesh)
      call state_unstr_init(state, mesh)
      call set_dam_break_ic(state, mesh)

      call run_dam_break(state, mesh)

      do k = 1, mesh%ncells_owned
         if (mesh%centroid_x(k) < 4.0_wp) then
            call check(error, abs(state%barotropic%stage_cv(k) - H_LEFT) < 1.0e-6_wp, &
                       "left undisturbed: stage = H_LEFT")
            if (allocated(error)) return
         else if (mesh%centroid_x(k) > 6.0_wp) then
            call check(error, abs(state%barotropic%stage_cv(k) - H_RIGHT) < 1.0e-6_wp, &
                       "right undisturbed: stage = H_RIGHT")
            if (allocated(error)) return
         end if
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dam_break_endpoints

   subroutine test_dam_break_convergence(error)
      !! Mesh-refinement convergence study. Run the same dam break at
      !! 50, 100, 200 cells along the channel, compute the relative
      !! L1 error vs the exact Ritter solution at each resolution,
      !! and check the L1 convergence rate.
      !!
      !! Theoretical expectation: with a second-order slope-limited
      !! central-upwind scheme on a problem that has a shock plus a
      !! smooth rarefaction, the L1 convergence rate is bounded by
      !!   * 2.0 in smooth regions (full second order)
      !!   * 1.0 at the shock (a TVD scheme is at most first order
      !!     at a shock)
      !! The whole-domain L1 norm averages between these two regimes.
      !! ANUGA-style schemes typically report L1 rates around 0.95 to
      !! 1.05 on the classical wet-bed dam break. We assert
      !! `rate >= 0.8` as a meaningful but not brittle threshold.
      !!
      !! We also assert that the error monotonically decreases with
      !! refinement -- if the coarse-to-fine error ever increased,
      !! that would be a strong sign of a bug, regardless of the rate.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N_LEVELS = 3
      integer, parameter :: NX_LEVELS(N_LEVELS) = [50, 100, 200]
      real(wp) :: errors(N_LEVELS)
      real(wp) :: rate_total
      integer :: level
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      do level = 1, N_LEVELS
         call build_channel_mesh(mesh, nx_in=NX_LEVELS(level), ny_in=NY)
         call state_unstr_init(state, mesh)
         call set_dam_break_ic(state, mesh)

         call run_dam_break(state, mesh)

         errors(level) = relative_l1_error_vs_ritter(state, mesh)

         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
      end do

      ! Errors must monotonically decrease as we refine.
      call check(error, errors(2) < errors(1), &
                 "L1 error must decrease from level 1 to level 2")
      if (allocated(error)) return
      call check(error, errors(3) < errors(2), &
                 "L1 error must decrease from level 2 to level 3")
      if (allocated(error)) return

      ! Convergence rate from the coarsest to the finest mesh:
      !   rate = log(error_coarse / error_fine) / log(nx_fine / nx_coarse)
      ! With NX_LEVELS = [50, 200] this is log2(error_50/error_200) / 2.
      rate_total = log(errors(1)/errors(N_LEVELS)) &
                   /log(real(NX_LEVELS(N_LEVELS), wp)/real(NX_LEVELS(1), wp))

      ! First-order in shocks + second-order in smooth regions gives an
      ! L1 rate in the ~0.85-1.05 range on the classical wet-bed dam break.
      ! With the ANUGA-matched CFL length scale (min centroid-to-edge-midpoint
      ! radius) the measured rate here is ~0.88. The threshold of 0.85 still
      ! catches stalled convergence (rate ~ 0.5) while not being brittle to
      ! the dt/temporal-error contribution at these resolutions.
      call check(error, rate_total >= 0.85_wp, &
                 "L1 convergence rate from level 1 to level 3 should be >= 0.85")
   end subroutine test_dam_break_convergence

end module test_dam_break_unstr
