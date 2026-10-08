!! Well-balanced C-property: lake at rest over irregular bathymetry
module test_c_property
   !! The **C-property** (Bermudez & Vazquez 1994): a fluid at rest over
   !! *arbitrary* bathymetry must stay EXACTLY at rest.  The discrete
   !! pressure-gradient and bed-slope source terms must cancel to round-off,
   !! which is the defining property of the Audusse et al. (2004) hydrostatic
   !! reconstruction used by `rki_kernel_flux` (structured) and the
   !! well-balanced KNP bed-source treatment used by `rki_kernel_flux_unstr`.
   !!
   !! This is an *algebraic identity*, not an approximation: the oracle is the
   !! initial condition itself and the only defensible tolerance is round-off
   !! accumulated over the step count.  A failure here is a well-balancedness
   !! bug in the kernel — never a tolerance to relax.
   !!
   !! ## Relationship to `test_lake_at_rest`
   !!
   !! `test_lake_at_rest` already asserts the C-property, but on a *smooth,
   !! fully-submerged* parabolic bump, structured-only, 2D-only, for 100 steps
   !! at a 1e-10 absolute tolerance.  This module strengthens it on five axes,
   !! each of which is a distinct failure mode:
   !!
   !!   1. **Irregular** bathymetry — two incommensurate sinusoids plus two
   !!      Gaussians, so no grid-aligned symmetry can mask an imbalance.
   !!   2. A **partially-emergent** bump (peak above the free surface), so the
   !!      wet/dry front sits ON the topography.  This is the case where
   !!      well-balancedness classically breaks.
   !!   3. **1000 steps** rather than 100 — a linearly-growing imbalance that
   !!      hides under 1e-10 at 100 steps shows at 1000.
   !!   4. **Both backends** — the unstructured KNP path places eta at the
   !!      triangle circumcentre, a different balance argument from the
   !!      structured A-grid reconstruction.
   !!   5. **3D multilayer** with a uniform-density column, which must also
   !!      stay at rest: this exercises the baroclinic pressure gradient
   !!      (which must vanish identically for flat rho) and the per-layer
   !!      thickness bookkeeping over the same irregular bed.
   !!
   !! ## Why the assertion is on MOMENTUM, not velocity
   !!
   !! The plan's metric is `max|u|, max|v|`.  Velocity `u = hu/h` is
   !! ill-conditioned at the wet/dry edge: a cell with `h ~ 1e-7 m` turns a
   !! perfectly acceptable `hu ~ 1e-15` into `u ~ 1e-8` and the test would
   !! report a false violation whose magnitude depends only on how close a
   !! cell sits to the dry threshold.  Momentum is the quantity the flux
   !! balance actually conserves, so `hu`/`hv` is the well-posed observable.
   !! We assert on momentum scaled by `H*sqrt(g*H)` (the natural momentum
   !! scale), which is *equivalent* to the velocity bound in every cell of
   !! O(H) depth and strictly better-behaved in the thin ones.  A depth-gated
   !! velocity check is kept as a second assertion so the plan's literal
   !! metric is still covered where it is meaningful.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, PI
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
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
   use rki_solver_unstr, only: solver_unstr_step
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   implicit none
   private

   public :: collect_c_property_tests

   ! ---- Domain (shared by both backends so the bathymetry is identical) ----
   integer, parameter :: NX = 100
   integer, parameter :: NY = 100
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = 1.0_wp
   real(wp), parameter :: DY = 1.0_wp
   real(wp), parameter :: CFL_NUM = 0.45_wp
   integer, parameter :: N_STEPS = 1000

   ! ---- Free surface ----
   real(wp), parameter :: ETA0 = 1.0_wp
      !! Flat lake surface (m).  Also the depth scale H: the deepest water is
      !! ETA0 (where b = 0) so H = ETA0 exactly.

   ! ---- Irregular bathymetry (see bed_elevation) ----
   real(wp), parameter :: B_OFFSET = 0.20_wp    !! lifts the ripple clear of b<0
   real(wp), parameter :: B_RIPPLE = 0.18_wp    !! amplitude of the sinusoid pair
   real(wp), parameter :: LAMBDA_X = 23.0_wp    !! incommensurate with LAMBDA_Y
   real(wp), parameter :: LAMBDA_Y = 17.0_wp    !! ... and with the grid spacing
   real(wp), parameter :: B_SUB = 0.45_wp       !! submerged Gaussian height
   real(wp), parameter :: X_SUB = 30.0_wp, Y_SUB = 30.0_wp
   real(wp), parameter :: R_SUB = 200.0_wp
   real(wp), parameter :: B_EMERGE = 1.30_wp    !! emergent Gaussian: peak > ETA0
   real(wp), parameter :: X_EMERGE = 70.0_wp, Y_EMERGE = 70.0_wp
   real(wp), parameter :: R_EMERGE = 60.0_wp

   ! ---- Tolerances: round-off accumulated over N_STEPS, nothing more ----
   real(wp), parameter :: MOM_SCALE = ETA0*sqrt(GRAVITY*ETA0)
      !! Natural momentum scale H*sqrt(gH) (m^2/s).
   real(wp), parameter :: VEL_SCALE = sqrt(GRAVITY*ETA0)
      !! Natural velocity scale sqrt(gH) (m/s).
   real(wp), parameter :: MOM_TOL = 1.0e-12_wp*MOM_SCALE
   real(wp), parameter :: VEL_TOL = 1.0e-12_wp*VEL_SCALE
   real(wp), parameter :: ETA_TOL = 1.0e-11_wp*ETA0

   real(wp), parameter :: H_VEL_MIN = 0.1_wp*ETA0
      !! Depth gate for the velocity assertion.  Cells shallower than this sit
      !! in the wet/dry transition where u = hu/h is ill-conditioned; the
      !! momentum assertion covers them.

   ! ---- Multilayer ----
   integer, parameter :: NZ = 5
   real(wp), parameter :: RHO_0 = 1000.0_wp

contains

   subroutine collect_c_property_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("c_property_structured_momentum", test_str_momentum), &
                  new_unittest("c_property_structured_velocity", test_str_velocity), &
                  new_unittest("c_property_structured_surface", test_str_surface), &
                  new_unittest("c_property_structured_bathy_is_emergent", test_str_emergent), &
                  new_unittest("c_property_unstr_momentum", test_unstr_momentum), &
                  new_unittest("c_property_unstr_surface", test_unstr_surface), &
                  new_unittest("c_property_multilayer_momentum", test_ml_momentum), &
                  new_unittest("c_property_multilayer_surface", test_ml_surface) &
                  ]
   end subroutine collect_c_property_tests

   ! ------------------------------------------------------------------
   ! Bathymetry — shared closed form, evaluated at a physical (x, y)
   ! ------------------------------------------------------------------

   pure function bed_elevation(x, y) result(b)
      !! Irregular bed: two incommensurate sinusoids (no grid-aligned
      !! symmetry) + a fully-submerged Gaussian + an EMERGENT Gaussian whose
      !! peak (B_OFFSET + B_EMERGE = 1.5 m) stands proud of the ETA0 = 1 m
      !! surface, so the lake contains a dry island with a wet/dry front
      !! running over sloping topography.
      !!
      !! Clamped at 0 from below so the bed never goes negative; the clamp
      !! also introduces a slope discontinuity, which a well-balanced scheme
      !! must handle just as exactly as a smooth one.
      real(wp), intent(in) :: x, y
      real(wp) :: b
      real(wp) :: r2_sub, r2_emerge

      r2_sub = (x - X_SUB)**2 + (y - Y_SUB)**2
      r2_emerge = (x - X_EMERGE)**2 + (y - Y_EMERGE)**2

      b = B_OFFSET &
          + B_RIPPLE*sin(2.0_wp*PI*x/LAMBDA_X)*cos(2.0_wp*PI*y/LAMBDA_Y) &
          + B_SUB*exp(-r2_sub/R_SUB) &
          + B_EMERGE*exp(-r2_emerge/R_EMERGE)

      b = max(b, 0.0_wp)
   end function bed_elevation

   pure function cell_x(i) result(x)
      !! Physical x of structured cell centre i (ghost-inclusive indexing).
      integer, intent(in) :: i
      real(wp) :: x
      x = (real(i - NGHOST, wp) - 0.5_wp)*DX
   end function cell_x

   pure function cell_y(j) result(y)
      integer, intent(in) :: j
      real(wp) :: y
      y = (real(j - NGHOST, wp) - 0.5_wp)*DY
   end function cell_y

   ! ------------------------------------------------------------------
   ! Structured 2D barotropic
   ! ------------------------------------------------------------------

   subroutine setup_structured(state)
      !! Flat surface at ETA0 over the irregular bed, at rest.
      !!
      !! NOTE the ghost rows are filled by the same formula as the interior.
      !! A formula bathymetry setter that leaves ghosts at zero produces a
      !! spurious bed step at the wall-adjacent face and destroys the balance
      !! (this is the `formula-bathy-ghost` failure class) — filling the FULL
      !! ghost-inclusive index range is load-bearing, not incidental.
      type(state_t), intent(inout) :: state

      integer :: i, j

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%b(i, j) = bed_elevation(cell_x(i), cell_y(j))
            state%barotropic%h(i, j) = max(ETA0 - state%barotropic%b(i, j), 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do
   end subroutine setup_structured

   subroutine run_structured(state)
      type(state_t), intent(inout) :: state
      integer :: n

      call solver_enter_data(state)
      do n = 1, N_STEPS
         call solver_step(state)
      end do
      call solver_exit_data(state)
   end subroutine run_structured

   subroutine test_str_emergent(error)
      !! Guard on the FIXTURE, not the solver: assert the bed really does
      !! pierce the surface and really does leave deep water elsewhere.  If a
      !! future edit to the bathymetry constants silently submerges the
      !! island, the C-property tests above would keep passing while no longer
      !! testing the wet/dry balance they exist to test.
      type(error_type), allocatable, intent(out) :: error

      integer :: i, j
      real(wp) :: b, b_max, h_max
      integer :: n_dry, n_wet

      b_max = 0.0_wp
      h_max = 0.0_wp
      n_dry = 0
      n_wet = 0
      do j = 1, NY
         do i = 1, NX
            b = bed_elevation(cell_x(i + NGHOST), cell_y(j + NGHOST))
            b_max = max(b_max, b)
            h_max = max(h_max, max(ETA0 - b, 0.0_wp))
            if (ETA0 - b <= 0.0_wp) then
               n_dry = n_dry + 1
            else
               n_wet = n_wet + 1
            end if
         end do
      end do

      call check(error, b_max > ETA0, &
                 "fixture: bathymetry never emerges above the free surface")
      if (allocated(error)) return
      call check(error, n_dry >= 10, &
                 "fixture: emergent island spans too few cells to be a real wet/dry front")
      if (allocated(error)) return
      call check(error, n_wet > 8000, &
                 "fixture: too little wet domain left to be a lake")
      if (allocated(error)) return
      call check(error, h_max > 0.5_wp*ETA0, &
                 "fixture: no deep water left, velocity gate would exclude everything")
   end subroutine test_str_emergent

   subroutine test_str_momentum(error)
      !! Primary C-property assertion: no spurious momentum anywhere.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: max_hu, max_hv

      call setup_structured(state)
      call run_structured(state)

      max_hu = 0.0_wp
      max_hv = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            max_hu = max(max_hu, abs(state%barotropic%hu(i, j)))
            max_hv = max(max_hv, abs(state%barotropic%hv(i, j)))
         end do
      end do

      call check(error, max_hu < MOM_TOL, &
                 "structured C-property: spurious x-momentum over irregular bed")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, max_hv < MOM_TOL, &
                 "structured C-property: spurious y-momentum over irregular bed")

      call state%destroy()
   end subroutine test_str_momentum

   subroutine test_str_velocity(error)
      !! The plan's literal metric (max|u|, max|v| < 1e-12*sqrt(gH)), applied
      !! where it is well-conditioned: cells at least H_VEL_MIN deep.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: max_u, max_v, h

      call setup_structured(state)
      call run_structured(state)

      max_u = 0.0_wp
      max_v = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            h = state%barotropic%h(i, j)
            if (h > H_VEL_MIN) then
               max_u = max(max_u, abs(state%barotropic%hu(i, j)/h))
               max_v = max(max_v, abs(state%barotropic%hv(i, j)/h))
            end if
         end do
      end do

      call check(error, max_u < VEL_TOL, &
                 "structured C-property: spurious x-velocity in the deep interior")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, max_v < VEL_TOL, &
                 "structured C-property: spurious y-velocity in the deep interior")

      call state%destroy()
   end subroutine test_str_velocity

   subroutine test_str_surface(error)
      !! Free surface stays flat at ETA0 across the WET domain.  Dry cells are
      !! excluded by construction: there eta = b > ETA0 and the "flat surface"
      !! statement does not apply.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: max_eta_err, eta_now

      call setup_structured(state)
      call run_structured(state)

      max_eta_err = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            if (state%barotropic%h(i, j) > DRY_TOLERANCE) then
               eta_now = state%barotropic%h(i, j) + state%barotropic%b(i, j)
               max_eta_err = max(max_eta_err, abs(eta_now - ETA0))
            end if
         end do
      end do

      call check(error, max_eta_err < ETA_TOL, &
                 "structured C-property: free surface departed from flat")

      call state%destroy()
   end subroutine test_str_surface

   ! ------------------------------------------------------------------
   ! Structured 3D multilayer (uniform density)
   ! ------------------------------------------------------------------

   subroutine setup_multilayer(state)
      !! Same lake, resolved into NZ sigma layers with a UNIFORM density
      !! column.  With rho constant the baroclinic pressure gradient must
      !! vanish identically, so the 3D problem reduces to the 2D one and the
      !! same exact-rest oracle applies.  Any imbalance here that the 2D case
      !! does not show is in the layer machinery (BPG, thickness
      !! distribution, per-layer bookkeeping), not in the flux.
      type(state_t), intent(inout) :: state

      integer :: i, j

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%rho_0 = RHO_0
      state%multilayer%beta_S = 0.0_wp   ! decouple S from density
      state%multilayer%alpha_T = 0.0_wp  ! decouple T from density

      call state%init(NX, NY, NGHOST, DX, DY, CFL_NUM)

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            state%barotropic%b(i, j) = bed_elevation(cell_x(i), cell_y(j))
            state%barotropic%h(i, j) = max(ETA0 - state%barotropic%b(i, j), 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Uniform density: with beta_S = alpha_T = 0 the EOS returns rho_0
      ! everywhere, so set it directly and skip the EOS call.
      state%multilayer%rho_layer(:, :, :) = RHO_0
   end subroutine setup_multilayer

   subroutine test_ml_momentum(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j, k
      real(wp) :: max_hu, max_hv

      call setup_multilayer(state)
      call run_structured(state)

      max_hu = 0.0_wp
      max_hv = 0.0_wp
      do k = 1, NZ
         do j = NGHOST + 1, state%grid%ny_total - NGHOST
            do i = NGHOST + 1, state%grid%nx_total - NGHOST
               max_hu = max(max_hu, abs(state%multilayer%hu_layer(i, j, k)))
               max_hv = max(max_hv, abs(state%multilayer%hv_layer(i, j, k)))
            end do
         end do
      end do

      call check(error, max_hu < MOM_TOL, &
                 "multilayer C-property: spurious per-layer x-momentum")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, max_hv < MOM_TOL, &
                 "multilayer C-property: spurious per-layer y-momentum")

      call state%destroy()
   end subroutine test_ml_momentum

   subroutine test_ml_surface(error)
      !! Layer thicknesses must still sum to the original depth and the
      !! surface must still be flat.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j, k
      real(wp) :: max_eta_err, h_sum

      call setup_multilayer(state)
      call run_structured(state)

      max_eta_err = 0.0_wp
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            h_sum = 0.0_wp
            do k = 1, NZ
               h_sum = h_sum + state%multilayer%h_layer(i, j, k)
            end do
            if (h_sum > DRY_TOLERANCE) then
               max_eta_err = max(max_eta_err, &
                                 abs(h_sum + state%barotropic%b(i, j) - ETA0))
            end if
         end do
      end do

      call check(error, max_eta_err < ETA_TOL, &
                 "multilayer C-property: summed layer thickness departed from flat surface")

      call state%destroy()
   end subroutine test_ml_surface

   ! ------------------------------------------------------------------
   ! Unstructured 2D barotropic
   ! ------------------------------------------------------------------

   subroutine build_lake_mesh(mesh, nx_in, ny_in)
      !! Triangulate [0, NX*DX] x [0, NY*DY], two CCW triangles per quad.
      !! Same convention as `test_dam_break_unstr`.
      type(mesh_unstr_t), intent(out) :: mesh
      integer, intent(in) :: nx_in, ny_in

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: i, j, nv, ntri, idx
      integer :: v00, v10, v11, v01, t
      real(wp) :: hx, hy

      hx = real(NX, wp)*DX/real(nx_in, wp)
      hy = real(NY, wp)*DY/real(ny_in, wp)

      nv = (nx_in + 1)*(ny_in + 1)
      allocate (vx(nv), vy(nv))
      do j = 0, ny_in
         do i = 0, nx_in
            idx = j*(nx_in + 1) + i + 1
            vx(idx) = real(i, wp)*hx
            vy(idx) = real(j, wp)*hy
         end do
      end do

      ntri = 2*nx_in*ny_in
      allocate (tv(3, ntri))
      do j = 0, ny_in - 1
         do i = 0, nx_in - 1
            v00 = j*(nx_in + 1) + i + 1
            v10 = j*(nx_in + 1) + (i + 1) + 1
            v11 = (j + 1)*(nx_in + 1) + (i + 1) + 1
            v01 = (j + 1)*(nx_in + 1) + i + 1
            t = 2*(j*nx_in + i) + 1
            tv(:, t) = [v00, v10, v11]
            tv(:, t + 1) = [v00, v11, v01]
         end do
      end do

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
      deallocate (vx, vy, tv)
   end subroutine build_lake_mesh

   subroutine setup_unstr(state, mesh)
      !! Same closed-form bed sampled at triangle centroids; stage flat at
      !! ETA0, at rest.  The unstructured state carries absolute `stage_cv`
      !! (= eta) rather than depth, so "flat lake" is literally
      !! `stage_cv = ETA0` — but only where that leaves positive depth; over
      !! the emergent island the stage must be the bed itself (zero depth),
      !! otherwise the IC starts with water above the island.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k
      real(wp) :: b

      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      do k = 1, mesh%ncells_total
         b = bed_elevation(mesh%centroid_x(k), mesh%centroid_y(k))
         state%barotropic%bed_cv(k) = b
         state%barotropic%stage_cv(k) = max(ETA0, b)
      end do
   end subroutine setup_unstr

   subroutine run_unstr(state, mesh, dt)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt
      integer :: n

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      do n = 1, N_STEPS
         call solver_unstr_step(state, mesh, dt)
      end do
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_unstr

   subroutine setup_and_run_unstr(state, mesh)
      !! Shared fixture.  A coarser triangulation than the structured case
      !! (50x50 quads = 5000 triangles) keeps 1000 steps affordable while
      !! still resolving both Gaussians and the wet/dry front.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      real(wp) :: dt

      call build_lake_mesh(mesh, 50, 50)
      call state_unstr_init(state, mesh)
      call setup_unstr(state, mesh)

      ! Fixed dt from the CFL condition.  The controlling length is
      ! `mesh%radius` — the min centroid-to-edge-midpoint distance, which is
      ! what `compute_timestep_unstr` divides by — NOT the quad edge: for the
      ! right-triangle split used here the radius is ~3.4x smaller, and sizing
      ! dt off the quad overshoots the stability limit and blows the run up
      ! (IEEE_OVERFLOW).  Take it from the mesh so the fixture cannot drift
      ! away from the solver's own criterion.  The lake is at rest, so the
      ! wave speed is sqrt(g*ETA0); a fixed step keeps the test deterministic.
      dt = CFL_NUM*minval(mesh%radius(1:mesh%ncells_owned))/sqrt(GRAVITY*ETA0)

      call run_unstr(state, mesh, dt)
   end subroutine setup_and_run_unstr

   subroutine test_unstr_momentum(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: k
      real(wp) :: max_xm, max_ym

      call setup_and_run_unstr(state, mesh)

      max_xm = 0.0_wp
      max_ym = 0.0_wp
      do k = 1, state%ncells_owned
         max_xm = max(max_xm, abs(state%barotropic%xmom_cv(k)))
         max_ym = max(max_ym, abs(state%barotropic%ymom_cv(k)))
      end do

      call check(error, max_xm < MOM_TOL, &
                 "unstructured C-property: spurious x-momentum over irregular bed")
      if (allocated(error)) then
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if
      call check(error, max_ym < MOM_TOL, &
                 "unstructured C-property: spurious y-momentum over irregular bed")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_momentum

   subroutine test_unstr_surface(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: k
      real(wp) :: max_eta_err, depth

      call setup_and_run_unstr(state, mesh)

      max_eta_err = 0.0_wp
      do k = 1, state%ncells_owned
         depth = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)
         if (depth > DRY_TOLERANCE) then
            max_eta_err = max(max_eta_err, &
                              abs(state%barotropic%stage_cv(k) - ETA0))
         end if
      end do

      call check(error, max_eta_err < ETA_TOL, &
                 "unstructured C-property: free surface departed from flat")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_surface

end module test_c_property
