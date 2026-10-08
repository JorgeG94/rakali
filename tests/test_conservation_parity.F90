!! Closed-budget conservation on the unstructured explicit path.
module test_conservation_parity
   !! `test_coastal_conservation` asserts mass / salt / heat closure on the
   !! **structured explicit** path.  Two of the Phase 1 conservation fixes (the
   !! jacobian BPG default and the clamp conservation guard) reach the
   !! unstructured path by shared kernels rather than by measurement — and "by
   !! construction" is precisely the reasoning that let the original defects
   !! sit undetected until someone measured.  This module measures.
   !!
   !! Path covered here:
   !!
   !!   * unstructured x explicit      (KNP on triangles)
   !!
   !! ## NOT covered here: unstructured x semi-implicit
   !!
   !! The fourth cell of the matrix is a real gap, stated rather than glossed.
   !! The unstructured SI path is TPFA and places eta at the triangle
   !! CIRCUMCENTRE, so it requires a Delaunay-quality mesh: `mesh_unstr_t`
   !! carries a fail-loud quality gate (`max_dev_angle`, `min_edge_delta`) that
   !! rejects coincident circumcentres, and the right-triangle split used by
   !! every quad-derived test mesh in this suite — including the one below —
   !! degenerates exactly that way.  Covering it needs the per-quad
   !! diagonal-choice construction that `test_si_unstr_ml_wetdry` uses
   !! (`delaunay_diag_bl_tr`), plus `state_unstr_enable_si`, whose `ok` flag
   !! MUST then be asserted — that existing test does `if (.not. ok) return`,
   !! so a mesh that failed the gate would make it silently vacuous.
   !!
   !! ## Why the unstructured cases are here rather than strengthening
   !! `test_multilayer_unstr`
   !!
   !! That module's `tracer_mass_conserved_in_closed_basin` is a **weak**
   !! conservation check and was never meant to be more: a four-triangle mesh,
   !! a flat bed, 25 steps, `beta_S = alpha_T = 0` so there is no baroclinic
   !! coupling at all, and it asserts on the tracer only — never on water mass.
   !! It is a useful smoke test for the advection split and is left alone.
   !! The cases here use a real mesh, a SLOPING bed, an active density field driving
   !! genuine flow, and assertions on mass as well as salt.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, PI
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_destroy, state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_solver_unstr, only: solver_unstr_step, solver_unstr_init_multilayer
   implicit none
   private

   public :: collect_conservation_parity_tests

   ! ---- Shared physical setup ----
   integer, parameter :: NZ = 4
   real(wp), parameter :: ETA0 = 10.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_BOTTOM = 35.0_wp, S_SURFACE = 30.0_wp
   real(wp), parameter :: T_BOTTOM = 8.0_wp, T_SURFACE = 18.0_wp
   real(wp), parameter :: S_REF = 35.0_wp, T_REF = 10.0_wp
   real(wp), parameter :: BETA_S = 0.8_wp, ALPHA_T = 0.2_wp
   real(wp), parameter :: B_RIPPLE = 1.0_wp
   real(wp), parameter :: LAMBDA = 37.0_wp

   ! ---- Unstructured ----
   integer, parameter :: NQ = 24
      !! Quads per side; 2 triangles each => 1152 cells.
   real(wp), parameter :: L_UNSTR = 160.0_wp
   integer, parameter :: N_STEPS_UNSTR = 200

   ! ---- Tolerances: identity tolerances, as in test_coastal_conservation ----
   real(wp), parameter :: MASS_TOL = 1.0e-12_wp
   real(wp), parameter :: TRACER_TOL = 1.0e-10_wp

contains

   subroutine collect_conservation_parity_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("unstr_explicit_conserves_mass", test_unstr_mass), &
                  new_unittest("unstr_explicit_conserves_salt", test_unstr_salt), &
                  new_unittest("unstr_fixture_is_stratified_and_sloping", test_unstr_fixture), &
                  new_unittest("unstr_fixture_actually_flows", test_unstr_flows) &
                  ]
   end subroutine collect_conservation_parity_tests

   pure function bed_at(x, y) result(b)
      !! Gentle sinusoidal bed. Shared by both fixtures so "sloping" means the
      !! same thing on triangles as on the Cartesian grid.
      real(wp), intent(in) :: x, y
      real(wp) :: b

      b = B_RIPPLE*(1.0_wp + sin(2.0_wp*PI*x/LAMBDA)*cos(2.0_wp*PI*y/LAMBDA))
      b = max(b, 0.0_wp)
   end function bed_at

   ! ------------------------------------------------------------------
   ! Unstructured x explicit
   ! ------------------------------------------------------------------

   subroutine build_square_mesh(mesh)
      !! [0, L]^2 triangulated as two CCW triangles per quad.
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: i, j, nv, ntri, idx, v00, v10, v11, v01, t
      real(wp) :: hq

      hq = L_UNSTR/real(NQ, wp)

      nv = (NQ + 1)*(NQ + 1)
      allocate (vx(nv), vy(nv))
      do j = 0, NQ
         do i = 0, NQ
            idx = j*(NQ + 1) + i + 1
            vx(idx) = real(i, wp)*hq
            vy(idx) = real(j, wp)*hq
         end do
      end do

      ntri = 2*NQ*NQ
      allocate (tv(3, ntri))
      do j = 0, NQ - 1
         do i = 0, NQ - 1
            v00 = j*(NQ + 1) + i + 1
            v10 = j*(NQ + 1) + (i + 1) + 1
            v11 = (j + 1)*(NQ + 1) + (i + 1) + 1
            v01 = (j + 1)*(NQ + 1) + i + 1
            t = 2*(j*NQ + i) + 1
            tv(:, t) = [v00, v10, v11]
            tv(:, t + 1) = [v00, v11, v01]
         end do
      end do

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
      deallocate (vx, vy, tv)
   end subroutine build_square_mesh

   subroutine setup_unstr(state, mesh)
      !! Sloping bed, stratified, active density, walls (the default when no
      !! boundary groups are tagged).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      integer :: cell, k
      real(wp) :: b, f, s_val, t_val

      call build_square_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      state%multilayer%rho_0 = RHO_0
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%S_ref = S_REF
      state%multilayer%T_ref = T_REF
      state%multilayer%initial_salinity = 0.5_wp*(S_BOTTOM + S_SURFACE)
      state%multilayer%initial_temperature = 0.5_wp*(T_BOTTOM + T_SURFACE)

      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      do cell = 1, mesh%ncells_total
         b = bed_at(mesh%centroid_x(cell), mesh%centroid_y(cell))
         state%barotropic%bed_cv(cell) = b
         state%barotropic%stage_cv(cell) = max(ETA0, b)
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call solver_unstr_init_multilayer(state, mesh)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      do k = 1, NZ
         f = real(k - 1, wp)/real(NZ - 1, wp)
         s_val = S_BOTTOM + (S_SURFACE - S_BOTTOM)*f
         t_val = T_BOTTOM + (T_SURFACE - T_BOTTOM)*f
         do cell = 1, mesh%ncells_total
            associate (ml => state%multilayer)
               ml%tracers(ml%idx_salinity)%hTr(k, cell) = ml%h_layer(k, cell)*s_val
               ml%tracers(ml%idx_temperature)%hTr(k, cell) = ml%h_layer(k, cell)*t_val
            end associate
         end do
      end do
   end subroutine setup_unstr

   subroutine run_unstr(state, mesh, dt)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt
      integer :: n

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      do n = 1, N_STEPS_UNSTR
         call solver_unstr_step(state, mesh, dt)
      end do
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_unstr

   subroutine unstr_totals(state, mesh, mass, salt)
      !! Area-weighted totals over the owned cells.  `mesh%area` varies per
      !! triangle, so this is a weighted sum, not a count.
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(out) :: mass, salt

      integer :: cell, k
      real(wp) :: h_k

      mass = 0.0_wp
      salt = 0.0_wp
      do cell = 1, state%ncells_owned
         do k = 1, NZ
            h_k = state%multilayer%h_layer(k, cell)
            if (h_k > 0.0_wp) then
               mass = mass + h_k*mesh%area(cell)
               salt = salt + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(k, cell) &
                      *mesh%area(cell)
            end if
         end do
      end do
   end subroutine unstr_totals

   subroutine unstr_dt(mesh, dt)
      !! CFL step from the mesh's own criterion — `mesh%radius` is the min
      !! centroid-to-edge-midpoint distance that `compute_timestep_unstr`
      !! divides by.  Sizing dt off the quad edge instead overshoots
      !! stability for a right-triangle split (see `test_c_property`).
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(out) :: dt
      real(wp), parameter :: GRAV = 9.80665_wp

      dt = 0.45_wp*minval(mesh%radius(1:mesh%ncells_owned))/sqrt(GRAV*ETA0)
   end subroutine unstr_dt

   subroutine test_unstr_mass(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: m0, s0, m1, s1, dt, rel

      call setup_unstr(state, mesh)
      call unstr_totals(state, mesh, m0, s0)
      call unstr_dt(mesh, dt)
      call run_unstr(state, mesh, dt)
      call unstr_totals(state, mesh, m1, s1)

      rel = abs(m1 - m0)/max(m0, tiny(1.0_wp))
      call check(error, rel < MASS_TOL, &
                 "unstructured closed basin does not conserve water mass")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_mass

   subroutine test_unstr_salt(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: m0, s0, m1, s1, dt, rel

      call setup_unstr(state, mesh)
      call unstr_totals(state, mesh, m0, s0)
      call unstr_dt(mesh, dt)
      call run_unstr(state, mesh, dt)
      call unstr_totals(state, mesh, m1, s1)

      rel = abs(s1 - s0)/max(abs(s0), tiny(1.0_wp))
      call check(error, rel < TRACER_TOL, &
                 "unstructured closed basin does not conserve salt")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_salt

   subroutine test_unstr_fixture(error)
      !! Guard: the unstructured fixture must really be stratified and really
      !! be sloping, or the two assertions above degrade into the weak
      !! flat/uniform check they exist to replace.
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: cell
      real(wp) :: s_bot, s_top, b_min, b_max

      call setup_unstr(state, mesh)

      associate (ml => state%multilayer)
         s_bot = ml%tracers(ml%idx_salinity)%hTr(1, 1)/max(ml%h_layer(1, 1), tiny(1.0_wp))
         s_top = ml%tracers(ml%idx_salinity)%hTr(NZ, 1)/max(ml%h_layer(NZ, 1), tiny(1.0_wp))
      end associate

      b_min = huge(1.0_wp)
      b_max = -huge(1.0_wp)
      do cell = 1, state%ncells_owned
         b_min = min(b_min, state%barotropic%bed_cv(cell))
         b_max = max(b_max, state%barotropic%bed_cv(cell))
      end do

      call check(error, abs(s_bot - s_top) > 1.0_wp, &
                 "fixture: unstructured column is not stratified")
      if (allocated(error)) then
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if
      call check(error, (b_max - b_min) > 0.5_wp, &
                 "fixture: unstructured bed is effectively flat")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_fixture

   subroutine test_unstr_flows(error)
      !! Same guard on the unstructured fixture.  Measured: max|xmom| ~ 1.1e-2.
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: cell
      real(wp) :: dt, mx

      call setup_unstr(state, mesh)
      call unstr_dt(mesh, dt)
      call run_unstr(state, mesh, dt)

      mx = 0.0_wp
      do cell = 1, state%ncells_owned
         mx = max(mx, abs(state%barotropic%xmom_cv(cell)))
      end do

      call check(error, mx > 1.0e-4_wp, &
                 "fixture: unstructured basin never developed flow, conservation test is vacuous")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstr_flows

end module test_conservation_parity
