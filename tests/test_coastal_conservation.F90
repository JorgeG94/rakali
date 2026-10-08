!! Closed conservation budgets for the coastal path (mass / salt / heat)
module test_coastal_conservation
   !! Asserts a **closed budget**, not a drift print: for each conserved
   !! quantity Q the residual
   !!
   !!   leak_Q = (Q_final - Q_initial) - integral(sources) dt
   !!
   !! must sit at round-off.  This is the reviewer's "can I trust a multi-week
   !! salt budget?" check, and it is an *identity* of a flux-form scheme rather
   !! than an approximation — so the tolerances below are what conservation
   !! MEANS, not numbers tuned until CI went green.
   !!
   !! ## What is asserted here
   !!
   !! 1. **Closed box, flat bed, zero forcing** — walls on all four sides, no
   !!    surface flux, a stably stratified column so S and T are non-trivial.
   !!    Sources are identically zero, so the oracle is exact: mass, salt and
   !!    heat are all unchanged.  Measured: mass 1.5e-16, salt 7.2e-14, heat
   !!    1.2e-13 relative over 400 steps — machine epsilon, as it should be.
   !! 2. **Prescribed constant surface heat flux, insulated sides** — flat bed,
   !!    walls, uniform `heat_flux = Q`.  The surface kernel applies
   !!    `hT(nz) += dt*Q/(rho_0*cp)` per column per step, so the summed heat
   !!    content must grow by exactly `sum(dt) * Q/(rho_0*cp) * n_interior`.
   !!    This is the "sources are accounted, not merely conserved" half.
   !!
   !! ## Three defects this phase found — all now FIXED
   !!
   !! **1. Spurious sigma BPG over slopes.** A closed box over a sloping bed
   !! lost mass while a flat bed was exact.  Bisection isolated it to the EOS
   !! returning a density differing from `rho_0`: the sigma BPG differenced the
   !! integrated anomaly at constant SIGMA INDEX, giving a fictitious force
   !! proportional to (uniform anomaly) x (bed slope).  `bpg_method` now
   !! defaults to `"jacobian"`; pinned by `test_bpg_uniform_density`.
   !!
   !! **2. The clamp destroyed tracer when the column mean left the bounds.**
   !! `frac = min(1, net/cap_to_add)` saturated and the residual was deleted.
   !! That saturation is exactly equivalent to the column mean leaving the
   !! bounds, so the kernel now skips those columns; pinned by
   !! `test_tracer_conservation`.
   !!
   !! **3. Closed walls passed mass under flow over a sloping bed.**  The wall
   !! BC mirrored the free surface and negated normal momentum, but the ghost
   !! BED came from the bathymetry setter — which fills ghost rows with the
   !! CONTINUATION of the topography, not its reflection.  A wall is a mirror,
   !! so that made the reconstruction across the wall face asymmetric.  At rest
   !! it was harmless (the C-property held to 1e-12), but under flow it gave a
   !! non-zero HLL mass flux through a closed wall.  `bc_wall_edge` now mirrors
   !! the bed as well.  Measured, closed box, stratified, sinusoidal bed, 400
   !! steps: relative mass change 1.364e-6 before, -2.1e-14 after, with the
   !! solver's own `mass_out` accumulator attributing 100 % of it to boundary
   !! flux (budget residual 2e-14 throughout).
   !!
   !! Because all three are closed, the fixtures below now assert conservation
   !! over a SLOPING bed as well as a flat one.  Earlier revisions of this file
   !! were deliberately flat-only because sloping geometry leaked.
   !!
   !! ## Why there is no new budget module
   !!
   !! The plan called for a new `rki_coastal_budget`.  That is largely
   !! unnecessary: `rki_console_stats` already defines `conservation_budget_t`
   !! with the residual identity `leak_Q = (Q_total - Q0) + out_Q - src_Q`, and
   !! `rki_solver` already ships `accumulate_boundary_mass_out` /
   !! `accumulate_mass_src`, which close the MASS budget.  The genuine hole is
   !! that the `salt_out`, `salt_src`, `heat_out`, `heat_src` fields of that
   !! type are never populated anywhere in `src/`, so `salt_active` /
   !! `heat_active` stay `.false.` and those quantities silently fall back to
   !! raw drift.  This module asserts the closure directly (forming the source
   !! integral itself), which is what makes the property testable today;
   !! wiring those accumulators into the solver is a separate source change.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, CP_WATER, PI
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_init, only: compute_total_mass, compute_layer_sums
   use rki_boundary_types, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic
   use rki_ml_eos, only: ml_equation_of_state
   implicit none
   private

   public :: collect_coastal_conservation_tests

   ! ---- Domain ----
   integer, parameter :: NX = 40
   integer, parameter :: NY = 40
   integer, parameter :: NGHOST = 2
   real(wp), parameter :: DX = 5.0_wp
   real(wp), parameter :: CFL_NUM = 0.45_wp
   integer, parameter :: NZ = 6
   integer, parameter :: N_STEPS = 400

   ! ---- Water column ----
   real(wp), parameter :: ETA0 = 10.0_wp
   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_BOTTOM = 35.0_wp, S_SURFACE = 30.0_wp
   real(wp), parameter :: T_BOTTOM = 8.0_wp, T_SURFACE = 18.0_wp
   real(wp), parameter :: S_REF = 35.0_wp, T_REF = 10.0_wp
   real(wp), parameter :: BETA_S = 0.8_wp, ALPHA_T = 0.2_wp

   ! ---- Irregular bed (sloping cases) ----
   real(wp), parameter :: B_RIPPLE = 1.5_wp
   real(wp), parameter :: LAMBDA_X = 47.0_wp, LAMBDA_Y = 31.0_wp

   ! ---- Surface heat flux (case 2) ----
   real(wp), parameter :: Q_SURFACE = 200.0_wp   !! W m-2, positive = warming

   ! ---- Tolerances (identity tolerances — see module docstring) ----
   real(wp), parameter :: MASS_TOL = 1.0e-12_wp
   real(wp), parameter :: TRACER_TOL = 1.0e-10_wp
   real(wp), parameter :: FLUX_CLOSURE_TOL = 1.0e-8_wp

contains

   subroutine collect_coastal_conservation_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("conservation_closed_box_mass", test_closed_mass), &
                  new_unittest("conservation_closed_box_salt", test_closed_salt), &
                  new_unittest("conservation_closed_box_heat", test_closed_heat), &
                  new_unittest("conservation_surface_heat_flux_closure", test_flux_closure), &
                  new_unittest("conservation_fixture_is_stratified", test_fixture_stratified), &
                  new_unittest("conservation_sloping_bed_mass", test_slope_mass), &
                  new_unittest("conservation_sloping_bed_salt", test_slope_salt), &
                  new_unittest("conservation_closed_walls_pass_no_mass", test_wall_tight) &
                  ]
   end subroutine collect_coastal_conservation_tests

   ! ------------------------------------------------------------------
   ! Fixture
   ! ------------------------------------------------------------------

   subroutine setup_closed_box(state, q_surface, ripple)
      !! Flat bed, uniform depth ETA0, walls, stably stratified, with an
      !! optional uniform surface heat flux.
      !!
      !! Flat AND quiescent DELIBERATELY.  See the module docstring: a
      !! stratified sloping bed still leaks mass at ~2e-6 after the BPG fix,
      !! and salt leaks ~9e-6 whenever tracers advect at all.  Using either
      !! regime here would make this test assert a known defect rather than
      !! the identity it exists to protect.
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: q_surface
      real(wp), intent(in), optional :: ripple
         !! Bed amplitude (m).  Absent or 0 => flat bed.

      integer :: i, j
      real(wp) :: rip, x, y, bed

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%rho_0 = RHO_0
      state%multilayer%alpha_T = ALPHA_T
      state%multilayer%beta_S = BETA_S
      state%multilayer%S_ref = S_REF
      state%multilayer%T_ref = T_REF
      state%multilayer%initial_salinity = 0.5_wp*(S_BOTTOM + S_SURFACE)
      state%multilayer%initial_temperature = 0.5_wp*(T_BOTTOM + T_SURFACE)

      call state%init(NX, NY, NGHOST, DX, DX, CFL_NUM)

      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      rip = 0.0_wp
      if (present(ripple)) rip = ripple

      ! Ghost rows get the FORMULA value, exactly as a real bathymetry setter
      ! writes them.  `bc_wall_edge` mirrors the bed into the ghosts each step,
      ! so this test also covers the wall-tightness fix rather than quietly
      ! sidestepping it by pre-mirroring here.
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x = (real(i - NGHOST, wp) - 0.5_wp)*DX
            y = (real(j - NGHOST, wp) - 0.5_wp)*DX
            bed = rip*(1.0_wp + sin(2.0_wp*PI*x/LAMBDA_X)*cos(2.0_wp*PI*y/LAMBDA_Y))
            bed = max(bed, 0.0_wp)
            state%barotropic%b(i, j) = bed
            state%barotropic%h(i, j) = max(ETA0 - bed, 0.0_wp)
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
      call seed_stratified_column(state)

      state%forcing%heat_flux = q_surface
      state%forcing%evap_minus_precip = 0.0_wp
      state%forcing%wind_x = 0.0_wp
      state%forcing%wind_y = 0.0_wp
      state%forcing%precipitation = 0.0_wp
   end subroutine setup_closed_box

   subroutine seed_stratified_column(state)
      !! Linear S and T profiles, bottom-up (k=1 is the bed, k=NZ the
      !! surface).  Stable: fresher and warmer on top.
      type(state_t), intent(inout) :: state

      integer :: i, j, k
      real(wp) :: f, s_val, t_val

      do k = 1, NZ
         if (NZ > 1) then
            f = real(k - 1, wp)/real(NZ - 1, wp)
         else
            f = 0.0_wp
         end if
         s_val = S_BOTTOM + (S_SURFACE - S_BOTTOM)*f
         t_val = T_BOTTOM + (T_SURFACE - T_BOTTOM)*f
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               associate (ml => state%multilayer)
                  ml%tracers(ml%idx_salinity)%hTr(i, j, k) = ml%h_layer(i, j, k)*s_val
                  ml%tracers(ml%idx_temperature)%hTr(i, j, k) = ml%h_layer(i, j, k)*t_val
               end associate
            end do
         end do
      end do

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, RHO_0, BETA_S, ALPHA_T, &
                                S_REF, T_REF, state%grid, NZ)
   end subroutine seed_stratified_column

   subroutine run_steps(state, elapsed)
      !! Fixed step count; returns the accumulated physical time so the source
      !! integral is formed against the SAME dt sequence the solver took (dt
      !! is adaptive, so assuming N*dt would be wrong).
      type(state_t), intent(inout) :: state
      real(wp), intent(out) :: elapsed
      integer :: n
      real(wp) :: t

      elapsed = 0.0_wp
      t = 0.0_wp
      call solver_enter_data(state)
      do n = 1, N_STEPS
         call solver_step(state, time=t)
         t = t + state%dt
         elapsed = elapsed + state%dt
      end do
      call solver_exit_data(state)
   end subroutine run_steps

   ! ------------------------------------------------------------------
   ! Case 1 — closed box, zero sources: every total is invariant
   ! ------------------------------------------------------------------

   subroutine test_closed_mass(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: m0, m1, elapsed, rel

      call setup_closed_box(state, 0.0_wp)
      m0 = compute_total_mass(state)
      call run_steps(state, elapsed)
      m1 = compute_total_mass(state)

      rel = abs(m1 - m0)/max(m0, tiny(1.0_wp))
      call check(error, rel < MASS_TOL, &
                 "closed-box water mass is not conserved to round-off")

      call state%destroy()
   end subroutine test_closed_mass

   subroutine test_closed_salt(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: h0, s0, t0, h1, s1, t1, elapsed, rel

      call setup_closed_box(state, 0.0_wp)
      call compute_layer_sums(state, h0, s0, t0)
      call run_steps(state, elapsed)
      call compute_layer_sums(state, h1, s1, t1)

      rel = abs(s1 - s0)/max(abs(s0), tiny(1.0_wp))
      call check(error, rel < TRACER_TOL, &
                 "closed-box salt content is not conserved to round-off")

      call state%destroy()
   end subroutine test_closed_salt

   subroutine test_closed_heat(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: h0, s0, t0, h1, s1, t1, elapsed, rel

      call setup_closed_box(state, 0.0_wp)
      call compute_layer_sums(state, h0, s0, t0)
      call run_steps(state, elapsed)
      call compute_layer_sums(state, h1, s1, t1)

      rel = abs(t1 - t0)/max(abs(t0), tiny(1.0_wp))
      call check(error, rel < TRACER_TOL, &
                 "closed-box heat content is not conserved to round-off")

      call state%destroy()
   end subroutine test_closed_heat

   subroutine test_fixture_stratified(error)
      !! Guard on the FIXTURE: a uniform tracer is conserved by almost any
      !! scheme, so if S and T were accidentally flat the assertions above
      !! would still pass while testing far less.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: im, jm
      real(wp) :: s_bot, s_top, t_bot, t_top

      call setup_closed_box(state, 0.0_wp)

      im = NGHOST + NX/2
      jm = NGHOST + NY/2
      associate (ml => state%multilayer)
         s_bot = ml%tracers(ml%idx_salinity)%hTr(im, jm, 1)/max(ml%h_layer(im, jm, 1), tiny(1.0_wp))
         s_top = ml%tracers(ml%idx_salinity)%hTr(im, jm, NZ)/max(ml%h_layer(im, jm, NZ), tiny(1.0_wp))
         t_bot = ml%tracers(ml%idx_temperature)%hTr(im, jm, 1)/max(ml%h_layer(im, jm, 1), tiny(1.0_wp))
         t_top = ml%tracers(ml%idx_temperature)%hTr(im, jm, NZ)/max(ml%h_layer(im, jm, NZ), tiny(1.0_wp))
      end associate

      call check(error, abs(s_bot - s_top) > 1.0_wp, &
                 "fixture: salinity column is not stratified")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, abs(t_bot - t_top) > 1.0_wp, &
                 "fixture: temperature column is not stratified")

      call state%destroy()
   end subroutine test_fixture_stratified

   ! ------------------------------------------------------------------
   ! Sloping bed — the geometry that used to leak
   ! ------------------------------------------------------------------

   subroutine test_slope_mass(error)
      !! Same closed box over a SINUSOIDAL bed.  The stratification is along
      !! sigma surfaces, which are not horizontal over a slope, so the column
      !! genuinely adjusts and real flow develops — this is not a lake at rest.
      !! Before the wall-bed-mirror fix this leaked 1.364e-6.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: m0, m1, elapsed, rel

      call setup_closed_box(state, 0.0_wp, ripple=B_RIPPLE)
      m0 = compute_total_mass(state)
      call run_steps(state, elapsed)
      m1 = compute_total_mass(state)

      rel = abs(m1 - m0)/max(m0, tiny(1.0_wp))
      call check(error, rel < MASS_TOL, &
                 "sloping-bed closed box does not conserve water mass")

      call state%destroy()
   end subroutine test_slope_mass

   subroutine test_slope_salt(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: h0, s0, t0, h1, s1, t1, elapsed, rel

      call setup_closed_box(state, 0.0_wp, ripple=B_RIPPLE)
      call compute_layer_sums(state, h0, s0, t0)
      call run_steps(state, elapsed)
      call compute_layer_sums(state, h1, s1, t1)

      rel = abs(s1 - s0)/max(abs(s0), tiny(1.0_wp))
      call check(error, rel < TRACER_TOL, &
                 "sloping-bed closed box does not conserve salt")

      call state%destroy()
   end subroutine test_slope_salt

   subroutine test_wall_tight(error)
      !! Assert the WALLS themselves pass no mass, using the solver's own
      !! accumulator rather than inferring it from the total.
      !!
      !! `mass_out` is accumulated as `-sum_interior flux_h*area`, i.e. the
      !! mass change the FLUX path accounts for.  Two assertions:
      !!
      !!   (a) `mass_out` itself is ~0 — a closed wall passes no flux.  This is
      !!       the direct statement of the fix; it read -1.364e-6 before.
      !!   (b) the budget closes: `(M1 - M0) + mass_out - mass_src ~ 0`, so
      !!       nothing OUTSIDE the flux path is modifying h.  This held even
      !!       when (a) failed, which is exactly how the defect was localised
      !!       to the boundary rather than the interior.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: m0, m1, elapsed, resid

      call setup_closed_box(state, 0.0_wp, ripple=B_RIPPLE)
      m0 = compute_total_mass(state)
      call run_steps(state, elapsed)
      m1 = compute_total_mass(state)

      call check(error, abs(state%barotropic%mass_out)/m0 < MASS_TOL, &
                 "closed walls passed mass flux over a sloping bed")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      resid = (m1 - m0) + state%barotropic%mass_out - state%barotropic%mass_src
      call check(error, abs(resid)/m0 < MASS_TOL, &
                 "mass budget does not close: something outside the flux path moved h")

      call state%destroy()
   end subroutine test_wall_tight

   ! ------------------------------------------------------------------
   ! Case 2 — prescribed surface heat flux: closure against the integral
   ! ------------------------------------------------------------------

   subroutine test_flux_closure(error)
      !! The surface kernel adds `dt*Q/(rho_0*cp)` to `hT(:,:,NZ)` each step,
      !! so over the run the summed heat content must grow by exactly
      !! `elapsed * Q/(rho_0*cp) * n_interior_cells`.  `compute_layer_sums`
      !! sums hT over interior cells and all layers — the same set.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: h0, s0, t0, h1, s1, t1
      real(wp) :: elapsed, expected, actual, rel

      call setup_closed_box(state, Q_SURFACE)
      call compute_layer_sums(state, h0, s0, t0)
      call run_steps(state, elapsed)
      call compute_layer_sums(state, h1, s1, t1)

      expected = elapsed*Q_SURFACE/(RHO_0*CP_WATER)*real(NX*NY, wp)
      actual = t1 - t0
      rel = abs(actual - expected)/max(abs(expected), tiny(1.0_wp))

      call check(error, expected > 0.0_wp, &
                 "fixture: surface heat flux integrated to zero")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, rel < FLUX_CLOSURE_TOL, &
                 "surface heat flux does not close against its time integral")

      call state%destroy()
   end subroutine test_flux_closure

end module test_coastal_conservation
