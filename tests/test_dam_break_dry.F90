!! Dry-bed dam break: the vacuum-front rarefaction against Ritter (1892).
module test_dam_break_dry
   !! Both existing dam-break tests are **wet-bed** (`test_dam_break`,
   !! `test_dam_break_unstr`, and the shock-limited `test_convergence`). The
   !! dry leading edge — where the depth goes to zero and the scheme has to
   !! hold positivity at a vacuum front — is only NaN-smoke-tested anywhere in
   !! the suite. This checks it against the exact solution.
   !!
   !! ## Oracle
   !!
   !! Ritter (1892), `h_R = 0`, in `ritter_dry_sample`. With
   !! `c0 = sqrt(g h_l)` and `xi = (x - x_dam)/t` the fan is a parabola in
   !! `xi` reaching zero exactly at the front, and the front advances at
   !! `2 c0` — twice the still-water wave speed. There is no shock and no star
   !! region, which is why the dry branch is a separate closed form rather
   !! than a limit of the wet Newton solve (as `h_R -> 0` the right shock
   !! degenerates into the rarefaction tip and the iteration loses its
   !! bracket).
   !!
   !! ## Locating the front honestly
   !!
   !! The exact profile approaches the front QUADRATICALLY: a distance `d`
   !! behind it, `h = (d/t)^2 / (9g)`. At `t = 10 s`, one cell back
   !! (`d = 0.25 m`) that is `h ~ 7e-6 m` — right on `DRY_TOLERANCE`. So
   !! "rightmost cell with `h > eps`" is a violently eps-sensitive way to
   !! locate the front, and comparing such a number against the analytic
   !! `x_f = x_dam + 2 c0 t` would mostly measure the threshold.
   !!
   !! Instead the front is located **consistently on both sides**: the same
   !! depth threshold `H_FRONT` is applied to the numerical profile and to the
   !! exact one. Inverting the parabola gives the exact threshold crossing in
   !! closed form,
   !!
   !!     x_thresh = x_f - 3 t sqrt(g H_FRONT)
   !!
   !! so the comparison is like-for-like and the tolerance measures the
   !! scheme, not the probe. `H_FRONT = 1e-3 * h_l` sits ~3 m behind the true
   !! front here — well clear of the round-off tail, still deep inside the
   !! region the parabola describes.
   !!
   !! ## FINDING (2026-08, updated after the thin-cell fix) — the tongue is
   !! now largely represented; a sub-vanished-scale deficit remains
   !!
   !! The original finding: the whole leading tongue below h ~ 0.03 m was
   !! missing and the tip advanced at 1.49 c0.  Cause, established by
   !! dose-response (varying both knobs in `rki_constants` and measuring
   !! the tip): TWO REDUNDANT momentum clamps in the update kernels — the
   !! hard `THIN_LAYER_THRESHOLD` momentum zero below 1 mm AND the
   !! `FROUDE_CAP = 10` rescale (the exact tongue has
   !! `Fr = 2(sqrt(h_l/h) - 1)`, which exceeds 10 exactly where
   !! h < h_l/36 ~ 0.028 m — precisely the missing depths; the leading wet
   !! cell measured at exactly Fr = 10).  Relaxing either alone changed
   !! little because the other kept binding — which is why each was once
   !! separately "eliminated" as a cause.  Relaxing both gave 1.888 c0 but
   !! also a thin-film velocity runaway and dt collapse on beaches.
   !!
   !! The fix (see `DESING_EPS` / `FROUDE_BACKSTOP` in `rki_constants`):
   !! Kurganov & Petrova (2007) desingularized velocity below
   !! eps = 1.5e-4 m plus a Fr = 50 blowup backstop.  Measured front
   !! position after the fix (structured, 800 cells, t = 10 s, 62.6 m of
   !! travel; baseline in parentheses):
   !!
   !!     thresh/h_l    x_num     x_ref    lag(m)     rel
   !!       1e-1       132.88    132.92     0.05    0.0008  (0.05)
   !!       1e-2       152.88    153.24     0.36    0.006   (7.11)
   !!       1e-3       153.62    159.66     6.04    0.096   (13.29)
   !!       1e-6       154.38    162.54     8.16    0.130   (15.91)
   !!
   !! The 1-cm front is now within 0.6 % of exact (was 11.4 %), so
   !! `H_FRONT_STRUCT` is tightened to 1e-2*h_l at the 1.5 % tolerance.
   !! The residual deficit lives below ~1 mm: the desingularization damps
   !! momentum below eps by design (that suppression is what prevents the
   !! measured beach runaway), so the tip advances at ~1.74 c0 rather than
   !! 2 c0 — vs a measured no-clamp ceiling of 1.888 c0 at this
   !! resolution (1st-order smearing accounts for the rest).
   !!
   !! The unstructured tip is unchanged by the same treatment of its
   !! protect kernel (1.586 c0 to three decimals with the fix in or out):
   !! its binding constraint is NOT the protect clamp but the KNP
   !! reconstruction's dry-cell handling (`centroid_pass` zeroes edge
   !! momenta below `minimum_allowed_height`) and coarser dx = 0.5 m.
   !! Its 1-cm front lags 4.40 m (7 %), so `H_FRONT_UNSTR` stays at
   !! 1e-1*h_l — an explicitly recorded residual, not a validated tongue.
   !!
   !! ## Tolerances
   !!
   !! HLL is 1st order at the vacuum-adjacent rarefaction tip, so these are
   !! order-honest bounds rather than aspirations: the wet-bed test achieves
   !! 2.5 % L1 and the dry front is strictly harder. Front position is
   !! normalised by the distance travelled, not by domain length, so it does
   !! not flatter itself as the run gets longer.
   !!
   !! ## Not run on the semi-implicit path
   !!
   !! Deliberate, per the campaign plan. Casulli's semi-implicit mode uses
   !! Eulerian-Lagrangian advection, which smears a sharp dry front by
   !! construction; it is not the intended tool for this regime and asserting
   !! a front position on it would be measuring the wrong thing. The explicit
   !! Godunov path is what buys faithful front physics, and that is what is
   !! validated here.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_init, only: compute_total_mass
   use rki_boundary_types, only: bc_type_from_string
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_destroy, state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_solver_unstr, only: solver_unstr_step
   use ritter_dambreak, only: ritter_dry_sample, ritter_dry_front_position
   implicit none
   private

   public :: collect_dam_break_dry_tests

   ! ---- Channel ----
   real(wp), parameter :: LX = 200.0_wp
   real(wp), parameter :: X_DAM = 100.0_wp
   real(wp), parameter :: H_LEFT = 1.0_wp
   real(wp), parameter :: T_END = 10.0_wp
      !! c0 = 3.13 m/s, so the front reaches x = 162.7 m and the rarefaction
      !! head retreats to x = 68.7 m — both comfortably clear of the walls, so
      !! the boundary treatment never enters the comparison.
   integer, parameter :: NX = 800
   integer, parameter :: NY = 4
   integer, parameter :: NGHOST = 2

   real(wp), parameter :: H_FRONT_STRUCT = 1.0e-2_wp*H_LEFT
      !! Structured front-detection depth, applied identically to the
      !! numerical and exact profiles.  Tightened from 1e-1*h_l after the
      !! thin-cell desingularization fix: at 1e-2*h_l the two now agree to
      !! 0.36 m out of 62.6 m travelled (0.6 %), inside the 1.5 %
      !! tolerance.  Below ~1e-3*h_l the desingularization band still lags
      !! (see the module docstring) — tightening further would assert a
      !! depth the scheme deliberately damps.
   real(wp), parameter :: H_FRONT_UNSTR = 1.0e-1_wp*H_LEFT
      !! Unstructured front-detection depth.  NOT tightened: the KNP
      !! path's 1-cm front still lags 7 % (its dry-cell reconstruction,
      !! not the protect clamp, is binding — see the module docstring), so
      !! 1e-1*h_l remains the deepest honestly-validated contour there
      !! (0.09 m lag, 0.14 %).

   ! ---- Tolerances ----
   real(wp), parameter :: FRONT_TOL = 0.015_wp    !! 1.5 % of distance travelled
   real(wp), parameter :: L1_TOL = 0.05_wp        !! 5 % relative L1 on depth
   real(wp), parameter :: MASS_TOL = 1.0e-10_wp

   ! ---- Unstructured channel ----
   integer, parameter :: NQX = 400, NQY = 2
   real(wp), parameter :: LY_UNSTR = 1.0_wp

contains

   subroutine collect_dam_break_dry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("dry_front_position_structured", test_front_structured), &
                  new_unittest("dry_profile_l1_structured", test_l1_structured), &
                  new_unittest("dry_no_negative_depth_structured", test_positivity_structured), &
                  new_unittest("dry_mass_conserved_structured", test_mass_structured), &
                  new_unittest("dry_front_position_unstructured", test_front_unstr), &
                  new_unittest("dry_profile_l1_unstructured", test_l1_unstr), &
                  new_unittest("dry_tip_celerity_is_known_deficient", test_tip_celerity), &
                  new_unittest("dry_tip_celerity_unstr_is_known_deficient", test_tip_celerity_unstr) &
                  ]
   end subroutine collect_dam_break_dry_tests

   pure function exact_threshold_x(t, h_front) result(x_thresh)
      !! Position where the EXACT profile passes through `h_front`.
      !! Inverting `h = (d/t)^2/(9g)` for `d` gives
      !! `d = 3 t sqrt(g h_front)` behind the front.
      real(wp), intent(in) :: t, h_front
      real(wp) :: x_thresh

      x_thresh = ritter_dry_front_position(H_LEFT, X_DAM, t) &
                 - 3.0_wp*t*sqrt(GRAVITY*h_front)
   end function exact_threshold_x

   ! ------------------------------------------------------------------
   ! Structured
   ! ------------------------------------------------------------------

   subroutine run_structured(state, t_final)
      !! Dry-bed dam break on the structured backend, CFL-limited stepping to
      !! `T_END`. `dt_max` clips the final step so every run lands exactly on
      !! `T_END` — the oracle is evaluated at that time and a partial step
      !! would bias the front.
      type(state_t), intent(inout) :: state
      real(wp), intent(out) :: t_final

      integer :: i, j
      real(wp) :: dx, x, t

      dx = LX/real(NX, wp)
      call state%init(NX, NY, NGHOST, dx, dx, 0.45_wp)
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x = (real(i - NGHOST, wp) - 0.5_wp)*dx
            state%barotropic%b(i, j) = 0.0_wp
            if (x < X_DAM) then
               state%barotropic%h(i, j) = H_LEFT
            else
               state%barotropic%h(i, j) = 0.0_wp
            end if
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      call solver_enter_data(state)
      t = 0.0_wp
      do while (t < T_END)
         call solver_step(state, dt_max=T_END - t, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)
      t_final = t
   end subroutine run_structured

   subroutine structured_front(state, x_front)
      !! Rightmost interior cell centre whose depth exceeds `H_FRONT_STRUCT`.
      type(state_t), intent(in) :: state
      real(wp), intent(out) :: x_front

      integer :: i, j_mid
      real(wp) :: dx

      dx = LX/real(NX, wp)
      j_mid = NGHOST + NY/2
      x_front = X_DAM
      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         if (state%barotropic%h(i, j_mid) > H_FRONT_STRUCT) then
            x_front = (real(i - NGHOST, wp) - 0.5_wp)*dx
         end if
      end do
   end subroutine structured_front

   subroutine test_front_structured(error)
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, x_num, x_ref, travelled, rel

      call run_structured(state, t_final)
      call structured_front(state, x_num)

      x_ref = exact_threshold_x(t_final, H_FRONT_STRUCT)
      travelled = ritter_dry_front_position(H_LEFT, X_DAM, t_final) - X_DAM
      rel = abs(x_num - x_ref)/travelled

      call check(error, rel < FRONT_TOL, &
                 "structured dry-front position error exceeds 1.5% of distance travelled")

      call state%destroy()
   end subroutine test_front_structured

   subroutine test_l1_structured(error)
      !! Relative L1 on depth across the whole channel, so the dry region is
      !! included and any spurious water ahead of the front counts against us.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j_mid
      real(wp) :: t_final, dx, x, h_ex, u_ex, err_sum, ref_sum, rel

      call run_structured(state, t_final)

      dx = LX/real(NX, wp)
      j_mid = NGHOST + NY/2
      err_sum = 0.0_wp
      ref_sum = 0.0_wp
      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         x = (real(i - NGHOST, wp) - 0.5_wp)*dx
         call ritter_dry_sample(H_LEFT, X_DAM, t_final, x, h_ex, u_ex)
         err_sum = err_sum + abs(state%barotropic%h(i, j_mid) - h_ex)*dx
         ref_sum = ref_sum + h_ex*dx
      end do
      rel = err_sum/max(ref_sum, tiny(1.0_wp))

      call check(error, rel < L1_TOL, &
                 "structured dry dam-break L1 depth error exceeds 5%")

      call state%destroy()
   end subroutine test_l1_structured

   subroutine test_positivity_structured(error)
      !! A vacuum front is where a Godunov scheme most easily produces
      !! negative depth. Nothing may go below zero anywhere, ever.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j
      real(wp) :: t_final, h_min

      call run_structured(state, t_final)

      h_min = huge(1.0_wp)
      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            h_min = min(h_min, state%barotropic%h(i, j))
         end do
      end do

      call check(error, h_min >= 0.0_wp, &
                 "negative depth produced at the dry front")

      call state%destroy()
   end subroutine test_positivity_structured

   subroutine test_mass_structured(error)
      !! The front advances into a closed domain, so total water mass is
      !! invariant. This is the wet/dry path's conservation check — the
      !! Phase 1 fixtures are all fully wet.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      real(wp) :: t_final, m0, m1, rel

      call run_structured(state, t_final)
      ! Re-derive the initial mass analytically: the IC is exactly H_LEFT over
      ! the cells left of the dam, so no separate pre-run measurement is
      ! needed and the comparison cannot drift with the fixture.
      m1 = compute_total_mass(state)
      m0 = H_LEFT*X_DAM*(real(NY, wp)*(LX/real(NX, wp)))
      rel = abs(m1 - m0)/m0

      call check(error, rel < MASS_TOL, &
                 "dry dam-break did not conserve water mass")

      call state%destroy()
   end subroutine test_mass_structured

   subroutine test_tip_celerity(error)
      !! CHARACTERISATION, not validation — and deliberately two-sided.
      !!
      !! The exact tip advances at `2 c0`.  After the thin-cell
      !! desingularization fix the model's leading edge (any depth above
      !! 1e-6 m) advances at ~1.74 c0 — up from the pre-fix 1.49 c0, vs a
      !! measured 1.888 c0 ceiling with all thin-cell control removed (a
      !! configuration that is unstable on beaches, see FROUDE_BACKSTOP).
      !! The residual deficit is the KP07 damping band below
      !! DESING_EPS = 1.5e-4 m plus 1st-order tip smearing.
      !!
      !! The LOWER bound is a regression guard: if the tip slows below
      !! 1.65, thin-cell mobility got worse and this fails.  The UPPER
      !! bound fires if someone closes the residual gap (e.g. a flux-side
      !! thin-tongue treatment): past 1.85, retire this guard and tighten
      !! `H_FRONT_STRUCT` / `FRONT_TOL` further.  A one-sided assertion
      !! would let such a fix land silently and leave this file
      !! permanently understating the model.
      type(error_type), allocatable, intent(out) :: error

      type(state_t) :: state
      integer :: i, j_mid
      real(wp) :: t_final, dx, x_tip, ratio, c0

      call run_structured(state, t_final)

      dx = LX/real(NX, wp)
      j_mid = NGHOST + NY/2
      x_tip = X_DAM
      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         if (state%barotropic%h(i, j_mid) > 1.0e-6_wp) then
            x_tip = (real(i - NGHOST, wp) - 0.5_wp)*dx
         end if
      end do

      c0 = sqrt(GRAVITY*H_LEFT)
      ratio = ((x_tip - X_DAM)/t_final)/c0

      call check(error, ratio > 1.65_wp, &
                 "dry-front tip celerity regressed below its post-fix 1.74 c0")
      if (allocated(error)) then
         call state%destroy()
         return
      end if
      call check(error, ratio < 1.85_wp, &
                 "dry-front tip celerity improved past 1.85 c0 — the residual "// &
                 "sub-mm deficit appears fixed; retire this guard and tighten "// &
                 "the front-position tests further")

      call state%destroy()
   end subroutine test_tip_celerity

   ! ------------------------------------------------------------------
   ! Unstructured
   ! ------------------------------------------------------------------

   subroutine build_channel(mesh)
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: i, j, nv, ntri, idx, v00, v10, v11, v01, t
      real(wp) :: hx, hy

      hx = LX/real(NQX, wp)
      hy = LY_UNSTR/real(NQY, wp)

      nv = (NQX + 1)*(NQY + 1)
      allocate (vx(nv), vy(nv))
      do j = 0, NQY
         do i = 0, NQX
            idx = j*(NQX + 1) + i + 1
            vx(idx) = real(i, wp)*hx
            vy(idx) = real(j, wp)*hy
         end do
      end do

      ntri = 2*NQX*NQY
      allocate (tv(3, ntri))
      do j = 0, NQY - 1
         do i = 0, NQX - 1
            v00 = j*(NQX + 1) + i + 1
            v10 = j*(NQX + 1) + (i + 1) + 1
            v11 = (j + 1)*(NQX + 1) + (i + 1) + 1
            v01 = (j + 1)*(NQX + 1) + i + 1
            t = 2*(j*NQX + i) + 1
            tv(:, t) = [v00, v10, v11]
            tv(:, t + 1) = [v00, v11, v01]
         end do
      end do

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
      deallocate (vx, vy, tv)
   end subroutine build_channel

   subroutine run_unstr(state, mesh, t_final)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(out) :: t_final

      integer :: k, n, n_steps
      real(wp) :: dt

      call build_channel(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp
      do k = 1, mesh%ncells_total
         if (mesh%centroid_x(k) < X_DAM) then
            state%barotropic%stage_cv(k) = H_LEFT
         else
            state%barotropic%stage_cv(k) = 0.0_wp
         end if
      end do

      ! Fixed dt from the mesh's own CFL length (min centroid-to-edge-midpoint
      ! distance), chosen so an integer number of steps lands exactly on T_END.
      dt = 0.45_wp*minval(mesh%radius(1:mesh%ncells_owned))/sqrt(GRAVITY*H_LEFT)
      n_steps = ceiling(T_END/dt)
      dt = T_END/real(n_steps, wp)

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      do n = 1, n_steps
         call solver_unstr_step(state, mesh, dt)
      end do
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      t_final = T_END
   end subroutine run_unstr

   subroutine test_front_unstr(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: k
      real(wp) :: t_final, x_num, x_ref, travelled, rel, depth

      call run_unstr(state, mesh, t_final)

      x_num = X_DAM
      do k = 1, state%ncells_owned
         depth = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)
         if (depth > H_FRONT_UNSTR) x_num = max(x_num, mesh%centroid_x(k))
      end do

      x_ref = exact_threshold_x(t_final, H_FRONT_UNSTR)
      travelled = ritter_dry_front_position(H_LEFT, X_DAM, t_final) - X_DAM
      rel = abs(x_num - x_ref)/travelled

      call check(error, rel < FRONT_TOL, &
                 "unstructured dry-front position error exceeds 1.5% of distance travelled")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_front_unstr

   subroutine test_l1_unstr(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: k
      real(wp) :: t_final, h_ex, u_ex, err_sum, ref_sum, rel, depth

      call run_unstr(state, mesh, t_final)

      err_sum = 0.0_wp
      ref_sum = 0.0_wp
      do k = 1, state%ncells_owned
         depth = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)
         call ritter_dry_sample(H_LEFT, X_DAM, t_final, mesh%centroid_x(k), h_ex, u_ex)
         err_sum = err_sum + abs(depth - h_ex)*mesh%area(k)
         ref_sum = ref_sum + h_ex*mesh%area(k)
      end do
      rel = err_sum/max(ref_sum, tiny(1.0_wp))

      call check(error, rel < L1_TOL, &
                 "unstructured dry dam-break L1 depth error exceeds 5%")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_l1_unstr

   subroutine test_tip_celerity_unstr(error)
      !! CHARACTERISATION, two-sided — the unstructured (KNP) tip deficit.
      !!
      !! Measured at 1.586 c0 (dx = 0.5 m) — and IDENTICAL to three
      !! decimals with the protect-kernel momentum clamp replaced or not,
      !! so unlike the structured path this deficit is NOT the update-side
      !! clamps: the binding constraint is the reconstruction's dry-cell
      !! handling (`centroid_pass` zeroes edge momenta below
      !! `minimum_allowed_height`) and the coarser mesh.  The LOWER bound
      !! guards regression; the UPPER bound fires when someone fixes the
      !! reconstruction — the prompt to retire this guard and tighten
      !! `H_FRONT_UNSTR` onto a deeper contour.
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      integer :: k
      real(wp) :: t_final, x_tip, ratio, c0, depth

      call run_unstr(state, mesh, t_final)

      x_tip = X_DAM
      do k = 1, state%ncells_owned
         depth = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)
         if (depth > 1.0e-6_wp) x_tip = max(x_tip, mesh%centroid_x(k))
      end do

      c0 = sqrt(GRAVITY*H_LEFT)
      ratio = ((x_tip - X_DAM)/t_final)/c0

      call check(error, ratio > 1.45_wp, &
                 "unstructured dry-front tip celerity regressed below its known 1.59 c0")
      if (allocated(error)) then
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if
      call check(error, ratio < 1.75_wp, &
                 "unstructured tip celerity improved past 1.75 c0 — the "// &
                 "reconstruction-side deficit appears fixed; retire this "// &
                 "guard and tighten H_FRONT_UNSTR onto a deeper contour")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_tip_celerity_unstr

end module test_dam_break_dry
