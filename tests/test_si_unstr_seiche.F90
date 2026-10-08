!! Closed-basin seiche test for the unstructured semi-implicit free surface.
module test_si_unstr_seiche
   !! Rectangular closed basin triangulated into two right triangles per
   !! structured quad, flat bed, still depth H.  IC: eta = A*cos(pi*x/L),
   !! U=0 (A << H, linear).  The fundamental standing-wave (seiche) period
   !! is T = 2*L/sqrt(g*H).
   !!
   !! Asserted:
   !!   * PERIOD within a few % of the analytic T at a moderate Courant
   !!     (the theta-method temporal dispersion + the coarse triangulated
   !!     antinode sampling set the bar; this checks the TPFA operator
   !!     reproduces the wave physics),
   !!   * MASS (total wetted volume sum_c area_c*h_c) conserved at round-off
   !!     — the SI edge-flux divergence telescopes over interior edges and
   !!     every wall edge carries exactly zero transport (Hf=0),
   !!   * STABILITY (non-growing amplitude) at large Courant.
   !!
   !! The area-weighted CG is the load-bearing detail: without the area
   !! weight the TPFA operator is non-symmetric on the variable-area
   !! triangulation and CG stalls / the period is wrong.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_enable_si, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_si_unstr_state, only: SI_PRECOND_NONE, SI_PRECOND_JACOBI
   use rki_timestepper, only: TIME_SEMI_IMPLICIT
   use rki_solver_unstr, only: solver_unstr_step
   use rki_mesh_unstr, only: mesh_unstr_si_quality_gate
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_unstr_seiche_tests

   real(wp), parameter :: LX = 4000.0_wp
      !! Basin length (m).
   real(wp), parameter :: LY = 800.0_wp
      !! Basin width (m); IC is y-uniform => 1-D dynamics.
   integer, parameter :: NXQ = 40
      !! Quads along x.
   integer, parameter :: NYQ = 8
      !! Quads along y.
   real(wp), parameter :: H0 = 10.0_wp
      !! Still water depth (m).
   real(wp), parameter :: AMP = 0.01_wp
      !! Seiche amplitude (m), << H0 (linear).
   real(wp), parameter :: MASS_DRIFT_REL = 1.0e-11_wp
      !! Max allowed |sum(area*h) - sum(area*h0)| / sum(area*h0).  Round-off
      !! scale by construction (edge-flux telescoping + exactly-zero wall
      !! edges); the only contributors are the CG relative tol (1e-10) and
      !! float round-off.

contains

   subroutine collect_si_unstr_seiche_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_seiche_period", test_period), &
                  new_unittest("si_unstr_seiche_stable_jacobi", test_stability_jacobi), &
                  new_unittest("si_unstr_seiche_mesh_quality", test_mesh_quality), &
                  new_unittest("si_unstr_bad_mesh_gate_rejects", test_bad_mesh_gate), &
                  new_unittest("si_unstr_flipped_pair_gate_rejects", test_flipped_pair_gate), &
                  new_unittest("si_unstr_seiche_period_converges", test_period_convergence) &
                  ]
   end subroutine collect_si_unstr_seiche_tests

   subroutine build_rect_mesh(mesh, lx, ly, nxq, nyq)
      !! Triangulate [0,lx] x [0,ly] into nxq*nyq quads, each split along the
      !! (0,0)->(1,1) diagonal into two CCW triangles.  Interior vertices are
      !! JITTERED so the triangles are acute Delaunay (an axis-aligned
      !! right-triangle split would put each circumcentre ON the shared
      !! hypotenuse => coincident circumcentres => edge_delta=0, which the
      !! TPFA quality gate rightly rejects).  Boundary vertices stay ON the
      !! box so the domain is exactly rectangular for the wall BC.  The jitter
      !! is deterministic (a fixed hash) so the test is reproducible.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp), intent(in) :: lx, ly
      integer, intent(in) :: nxq, nyq

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: nvx, nvy, nv, ntri, i, j, vid, t
      real(wp) :: hx, hy, jx, jy
      real(wp), parameter :: JIT = 0.28_wp

      nvx = nxq + 1
      nvy = nyq + 1
      nv = nvx*nvy
      ntri = 2*nxq*nyq
      hx = lx/real(nxq, wp)
      hy = ly/real(nyq, wp)

      allocate (vx(nv), vy(nv), tv(3, ntri))
      do j = 0, nyq
         do i = 0, nxq
            vid = j*nvx + i + 1
            jx = 0.0_wp
            jy = 0.0_wp
            ! Jitter interior vertices only (keep the boundary on the box).
            if (i > 0 .and. i < nxq) jx = JIT*hx*(hash01(i, j, 1) - 0.5_wp)
            if (j > 0 .and. j < nyq) jy = JIT*hy*(hash01(i, j, 2) - 0.5_wp)
            vx(vid) = real(i, wp)*hx + jx
            vy(vid) = real(j, wp)*hy + jy
         end do
      end do

      t = 0
      do j = 0, nyq - 1
         do i = 0, nxq - 1
            ! Quad corners (1-based vertex ids).
            block
               integer :: bl, br, tl, tr
               logical :: diag_bl_tr
               bl = j*nvx + i + 1
               br = j*nvx + i + 2
               tl = (j + 1)*nvx + i + 1
               tr = (j + 1)*nvx + i + 2
               ! Pick the Delaunay diagonal (the in-circle test) so the two
               ! triangles are the locally-equiangular pair — circumcentres
               ! sit strictly inside acute triangles (edge_delta > 0, small
               ! deviation angle).  A fixed diagonal on a jittered grid would
               ! occasionally give an obtuse triangle the gate rejects.
               diag_bl_tr = delaunay_diag_bl_tr( &
                            vx(bl), vy(bl), vx(br), vy(br), &
                            vx(tr), vy(tr), vx(tl), vy(tl))
               if (diag_bl_tr) then
                  t = t + 1
                  tv(:, t) = [bl, br, tr]
                  t = t + 1
                  tv(:, t) = [bl, tr, tl]
               else
                  t = t + 1
                  tv(:, t) = [bl, br, tl]
                  t = t + 1
                  tv(:, t) = [br, tr, tl]
               end if
            end block
         end do
      end do

      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_rect_mesh

   pure function hash01(i, j, s) result(r)
      !! Deterministic pseudo-random value in [0,1) from a 3-integer key.
      integer, intent(in) :: i, j, s
      real(wp) :: r
      integer(kind=8) :: h
      h = int(i, 8)*73856093_8 + int(j, 8)*19349663_8 + int(s, 8)*83492791_8
      h = iand(h, 2147483647_8)
      r = real(modulo(h, 100003_8), wp)/100003.0_wp
   end function hash01

   pure function delaunay_diag_bl_tr(xbl, ybl, xbr, ybr, xtr, ytr, xtl, ytl) &
      result(use_bl_tr)
      !! In-circle test choosing the Delaunay diagonal of the quad
      !! (bl, br, tr, tl) CCW.  Returns .true. to split along bl--tr (the
      !! diagonal that keeps both triangles Delaunay), .false. for br--tl.
      !! Uses the standard in-circle determinant on triangle (bl, br, tr):
      !! if tl is INSIDE its circumcircle the bl--tr diagonal is NOT
      !! Delaunay, so flip.
      real(wp), intent(in) :: xbl, ybl, xbr, ybr, xtr, ytr, xtl, ytl
      logical :: use_bl_tr
      real(wp) :: ax, ay, bx, by, cx, cy, incirc

      ax = xbl - xtl; ay = ybl - ytl
      bx = xbr - xtl; by = ybr - ytl
      cx = xtr - xtl; cy = ytr - ytl
      incirc = (ax*ax + ay*ay)*(bx*cy - cx*by) &
               - (bx*bx + by*by)*(ax*cy - cx*ay) &
               + (cx*cx + cy*cy)*(ax*by - bx*ay)
      ! (bl,br,tr) CCW => incirc > 0 means tl is inside the circumcircle =>
      ! flip to br--tl; else keep bl--tr.
      use_bl_tr = (incirc <= 0.0_wp)
   end function delaunay_diag_bl_tr

   subroutine setup_seiche(state, mesh, dt, precond)
      !! Flat bed, still depth H0 + seiche IC at the circumcentres, SI mode.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt
      integer, intent(in) :: precond

      integer :: k
      real(wp) :: xcc, eta_ic
      logical :: ok

      call state_unstr_init(state, mesh)
      state%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state, mesh, precond, 1.0e-10_wp, dt, ok)
      if (.not. ok) then
         call global_logger%error("setup_seiche: SI quality gate failed on a "// &
                                  "clean rectangle mesh (unexpected)")
         error stop "setup_seiche quality gate"
      end if

      ! Flat bed b=0; stage = h = H0 + A*cos(pi*x/L) at the circumcentre.
      do k = 1, mesh%ncells_total
         xcc = mesh%circumcentre_x(k)
         eta_ic = AMP*cos(4.0_wp*atan(1.0_wp)*xcc/LX)
         state%barotropic%bed_cv(k) = 0.0_wp
         state%barotropic%stage_cv(k) = H0 + eta_ic
         state%barotropic%height_cv(k) = H0 + eta_ic
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do
   end subroutine setup_seiche

   function mass_sum(state, mesh) result(s)
      !! Total wetted volume sum_c area_c * h_c (host copy).
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp) :: s
      integer :: k
      s = 0.0_wp
      do k = 1, mesh%ncells_owned
         s = s + mesh%area(k)*max(state%barotropic%stage_cv(k) &
                                  - state%barotropic%bed_cv(k), 0.0_wp)
      end do
   end function mass_sum

   subroutine run_seiche(state, mesh, n_steps, dt, eta_series, xprobe, &
                         hsum0, hdrift, iters_mean)
      !! Advance n_steps, sampling eta at a probe cell near the west wall.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt, xprobe
      real(wp), intent(out) :: eta_series(n_steps)
      real(wp), intent(out) :: hsum0, hdrift, iters_mean

      integer :: n, k, kprobe, iters_total
      real(wp) :: best

      ! Probe = the owned cell whose circumcentre is closest to (xprobe, LY/2).
      kprobe = 1
      best = huge(1.0_wp)
      do k = 1, mesh%ncells_owned
         if (abs(mesh%circumcentre_x(k) - xprobe) &
             + abs(mesh%circumcentre_y(k) - 0.5_wp*LY) < best) then
            best = abs(mesh%circumcentre_x(k) - xprobe) &
                   + abs(mesh%circumcentre_y(k) - 0.5_wp*LY)
            kprobe = k
         end if
      end do

      iters_total = 0
      hdrift = 0.0_wp

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      hsum0 = mass_sum(state, mesh)

      do n = 1, n_steps
         call solver_unstr_step(state, mesh, dt)
         iters_total = iters_total + state%si%cg_iters
         !$omp target update from(state%barotropic%stage_cv)
         eta_series(n) = state%barotropic%stage_cv(kprobe) - H0
         hdrift = max(hdrift, abs(mass_sum(state, mesh) - hsum0))
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      iters_mean = real(iters_total, wp)/real(n_steps, wp)
   end subroutine run_seiche

   function measure_period(eta_series, n, dt) result(t_meas)
      !! Period from up-going zero crossings of (eta - mean).
      integer, intent(in) :: n
      real(wp), intent(in) :: eta_series(n)
      real(wp), intent(in) :: dt
      real(wp) :: t_meas

      integer :: k, n_cross, first_cross, last_cross
      real(wp) :: mean_eta, prev, cur, frac, t_cross, t_first, t_last

      mean_eta = sum(eta_series)/real(n, wp)
      n_cross = 0; first_cross = 0; last_cross = 0
      t_first = 0.0_wp; t_last = 0.0_wp
      do k = 2, n
         prev = eta_series(k - 1) - mean_eta
         cur = eta_series(k) - mean_eta
         if (prev <= 0.0_wp .and. cur > 0.0_wp) then
            frac = -prev/(cur - prev)
            t_cross = (real(k - 1, wp) + frac)*dt
            n_cross = n_cross + 1
            if (n_cross == 1) then
               first_cross = n_cross; t_first = t_cross
            end if
            last_cross = n_cross; t_last = t_cross
         end if
      end do
      if (last_cross > first_cross) then
         t_meas = (t_last - t_first)/real(last_cross - first_cross, wp)
      else
         t_meas = -1.0_wp
      end if
   end function measure_period

   subroutine test_period(error)
      !! Period within 4% of analytic T at Courant ~5, plain CG.
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: t_analytic, dt, courant, t_meas, dx_eff
      real(wp) :: hsum0, hdrift, iters_mean, err_rel
      real(wp), allocatable :: eta_series(:)
      integer :: n_steps, n_periods

      call build_rect_mesh(mesh, LX, LY, NXQ, NYQ)
      dx_eff = LX/real(NXQ, wp)
      t_analytic = 2.0_wp*LX/sqrt(GRAVITY*H0)
      courant = 5.0_wp
      dt = courant*dx_eff/sqrt(GRAVITY*H0)
      n_periods = 5
      n_steps = ceiling(real(n_periods, wp)*t_analytic/dt)
      allocate (eta_series(n_steps))

      call setup_seiche(state, mesh, dt, SI_PRECOND_NONE)
      call run_seiche(state, mesh, n_steps, dt, eta_series, 0.05_wp*LX, &
                      hsum0, hdrift, iters_mean)
      t_meas = measure_period(eta_series, n_steps, dt)

      call global_logger%info("[si-unstr seiche C5] T_analytic="//rstr(t_analytic)// &
                              " T_meas="//rstr(t_meas)//" dt="//rstr(dt)// &
                              " CG_iters/step="//rstr(iters_mean)// &
                              " mass_drift="//rstr(hdrift))

      call check(error, t_meas > 0.0_wp, "si-unstr seiche: no crossings")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      err_rel = abs(t_meas - t_analytic)/t_analytic
      call check(error, err_rel < 0.04_wp, &
                 "si-unstr seiche period error "//rstr(err_rel*100.0_wp)//"% > 4%")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, hdrift < MASS_DRIFT_REL*hsum0, &
                 "si-unstr seiche: mass drift "//rstr(hdrift)//" > bound "// &
                 rstr(MASS_DRIFT_REL*hsum0))
      call cleanup(state, mesh)
   end subroutine test_period

   function seiche_period_error(nxq, nyq, courant) result(err_rel)
      !! Relative error of the measured seiche period vs the analytic
      !! `2*LX/sqrt(g*H0)` on an `nxq x nyq` triangulated rectangle at a
      !! fixed advective-free Courant number `courant`.
      integer, intent(in) :: nxq, nyq
      real(wp), intent(in) :: courant
      real(wp) :: err_rel

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: t_analytic, dt, t_meas, dx_eff
      real(wp) :: hsum0, hdrift, iters_mean
      real(wp), allocatable :: eta_series(:)
      integer :: n_steps

      call build_rect_mesh(mesh, LX, LY, nxq, nyq)
      dx_eff = LX/real(nxq, wp)
      t_analytic = 2.0_wp*LX/sqrt(GRAVITY*H0)
      dt = courant*dx_eff/sqrt(GRAVITY*H0)
      n_steps = ceiling(5.0_wp*t_analytic/dt)
      allocate (eta_series(n_steps))

      call setup_seiche(state, mesh, dt, SI_PRECOND_NONE)
      call run_seiche(state, mesh, n_steps, dt, eta_series, 0.05_wp*LX, &
                      hsum0, hdrift, iters_mean)
      t_meas = measure_period(eta_series, n_steps, dt)
      call cleanup(state, mesh)

      if (t_meas > 0.0_wp) then
         err_rel = abs(t_meas - t_analytic)/t_analytic
      else
         err_rel = huge(1.0_wp)
      end if
   end function seiche_period_error

   subroutine test_period_convergence(error)
      !! Resolution study against the analytic period at a fixed Courant
      !! number: the period error must converge at ~2nd order as the
      !! triangulation is refined.  This ties the TPFA operator to the analytic solution
      !! independently of any other discretisation.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: e_coarse, e_mid, e_fine

      e_coarse = seiche_period_error(20, 4, 2.0_wp)
      e_mid = seiche_period_error(40, 8, 2.0_wp)
      e_fine = seiche_period_error(80, 16, 2.0_wp)

      call global_logger%info("[si-unstr seiche convergence] period error % "// &
                              "20x4="//rstr(e_coarse*100.0_wp)// &
                              " 40x8="//rstr(e_mid*100.0_wp)// &
                              " 80x16="//rstr(e_fine*100.0_wp))
      ! Second-order spatial convergence is a factor ~4 per halving of dx
      ! (measured 3.9x and 3.9x); require > 3x so a first-order or stalled
      ! operator fails.
      call check(error, e_coarse > 3.0_wp*e_mid .and. e_mid > 3.0_wp*e_fine, &
                 "si-unstr seiche period error is not converging at ~2nd order")
      if (allocated(error)) return
      call check(error, e_fine < 0.002_wp, &
                 "si-unstr seiche fine-mesh period error "// &
                 rstr(e_fine*100.0_wp)//"% > 0.2%")
   end subroutine test_period_convergence

   subroutine test_stability_jacobi(error)
      !! Non-growing amplitude at Courant 20 with the Jacobi preconditioner
      !! (which does real work on the variable-area triangulation).
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: t_analytic, dt, courant, dx_eff
      real(wp) :: hsum0, hdrift, iters_mean
      real(wp), allocatable :: eta_series(:)
      integer :: n_steps, n_periods, steps_per_period, k, k0
      real(wp) :: amp_first, amp_last

      call build_rect_mesh(mesh, LX, LY, NXQ, NYQ)
      dx_eff = LX/real(NXQ, wp)
      t_analytic = 2.0_wp*LX/sqrt(GRAVITY*H0)
      courant = 20.0_wp
      dt = courant*dx_eff/sqrt(GRAVITY*H0)
      n_periods = 6
      steps_per_period = max(1, nint(t_analytic/dt))
      n_steps = n_periods*steps_per_period
      allocate (eta_series(n_steps))

      call setup_seiche(state, mesh, dt, SI_PRECOND_JACOBI)
      call run_seiche(state, mesh, n_steps, dt, eta_series, 0.05_wp*LX, &
                      hsum0, hdrift, iters_mean)

      call global_logger%info("[si-unstr seiche C20 jacobi] dt="//rstr(dt)// &
                              " CG_iters/step="//rstr(iters_mean)// &
                              " mass_drift="//rstr(hdrift))

      amp_first = 0.0_wp
      do k = 1, steps_per_period
         amp_first = max(amp_first, abs(eta_series(k)))
      end do
      amp_last = 0.0_wp
      k0 = (n_periods - 1)*steps_per_period + 1
      do k = k0, n_steps
         amp_last = max(amp_last, abs(eta_series(k)))
      end do

      call check(error, amp_last <= amp_first*1.02_wp, &
                 "si-unstr seiche amplitude GREW at C20 (unstable): first="// &
                 rstr(amp_first)//" last="//rstr(amp_last))
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, amp_last > 0.05_wp*amp_first, &
                 "si-unstr seiche amplitude COLLAPSED at C20 (over-damped)")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, hdrift < MASS_DRIFT_REL*hsum0, &
                 "si-unstr seiche C20: mass drift "//rstr(hdrift)//" > bound")
      call cleanup(state, mesh)
   end subroutine test_stability_jacobi

   subroutine test_mesh_quality(error)
      !! The clean structured-Delaunay rectangle passes the SI quality gate
      !! with a tiny max deviation angle (a correctness assertion on the
      !! circumcentre-line orthogonality).
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh

      call build_rect_mesh(mesh, LX, LY, NXQ, NYQ)
      call global_logger%info("[si-unstr mesh] max_dev_angle="// &
                              rstr(mesh%max_dev_angle)//" deg, min_edge_delta="// &
                              rstr(mesh%min_edge_delta)//" m")
      call check(error, mesh%max_dev_angle < 1.0_wp, &
                 "si-unstr: clean Delaunay rectangle has max_dev_angle "// &
                 rstr(mesh%max_dev_angle)//" deg (expected < 1 deg)")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if
      call check(error, mesh%min_edge_delta > 0.0_wp, &
                 "si-unstr: non-positive min_edge_delta on a clean mesh")
      call mesh_unstr_destroy(mesh)
   end subroutine test_mesh_quality

   subroutine test_bad_mesh_gate(error)
      !! A deliberately BAD mesh — the unjittered axis-aligned right-triangle
      !! split — puts each circumcentre ON the shared hypotenuse, so two
      !! triangles sharing a hypotenuse have coincident circumcentres
      !! (edge_delta = 0).  The fail-loud TPFA quality gate must REJECT it
      !! (ok=.false.), the house-style guard against silently solving a
      !! subtly-wrong (or singular) eta on a non-Delaunay mesh.
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: nvx, i, j, vid, t, nxq, nyq, nv, ntri
      real(wp) :: hx, hy
      logical :: ok

      ! Small unjittered right-triangle grid: circumcentre on the hypotenuse.
      nxq = 4; nyq = 4
      nvx = nxq + 1
      nv = nvx*(nyq + 1)
      ntri = 2*nxq*nyq
      hx = 100.0_wp; hy = 100.0_wp
      allocate (vx(nv), vy(nv), tv(3, ntri))
      do j = 0, nyq
         do i = 0, nxq
            vid = j*nvx + i + 1
            vx(vid) = real(i, wp)*hx
            vy(vid) = real(j, wp)*hy
         end do
      end do
      t = 0
      do j = 0, nyq - 1
         do i = 0, nxq - 1
            block
               integer :: bl, br, tl, tr
               bl = j*nvx + i + 1
               br = j*nvx + i + 2
               tl = (j + 1)*nvx + i + 1
               tr = (j + 1)*nvx + i + 2
               t = t + 1
               tv(:, t) = [bl, br, tr]
               t = t + 1
               tv(:, t) = [bl, tr, tl]
            end block
         end do
      end do
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)

      call mesh_unstr_si_quality_gate(mesh, ok)
      call global_logger%info("[si-unstr bad mesh] min_edge_delta="// &
                              rstr(mesh%min_edge_delta)//" m, gate ok="// &
                              merge("T", "F", ok))
      call check(error,.not. ok, &
                 "si-unstr: the quality gate ACCEPTED a degenerate "// &
                 "right-triangle mesh (min_edge_delta="//rstr(mesh%min_edge_delta)// &
                 " m) — it must fail loud")
      call mesh_unstr_destroy(mesh)
   end subroutine test_bad_mesh_gate

   subroutine test_flipped_pair_gate(error)
      !! A LOCALLY NON-DELAUNAY (flipped) pair must be rejected with a
      !! STRICTLY NEGATIVE signed projection — the case the deviation angle
      !! can never flag (both circumcentres always sit on the shared edge's
      !! perpendicular bisector) and the coincident-circumcentre (=0) test
      !! above does not exercise.  Construction: shared edge (0,0)-(100,0);
      !! T1 apex far below at (50,-300) => a huge circumcircle whose top
      !! reaches y ~ +8.3; T2 apex at (50,5) lies INSIDE it (non-Delaunay).
      !! Both circumcentres land BELOW the shared edge (cc1 = (50,-145.8),
      !! cc2 = (50,-247.5)), so the neighbour circumcentre projects
      !! NEGATIVELY on T1's outward (+y) normal: -101.7 m.  The test first
      !! VERIFIES the construction (min_edge_delta strictly < 0), then
      !! asserts the gate rejects.
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)
      logical :: ok

      vx = [0.0_wp, 100.0_wp, 50.0_wp, 50.0_wp]
      vy = [0.0_wp, 0.0_wp, -300.0_wp, 5.0_wp]
      ! T1 = (V1, V3, V2) CCW (apex below), T2 = (V1, V2, V4) CCW (apex
      ! above); shared edge V1-V2.
      tv(:, 1) = [1, 3, 2]
      tv(:, 2) = [1, 2, 4]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)

      call mesh_unstr_si_quality_gate(mesh, ok)
      call global_logger%info("[si-unstr flipped pair] min_edge_delta="// &
                              rstr(mesh%min_edge_delta)//" m (analytic -101.7), "// &
                              "gate ok="//merge("T", "F", ok))
      ! Verify the construction actually produced a strictly negative
      ! projection (not the coincident =0 degenerate case).
      call check(error, mesh%min_edge_delta < -1.0_wp, &
                 "flipped pair: construction failed — min_edge_delta="// &
                 rstr(mesh%min_edge_delta)//" is not strictly negative")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if
      call check(error,.not. ok, &
                 "flipped pair: the quality gate ACCEPTED a locally "// &
                 "non-Delaunay pair (signed projection "// &
                 rstr(mesh%min_edge_delta)//" m) — it must fail loud")
      call mesh_unstr_destroy(mesh)
   end subroutine test_flipped_pair_gate

   subroutine cleanup(state, mesh)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine cleanup

   pure function rstr(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(ES12.5)') x
      s = trim(adjustl(buf))
   end function rstr

end module test_si_unstr_seiche
