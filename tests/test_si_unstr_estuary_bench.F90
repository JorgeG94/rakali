!! Acceptance bench: unstructured estuary — SI vs explicit KNP head-to-head.
module test_si_unstr_estuary_bench
   !! The P3.3 acceptance case: a FUNNEL-mesh estuary (width tapering
   !! 2 km -> 600 m over 12 km, jittered Delaunay-diagonal triangulation),
   !! deepened flats (bed -8 m at the mouth -> -2 m at the head — all-wet,
   !! the bench isolates the solver tier, not wet/dry), wall + clamped
   !! tidal west edge (M2, 0.5 m), nz=4 sigma layers with S/T active, and
   !! Manning friction — run THREE ways on the identical mesh + physics:
   !!
   !!   A. explicit KNP multilayer (RK2, gravity-wave-CFL dt — the
   !!      production explicit path),
   !!   B. SI plain CG, rtol = 1e-10 (the reference SI tier),
   !!   C. SI Jacobi + rtol = 1e-8 (the tuned tier from the
   !!      preconditioner study).
   !!
   !! ASSERTED (the trust story): tidal amplitude at three probes
   !! (mouth / mid / head) agrees between the explicit and both SI runs
   !! within the theta-damping budget, and the phase (peak time) agrees
   !! within a small fraction of the period.  LOGGED (the cost story):
   !! wall-clock per run + steps + CG iters/step — the honest speedup is
   !! computed OUTSIDE ctest by the median-of-3 clean-run protocol (this
   !! binary is deterministic; run it three times and take the medians).
   !! NOTE this is the plain/Jacobi CG tier — the BASELINE record: expect
   !! the ~4-8x band the plan budgeted for triangles, NOT the structured
   !! MG 25x.  The shipped agglomeration MG (`precond="mg"`,
   !! rki_si_unstr_mg) is the plateau-breaker for this tier (1526 -> 16
   !! iterations on the 480-tri funnel; see test_si_unstr_mg).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, TWO_PI
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_enable_si, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_si_unstr_state, only: SI_PRECOND_NONE, SI_PRECOND_JACOBI
   use rki_boundary_unstr, only: BC_WALL, BC_CLAMPED, BC_TIDAL
   use rki_timestepper, only: TIME_SEMI_IMPLICIT
   use rki_solver_unstr, only: solver_unstr_step, solver_unstr_advance
   use rki_si_unstr_step, only: SI_UNSTR_THETA
   use rki_ml_dynamics_unstr, only: ml_init_layers_unstr
   use rki_ml_tracers_unstr, only: ml_seed_standard_tracers_unstr
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_unstr_estuary_bench_tests

   ! Funnel geometry + physics.
   real(wp), parameter :: LX = 12000.0_wp
      !! Estuary length (m).
   real(wp), parameter :: W_MOUTH = 2000.0_wp
      !! Width at the mouth (m).
   real(wp), parameter :: W_HEAD = 600.0_wp
      !! Width at the head (m).
   integer, parameter :: NXQ = 40, NYQ = 6
      !! CI size (480 tris, ~30 s gfortran / ~50 s GPU for the 3 runs).
      !! GPU SCALING STORY (V100 sweep, 2026-07-03): SMALL meshes are
      !! LAUNCH-LATENCY-BOUND — neither solver saturates the GPU below
      !! ~10k tris (explicit throughput grows 0.38 -> 3.73 Mcells/s over
      !! the 480 -> 86,400-tri sweep), which is why this CI size measures
      !! only 1.5x (plain) / 3.1x (jacobi).  At FIXED dx the CG iteration
      !! count PLATEAUS (~2500/step, tidal-wavelength-bounded, and the
      !! SAME count on CPU and GPU — conditioning, not a GPU effect), so
      !! the speedup GROWS with scale: at 86,400 tris SI jacobi/rtol=1e-8
      !! measures 9.2x vs the explicit KNP (median-of-3).  plain CG at
      !! rtol=1e-10 hits the 4000-iteration cap on funnel meshes >= 7680
      !! tris, so the practical envelope for THIS no-MG tier is
      !! jacobi/1e-8 — or switch to the shipped `precond="mg"`
      !! (agglomeration multigrid), which collapses the iteration count.
      !! Follow-up: async launch batching is the small-mesh lever.
   integer, parameter :: NZ = 4
   real(wp), parameter :: DEPTH_MOUTH = 8.0_wp
      !! Still depth at the mouth (m) — "deepened flats", all-wet.
   real(wp), parameter :: DEPTH_HEAD = 2.0_wp
      !! Still depth at the head (m).
   real(wp), parameter :: TIDE_AMP = 0.5_wp
   real(wp), parameter :: TIDE_PERIOD = 44712.0_wp
      !! M2 (s).
   real(wp), parameter :: MANNING = 0.02_wp
   real(wp), parameter :: SI_COURANT = 20.0_wp
      !! Barotropic Courant of the SI runs (dt ~ minutes).
   integer, parameter :: WEST_GROUP = 1
   integer, parameter :: WALL_GROUP = 2
   real(wp), parameter :: AMP_REL_TOL = 0.12_wp
      !! Probe-amplitude agreement gate SI vs explicit.  Budget: the
      !! theta=0.6 SI damping at Courant ~20 plus the first-order upwind
      !! explicit dissipation differ by a few % on a damped co-oscillating
      !! estuary; 12% is the honest envelope (measured values logged).
   real(wp), parameter :: PHASE_TOL_FRAC = 0.05_wp
      !! Peak-time agreement gate as a fraction of the M2 period.

contains

   subroutine collect_si_unstr_estuary_bench_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_estuary_bench", test_estuary_bench) &
                  ]
   end subroutine collect_si_unstr_estuary_bench_tests

   subroutine build_funnel_mesh(mesh)
      !! Funnel triangulation: mapped quads with width W(x) tapering
      !! linearly mouth -> head, gentle jitter + per-quad Delaunay
      !! diagonals (the acute construction; the hardened SIGNED-delta
      !! quality gate verifies the taper kept every pair locally Delaunay).
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :), bc_tags(:, :)
      integer :: nvx, nv, ntri, i, j, vid, t
      real(wp) :: hx, wx, jx, jy, yfrac
      real(wp), parameter :: JIT = 0.18_wp

      nvx = NXQ + 1
      nv = nvx*(NYQ + 1)
      ntri = 2*NXQ*NYQ
      hx = LX/real(NXQ, wp)

      allocate (vx(nv), vy(nv), tv(3, ntri))
      do j = 0, NYQ
         do i = 0, NXQ
            vid = j*nvx + i + 1
            wx = W_MOUTH + (W_HEAD - W_MOUTH)*(real(i, wp)/real(NXQ, wp))
            yfrac = real(j, wp)/real(NYQ, wp)
            jx = 0.0_wp
            jy = 0.0_wp
            if (i > 0 .and. i < NXQ) jx = JIT*hx*(hash01(i, j, 1) - 0.5_wp)
            if (j > 0 .and. j < NYQ) jy = JIT*(wx/real(NYQ, wp)) &
                                          *(hash01(i, j, 2) - 0.5_wp)
            vx(vid) = real(i, wp)*hx + jx
            ! Centre the funnel on y = 0 (banks at +-W(x)/2).
            vy(vid) = (yfrac - 0.5_wp)*wx + jy
         end do
      end do

      t = 0
      do j = 0, NYQ - 1
         do i = 0, NXQ - 1
            block
               integer :: bl, br, tl, tr
               logical :: diag_bl_tr
               bl = j*nvx + i + 1
               br = j*nvx + i + 2
               tl = (j + 1)*nvx + i + 1
               tr = (j + 1)*nvx + i + 2
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
      allocate (bc_tags(3, ntri), source=0)
      do t = 1, mesh%ncells_owned
         do i = 1, 3
            if (mesh%neighbours(i, t) < 0) then
               if (mesh%edge_mid_x(i, t) < 0.25_wp*hx) then
                  bc_tags(i, t) = WEST_GROUP
               else
                  bc_tags(i, t) = WALL_GROUP
               end if
            end if
         end do
      end do
      call mesh_unstr_destroy(mesh)
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv, bc_tags=bc_tags)
   end subroutine build_funnel_mesh

   pure function hash01(i, j, s) result(r)
      integer, intent(in) :: i, j, s
      real(wp) :: r
      integer(kind=8) :: h
      h = int(i, 8)*73856093_8 + int(j, 8)*19349663_8 + int(s, 8)*83492791_8
      h = iand(h, 2147483647_8)
      r = real(modulo(h, 100003_8), wp)/100003.0_wp
   end function hash01

   pure function delaunay_diag_bl_tr(xbl, ybl, xbr, ybr, xtr, ytr, xtl, ytl) &
      result(use_bl_tr)
      real(wp), intent(in) :: xbl, ybl, xbr, ybr, xtr, ytr, xtl, ytl
      logical :: use_bl_tr
      real(wp) :: ax, ay, bx, by, cx, cy, incirc
      ax = xbl - xtl; ay = ybl - ytl
      bx = xbr - xtl; by = ybr - ytl
      cx = xtr - xtl; cy = ytr - ytl
      incirc = (ax*ax + ay*ay)*(bx*cy - cx*by) &
               - (bx*bx + by*by)*(ax*cy - cx*ay) &
               + (cx*cx + cy*cy)*(ax*by - bx*ay)
      use_bl_tr = (incirc <= 0.0_wp)
   end function delaunay_diag_bl_tr

   subroutine setup_common(state, mesh, use_si, precond, rtol, dt_si)
      !! Shared physics: funnel bathymetry, still-water IC, S/T, Manning,
      !! clamped tidal west edge.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      logical, intent(in) :: use_si
      integer, intent(in) :: precond
      real(wp), intent(in) :: rtol, dt_si

      integer :: k
      logical :: ok
      real(wp) :: xfrac

      call state_unstr_init(state, mesh, nz_ml=NZ)
      state%barotropic%manning_n = MANNING
      state%multilayer%initial_salinity = 35.0_wp
      state%multilayer%initial_temperature = 15.0_wp
      state%bc_type(WEST_GROUP) = merge(BC_CLAMPED, BC_TIDAL, use_si)
      state%bc_type(WALL_GROUP) = BC_WALL
      if (use_si) then
         state%time_scheme = TIME_SEMI_IMPLICIT
         call state_unstr_enable_si(state, mesh, precond, rtol, dt_si, ok)
         if (.not. ok) error stop "estuary_bench: quality gate failed"
      end if

      do k = 1, mesh%ncells_total
         xfrac = mesh%circumcentre_x(k)/LX
         state%barotropic%bed_cv(k) = -(DEPTH_MOUTH &
                                        + (DEPTH_HEAD - DEPTH_MOUTH)*xfrac)
         state%barotropic%stage_cv(k) = 0.0_wp
         state%barotropic%height_cv(k) = -state%barotropic%bed_cv(k)
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do
   end subroutine setup_common

   subroutine find_probes(mesh, kp)
      !! Probe cells nearest (0.1, 0.5, 0.9)*LX on the channel axis y=0.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(out) :: kp(3)

      real(wp), parameter :: XF(3) = [0.1_wp, 0.5_wp, 0.9_wp]
      integer :: p, k
      real(wp) :: best, d

      do p = 1, 3
         kp(p) = 1
         best = huge(1.0_wp)
         do k = 1, mesh%ncells_owned
            d = abs(mesh%circumcentre_x(k) - XF(p)*LX) &
                + abs(mesh%circumcentre_y(k))
            if (d < best) then
               best = d
               kp(p) = k
            end if
         end do
      end do
   end subroutine find_probes

   subroutine run_case(mesh, use_si, precond, rtol, dt_si, n_samp, kp, &
                       eta_ser, wall_s, iters_mean, n_steps_out)
      !! One full bench run: 1.5 tidal cycles, sampling the probes on the
      !! SI-step cadence (the explicit run advances to the same sample
      !! times).  Returns wall-clock (s) around the time loop only.
      type(mesh_unstr_t), intent(inout) :: mesh
      logical, intent(in) :: use_si
      integer, intent(in) :: precond
      real(wp), intent(in) :: rtol, dt_si
      integer, intent(in) :: n_samp
      integer, intent(in) :: kp(3)
      real(wp), intent(out) :: eta_ser(n_samp, 3)
      real(wp), intent(out) :: wall_s, iters_mean
      integer, intent(out) :: n_steps_out

      type(unstr_state_t) :: state
      integer :: n, p
      integer(kind=8) :: c0, c1, crate
      real(wp) :: t_bc, eta_bc, om
      integer :: iters_total

      om = TWO_PI/TIDE_PERIOD
      call setup_common(state, mesh, use_si, precond, rtol, dt_si)
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_init_layers_unstr(state)
      call ml_seed_standard_tracers_unstr(state)

      iters_total = 0
      call system_clock(c0, crate)
      do n = 1, n_samp
         if (use_si) then
            t_bc = state%time + SI_UNSTR_THETA*dt_si
            eta_bc = TIDE_AMP*sin(om*t_bc)
            state%bc_eta_ext(WEST_GROUP) = eta_bc
            !$omp target update to(state%bc_eta_ext)
            call solver_unstr_step(state, mesh, dt_si)
            iters_total = iters_total + state%si%cg_iters
         else
            ! Explicit KNP: the Flather/tidal ghost reads bc_eta_ext each
            ! flux pass; advance owns the CFL dt.
            eta_bc = TIDE_AMP*sin(om*(state%time + 0.5_wp*dt_si))
            state%bc_eta_ext(WEST_GROUP) = eta_bc
            !$omp target update to(state%bc_eta_ext)
            call solver_unstr_advance(state, mesh, real(n, wp)*dt_si)
         end if
         !$omp target update from(state%barotropic%stage_cv)
         do p = 1, 3
            eta_ser(n, p) = state%barotropic%stage_cv(kp(p))
         end do
      end do
      call system_clock(c1)
      wall_s = real(c1 - c0, wp)/real(crate, wp)
      n_steps_out = state%step_count
      iters_mean = real(iters_total, wp)/real(max(n_samp, 1), wp)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
   end subroutine run_case

   subroutine amp_phase(eta, n_samp, n_lastcycle, dt, amp, t_peak)
      !! Amplitude ((max-min)/2) and peak time over the LAST full cycle.
      integer, intent(in) :: n_samp, n_lastcycle
      real(wp), intent(in) :: eta(n_samp)
      real(wp), intent(in) :: dt
      real(wp), intent(out) :: amp, t_peak

      integer :: n, n0
      real(wp) :: emax, emin

      n0 = n_samp - n_lastcycle + 1
      emax = -huge(1.0_wp)
      emin = huge(1.0_wp)
      t_peak = 0.0_wp
      do n = n0, n_samp
         if (eta(n) > emax) then
            emax = eta(n)
            t_peak = real(n, wp)*dt
         end if
         emin = min(emin, eta(n))
      end do
      amp = 0.5_wp*(emax - emin)
   end subroutine amp_phase

   subroutine test_estuary_bench(error)
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      real(wp) :: dx_eff, dt_si, om
      integer :: kp(3), n_samp, n_cyc_samp, p
      real(wp), allocatable :: eta_exp(:, :), eta_sip(:, :), eta_sij(:, :)
      real(wp) :: wall_exp, wall_sip, wall_sij, im_p, im_j, im_e
      integer :: steps_exp, steps_sip, steps_sij
      real(wp) :: amp_e(3), tp_e(3), amp_p(3), tp_p(3), amp_j(3), tp_j(3)
      real(wp) :: arel, prel
      character(len=1), parameter :: nl = char(10)

      call build_funnel_mesh(mesh)
      dx_eff = LX/real(NXQ, wp)
      dt_si = SI_COURANT*dx_eff/sqrt(GRAVITY*DEPTH_MOUTH)
      om = TWO_PI/TIDE_PERIOD
      n_cyc_samp = nint(TIDE_PERIOD/dt_si)
      n_samp = nint(1.5_wp*TIDE_PERIOD/dt_si)
      allocate (eta_exp(n_samp, 3), eta_sip(n_samp, 3), eta_sij(n_samp, 3))

      call find_probes(mesh, kp)

      ! A. explicit KNP (RK2, CFL dt).
      call run_case(mesh, .false., SI_PRECOND_NONE, 1.0e-10_wp, dt_si, &
                    n_samp, kp, eta_exp, wall_exp, im_e, steps_exp)
      ! B. SI plain CG, rtol 1e-10.
      call run_case(mesh, .true., SI_PRECOND_NONE, 1.0e-10_wp, dt_si, &
                    n_samp, kp, eta_sip, wall_sip, im_p, steps_sip)
      ! C. SI Jacobi, rtol 1e-8 (the tuned tier).
      call run_case(mesh, .true., SI_PRECOND_JACOBI, 1.0e-8_wp, dt_si, &
                    n_samp, kp, eta_sij, wall_sij, im_j, steps_sij)

      do p = 1, 3
         call amp_phase(eta_exp(:, p), n_samp, n_cyc_samp, dt_si, amp_e(p), tp_e(p))
         call amp_phase(eta_sip(:, p), n_samp, n_cyc_samp, dt_si, amp_p(p), tp_p(p))
         call amp_phase(eta_sij(:, p), n_samp, n_cyc_samp, dt_si, amp_j(p), tp_j(p))
      end do

      call global_logger%info("[si-unstr estuary BENCH] ncell="// &
                              rstr(real(mesh%ncells_owned, wp))// &
                              " nz="//rstr(real(NZ, wp))// &
                              " dt_si="//rstr(dt_si)//" s"//nl// &
                              "  explicit KNP : wall="//rstr(wall_exp)//" s, steps="// &
                              rstr(real(steps_exp, wp))//nl// &
                              "  SI plain 1e-10: wall="//rstr(wall_sip)//" s, steps="// &
                              rstr(real(steps_sip, wp))//", CG/step="//rstr(im_p)//nl// &
                              "  SI jacobi 1e-8: wall="//rstr(wall_sij)//" s, steps="// &
                              rstr(real(steps_sij, wp))//", CG/step="//rstr(im_j)//nl// &
                              "  speedup(this run): plain "// &
                              rstr(wall_exp/max(wall_sip, 1.0e-9_wp))// &
                              "x, jacobi "//rstr(wall_exp/max(wall_sij, 1.0e-9_wp))//"x")
      do p = 1, 3
         call global_logger%info("[si-unstr estuary probe"//rstr(real(p, wp))// &
                                 "] amp exp/plain/jacobi = "//rstr(amp_e(p))//" / "// &
                                 rstr(amp_p(p))//" / "//rstr(amp_j(p))// &
                                 "; peak-t = "//rstr(tp_e(p))//" / "//rstr(tp_p(p))// &
                                 " / "//rstr(tp_j(p)))
      end do

      ! Physics agreement gates (both SI tiers vs explicit, all probes).
      do p = 1, 3
         arel = abs(amp_p(p) - amp_e(p))/max(amp_e(p), 1.0e-12_wp)
         call check(error, arel < AMP_REL_TOL, &
                    "estuary: SI-plain amplitude off by "//rstr(arel*100.0_wp)// &
                    "% at probe "//rstr(real(p, wp)))
         if (allocated(error)) exit
         arel = abs(amp_j(p) - amp_e(p))/max(amp_e(p), 1.0e-12_wp)
         call check(error, arel < AMP_REL_TOL, &
                    "estuary: SI-jacobi amplitude off by "//rstr(arel*100.0_wp)// &
                    "% at probe "//rstr(real(p, wp)))
         if (allocated(error)) exit
         prel = abs(tp_p(p) - tp_e(p))/TIDE_PERIOD
         call check(error, prel < PHASE_TOL_FRAC, &
                    "estuary: SI-plain peak time off by "//rstr(prel*100.0_wp)// &
                    "% of T at probe "//rstr(real(p, wp)))
         if (allocated(error)) exit
         prel = abs(tp_j(p) - tp_e(p))/TIDE_PERIOD
         call check(error, prel < PHASE_TOL_FRAC, &
                    "estuary: SI-jacobi peak time off by "//rstr(prel*100.0_wp)// &
                    "% of T at probe "//rstr(real(p, wp)))
         if (allocated(error)) exit
      end do

      call mesh_unstr_destroy(mesh)
   end subroutine test_estuary_bench

   pure function rstr(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(ES12.5)') x
      s = trim(adjustl(buf))
   end function rstr

end module test_si_unstr_estuary_bench
