!! Barotropic-limit gate for the unstructured 3D multilayer SI step.
module test_si_unstr_ml_baro_limit
   !! A closed-basin seiche with nz=4 UNIFORM density, zero vertical
   !! viscosity and zero bottom drag must reproduce the 2D unstructured SI
   !! seiche: the column solve degenerates to the identity (G_k = u^n_k,
   !! R_k = -1, Hhat = H, Gbar = column transport), so the 3D step IS the
   !! 2D step up to the per-layer reconstruction round-trip.
   !!
   !! The equivalence is essentially EXACT (measured pointwise eta
   !! difference ~1.4e-10 m on a 0.01 m seiche over 5 periods) because the
   !! unstructured ML predictor/recovery use the ANCHOR + INCREMENT form:
   !! the persisted per-layer edge transports are the anchor and only the
   !! implicit-mixing / explicit-tendency increment (h_k*G_k - hu0_k)
   !! routes through the cell average — in the identity-column limit the
   !! increment vanishes and every layer is exactly dsig_k times the 2D
   !! edge update.  (A pure cell-averaged Gbar predictor — the literal
   !! structured shape — would differ from the edge-persistent unstructured
   !! 2D path by a reconstruction filter; the anchor form removes that.)
   !! Residual difference: FP association + the Hhat-vs-h Hf fill
   !! round-off + CG-tolerance accumulation.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_enable_si, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_si_unstr_state, only: SI_PRECOND_NONE
   use rki_timestepper, only: TIME_SEMI_IMPLICIT
   use rki_solver_unstr, only: solver_unstr_step
   use rki_ml_dynamics_unstr, only: ml_init_layers_unstr
   use rki_ml_tracers_unstr, only: ml_seed_standard_tracers_unstr
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_unstr_ml_baro_limit_tests

   real(wp), parameter :: LX = 4000.0_wp
   real(wp), parameter :: LY = 800.0_wp
   integer, parameter :: NXQ = 40
   integer, parameter :: NYQ = 8
   integer, parameter :: NZ = 4
   real(wp), parameter :: H0 = 10.0_wp
   real(wp), parameter :: AMP = 0.01_wp
   real(wp), parameter :: PERIOD_REL_TOL = 1.0e-4_wp
      !! 2D vs 3D measured-period agreement (measured identical to 6
      !! digits; both runs share dt + sampling so quantisation cancels).
   real(wp), parameter :: ETA_POINTWISE_TOL = 1.0e-6_wp
      !! Max |eta_3d - eta_2d| over the run, as a FRACTION of the seiche
      !! amplitude.  Measured 1.4e-8 of amplitude (1.4e-10 m) — the
      !! anchor+increment exact limit; the gate leaves ~70x for
      !! reduction-order variation across compilers/GPU.
   real(wp), parameter :: MASS_DRIFT_REL = 1.0e-11_wp
      !! Same round-off mass gate as the 2D test — the ML eta update is the
      !! identical telescoping edge-flux continuity.

contains

   subroutine collect_si_unstr_ml_baro_limit_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_ml_baro_limit", test_baro_limit) &
                  ]
   end subroutine collect_si_unstr_ml_baro_limit_tests

   subroutine build_rect_mesh(mesh)
      !! Jittered Delaunay-diagonal rectangle triangulation (the P3.1 seiche
      !! construction — acute triangles, passes the TPFA quality gate).
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: nvx, nv, ntri, i, j, vid, t
      real(wp) :: hx, hy, jx, jy
      real(wp), parameter :: JIT = 0.28_wp

      nvx = NXQ + 1
      nv = nvx*(NYQ + 1)
      ntri = 2*NXQ*NYQ
      hx = LX/real(NXQ, wp)
      hy = LY/real(NYQ, wp)

      allocate (vx(nv), vy(nv), tv(3, ntri))
      do j = 0, NYQ
         do i = 0, NXQ
            vid = j*nvx + i + 1
            jx = 0.0_wp
            jy = 0.0_wp
            if (i > 0 .and. i < NXQ) jx = JIT*hx*(hash01(i, j, 1) - 0.5_wp)
            if (j > 0 .and. j < NYQ) jy = JIT*hy*(hash01(i, j, 2) - 0.5_wp)
            vx(vid) = real(i, wp)*hx + jx
            vy(vid) = real(j, wp)*hy + jy
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
   end subroutine build_rect_mesh

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

   subroutine setup_state(state, mesh, dt, with_ml)
      !! Common seiche IC; `with_ml` selects the nz=4 uniform-density
      !! inviscid multilayer configuration (the identity-column limit).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt
      logical, intent(in) :: with_ml

      integer :: k
      real(wp) :: xcc, eta_ic
      logical :: ok

      if (with_ml) then
         call state_unstr_init(state, mesh, nz_ml=NZ)
         ! Identity-column limit: uniform density (EOS off), zero vertical
         ! viscosity, zero bottom drag, no Smagorinsky, no wind, f=0.
         state%multilayer%alpha_T = 0.0_wp
         state%multilayer%beta_S = 0.0_wp
         state%multilayer%vmix_nu_max = 0.0_wp
         state%multilayer%vmix_nu_bg = 0.0_wp
         state%multilayer%smagorinsky_cs = 0.0_wp
      else
         call state_unstr_init(state, mesh)
      end if
      state%barotropic%manning_n = 0.0_wp
      state%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state, mesh, SI_PRECOND_NONE, 1.0e-10_wp, dt, ok)
      if (.not. ok) error stop "si_unstr_ml_baro_limit: quality gate failed"

      do k = 1, mesh%ncells_total
         xcc = mesh%circumcentre_x(k)
         eta_ic = AMP*cos(4.0_wp*atan(1.0_wp)*xcc/LX)
         state%barotropic%bed_cv(k) = 0.0_wp
         state%barotropic%stage_cv(k) = H0 + eta_ic
         state%barotropic%height_cv(k) = H0 + eta_ic
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do
   end subroutine setup_state

   function mass_sum(state, mesh) result(s)
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

   subroutine run_series(state, mesh, n_steps, dt, with_ml, eta_series, &
                         hsum0, hdrift)
      !! Advance n_steps sampling the probe cell; ML runs init layers +
      !! tracers on the device first.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      integer, intent(in) :: n_steps
      real(wp), intent(in) :: dt
      logical, intent(in) :: with_ml
      real(wp), intent(out) :: eta_series(n_steps)
      real(wp), intent(out) :: hsum0, hdrift

      integer :: n, k, kprobe
      real(wp) :: best, xprobe

      xprobe = 0.05_wp*LX
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

      hdrift = 0.0_wp

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      if (with_ml) then
         call ml_init_layers_unstr(state)
         call ml_seed_standard_tracers_unstr(state)
      end if

      hsum0 = mass_sum(state, mesh)
      do n = 1, n_steps
         call solver_unstr_step(state, mesh, dt)
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         eta_series(n) = state%barotropic%stage_cv(kprobe) - H0
         hdrift = max(hdrift, abs(mass_sum(state, mesh) - hsum0))
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_series

   function measure_period(eta_series, n, dt) result(t_meas)
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

   subroutine test_baro_limit(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state_2d, state_3d
      type(mesh_unstr_t) :: mesh
      real(wp) :: t_analytic, dt, courant, dx_eff
      real(wp) :: hsum0_2d, hdrift_2d, hsum0_3d, hdrift_3d
      real(wp) :: t_2d, t_3d, period_rel, ptw_max
      real(wp), allocatable :: eta_2d(:), eta_3d(:)
      integer :: n_steps, n_periods, n

      dx_eff = LX/real(NXQ, wp)
      t_analytic = 2.0_wp*LX/sqrt(GRAVITY*H0)
      courant = 5.0_wp
      dt = courant*dx_eff/sqrt(GRAVITY*H0)
      n_periods = 5
      n_steps = ceiling(real(n_periods, wp)*t_analytic/dt)
      allocate (eta_2d(n_steps), eta_3d(n_steps))

      ! 2D reference run.
      call build_rect_mesh(mesh)
      call setup_state(state_2d, mesh, dt, with_ml=.false.)
      call run_series(state_2d, mesh, n_steps, dt, .false., eta_2d, &
                      hsum0_2d, hdrift_2d)
      call state_unstr_destroy(state_2d)
      call mesh_unstr_destroy(mesh)

      ! nz=4 uniform-density inviscid ML run on the identical mesh.
      call build_rect_mesh(mesh)
      call setup_state(state_3d, mesh, dt, with_ml=.true.)
      call run_series(state_3d, mesh, n_steps, dt, .true., eta_3d, &
                      hsum0_3d, hdrift_3d)
      call state_unstr_destroy(state_3d)
      call mesh_unstr_destroy(mesh)

      t_2d = measure_period(eta_2d, n_steps, dt)
      t_3d = measure_period(eta_3d, n_steps, dt)
      ptw_max = 0.0_wp
      do n = 1, n_steps
         ptw_max = max(ptw_max, abs(eta_3d(n) - eta_2d(n)))
      end do

      call global_logger%info("[si-unstr ml-baro] T_2d="//rstr(t_2d)// &
                              " T_3d="//rstr(t_3d)//" ptw_max="//rstr(ptw_max)// &
                              " ("//rstr(100.0_wp*ptw_max/AMP)//"% of amp)"// &
                              " mass2d="//rstr(hdrift_2d)// &
                              " mass3d="//rstr(hdrift_3d))

      call check(error, t_2d > 0.0_wp .and. t_3d > 0.0_wp, &
                 "si-unstr ml-baro: period measurement failed")
      if (allocated(error)) return
      period_rel = abs(t_3d - t_2d)/t_2d
      call check(error, period_rel < PERIOD_REL_TOL, &
                 "si-unstr ml-baro: 3D vs 2D period disagree by "// &
                 rstr(period_rel*100.0_wp)//"%")
      if (allocated(error)) return
      call check(error, ptw_max < ETA_POINTWISE_TOL*AMP, &
                 "si-unstr ml-baro: pointwise eta diff "//rstr(ptw_max)// &
                 " exceeds "//rstr(ETA_POINTWISE_TOL*AMP)// &
                 " (reconstruction-filter budget)")
      if (allocated(error)) return
      call check(error, hdrift_3d < MASS_DRIFT_REL*hsum0_3d, &
                 "si-unstr ml-baro: 3D mass drift "//rstr(hdrift_3d)// &
                 " above round-off bound")
   end subroutine test_baro_limit

   pure function rstr(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(ES12.5)') x
      s = trim(adjustl(buf))
   end function rstr

end module test_si_unstr_ml_baro_limit
