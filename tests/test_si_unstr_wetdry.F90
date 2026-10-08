!! Casulli wetting/drying tests for the unstructured semi-implicit surface.
module test_si_unstr_wetdry
   !! Three gates for the P3.3 unstructured SI wet/dry:
   !!
   !! 1. THACKER planar oscillation on a triangulated paraboloid basin
   !!    (Thacker 1981; the prototype's disc case on a jittered-Delaunay
   !!    rectangle — the dry corners never wet, so the domain shape is
   !!    immaterial).  u is spatially uniform => advection == 0 => the
   !!    linear-momentum SI core is exact: a legitimate finite-amplitude
   !!    MOVING-SHORELINE benchmark.  Gates vs the prototype disc numbers
   !!    (period err ~0%, shore excursion ~= gamma0*R^2/D0 = 200 m to mesh
   !!    quantisation, minV = 0 exactly, mass at round-off, Newton 1-3
   !!    outer / max 5, nonconv 0).
   !!
   !! 2. TIDAL-FLAT cycling (sloping bed, clamped tidal west edge): wet
   !!    fraction breathes, min cell volume >= 0, and the DOMAIN WATER
   !!    BUDGET closes against the integrated open-edge flux (LOGGED — the
   !!    open-boundary closure diagnostic).  This is ALSO the R3 dead-edge
   !!    guard: interior edges dry and rewet every cycle, and a phantom
   !!    (1-theta)*U^n flux on a dead interior edge would break the budget
   !!    and the mass round-off — tested, not assumed.
   !!
   !! 3. ALL-WET agreement: on a deep basin that never dries, the wet/dry
   !!    solver must agree with the linear solver (~1e-9 — bounded by the
   !!    Newton-vs-direct stopping criteria, not bitwise; TRUE byte
   !!    identity is the wetdry=.false. gate, where the wet/dry code is
   !!    never called — anchored by the untouched `test_si_unstr_seiche`
   !!    numbers).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, TWO_PI
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_enable_si, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_si_unstr_state, only: SI_PRECOND_JACOBI
   use rki_boundary_unstr, only: BC_WALL, BC_CLAMPED
   use rki_timestepper, only: TIME_SEMI_IMPLICIT
   use rki_solver_unstr, only: solver_unstr_step
   use rki_si_unstr_step, only: SI_UNSTR_THETA
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_unstr_wetdry_tests

   integer, parameter :: WEST_GROUP = 1
   integer, parameter :: WALL_GROUP = 2

contains

   subroutine collect_si_unstr_wetdry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_wetdry_thacker", test_thacker), &
                  new_unittest("si_unstr_wetdry_tidal_flat", test_tidal_flat), &
                  new_unittest("si_unstr_wetdry_all_wet", test_all_wet) &
                  ]
   end subroutine collect_si_unstr_wetdry_tests

   ! ---------------------------------------------------------------------
   ! Jittered Delaunay-diagonal rectangle triangulation (the P3.1 acute
   ! construction), with optional per-edge BC tags for a clamped west edge.
   ! ---------------------------------------------------------------------
   subroutine build_rect_mesh(mesh, x0, lx, ly, nxq, nyq, tag_west)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp), intent(in) :: x0, lx, ly
      integer, intent(in) :: nxq, nyq
      logical, intent(in) :: tag_west

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :), bc_tags(:, :)
      integer :: nvx, nv, ntri, i, j, vid, t
      real(wp) :: hx, hy, jx, jy
      real(wp), parameter :: JIT = 0.28_wp

      nvx = nxq + 1
      nv = nvx*(nyq + 1)
      ntri = 2*nxq*nyq
      hx = lx/real(nxq, wp)
      hy = ly/real(nyq, wp)

      allocate (vx(nv), vy(nv), tv(3, ntri))
      do j = 0, nyq
         do i = 0, nxq
            vid = j*nvx + i + 1
            jx = 0.0_wp
            jy = 0.0_wp
            if (i > 0 .and. i < nxq) jx = JIT*hx*(hash01(i, j, 1) - 0.5_wp)
            if (j > 0 .and. j < nyq) jy = JIT*hy*(hash01(i, j, 2) - 0.5_wp)
            vx(vid) = x0 + real(i, wp)*hx + jx
            vy(vid) = real(j, wp)*hy + jy
         end do
      end do

      t = 0
      do j = 0, nyq - 1
         do i = 0, nxq - 1
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
      if (tag_west) then
         allocate (bc_tags(3, ntri), source=0)
         do t = 1, mesh%ncells_owned
            do i = 1, 3
               if (mesh%neighbours(i, t) < 0) then
                  if (mesh%edge_mid_x(i, t) < x0 + 0.25_wp*hx) then
                     bc_tags(i, t) = WEST_GROUP
                  else
                     bc_tags(i, t) = WALL_GROUP
                  end if
               end if
            end do
         end do
         call mesh_unstr_destroy(mesh)
         call mesh_unstr_init_from_triangles(mesh, vx, vy, tv, bc_tags=bc_tags)
      end if
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

   function wet_volume(state, mesh) result(v)
      !! Total wetted volume sum_c area_c*max(stage-b, 0) (host).
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp) :: v
      integer :: k
      v = 0.0_wp
      do k = 1, mesh%ncells_owned
         v = v + mesh%area(k)*max(state%barotropic%stage_cv(k) &
                                  - state%barotropic%bed_cv(k), 0.0_wp)
      end do
   end function wet_volume

   subroutine test_thacker(error)
      !! Thacker planar oscillation: paraboloid bed z_b = -D0*(1 - r^2/R^2),
      !! planar-tilt IC eta = gamma0*x; omega = sqrt(2*g*D0)/R.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: R_BASIN = 4000.0_wp
      real(wp), parameter :: D0 = 10.0_wp
      real(wp), parameter :: ETA0 = 0.5_wp
      integer, parameter :: NXQ = 40, NYQ = 40

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: lhalf, gamma0, omega, t_exact, dt, c_wave, dx_eff
      real(wp) :: r2, mass0, mdrift, min_v, shore_min, shore_max, shore_pp
      real(wp) :: t_meas, err_rel, xcc, ycc
      real(wp), allocatable :: eta_series(:)
      integer :: n_steps, n_periods, k, n, kprobe, newton_max
      logical :: ok
      real(wp) :: best, vk

      lhalf = 1.15_wp*R_BASIN
      gamma0 = ETA0/R_BASIN
      omega = sqrt(2.0_wp*GRAVITY*D0)/R_BASIN
      t_exact = TWO_PI/omega
      call build_rect_mesh(mesh, -lhalf, 2.0_wp*lhalf, 2.0_wp*lhalf, &
                           NXQ, NYQ, tag_west=.false.)

      c_wave = sqrt(GRAVITY*D0)
      dx_eff = 2.0_wp*lhalf/real(NXQ, wp)
      dt = 5.0_wp*dx_eff/c_wave
      n_periods = 3
      n_steps = ceiling(real(n_periods, wp)*t_exact/dt)
      dt = real(n_periods, wp)*t_exact/real(n_steps, wp)
      allocate (eta_series(n_steps))

      call state_unstr_init(state, mesh)
      state%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state, mesh, SI_PRECOND_JACOBI, 1.0e-10_wp, &
                                 dt, ok, wetdry=.true.)
      call check(error, ok, "thacker: quality gate failed")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if

      ! Paraboloid bed (y-centred at lhalf) + planar-tilt IC; dry cells sit
      ! at stage = bed.
      do k = 1, mesh%ncells_total
         xcc = mesh%circumcentre_x(k)
         ycc = mesh%circumcentre_y(k) - lhalf
         r2 = xcc*xcc + ycc*ycc
         state%barotropic%bed_cv(k) = -D0*(1.0_wp - r2/(R_BASIN*R_BASIN))
         state%barotropic%stage_cv(k) = max(gamma0*xcc, state%barotropic%bed_cv(k))
         state%barotropic%height_cv(k) = state%barotropic%stage_cv(k) &
                                         - state%barotropic%bed_cv(k)
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do

      ! Probe near (0.6R, centre).
      kprobe = 1
      best = huge(1.0_wp)
      do k = 1, mesh%ncells_owned
         if (abs(mesh%circumcentre_x(k) - 0.6_wp*R_BASIN) &
             + abs(mesh%circumcentre_y(k) - lhalf) < best) then
            best = abs(mesh%circumcentre_x(k) - 0.6_wp*R_BASIN) &
                   + abs(mesh%circumcentre_y(k) - lhalf)
            kprobe = k
         end if
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      mass0 = wet_volume(state, mesh)
      mdrift = 0.0_wp
      min_v = huge(1.0_wp)
      shore_min = huge(1.0_wp)
      shore_max = -huge(1.0_wp)
      newton_max = 0

      do n = 1, n_steps
         call solver_unstr_step(state, mesh, dt)
         newton_max = max(newton_max, state%si%newton_iters)
         !$omp target update from(state%barotropic%stage_cv)
         eta_series(n) = state%barotropic%stage_cv(kprobe)
         mdrift = max(mdrift, abs(wet_volume(state, mesh) - mass0))
         ! Shoreline x = max circumcentre-x over wet cells; min cell v.
         block
            real(wp) :: xs
            xs = -huge(1.0_wp)
            do k = 1, mesh%ncells_owned
               vk = max(state%barotropic%stage_cv(k) &
                        - state%barotropic%bed_cv(k), 0.0_wp)
               min_v = min(min_v, vk)
               if (vk > 1.0e-6_wp) xs = max(xs, mesh%circumcentre_x(k))
            end do
            shore_min = min(shore_min, xs)
            shore_max = max(shore_max, xs)
         end block
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      shore_pp = shore_max - shore_min
      t_meas = measure_period(eta_series, n_steps, dt)
      err_rel = abs(t_meas - t_exact)/t_exact

      call global_logger%info("[si-unstr thacker] T_exact="//rstr(t_exact)// &
                              " T_meas="//rstr(t_meas)//" (err "// &
                              rstr(err_rel*100.0_wp)//"%) shore_pp="// &
                              rstr(shore_pp)//" (exact "// &
                              rstr(gamma0*R_BASIN*R_BASIN/D0)//") minV="// &
                              rstr(min_v)//" mass_drift_rel="// &
                              rstr(mdrift/mass0)//" newton_max="// &
                              rstr(real(newton_max, wp)))

      call check(error, t_meas > 0.0_wp, "thacker: period measurement failed")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, err_rel < 0.05_wp, &
                 "thacker period err "//rstr(err_rel*100.0_wp)//"% > 5%")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      ! Shore excursion ~ gamma0*R^2/D0 = 200 m, quantised by dx ~ 230 m
      ! (prototype disc: 187-208 m at similar resolution).
      call check(error, shore_pp > 100.0_wp .and. shore_pp < 350.0_wp, &
                 "thacker shoreline excursion "//rstr(shore_pp)// &
                 " outside [100,350] (exact 200, mesh-quantised)")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, min_v >= 0.0_wp, &
                 "thacker: negative cell water depth "//rstr(min_v))
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, mdrift < 1.0e-11_wp*mass0, &
                 "thacker: mass drift "//rstr(mdrift/mass0)//" rel > 1e-11")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, newton_max <= 6, &
                 "thacker: Newton needed "//rstr(real(newton_max, wp))// &
                 " outers (prototype max 5)")
      call cleanup(state, mesh)
   end subroutine test_thacker

   subroutine test_tidal_flat(error)
      !! Sloping-bed channel, clamped tidal west edge: wet/dry breathing +
      !! the domain WATER BUDGET closed against the integrated open-edge
      !! flux (logged) — the R3 dead-interior-edge guard.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: LX = 12000.0_wp, LY = 1500.0_wp
      integer, parameter :: NXQ = 40, NYQ = 6
      real(wp), parameter :: TIDE_AMP = 1.5_wp
      real(wp), parameter :: TIDE_PERIOD = 44712.0_wp

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: dt, om, t_bc, eta_bc, theta
      real(wp) :: v_start, water_in, water_err, min_v, wet_lo, wet_hi
      integer :: n, n_steps, k, i, nwet, newton_max
      logical :: ok
      real(wp) :: vk, q_step, frac

      call build_rect_mesh(mesh, 0.0_wp, LX, LY, NXQ, NYQ, tag_west=.true.)
      dt = 12.0_wp*(LX/real(NXQ, wp))/sqrt(GRAVITY*6.0_wp)
      n_steps = ceiling(TIDE_PERIOD/dt)
      om = TWO_PI/TIDE_PERIOD
      theta = SI_UNSTR_THETA

      call state_unstr_init(state, mesh)
      state%time_scheme = TIME_SEMI_IMPLICIT
      state%bc_type(WEST_GROUP) = BC_CLAMPED
      state%bc_type(WALL_GROUP) = BC_WALL
      call state_unstr_enable_si(state, mesh, SI_PRECOND_JACOBI, 1.0e-10_wp, &
                                 dt, ok, wetdry=.true.)
      call check(error, ok, "tidal_flat: quality gate failed")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if

      ! Bed -6 m (mouth) -> +1.5 m (inner flat); IC still water at 0.
      do k = 1, mesh%ncells_total
         state%barotropic%bed_cv(k) = -6.0_wp &
                                      + 7.5_wp*(mesh%circumcentre_x(k)/LX)
         state%barotropic%stage_cv(k) = max(0.0_wp, state%barotropic%bed_cv(k))
         state%barotropic%height_cv(k) = state%barotropic%stage_cv(k) &
                                         - state%barotropic%bed_cv(k)
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      v_start = wet_volume(state, mesh)
      water_in = 0.0_wp
      water_err = 0.0_wp
      min_v = huge(1.0_wp)
      wet_lo = huge(1.0_wp)
      wet_hi = -huge(1.0_wp)
      newton_max = 0

      do n = 1, n_steps
         t_bc = state%time + theta*dt
         eta_bc = TIDE_AMP*sin(om*t_bc)
         state%bc_eta_ext(WEST_GROUP) = eta_bc
         !$omp target update to(state%bc_eta_ext)

         call solver_unstr_step(state, mesh, dt)
         newton_max = max(newton_max, state%si%newton_iters)

         !$omp target update from(state%barotropic%stage_cv)
         !$omp target update from(state%si%u_edge)
         !$omp target update from(state%si%u_n_edge)

         ! Open-edge inflow this step: the EXACT blended flux the eta
         ! update used, theta*U^{n+1} + (1-theta)*U^n, + along the OUTWARD
         ! normal => inflow is the negative.
         q_step = 0.0_wp
         do k = 1, mesh%ncells_owned
            do i = 1, 3
               if (mesh%neighbours(i, k) < 0) then
                  if (mesh%bc_group(i, k) == WEST_GROUP) then
                     q_step = q_step - (theta*state%si%u_edge(i, k) &
                                        + (1.0_wp - theta)*state%si%u_n_edge(i, k)) &
                              *mesh%edge_len(i, k)*dt
                  end if
               end if
            end do
         end do
         water_in = water_in + q_step
         water_err = max(water_err, &
                         abs(wet_volume(state, mesh) - v_start - water_in))

         nwet = 0
         do k = 1, mesh%ncells_owned
            vk = max(state%barotropic%stage_cv(k) &
                     - state%barotropic%bed_cv(k), 0.0_wp)
            min_v = min(min_v, vk)
            if (vk > 1.0e-6_wp) nwet = nwet + 1
         end do
         frac = real(nwet, wp)/real(mesh%ncells_owned, wp)
         wet_lo = min(wet_lo, frac)
         wet_hi = max(wet_hi, frac)
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      call global_logger%info("[si-unstr tidal-flat] WATER BUDGET |dV-inflow|/V0="// &
                              rstr(water_err/v_start)//" (V0="//rstr(v_start)// &
                              " m^3); wet fraction ["//rstr(wet_lo)//","// &
                              rstr(wet_hi)//"]; minV="//rstr(min_v)// &
                              "; newton_max="//rstr(real(newton_max, wp)))

      ! Breathing: the flats must actually dry and rewet (non-vacuous R3).
      call check(error, (wet_hi - wet_lo) > 0.05_wp, &
                 "tidal_flat: wet fraction did not breathe ("//rstr(wet_lo)// &
                 ".."//rstr(wet_hi)//") — R3 test is vacuous")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, min_v >= 0.0_wp, &
                 "tidal_flat: negative water depth "//rstr(min_v))
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      ! The budget closes to CG-tolerance scale: a phantom flux on a dead
      ! (dried) INTERIOR edge, or a frozen-Hf back-sub mismatch (R2), shows
      ! up here first (prototype closed to 6.7e-13 rel).
      call check(error, water_err < 1.0e-9_wp*v_start, &
                 "tidal_flat: water budget error "//rstr(water_err/v_start)// &
                 " rel > 1e-9 (dead-edge phantom flux / R2 breach)")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, newton_max <= 8, &
                 "tidal_flat: Newton needed "//rstr(real(newton_max, wp))//" outers")
      call cleanup(state, mesh)
   end subroutine test_tidal_flat

   subroutine test_all_wet(error)
      !! Deep flat basin (never dries): wetdry=.true. must agree with the
      !! linear solver to ~1e-9 (stopping-criteria bound, not bitwise).
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: LX = 4000.0_wp, LY = 800.0_wp
      integer, parameter :: NXQ = 20, NYQ = 4
      real(wp), parameter :: H0 = 50.0_wp
      real(wp), parameter :: AMP = 0.05_wp
      integer, parameter :: NSTEPS = 12

      type(unstr_state_t) :: state_wd, state_lin
      type(mesh_unstr_t) :: mesh
      real(wp) :: dt, diff_max
      real(wp), allocatable :: stage_wd(:)
      integer :: k, n
      logical :: ok

      call build_rect_mesh(mesh, 0.0_wp, LX, LY, NXQ, NYQ, tag_west=.false.)
      dt = 10.0_wp*(LX/real(NXQ, wp))/sqrt(GRAVITY*H0)

      ! Wet/dry run.
      call state_unstr_init(state_wd, mesh)
      state_wd%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state_wd, mesh, SI_PRECOND_JACOBI, 1.0e-12_wp, &
                                 dt, ok, wetdry=.true.)
      call check(error, ok, "all_wet: gate failed (wetdry)")
      if (allocated(error)) then
         call cleanup(state_wd, mesh); return
      end if
      call seiche_ic(state_wd, mesh, LX, H0, AMP)
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state_wd)
      do n = 1, NSTEPS
         call solver_unstr_step(state_wd, mesh, dt)
      end do
      !$omp target update from(state_wd%barotropic%stage_cv)
      allocate (stage_wd(mesh%ncells_owned))
      do k = 1, mesh%ncells_owned
         stage_wd(k) = state_wd%barotropic%stage_cv(k)
      end do
      call state_unstr_exit_data(state_wd)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state_wd)

      ! Linear run on the same mesh/IC (same tight rtol so the comparison
      ! is bounded by the Newton-vs-direct stopping difference).
      call state_unstr_init(state_lin, mesh)
      state_lin%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state_lin, mesh, SI_PRECOND_JACOBI, 1.0e-12_wp, &
                                 dt, ok)
      call check(error, ok, "all_wet: gate failed (linear)")
      if (allocated(error)) then
         call cleanup(state_lin, mesh); return
      end if
      call seiche_ic(state_lin, mesh, LX, H0, AMP)
      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state_lin)
      do n = 1, NSTEPS
         call solver_unstr_step(state_lin, mesh, dt)
      end do
      !$omp target update from(state_lin%barotropic%stage_cv)
      diff_max = 0.0_wp
      do k = 1, mesh%ncells_owned
         diff_max = max(diff_max, &
                        abs(stage_wd(k) - state_lin%barotropic%stage_cv(k)))
      end do
      call state_unstr_exit_data(state_lin)
      call mesh_unstr_exit_data(mesh)

      call global_logger%info("[si-unstr all-wet] max|stage_wd - stage_lin|="// &
                              rstr(diff_max)//" m over "// &
                              rstr(real(NSTEPS, wp))//" steps")

      call check(error, diff_max < 1.0e-9_wp, &
                 "all_wet: wet/dry vs linear disagree by "//rstr(diff_max)// &
                 " m (> 1e-9; the all-wet gated-mean Hf should reduce "// &
                 "exactly to the linear edge depth)")
      call cleanup(state_lin, mesh)
   end subroutine test_all_wet

   subroutine seiche_ic(state, mesh, lx, h0, amp)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: lx, h0, amp
      integer :: k
      do k = 1, mesh%ncells_total
         state%barotropic%bed_cv(k) = -h0
         state%barotropic%stage_cv(k) = &
            amp*cos(4.0_wp*atan(1.0_wp)*mesh%circumcentre_x(k)/lx)
         state%barotropic%height_cv(k) = state%barotropic%stage_cv(k) + h0
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do
   end subroutine seiche_ic

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

end module test_si_unstr_wetdry
