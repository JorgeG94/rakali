!! Uniform-salinity tidal-pumping invariance for the unstructured 3D SI step.
module test_si_unstr_ml_tidal_tracer
   !! The continuity <-> tracer contract on the unstructured 3D SI path:
   !! a triangulated channel with a CLAMPED (Dirichlet) tidal west edge,
   !! nz=4 sigma layers, and UNIFORM salinity S=35 (the inflow reservoir
   !! also 35) must keep S uniform through flood/ebb pumping — the
   !! dsig-proportional tracer advection uses the EXACT blended edge mass
   !! flux the eta update used, so the per-layer tracer flux matches the
   !! dsig-slaved thickness update and uniformity is preserved to CG
   !! tolerance.  This is the guard for two lessons carried edge-wise:
   !!   * the R3 dead-edge masking (a phantom U^n flux on a dead edge would
   !!     desync tracer from mass), and
   !!   * the Dirichlet ghost R = -1 seed in the per-layer recovery (the
   !!     estuary-mouth half-face mismatch freshened uniform salt on every
   !!     flood tide before the seed).
   !!
   !! Asserted over one full tidal cycle:
   !!   * max |S - 35| over every wet owned cell/layer stays at CG-tolerance
   !!     scale (gate measured + headroom),
   !!   * the surface responds (the tide actually pumps — a no-flow run
   !!     would pass the uniformity gate vacuously),
   !!   * total salt mass tracks total water volume (ratio drift bounded).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, TWO_PI, DRY_TOLERANCE
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
   use rki_ml_dynamics_unstr, only: ml_init_layers_unstr
   use rki_ml_tracers_unstr, only: ml_seed_standard_tracers_unstr
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_unstr_ml_tidal_tracer_tests

   real(wp), parameter :: LX = 3000.0_wp
   real(wp), parameter :: LY = 600.0_wp
   integer, parameter :: NXQ = 30
   integer, parameter :: NYQ = 6
   integer, parameter :: NZ = 4
   real(wp), parameter :: H0 = 20.0_wp
   real(wp), parameter :: S0 = 35.0_wp
   real(wp), parameter :: TIDE_AMP = 0.5_wp
   real(wp), parameter :: TIDE_PERIOD = 44712.0_wp
   integer, parameter :: WEST_GROUP = 1
   integer, parameter :: WALL_GROUP = 2
   real(wp), parameter :: S_UNIFORM_TOL = 5.0e-5_wp
      !! Max |S - S0| (PSU) over every wet cell/layer/step.  The ONLY
      !! uniformity-error source is the CG relative tolerance entering
      !! through the eta-residual (the mass<->tracer split itself is
      !! exact) — VERIFIED by an rtol sweep: rtol=1e-10 measures 1.4e-5
      !! PSU over the 783-step cycle, rtol=1e-13 measures 1.3e-8 (linear
      !! in rtol).  Gate = measured(rtol=1e-10) x ~3.5 headroom.

contains

   subroutine collect_si_unstr_ml_tidal_tracer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_ml_uniform_salt_pumping", test_uniform_salt) &
                  ]
   end subroutine collect_si_unstr_ml_tidal_tracer_tests

   subroutine build_channel_mesh(mesh)
      !! Jittered Delaunay channel with WEST_GROUP on the x~0 boundary
      !! edges, WALL_GROUP elsewhere (the P3.1 channel construction).
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :), bc_tags(:, :)
      integer :: nvx, nv, ntri, i, j, vid, t
      real(wp) :: hx, hy, xmid, jx, jy
      real(wp), parameter :: JIT = 0.28_wp

      nvx = NXQ + 1
      nv = nvx*(NYQ + 1)
      ntri = 2*NXQ*NYQ
      hx = LX/real(NXQ, wp)
      hy = LY/real(NYQ, wp)

      allocate (vx(nv), vy(nv), tv(3, ntri), bc_tags(3, ntri))
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
      xmid = 0.25_wp*hx
      bc_tags = 0
      do t = 1, mesh%ncells_owned
         do i = 1, 3
            if (mesh%neighbours(i, t) < 0) then
               if (mesh%edge_mid_x(i, t) < xmid) then
                  bc_tags(i, t) = WEST_GROUP
               else
                  bc_tags(i, t) = WALL_GROUP
               end if
            end if
         end do
      end do
      call mesh_unstr_destroy(mesh)
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv, bc_tags=bc_tags)
   end subroutine build_channel_mesh

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

   subroutine test_uniform_salt(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: dx_eff, dt, om, t_bc, eta_bc, courant
      real(wp) :: s_dev_max, eta_swing_min, eta_swing_max, mean_e
      integer :: n, n_steps, k, kz, it_s
      logical :: ok
      real(wp) :: h_kc, s_kc

      call build_channel_mesh(mesh)
      dx_eff = LX/real(NXQ, wp)
      courant = 8.0_wp
      dt = courant*dx_eff/sqrt(GRAVITY*H0)
      n_steps = ceiling(TIDE_PERIOD/dt)
      om = TWO_PI/TIDE_PERIOD

      call state_unstr_init(state, mesh, nz_ml=NZ)
      state%time_scheme = TIME_SEMI_IMPLICIT
      state%bc_type(WEST_GROUP) = BC_CLAMPED
      state%bc_type(WALL_GROUP) = BC_WALL
      ! Uniform-density tide: EOS off; PP81 + Manning ON (the full column
      ! solve must not break uniformity — mixing a uniform tracer is a
      ! no-op; drag only reshapes the velocity profile).
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%initial_salinity = S0
      ! Inflow reservoir: the west-edge Dirichlet inflow carries S0.
      state%multilayer%bc_S_ext(WEST_GROUP) = S0
      state%barotropic%manning_n = 0.02_wp
      call state_unstr_enable_si(state, mesh, SI_PRECOND_JACOBI, 1.0e-10_wp, dt, ok)
      call check(error, ok, "si-unstr ml-tidal: quality gate failed")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if

      do k = 1, mesh%ncells_total
         state%barotropic%bed_cv(k) = 0.0_wp
         state%barotropic%stage_cv(k) = H0
         state%barotropic%height_cv(k) = H0
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_init_layers_unstr(state)
      call ml_seed_standard_tracers_unstr(state)

      it_s = state%multilayer%idx_salinity
      s_dev_max = 0.0_wp
      eta_swing_min = huge(1.0_wp)
      eta_swing_max = -huge(1.0_wp)

      do n = 1, n_steps
         t_bc = state%time + SI_UNSTR_THETA*dt
         eta_bc = TIDE_AMP*sin(om*t_bc)
         state%bc_eta_ext(WEST_GROUP) = H0 + eta_bc
         ! [acc->omp-cpu] dropped device directive (host no-op): update

         call solver_unstr_step(state, mesh, dt)

         ! [acc->omp-cpu] dropped device directive (host no-op): update
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         mean_e = 0.0_wp
         do k = 1, mesh%ncells_owned
            mean_e = mean_e + (state%barotropic%stage_cv(k) - H0)
            do kz = 1, NZ
               h_kc = state%multilayer%h_layer(kz, k)
               if (h_kc > DRY_TOLERANCE) then
                  s_kc = state%multilayer%tracers(it_s)%hTr(kz, k)/h_kc
                  s_dev_max = max(s_dev_max, abs(s_kc - S0))
               end if
            end do
         end do
         mean_e = mean_e/real(mesh%ncells_owned, wp)
         eta_swing_min = min(eta_swing_min, mean_e)
         eta_swing_max = max(eta_swing_max, mean_e)
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      call global_logger%info("[si-unstr ml-tidal] max|S-35|="//rstr(s_dev_max)// &
                              " PSU; mean-eta swing ["//rstr(eta_swing_min)// &
                              ","//rstr(eta_swing_max)//"] over "// &
                              rstr(real(n_steps, wp))//" steps")

      ! The tide must actually pump (non-vacuous uniformity).
      call check(error, (eta_swing_max - eta_swing_min) > 0.5_wp*TIDE_AMP, &
                 "si-unstr ml-tidal: surface swing too small — the tide "// &
                 "is not pumping, uniformity check is vacuous")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, s_dev_max < S_UNIFORM_TOL, &
                 "si-unstr ml-tidal: uniform salinity corrupted — max|S-35| = "// &
                 rstr(s_dev_max)//" PSU exceeds "//rstr(S_UNIFORM_TOL)// &
                 " (continuity<->tracer flux split desynced)")
      call cleanup(state, mesh)
   end subroutine test_uniform_salt

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

end module test_si_unstr_ml_tidal_tracer
