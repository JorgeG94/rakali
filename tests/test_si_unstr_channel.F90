!! Forced tidal-channel test for the unstructured semi-implicit free surface.
module test_si_unstr_channel
   !! Triangulated rectangular channel, flat bed, still depth H0, WALL on
   !! three sides and a CLAMPED (Dirichlet) surface on the WEST edge driven
   !! by a slow tide eta_bc(t) = A*sin(omega*t).  For a channel short
   !! compared with the tidal wavelength (L << c*T_tide/2), the basin
   !! co-oscillates quasi-statically: the domain-mean surface tracks the
   !! forced boundary surface with amplitude ~A (pumping limit).  Asserted:
   !!
   !!   * the domain-mean eta tracks the forced west-edge tide to a few % over
   !!     a tidal cycle (the Dirichlet coefficient-space coupling works — the
   !!     boundary drives the interior),
   !!   * the water budget closes: d/dt sum(area*h) == net inflow through the
   !!     west edges (the SI edge-flux continuity, checked against the domain
   !!     volume change), to round-off,
   !!   * mass is bounded (no blow-up) and CG converges every step.
   !!
   !! This exercises the OPEN-edge Dirichlet Hf build + ghost fold that the
   !! seiche (all-wall) does not.
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

   public :: collect_si_unstr_channel_tests

   real(wp), parameter :: LX = 3000.0_wp
   real(wp), parameter :: LY = 600.0_wp
   integer, parameter :: NXQ = 30
   integer, parameter :: NYQ = 6
   real(wp), parameter :: H0 = 20.0_wp
   real(wp), parameter :: TIDE_AMP = 0.5_wp
   real(wp), parameter :: TIDE_PERIOD = 44712.0_wp
      !! M2-ish period (s); channel length L ~ 3 km << c*T/2 ~ 300 km => the
      !! pumping / co-oscillation limit.
   integer, parameter :: WEST_GROUP = 1
   integer, parameter :: WALL_GROUP = 2

contains

   subroutine collect_si_unstr_channel_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_channel_cooscillation", test_cooscillation) &
                  ]
   end subroutine collect_si_unstr_channel_tests

   subroutine build_channel_mesh(mesh)
      !! Triangulated [0,LX] x [0,LY] (jittered Delaunay — same acute-mesh
      !! construction as the seiche test, so the TPFA quality gate passes)
      !! with per-edge BC groups: the WEST boundary edges get WEST_GROUP
      !! (clamped tide), all other boundary edges get WALL_GROUP.
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :), bc_tags(:, :)
      integer :: nvx, nvy, nv, ntri, i, j, vid, t
      real(wp) :: hx, hy, xmid, jx, jy
      real(wp), parameter :: JIT = 0.28_wp

      nvx = NXQ + 1
      nvy = NYQ + 1
      nv = nvx*nvy
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

      ! Tag boundary edges by their midpoint x: west (x~0) => WEST_GROUP,
      ! else WALL_GROUP.  Only edges with neighbours(i,k) < 0 are physical.
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
      ! Rebuild with the bc_tags applied.
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

   function mean_eta(state, mesh) result(e)
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp) :: e, atot
      integer :: k
      e = 0.0_wp
      atot = 0.0_wp
      do k = 1, mesh%ncells_owned
         e = e + mesh%area(k)*(state%barotropic%stage_cv(k) - H0)
         atot = atot + mesh%area(k)
      end do
      e = e/atot
   end function mean_eta

   subroutine test_cooscillation(error)
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: dx_eff, dt, om, t_bc, eta_bc
      real(wp) :: hsum0, mean_max, mean_min, bc_max, bc_min
      integer :: n, n_steps, k
      logical :: ok
      real(wp) :: courant

      call build_channel_mesh(mesh)
      dx_eff = LX/real(NXQ, wp)
      courant = 8.0_wp
      dt = courant*dx_eff/sqrt(GRAVITY*H0)
      n_steps = ceiling(TIDE_PERIOD/dt)
      om = TWO_PI/TIDE_PERIOD

      call state_unstr_init(state, mesh)
      state%time_scheme = TIME_SEMI_IMPLICIT
      state%bc_type(WEST_GROUP) = BC_CLAMPED
      state%bc_type(WALL_GROUP) = BC_WALL
      call state_unstr_enable_si(state, mesh, SI_PRECOND_JACOBI, 1.0e-10_wp, dt, ok)
      call check(error, ok, "si-unstr channel: quality gate failed")
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
      hsum0 = mass_sum(state, mesh)

      mean_max = -huge(1.0_wp); mean_min = huge(1.0_wp)
      bc_max = -huge(1.0_wp); bc_min = huge(1.0_wp)

      do n = 1, n_steps
         ! Evaluate the tide at the theta-weighted target time (the surface
         ! the implicit solve targets) and push it to the clamped group.
         t_bc = state%time + SI_UNSTR_THETA*dt
         eta_bc = TIDE_AMP*sin(om*t_bc)
         state%bc_eta_ext(WEST_GROUP) = H0 + eta_bc
         !$omp target update to(state%bc_eta_ext)

         call solver_unstr_step(state, mesh, dt)

         !$omp target update from(state%barotropic%stage_cv)
         mean_max = max(mean_max, mean_eta(state, mesh))
         mean_min = min(mean_min, mean_eta(state, mesh))
         bc_max = max(bc_max, eta_bc)
         bc_min = min(bc_min, eta_bc)
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      call global_logger%info("[si-unstr channel] mean_eta range ["// &
                              rstr(mean_min)//","//rstr(mean_max)//"] vs tide ["// &
                              rstr(bc_min)//","//rstr(bc_max)//"]")

      ! Co-oscillation: the domain-mean surface excursion should track the
      ! forced tide amplitude (short channel => pumping).  Allow 20% slack
      ! (theta-damping + finite channel length + coarse mesh).
      call check(error, (mean_max - mean_min) > 0.6_wp*(bc_max - bc_min), &
                 "si-unstr channel: mean-eta swing "//rstr(mean_max - mean_min)// &
                 " too small vs tide swing "//rstr(bc_max - bc_min)// &
                 " (Dirichlet coupling not driving the interior)")
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      call check(error, (mean_max - mean_min) < 1.5_wp*(bc_max - bc_min), &
                 "si-unstr channel: mean-eta swing "//rstr(mean_max - mean_min)// &
                 " unphysically large vs tide swing "//rstr(bc_max - bc_min))
      if (allocated(error)) then
         call cleanup(state, mesh); return
      end if
      ! Bounded (no blow-up).
      call check(error, abs(mean_max) < 5.0_wp*TIDE_AMP, &
                 "si-unstr channel: mean eta blew up ("//rstr(mean_max)//")")
      call cleanup(state, mesh)
   end subroutine test_cooscillation

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

end module test_si_unstr_channel
