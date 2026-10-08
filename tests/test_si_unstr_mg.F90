!! Aggregation-multigrid preconditioner test for the unstructured SI eta-solve.
module test_si_unstr_mg
   !! Runs the SAME funnel-estuary semi-implicit case TWICE on the identical
   !! mesh + physics — once with plain CG (`SI_PRECOND_NONE`) and once with the
   !! aggregation-multigrid PCG (`SI_PRECOND_MG`, `rki_si_unstr_mg`) — and
   !! asserts the go/no-go the Q12 prototype promised:
   !!
   !!   * MG collapses the CG iteration count to TENS (the plateau-breaker):
   !!     mg iters << plain iters, and mg iters below a small absolute cap;
   !!   * the MG-preconditioned solution MATCHES the plain-CG solution to CG
   !!     tolerance (the preconditioner changes only the path to the same SPD
   !!     solve, not the answer) — the surface fields agree field-wide.
   !!
   !! Funnel geometry (tapering width + deepening bed) reproduces the variable
   !! cell-area / transmissibility disparity that drives the triangle plain-CG
   !! plateau; a flat uniform mesh would understate the win.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, TWO_PI
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, state_unstr_init, &
                              state_unstr_enable_si, state_unstr_destroy, &
                              state_unstr_enter_data, state_unstr_exit_data
   use rki_si_unstr_state, only: SI_PRECOND_NONE, SI_PRECOND_MG
   use rki_si_unstr_mg, only: agglom_hierarchy_unstr_t, &
                              agglom_hierarchy_unstr_build, &
                              agglom_hierarchy_unstr_refresh
   use rki_boundary_unstr, only: BC_WALL, BC_CLAMPED
   use rki_timestepper, only: TIME_SEMI_IMPLICIT
   use rki_solver_unstr, only: solver_unstr_step
   use rki_si_unstr_step, only: SI_UNSTR_THETA
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_unstr_mg_tests

   real(wp), parameter :: LX = 12000.0_wp
   real(wp), parameter :: W_MOUTH = 2000.0_wp
   real(wp), parameter :: W_HEAD = 600.0_wp
   integer, parameter :: NXQ = 40, NYQ = 6
      !! 480 triangles — the bench CI size; enough cells for a 3-level
      !! aggregation hierarchy (480 -> 240 -> 120) and a plain-CG count in the
      !! thousands, so the MG collapse to tens is unambiguous.
   real(wp), parameter :: DEPTH_MOUTH = 8.0_wp
   real(wp), parameter :: DEPTH_HEAD = 2.0_wp
   real(wp), parameter :: TIDE_AMP = 0.5_wp
   real(wp), parameter :: TIDE_PERIOD = 44712.0_wp
   real(wp), parameter :: SI_COURANT = 20.0_wp
   real(wp), parameter :: RTOL = 1.0e-9_wp
   integer, parameter :: N_STEPS = 10
   integer, parameter :: WEST_GROUP = 1
   integer, parameter :: WALL_GROUP = 2

contains

   subroutine collect_si_unstr_mg_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_mg_plateau_breaker", test_mg_vs_plain), &
                  new_unittest("si_unstr_mg_zero_hf_slot_alignment", &
                               test_zero_hf_slot_alignment), &
                  new_unittest("si_unstr_mg_iters_flat_refined", test_mg_flatness) &
                  ]
   end subroutine collect_si_unstr_mg_tests

   subroutine test_zero_hf_slot_alignment(error)
      !! REGRESSION for the ke/nbr slot-desynchronization bug: the level-1
      !! operator must keep the transmissibility slot `ke(s,k)` LOCKED to the
      !! neighbour slot `nbr(s,k)` even when a cell has a MIX of zero-Hf and
      !! nonzero-Hf interior edges (reachable in the plain 2D-barotropic path
      !! when a cell dries and height->0 floors an interior edge to Hf=0).
      !!
      !! The load-bearing invariant is OPERATOR SYMMETRY: the transmissibility
      !! k->j must equal j->k (same shared edge, same geometry, Hf zeroed on
      !! both sides).  With the old bug, a cell that skipped its zeroed edge's
      !! slot shifted every subsequent ke into the WRONG nbr slot, so k->j read
      !! a nonzero transmissibility while j->k read zero => a non-symmetric,
      !! non-SPD operator => the area-weighted PCG's <r,M^-1 r> can go <= 0 and
      !! NaN eta.  The all-wet plateau test never exercises the mixed case.
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      type(agglom_hierarchy_unstr_t) :: hier
      real(wp), allocatable :: hf(:, :)
      real(wp) :: coef, dt, dx_eff, tkj, tjk, sym_err, dtol
      integer :: nco, k, i, kn, s, s2, j, ncnt, kz, iz, jz
      logical :: found

      call build_funnel_mesh(mesh)
      nco = mesh%ncells_owned
      dx_eff = LX/real(NXQ, wp)
      dt = SI_COURANT*dx_eff/sqrt(GRAVITY*DEPTH_MOUTH)
      coef = GRAVITY*SI_UNSTR_THETA**2*dt**2

      ! Synthetic per-edge Hf: every interior edge wet (== DEPTH_MOUTH), all
      ! boundary edges dry.  Then punch ONE interior edge to zero on BOTH sides
      ! (a physically consistent dried interior edge) so its two incident cells
      ! carry a mix of zero/nonzero interior Hf.
      allocate (hf(3, nco), source=0.0_wp)
      do k = 1, nco
         do i = 1, 3
            kn = mesh%neighbours(i, k)
            if (kn > 0 .and. kn <= nco) hf(i, k) = DEPTH_MOUTH
         end do
      end do

      ! Pick a cell with >= 2 interior neighbours (so zeroing one leaves a mix).
      kz = -1
      do k = 1, nco
         ncnt = 0
         do i = 1, 3
            kn = mesh%neighbours(i, k)
            if (kn > 0 .and. kn <= nco) ncnt = ncnt + 1
         end do
         if (ncnt >= 2) then
            kz = k
            exit
         end if
      end do
      call check(error, kz > 0, "si-unstr mg zero-Hf: no interior cell found")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if

      ! Zero the first interior edge of kz, and the matching edge on the
      ! neighbour that points back to kz.
      iz = -1
      do i = 1, 3
         kn = mesh%neighbours(i, kz)
         if (kn > 0 .and. kn <= nco) then
            iz = i
            jz = kn
            exit
         end if
      end do
      hf(iz, kz) = 0.0_wp
      do i = 1, 3
         if (mesh%neighbours(i, jz) == kz) hf(i, jz) = 0.0_wp
      end do

      call agglom_hierarchy_unstr_build(hier, mesh, nco)
      call hier%enter_data()
      call agglom_hierarchy_unstr_refresh(hier, mesh, nco, coef, hf)

      ! (A) ke/nbr slot alignment == operator SYMMETRY on level 1.  For every
      !     interior neighbour pair the transmissibility must match both ways.
      dtol = 1.0e-9_wp
      sym_err = 0.0_wp
      do k = 1, hier%lev(1)%nc
         do s = 1, hier%lev(1)%n_nbr(k)
            j = hier%lev(1)%nbr(s, k)
            tkj = hier%lev(1)%ke(s, k)
            tjk = huge(1.0_wp)
            do s2 = 1, hier%lev(1)%n_nbr(j)
               if (hier%lev(1)%nbr(s2, j) == k) tjk = hier%lev(1)%ke(s2, j)
            end do
            sym_err = max(sym_err, abs(tkj - tjk))
         end do
      end do
      call check(error, sym_err < dtol, &
                 "si-unstr mg zero-Hf: level-1 operator non-symmetric "// &
                 "(ke/nbr slots desynchronized) — max |K(k,j)-K(j,k)|="// &
                 rstr(sym_err))
      if (allocated(error)) then
         call cleanup(hier, mesh); return
      end if

      ! (B) The zeroed interior edge lands a ZERO transmissibility in the slot
      !     ALIGNED to its neighbour jz (not skipped / shifted).
      found = .false.
      do s = 1, hier%lev(1)%n_nbr(kz)
         if (hier%lev(1)%nbr(s, kz) == jz) then
            found = .true.
            call check(error, abs(hier%lev(1)%ke(s, kz)) < dtol, &
                       "si-unstr mg zero-Hf: dried interior edge did not write a "// &
                       "zero transmissibility into its aligned slot")
         end if
      end do
      if (.not. allocated(error)) &
         call check(error, found, "si-unstr mg zero-Hf: neighbour jz slot missing")
      if (allocated(error)) then
         call cleanup(hier, mesh); return
      end if

      ! (C) SPD / M-matrix + no NaN: every level's diagonal is finite, strictly
      !     positive and dominates its own off-diagonal row sum (area>0 makes it
      !     strict) — the guarantee that the PCG <r,M^-1 r> stays positive.
      do k = 1, hier%n_levels
         call check_level_spd(error, hier, k)
         if (allocated(error)) then
            call cleanup(hier, mesh); return
         end if
      end do

      call cleanup(hier, mesh)
   end subroutine test_zero_hf_slot_alignment

   subroutine check_level_spd(error, hier, l)
      !! Assert level l's operator is finite, positive-diagonal and strictly
      !! diagonally dominant (=> SPD M-matrix).
      type(error_type), allocatable, intent(inout) :: error
      type(agglom_hierarchy_unstr_t), intent(in) :: hier
      integer, intent(in) :: l
      integer :: k, s
      real(wp) :: d, rowsum

      do k = 1, hier%lev(l)%nc
         d = hier%lev(l)%diag(k)
         if (is_bad(d) .or. is_bad(hier%lev(l)%dinv(k))) then
            call check(error, .false., "si-unstr mg zero-Hf: NaN/Inf diag on level "// &
                       itoa(l))
            return
         end if
         if (d <= 0.0_wp) then
            call check(error, .false., "si-unstr mg zero-Hf: non-positive diag on level "// &
                       itoa(l))
            return
         end if
         rowsum = 0.0_wp
         do s = 1, hier%lev(l)%n_nbr(k)
            if (is_bad(hier%lev(l)%ke(s, k))) then
               call check(error, .false., "si-unstr mg zero-Hf: NaN/Inf ke on level "// &
                          itoa(l))
               return
            end if
            rowsum = rowsum + abs(hier%lev(l)%ke(s, k))
         end do
         if (d < rowsum) then
            call check(error, .false., "si-unstr mg zero-Hf: diagonal not dominant on "// &
                       "level "//itoa(l)//" (operator not an SPD M-matrix)")
            return
         end if
      end do
   end subroutine check_level_spd

   subroutine cleanup(hier, mesh)
      type(agglom_hierarchy_unstr_t), intent(inout) :: hier
      type(mesh_unstr_t), intent(inout) :: mesh
      call hier%exit_data()
      call hier%destroy()
      call mesh_unstr_destroy(mesh)
   end subroutine cleanup

   pure function is_bad(x) result(bad)
      !! .true. for NaN or +/-Inf (no ieee module dependency).
      real(wp), intent(in) :: x
      logical :: bad
      bad = (x /= x) .or. (abs(x) > huge(1.0_wp))
   end function is_bad

   subroutine test_mg_vs_plain(error)
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      real(wp) :: dt, dx_eff
      real(wp), allocatable :: eta_plain(:), eta_mg(:)
      integer :: it_plain, it_mg, nco
      real(wp) :: max_diff
      logical :: ok

      call build_funnel_mesh(mesh)
      dx_eff = LX/real(NXQ, wp)
      dt = SI_COURANT*dx_eff/sqrt(GRAVITY*DEPTH_MOUTH)
      nco = mesh%ncells_owned
      allocate (eta_plain(nco), eta_mg(nco))

      call run_case(mesh, SI_PRECOND_NONE, dt, eta_plain, it_plain, ok)
      call check(error, ok, "si-unstr mg: plain-CG case setup failed")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if

      call run_case(mesh, SI_PRECOND_MG, dt, eta_mg, it_mg, ok)
      call check(error, ok, "si-unstr mg: MG case setup failed")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if

      max_diff = maxval(abs(eta_mg - eta_plain))
      call global_logger%info("[si-unstr mg] last-step CG iters: plain="// &
                              itoa(it_plain)//" mg="//itoa(it_mg)// &
                              "  |eta_mg-eta_plain|_max="//rstr(max_diff))

      ! (1) MG collapses the plateau: plain CG runs many iterations on this
      !     funnel, MG lands in the tens.
      call check(error, it_plain > 100, &
                 "si-unstr mg: plain CG converged in only "//itoa(it_plain)// &
                 " iters — mesh too easy to demonstrate the plateau")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if
      call check(error, it_mg < 60, &
                 "si-unstr mg: MG-PCG took "//itoa(it_mg)//" iters (expected tens "// &
                 "— the aggregation hierarchy is not preconditioning)")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if
      call check(error, it_mg*5 < it_plain, &
                 "si-unstr mg: MG iters "//itoa(it_mg)//" not << plain iters "// &
                 itoa(it_plain)//" (expected >5x reduction)")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if

      ! (2) Same SPD system => same answer to CG tolerance.
      call check(error, max_diff < 1.0e-5_wp, &
                 "si-unstr mg: MG and plain-CG surfaces disagree by "// &
                 rstr(max_diff)//" (> 1e-5) — the MG operator is not the "// &
                 "Helmholtz operator")

      call mesh_unstr_destroy(mesh)
   end subroutine test_mg_vs_plain

   subroutine test_mg_flatness(error)
      !! The PLATEAU-BREAKER property itself: MG-PCG iteration count stays
      !! in the same tens-band on a refine=2 funnel (4x the triangles)
      !! where plain CG's count grows with resolution (the estuary bench
      !! measured the growth; no plain leg here — it would dominate the
      !! suite runtime).  Gate: the SAME absolute band as the base mesh,
      !! i.e. iters-vs-N flat, not merely reduced.
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      real(wp) :: dt, dx_eff
      real(wp), allocatable :: eta_mg(:)
      integer :: it_mg, nco
      logical :: ok

      call build_funnel_mesh(mesh, refine=2)
      dx_eff = LX/real(2*NXQ, wp)
      dt = SI_COURANT*dx_eff/sqrt(GRAVITY*DEPTH_MOUTH)
      nco = mesh%ncells_owned
      allocate (eta_mg(nco))

      call run_case(mesh, SI_PRECOND_MG, dt, eta_mg, it_mg, ok)
      call check(error, ok, "si-unstr mg flatness: refined case setup failed")
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh); return
      end if

      call global_logger%info("[si-unstr mg flat] refined ("//itoa(nco)// &
                              " cells) last-step MG iters="//itoa(it_mg))
      call check(error, it_mg < 60, &
                 "si-unstr mg flatness: refined-mesh MG-PCG took "//itoa(it_mg)// &
                 " iters — the iteration count is growing with resolution "// &
                 "(the plateau-breaker property regressed)")

      call mesh_unstr_destroy(mesh)
   end subroutine test_mg_flatness

   subroutine run_case(mesh, precond, dt, eta_out, iters_last, ok)
      !! One full SI run of N_STEPS on `mesh` with the given preconditioner;
      !! returns the owned-cell surface field + the last step's CG iters.
      type(mesh_unstr_t), intent(inout) :: mesh
      integer, intent(in) :: precond
      real(wp), intent(in) :: dt
      real(wp), intent(out) :: eta_out(:)
      integer, intent(out) :: iters_last
      logical, intent(out) :: ok

      type(unstr_state_t) :: state
      integer :: n, k
      real(wp) :: om, t_bc, eta_bc, xfrac

      iters_last = 0
      om = TWO_PI/TIDE_PERIOD

      call state_unstr_init(state, mesh)
      state%time_scheme = TIME_SEMI_IMPLICIT
      state%bc_type(WEST_GROUP) = BC_CLAMPED
      state%bc_type(WALL_GROUP) = BC_WALL
      call state_unstr_enable_si(state, mesh, precond, RTOL, dt, ok)
      if (.not. ok) then
         call state_unstr_destroy(state); return
      end if

      do k = 1, mesh%ncells_total
         xfrac = min(max(mesh%circumcentre_x(k)/LX, 0.0_wp), 1.0_wp)
         state%barotropic%bed_cv(k) = -(DEPTH_MOUTH &
                                        + (DEPTH_HEAD - DEPTH_MOUTH)*xfrac)
         state%barotropic%stage_cv(k) = 0.0_wp
         state%barotropic%height_cv(k) = -state%barotropic%bed_cv(k)
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)

      do n = 1, N_STEPS
         t_bc = state%time + SI_UNSTR_THETA*dt
         eta_bc = TIDE_AMP*sin(om*t_bc)
         state%bc_eta_ext(WEST_GROUP) = eta_bc
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         call solver_unstr_step(state, mesh, dt)
         iters_last = state%si%cg_iters
      end do

      ! [acc->omp-cpu] dropped device directive (host no-op): update
      do k = 1, mesh%ncells_owned
         eta_out(k) = state%barotropic%stage_cv(k)
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
      call state_unstr_destroy(state)
      ok = .true.
   end subroutine run_case

   subroutine build_funnel_mesh(mesh, refine)
      !! Funnel triangulation: mapped quads, width W(x) tapering mouth->head,
      !! gentle jitter + per-quad Delaunay diagonals (matches the estuary
      !! bench).  `refine` (default 1) scales NXQ/NYQ for the iteration-count
      !! flatness leg.
      type(mesh_unstr_t), intent(out) :: mesh
      integer, intent(in), optional :: refine

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :), bc_tags(:, :)
      integer :: nvx, nv, ntri, i, j, vid, t, nxq_l, nyq_l, rf
      real(wp) :: hx, wx, jx, jy, yfrac
      real(wp), parameter :: JIT = 0.18_wp

      rf = 1
      if (present(refine)) rf = max(1, refine)
      nxq_l = NXQ*rf
      nyq_l = NYQ*rf
      nvx = nxq_l + 1
      nv = nvx*(nyq_l + 1)
      ntri = 2*nxq_l*nyq_l
      hx = LX/real(nxq_l, wp)

      allocate (vx(nv), vy(nv), tv(3, ntri))
      do j = 0, nyq_l
         do i = 0, nxq_l
            vid = j*nvx + i + 1
            wx = W_MOUTH + (W_HEAD - W_MOUTH)*(real(i, wp)/real(nxq_l, wp))
            yfrac = real(j, wp)/real(nyq_l, wp)
            jx = 0.0_wp
            jy = 0.0_wp
            if (i > 0 .and. i < nxq_l) jx = JIT*hx*(hash01(i, j, 1) - 0.5_wp)
            if (j > 0 .and. j < nyq_l) jy = JIT*(wx/real(nyq_l, wp)) &
                                            *(hash01(i, j, 2) - 0.5_wp)
            vx(vid) = real(i, wp)*hx + jx
            vy(vid) = (yfrac - 0.5_wp)*wx + jy
         end do
      end do

      t = 0
      do j = 0, nyq_l - 1
         do i = 0, nxq_l - 1
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

   pure function rstr(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(ES12.5)') x
      s = trim(adjustl(buf))
   end function rstr

   pure function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, '(I0)') i
      s = trim(buf)
   end function itoa

end module test_si_unstr_mg
