!! Casulli wet/dry for the unstructured 3D multilayer semi-implicit step.
module test_si_unstr_ml_wetdry
   !! Q8 gates for the unstructured ML+wet/dry port (the triangular twin of
   !! the earlier structured tidal_flat_cycle):
   !!
   !! 1. TIDAL-FLAT cycling — a triangulated sloping-bed channel (deep west
   !!    mouth -> east intertidal flat), clamped tidal west edge, nz=4 sigma
   !!    layers, uniform S=S0.  Over 2+ tidal cycles: wet fraction BREATHES
   !!    (the flat actually dries and re-wets), no negative column depth, no
   !!    NaN, uniform S preserved in the WELL-WET interior (the mass<->tracer
   !!    flux consistency — the dsig-proportional advection + dead-edge R3
   !!    zeroing carried edge-wise), Newton outers under the cap.  This is
   !!    the tidal-flat cycling gate.
   !!
   !! 2. DRY-COLUMN SANITY — a column that fully dries has ALL layers at
   !!    h_layer ~ 0 (dsig scale-together => the H_VANISHED regime, no
   !!    partial-layer artefact) and HOLDS its last-wet hTr (mass ~ 0); on
   !!    re-wetting the concentration recovers ~S0 (hold-last-wet).
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

   public :: collect_si_unstr_ml_wetdry_tests

   real(wp), parameter :: LX = 12000.0_wp, LY = 1500.0_wp
   integer, parameter :: NXQ = 40, NYQ = 6
   integer, parameter :: NZ = 4
   real(wp), parameter :: S0 = 35.0_wp
   real(wp), parameter :: BED_W = -6.0_wp
      !! Deep west channel bed (m); dt uses this depth for the wave CFL.
   real(wp), parameter :: BED_E = 1.5_wp
      !! East intertidal-flat bed (m) — dries at low water.
   real(wp), parameter :: TIDE_AMP = 1.5_wp
   real(wp), parameter :: TIDE_PERIOD = 44712.0_wp
   integer, parameter :: WEST_GROUP = 1
   integer, parameter :: WALL_GROUP = 2
   real(wp), parameter :: S_INV_TOL = 1.0e-4_wp
      !! Max |S - S0| (PSU) over the WELL-WET interior (h>1 m).  With the
      !! consistent (raw-h) column integrals + the wet-side front fallback
      !! the mass<->tracer chain is exact to solver tolerance even ACROSS
      !! the wet/dry front: measured 3.9e-8 PSU over 2.2 tidal cycles at
      !! rtol=1e-11 (gfortran CPU).  Gate = measured x ~2500 headroom —
      !! far below the O(1) PSU signal of either defect this test caught
      !! (the floored-Hhat conveyance mismatch and the empty-upwind
      !! tracer-free front flux).

contains

   subroutine collect_si_unstr_ml_wetdry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("si_unstr_ml_wetdry_tidal_flat", test_tidal_flat), &
                  new_unittest("si_unstr_ml_wetdry_dry_column", test_dry_column) &
                  ]
   end subroutine collect_si_unstr_ml_wetdry_tests

   ! ---------------------------------------------------------------------
   ! Jittered Delaunay-diagonal rectangle triangulation with a clamped west
   ! edge (the P3.1 channel construction, shared by the 2D wet/dry test).
   ! ---------------------------------------------------------------------
   subroutine build_flat_mesh(mesh)
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :), bc_tags(:, :)
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
   end subroutine build_flat_mesh

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
      ax = xbl - xtl
      ay = ybl - ytl
      bx = xbr - xtl
      by = ybr - ytl
      cx = xtr - xtl
      cy = ytr - ytl
      incirc = (ax*ax + ay*ay)*(bx*cy - cx*by) &
               - (bx*bx + by*by)*(ax*cy - cx*ay) &
               + (cx*cx + cy*cy)*(ax*by - bx*ay)
      use_bl_tr = (incirc <= 0.0_wp)
   end function delaunay_diag_bl_tr

   subroutine setup_flat(state, mesh, dt, om)
      !! Sloping bed, still water at eta=0 (east flat starts dry), nz=4
      !! sigma, uniform S=S0.  Clamped tidal west edge.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(out) :: dt, om
      integer :: k
      logical :: ok
      real(wp) :: frac

      call build_flat_mesh(mesh)
      dt = 12.0_wp*(LX/real(NXQ, wp))/sqrt(GRAVITY*abs(BED_W))
      om = TWO_PI/TIDE_PERIOD

      call state_unstr_init(state, mesh, nz_ml=NZ)
      state%time_scheme = TIME_SEMI_IMPLICIT
      state%bc_type(WEST_GROUP) = BC_CLAMPED
      state%bc_type(WALL_GROUP) = BC_WALL
      ! Uniform-density tide: EOS off; the inflow reservoir carries S0.
      ! Vertical mixing off.
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%initial_salinity = S0
      state%multilayer%bc_S_ext(WEST_GROUP) = S0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%kappa_S_bg = 0.0_wp
      state%multilayer%kappa_T_bg = 0.0_wp
      call state_unstr_enable_si(state, mesh, SI_PRECOND_JACOBI, 1.0e-11_wp, &
                                 dt, ok, wetdry=.true.)
      if (.not. ok) return

      do k = 1, mesh%ncells_total
         frac = max(0.0_wp, min(1.0_wp, mesh%circumcentre_x(k)/LX))
         state%barotropic%bed_cv(k) = BED_W + (BED_E - BED_W)*frac
         state%barotropic%stage_cv(k) = max(0.0_wp, state%barotropic%bed_cv(k))
         state%barotropic%height_cv(k) = state%barotropic%stage_cv(k) &
                                         - state%barotropic%bed_cv(k)
         state%barotropic%xmom_cv(k) = 0.0_wp
         state%barotropic%ymom_cv(k) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_init_layers_unstr(state)
      call ml_seed_standard_tracers_unstr(state)
   end subroutine setup_flat

   subroutine test_tidal_flat(error)
      !! Wet-dry-wet cycle on triangles: wet fraction breathes, min column
      !! depth >= 0, uniform S preserved in the well-wet interior.
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: dt, om, t_bc, eta_bc
      real(wp) :: min_h, s_err_deep, wet_lo, wet_hi
      integer :: n, n_steps, k, kz, it_s, nwet, newton_max
      logical :: bad_run, ok
      real(wp) :: h_tot, hk, sc, frac

      call setup_flat(state, mesh, dt, om)
      ok = state%use_semi_implicit
      call check(error, ok, "ml_wetdry tidal_flat: SI enable / quality gate failed")
      if (allocated(error)) then
         call cleanup(state, mesh)
         return
      end if

      it_s = state%multilayer%idx_salinity
      n_steps = ceiling(2.2_wp*TIDE_PERIOD/dt)
      min_h = huge(1.0_wp)
      s_err_deep = 0.0_wp
      wet_lo = huge(1.0_wp)
      wet_hi = -huge(1.0_wp)
      newton_max = 0
      bad_run = .false.

      ! Per-step column CLOSURE audit (the two invariants the port hardens):
      !   water: dh_c == -dt/area * sum_e mass_flux_e*len_e   (row-exact —
      !          catches any fill/recovery Hf inconsistency, e.g. the
      !          floored-Hhat defect this test found), and
      !   salt:  dsalt_c == S0 * dh_c under uniform S (catches any
      !          water-without-salt path, e.g. the empty-upwind front flux).
      block
         real(wp), allocatable :: h_before(:), s_before(:)
         real(wp) :: mass_close_max, salt_close_max
         integer :: i2
         real(wp) :: mflux_div, dh, dsalt
         allocate (h_before(mesh%ncells_owned), s_before(mesh%ncells_owned))
         mass_close_max = 0.0_wp
         salt_close_max = 0.0_wp
         ! Seed the host copies (layer init + tracer seed ran on the device).
         !$acc update self(state%barotropic%stage_cv)
         !$acc update self(state%multilayer%tracers(it_s)%hTr)
         do n = 1, n_steps
            t_bc = state%time + SI_UNSTR_THETA*dt
            eta_bc = TIDE_AMP*sin(om*t_bc)
            state%bc_eta_ext(WEST_GROUP) = eta_bc
            !$acc update device(state%bc_eta_ext)

            ! Pre-step column water + salt (host copies are current: step 1
            ! reads the IC, later steps the previous iteration's update self).
            do k = 1, mesh%ncells_owned
               h_before(k) = max(state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k), 0.0_wp)
               s_before(k) = 0.0_wp
               do kz = 1, NZ
                  s_before(k) = s_before(k) + state%multilayer%tracers(it_s)%hTr(kz, k)
               end do
            end do

            call solver_unstr_step(state, mesh, dt)
            newton_max = max(newton_max, state%si%newton_iters)

            !$acc update self(state%barotropic%stage_cv)
            !$acc update self(state%multilayer%h_layer)
            !$acc update self(state%multilayer%tracers(it_s)%hTr)
            !$acc update self(state%si%mass_flux_edge)

            nwet = 0
            do k = 1, mesh%ncells_owned
               h_tot = max(state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k), 0.0_wp)
               if (.not. (h_tot == h_tot) .or. abs(h_tot) > 1.0e12_wp) bad_run = .true.
               min_h = min(min_h, h_tot)
               if (h_tot > 1.0e-6_wp) nwet = nwet + 1
               ! Column closure vs the blended edge mass flux the eta solve used.
               mflux_div = 0.0_wp
               do i2 = 1, 3
                  mflux_div = mflux_div + state%si%mass_flux_edge(i2, k)*mesh%edge_len(i2, k)
               end do
               dh = h_tot - h_before(k)
               dsalt = -s_before(k)
               do kz = 1, NZ
                  dsalt = dsalt + state%multilayer%tracers(it_s)%hTr(kz, k)
               end do
               mass_close_max = max(mass_close_max, abs(dh + dt*mflux_div/mesh%area(k)))
               salt_close_max = max(salt_close_max, abs(dsalt + S0*dt*mflux_div/mesh%area(k)))
               if (h_tot > 1.0_wp) then
                  do kz = 1, NZ
                     hk = state%multilayer%h_layer(kz, k)
                     if (hk > 1.0e-6_wp) then
                        sc = state%multilayer%tracers(it_s)%hTr(kz, k)/hk
                        s_err_deep = max(s_err_deep, abs(sc - S0))
                     end if
                  end do
               end if
            end do
            frac = real(nwet, wp)/real(mesh%ncells_owned, wp)
            wet_lo = min(wet_lo, frac)
            wet_hi = max(wet_hi, frac)
         end do

         call global_logger%info("[si-unstr ml-wetdry tidal-flat] column closure "// &
                                 "max: water="//rstr(mass_close_max)//" m, salt="// &
                                 rstr(salt_close_max)//" m*PSU")
         call check(error, mass_close_max < 1.0e-7_wp, &
                    "ml_wetdry tidal_flat: water column closure "//rstr(mass_close_max)// &
                    " m > 1e-7 (eta-solve vs recovered-flux inconsistency)")
         if (allocated(error)) then
            call state_unstr_exit_data(state)
            call mesh_unstr_exit_data(mesh)
            call cleanup2(state, mesh)
            return
         end if
         call check(error, salt_close_max < 1.0e-4_wp, &
                    "ml_wetdry tidal_flat: salt column closure "//rstr(salt_close_max)// &
                    " m*PSU > 1e-4 (water moving without salt at the front)")
         if (allocated(error)) then
            call state_unstr_exit_data(state)
            call mesh_unstr_exit_data(mesh)
            call cleanup2(state, mesh)
            return
         end if
      end block

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      call global_logger%info("[si-unstr ml-wetdry tidal-flat] steps="// &
                              istr(n_steps)//" wet_frac=["//rstr(wet_lo)//","// &
                              rstr(wet_hi)//"] min_h="//rstr(min_h)// &
                              " max|S-35|@1m="//rstr(s_err_deep)// &
                              " newton_max="//istr(newton_max))

      call check(error,.not. bad_run, "ml_wetdry tidal_flat: non-finite / blow-up")
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      call check(error, min_h >= -1.0e-9_wp, &
                 "ml_wetdry tidal_flat: negative column depth "//rstr(min_h))
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      call check(error, (wet_hi - wet_lo) > 0.05_wp, &
                 "ml_wetdry tidal_flat: wet fraction did not breathe (["// &
                 rstr(wet_lo)//","//rstr(wet_hi)//"]) — the flat never wet/dried")
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      call check(error, s_err_deep < S_INV_TOL, &
                 "ml_wetdry tidal_flat: uniform S broke in the well-wet interior (h>1m): "// &
                 "max|S-35|="//rstr(s_err_deep)//" >= "//rstr(S_INV_TOL)// &
                 " (mass<->tracer flux inconsistency)")
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      call check(error, newton_max <= 10, &
                 "ml_wetdry tidal_flat: Newton needed "//istr(newton_max)//" outers")
      call cleanup2(state, mesh)
   end subroutine test_tidal_flat

   subroutine test_dry_column(error)
      !! Force the east flat fully dry at extreme low water, then re-flood:
      !! a dried column has all layers ~0 (dsig scale-together, H_VANISHED
      !! regime) and holds its last-wet salt; on re-wetting the
      !! concentration recovers ~S0.
      type(error_type), allocatable, intent(out) :: error

      type(unstr_state_t) :: state
      type(mesh_unstr_t) :: mesh
      real(wp) :: dt, om, t_bc, eta_bc
      integer :: n, n_low, k, kz, it_s, kdry, kwet
      logical :: ok, saw_dry, saw_rewet
      real(wp) :: h_tot, dry_layer_max, rewet_s_err, hk, sc, xmax

      call setup_flat(state, mesh, dt, om)
      ok = state%use_semi_implicit
      call check(error, ok, "ml_wetdry dry_column: SI enable failed")
      if (allocated(error)) then
         call cleanup(state, mesh)
         return
      end if
      it_s = state%multilayer%idx_salinity

      ! Probe the MID-flat cell nearest x = 0.85*LX (bed ~ +0.4 m): dry at
      ! the IC (eta=0 < bed), yet shallow enough that the amplified up-sweep
      ! re-floods it past the h>1 m re-wet gate (the east-most cell at
      ! bed ~ +1.5 m never gets that deep — a vacuous probe).
      kdry = 1
      xmax = huge(1.0_wp)
      do k = 1, mesh%ncells_owned
         if (abs(mesh%circumcentre_x(k) - 0.85_wp*LX) < xmax) then
            xmax = abs(mesh%circumcentre_x(k) - 0.85_wp*LX)
            kdry = k
         end if
      end do

      ! Drive the tide down hard (large negative forcing) so the flat drains,
      ! then back up so it re-floods.
      n_low = ceiling(0.5_wp*TIDE_PERIOD/dt)
      saw_dry = .false.
      saw_rewet = .false.
      dry_layer_max = 0.0_wp
      rewet_s_err = huge(1.0_wp)

      do n = 1, 2*n_low
         t_bc = state%time + SI_UNSTR_THETA*dt
         ! One full down-then-up sweep centred to expose dry + re-wet.
         eta_bc = -2.0_wp*TIDE_AMP*cos(om*t_bc)
         state%bc_eta_ext(WEST_GROUP) = eta_bc
         !$acc update device(state%bc_eta_ext)

         call solver_unstr_step(state, mesh, dt)

         !$acc update self(state%barotropic%stage_cv)
         !$acc update self(state%multilayer%h_layer)
         !$acc update self(state%multilayer%tracers(it_s)%hTr)

         h_tot = max(state%barotropic%stage_cv(kdry) - state%barotropic%bed_cv(kdry), 0.0_wp)
         if (h_tot <= DRY_TOLERANCE) then
            saw_dry = .true.
            ! All layers must have scaled together to ~0 (no partial layer).
            do kz = 1, NZ
               dry_layer_max = max(dry_layer_max, state%multilayer%h_layer(kz, kdry))
            end do
         else if (saw_dry .and. h_tot > 1.0_wp) then
            ! Re-wetted after a dry spell: the held salt must reconstitute ~S0.
            saw_rewet = .true.
            rewet_s_err = 0.0_wp
            do kz = 1, NZ
               hk = state%multilayer%h_layer(kz, kdry)
               if (hk > 1.0e-6_wp) then
                  sc = state%multilayer%tracers(it_s)%hTr(kz, kdry)/hk
                  rewet_s_err = max(rewet_s_err, abs(sc - S0))
               end if
            end do
         end if
      end do

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      call global_logger%info("[si-unstr ml-wetdry dry-column] saw_dry="// &
                              lstr(saw_dry)//" max_dry_layer_h="// &
                              rstr(dry_layer_max)//" saw_rewet="// &
                              lstr(saw_rewet)//" rewet|S-35|="//rstr(rewet_s_err))

      call check(error, saw_dry, &
                 "ml_wetdry dry_column: the flat cell never fully dried")
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      ! dsig scale-together: every layer collapses to the H_VANISHED regime.
      call check(error, dry_layer_max <= DRY_TOLERANCE, &
                 "ml_wetdry dry_column: a layer held finite thickness "// &
                 rstr(dry_layer_max)//" in a dried column (dsig scale-together "// &
                 "broken — partial-layer artefact)")
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      call check(error, saw_rewet, &
                 "ml_wetdry dry_column: the dried cell never re-flooded past "// &
                 "the h>1m gate (vacuous re-wet probe)")
      if (allocated(error)) then
         call cleanup2(state, mesh)
         return
      end if
      call check(error, rewet_s_err < 5.0e-2_wp, &
                 "ml_wetdry dry_column: re-wetted concentration "// &
                 rstr(rewet_s_err)//" off S0 (hold-last-wet broken)")
      call cleanup2(state, mesh)
   end subroutine test_dry_column

   subroutine cleanup(state, mesh)
      !! Pre-enter_data cleanup (enable failed before enter_data).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine cleanup

   subroutine cleanup2(state, mesh)
      !! Post-run cleanup (enter_data already exited by the caller).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine cleanup2

   pure function rstr(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, "(ES12.5)") x
      s = trim(adjustl(buf))
   end function rstr

   pure function istr(k) result(s)
      integer, intent(in) :: k
      character(len=:), allocatable :: s
      character(len=16) :: buf
      write (buf, "(I0)") k
      s = trim(adjustl(buf))
   end function istr

   pure function lstr(b) result(s)
      logical, intent(in) :: b
      character(len=:), allocatable :: s
      if (b) then
         s = "T"
      else
         s = "F"
      end if
   end function lstr

end module test_si_unstr_ml_wetdry
