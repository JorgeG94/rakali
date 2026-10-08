!! Tests for the semi-implicit vertical-mixing closures (k-epsilon + KPP)
!! on the unstructured backend (Q1: closure -> si%nu_iface seam).
module test_si_closures_unstr
   !! Five tests mirroring the Q1 prototype certification on the unstructured
   !! `(nz, ncells)` layout:
   !!
   !! 1. keps_si_dt_invariance_unstr — the k-eps column integrator's STEADY
   !!    state is dt-independent: frozen sheared column at dt=10 ref vs
   !!    dt=240/480, interior nu_t agrees to 2% (kernel-level on flat arrays).
   !! 2. keps_si_positivity_decay_unstr — k, eps >= floors under shock-off
   !!    decay at dt=480 (spin up with strong walls, then kill the forcing).
   !! 3. si_keps_wind_channel_unstr — full-state: SI+keps closed basin, assert
   !!    stability + tke/eps >= floors + wind-setup tilt sign.
   !! 4. si_kpp_stratified_unstr — full-state: SI+KPP stratified basin, assert
   !!    stability + T within seeded envelope + column salt conservation.
   !! 5. si_kpp_nonlocal_conservation_unstr — kernel-level: the implicit KPP
   !!    tracer kernel conserves column S/T content to 1e-12.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, CP_WATER, KEPS_STAB_CONSTANT
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
   use rki_ml_keps_unstr, only: ml_keps_update_unstr_impl
   use rki_si_unstr_closures, only: si_unstr_kpp_diffusion_tracers
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_closures_unstr_tests

   ! Kernel-level cell count for tests 1, 2, 5.
   integer, parameter :: N = 8
   real(wp), parameter :: K_MIN = 1.0e-8_wp, EPS_MIN = 1.0e-12_wp
      !! Floors mirrored from rki_ml_keps (positivity assertions).

   ! Full-state mesh parameters (tests 3, 4).
   real(wp), parameter :: LX = 4000.0_wp
   real(wp), parameter :: LY = 800.0_wp
   integer, parameter :: NXQ = 20
   integer, parameter :: NYQ = 4
   real(wp), parameter :: H0 = 20.0_wp
   real(wp), parameter :: JIT = 0.28_wp

contains

   subroutine collect_si_closures_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("keps_si_dt_invariance_unstr", test_keps_dt_invariance), &
                  new_unittest("keps_si_positivity_decay_unstr", test_keps_positivity_decay), &
                  new_unittest("si_keps_wind_channel_unstr", test_si_keps_wind_channel), &
                  new_unittest("si_keps_picard2_unstr", test_si_keps_picard2), &
                  new_unittest("si_kpp_stratified_unstr", test_si_kpp_stratified), &
                  new_unittest("si_kpp_nonlocal_conservation_unstr", test_si_kpp_nonlocal) &
                  ]
   end subroutine collect_si_closures_unstr_tests

   ! ---------------------------------------------------------------------------
   ! Test 1: frozen sheared column — the k-eps steady state is dt-independent.
   ! (nz, ncells) layout, explicit-shape flat arrays — kernel-level only.
   ! ---------------------------------------------------------------------------
   subroutine test_keps_dt_invariance(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 10
      real(wp), parameter :: T_END = 5.0e4_wp, USTAR = 0.01_wp
      real(wp) :: nut_ref(NZ + 1, N)
      real(wp) :: nut_si(NZ + 1, N)
      real(wp) :: rel, rel_max
      integer :: m, idt
      real(wp), parameter :: dts(2) = [240.0_wp, 480.0_wp]

      call run_frozen_column(NZ, 10.0_wp, T_END, USTAR, nut_ref)
      do idt = 1, 2
         call run_frozen_column(NZ, dts(idt), T_END, USTAR, nut_si)
         rel_max = 0.0_wp
         ! Interior interfaces only: 3..NZ-1 (skip the wall-pinned rows 2, NZ
         ! and the copy rows 1, NZ+1).
         do m = 3, NZ - 1
            rel = abs(nut_si(m, 4) - nut_ref(m, 4)) &
                  /max(nut_ref(m, 4), 1.0e-6_wp)
            if (rel > rel_max) rel_max = rel
         end do
         call global_logger%info("[si_closures_unstr dt_invariance] dt=" &
                                 //rstr(dts(idt))//" max rel dnu = "//rstr(rel_max))
         call check(error, rel_max < 0.02_wp, &
                    "k-eps steady nu_t (unstr) at dt="//rstr(dts(idt))// &
                    " deviates from dt=10 reference: "//rstr(rel_max))
         if (allocated(error)) return
      end do
   end subroutine test_keps_dt_invariance

   subroutine run_frozen_column(nz, dt, t_end, ustar, nut_out)
      !! Advance ml_keps_update_unstr_impl on a FROZEN sheared column to
      !! t_end at dt; return the nu_t field.  (nz, N) layout.
      integer, intent(in) :: nz
      real(wp), intent(in) :: dt, t_end, ustar
      real(wp), intent(out) :: nut_out(nz + 1, N)

      real(wp) :: h(nz, N), hu(nz, N), hv(nz, N), rho(nz, N)
      real(wp) :: tke(nz + 1, N), eps(nz + 1, N)
      real(wp) :: nut(nz + 1, N), kt(nz + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      integer :: k, istep, nsteps

      h = 2.0_wp                        ! H = 2*nz metres
      do k = 1, nz
         hu(k, :) = h(k, :)*(0.02_wp*real(k, wp))   ! mild linear shear
      end do
      hv = 0.0_wp
      rho = 1025.0_wp                   ! neutral
      tke = 1.0e-6_wp
      eps = 1.0e-8_wp
      nut = 1.0e-4_wp
      kt = 1.0e-4_wp
      ustar_s = ustar
      ustar_b = ustar
      height_cv = real(nz, wp)*2.0_wp
      nsteps = nint(t_end/dt)

      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      do istep = 1, nsteps
         call ml_keps_update_unstr_impl(nz, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                        ustar_s, ustar_b, height_cv, 1025.0_wp, dt, &
                                        .false., KEPS_STAB_CONSTANT)
      end do
      !$omp end target data
      nut_out = nut
   end subroutine run_frozen_column

   ! ---------------------------------------------------------------------------
   ! Test 2: shock-off decay at dt=480 — positivity floors hold.
   ! ---------------------------------------------------------------------------
   subroutine test_keps_positivity_decay(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 10
      real(wp), parameter :: DT = 480.0_wp
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), rho(NZ, N)
      real(wp) :: tke(NZ + 1, N), eps(NZ + 1, N)
      real(wp) :: nut(NZ + 1, N), kt(NZ + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      real(wp) :: tke_spun
      integer :: k, istep

      h = 2.0_wp
      do k = 1, NZ
         hu(k, :) = h(k, :)*(0.05_wp*real(k, wp))
      end do
      hv = 0.0_wp
      rho = 1025.0_wp
      tke = 1.0e-6_wp
      eps = 1.0e-8_wp
      nut = 1.0e-4_wp
      kt = 1.0e-4_wp
      height_cv = real(NZ, wp)*2.0_wp

      !$omp target data map(to: h, rho, hv, height_cv) map(tofrom: hu, tke, eps, nut, kt, ustar_s, ustar_b)
      ! Spin up 20 large steps with strong wall forcing.
      ustar_s = 0.05_wp
      ustar_b = 0.05_wp
      !$omp target update to(ustar_s, ustar_b)
      do istep = 1, 20
         call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                        ustar_s, ustar_b, height_cv, 1025.0_wp, DT, &
                                        .false., KEPS_STAB_CONSTANT)
      end do
      !$omp target update from(tke)
      tke_spun = tke(NZ/2, 4)
      ! Shock off: kill the walls AND the shear; decay 100 large steps.
      ustar_s = 0.0_wp
      ustar_b = 0.0_wp
      hu = 0.0_wp
      !$omp target update to(ustar_s, ustar_b, hu)
      do istep = 1, 100
         call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                        ustar_s, ustar_b, height_cv, 1025.0_wp, DT, &
                                        .false., KEPS_STAB_CONSTANT)
      end do
      !$omp end target data

      call check(error, all(tke >= K_MIN) .and. all(eps >= EPS_MIN), &
                 "k / eps (unstr) dropped below positivity floors under shock-off decay")
      if (allocated(error)) return
      call check(error, tke(NZ/2, 4) < tke_spun, &
                 "TKE (unstr) failed to decay after the forcing shock-off")
   end subroutine test_keps_positivity_decay

   ! ---------------------------------------------------------------------------
   ! Shared mesh builder (tests 3, 4).
   ! ---------------------------------------------------------------------------
   subroutine build_rect_mesh(mesh)
      !! Jittered Delaunay-diagonal rectangle triangulation, same construction
      !! as test_si_unstr_ml_baro_limit.
      type(mesh_unstr_t), intent(out) :: mesh

      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      integer :: nvx, nv, ntri, i, j, vid, t
      real(wp) :: hx, hy, jx, jy

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

   ! ---------------------------------------------------------------------------
   ! Test 3: SI + k-eps wind channel: stability + tke/eps >= floors + tilt sign.
   ! ---------------------------------------------------------------------------
   subroutine test_si_keps_wind_channel(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 6
      real(wp), parameter :: TAU = 0.1_wp   ! wind stress [Pa]
      real(wp), parameter :: DT = 120.0_wp  ! SI Courant ~ 8
      integer, parameter :: N_STEPS = 30    ! 60 minutes
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, nco
      logical :: ok
      real(wp) :: tke_floor_ok, eps_floor_ok
      real(wp) :: eta_west, eta_east, xcc_west, xcc_east, xcc

      call build_rect_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      state%multilayer%vmix_use_keps = .true.
      state%multilayer%keps_stability = KEPS_STAB_CONSTANT
      state%multilayer%vmix_nu_max = 0.01_wp
      state%multilayer%vmix_nu_bg = 1.0e-4_wp
      ! EOS off (uniform density): neutral column, clean wind-setup test.
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp
      ! Small Manning drag so the wind-setup signal saturates quickly.
      state%barotropic%manning_n = 0.01_wp
      ! Wind forcing (uniform scalar).
      state%barotropic%wind_tau_x = TAU
      state%barotropic%wind_tau_y = 0.0_wp

      state%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state, mesh, SI_PRECOND_NONE, 1.0e-10_wp, DT, ok)
      if (.not. ok) then
         call check(error, .false., "si_keps_wind_channel_unstr: quality gate failed")
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if

      ! Flat bed + rest IC.
      do c = 1, mesh%ncells_total
         state%barotropic%bed_cv(c) = 0.0_wp
         state%barotropic%stage_cv(c) = H0
         state%barotropic%height_cv(c) = H0
         state%barotropic%xmom_cv(c) = 0.0_wp
         state%barotropic%ymom_cv(c) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_init_layers_unstr(state)
      call ml_seed_standard_tracers_unstr(state)

      block
         integer :: n
         do n = 1, N_STEPS
            call solver_unstr_step(state, mesh, DT)
         end do
      end block

      !$omp target update from(state%barotropic%stage_cv, &
      !$omp&            state%multilayer%keps_tke, state%multilayer%keps_eps)

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      nco = mesh%ncells_owned

      ! Positivity check: all owned cells, all interfaces.
      tke_floor_ok = 1.0_wp
      eps_floor_ok = 1.0_wp
      do c = 1, nco
         if (any(state%multilayer%keps_tke(:, c) < K_MIN)) tke_floor_ok = 0.0_wp
         if (any(state%multilayer%keps_eps(:, c) < EPS_MIN)) eps_floor_ok = 0.0_wp
      end do

      call global_logger%info("[si_keps_wind_channel_unstr] tke_ok=" &
                              //rstr(tke_floor_ok)//" eps_ok="//rstr(eps_floor_ok))
      call check(error, tke_floor_ok > 0.5_wp, &
                 "SI k-eps (unstr): tke below K_MIN floor after wind-channel run")
      if (allocated(error)) then
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if
      call check(error, eps_floor_ok > 0.5_wp, &
                 "SI k-eps (unstr): eps below EPS_MIN floor after wind-channel run")
      if (allocated(error)) then
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if

      ! Wind-setup tilt sign: eastward wind => stage piles up in the east
      ! (positive tilt: east > west).
      eta_west = 0.0_wp
      eta_east = 0.0_wp
      xcc_west = huge(1.0_wp)
      xcc_east = -huge(1.0_wp)
      do c = 1, nco
         xcc = mesh%circumcentre_x(c)
         if (xcc < xcc_west) then
            xcc_west = xcc
            eta_west = state%barotropic%stage_cv(c)
         end if
         if (xcc > xcc_east) then
            xcc_east = xcc
            eta_east = state%barotropic%stage_cv(c)
         end if
      end do
      call global_logger%info("[si_keps_wind_channel_unstr] tilt=" &
                              //rstr(eta_east - eta_west)//" (east-west)")
      call check(error, eta_east > eta_west, &
                 "SI k-eps (unstr): wind-setup tilt has wrong sign (east <= west)")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_si_keps_wind_channel

   ! ---------------------------------------------------------------------------
   ! Test 3b: keps_si_picard = 2 (unstr) — the Picard save slots allocate,
   ! map, and run; the two-pass steady state stays close to single-pass.
   ! IMPORTANT ORDERING (the driver contract this test guards): vmix_use_keps
   ! + keps_si_picard are set BEFORE state_unstr_init, because the save-slot
   ! allocation is gated on them at init (rki_driver syncs them pre-init for
   ! the same reason).
   ! ---------------------------------------------------------------------------
   subroutine test_si_keps_picard2(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 6
      real(wp), parameter :: TAU = 0.1_wp
      real(wp), parameter :: DT = 120.0_wp
      integer, parameter :: N_STEPS = 30
      real(wp) :: xm_p1, xm_p2, rel
      logical :: ok1, ok2, pos_ok

      call run_picard_channel(1, NZ, TAU, DT, N_STEPS, xm_p1, pos_ok, ok1)
      if (.not. ok1) then
         call check(error, .false., "si_keps_picard2_unstr: p1 quality gate failed")
         return
      end if
      call run_picard_channel(2, NZ, TAU, DT, N_STEPS, xm_p2, pos_ok, ok2)
      if (.not. ok2) then
         call check(error, .false., "si_keps_picard2_unstr: p2 quality gate failed")
         return
      end if

      call check(error, pos_ok, &
                 "picard=2 (unstr): tke/eps below floors after wind-channel run")
      if (allocated(error)) return

      rel = abs(xm_p2 - xm_p1)/max(abs(xm_p1), 1.0e-10_wp)
      call global_logger%info("[si_keps_picard2_unstr] xmom_p1="//rstr(xm_p1)// &
                              " xmom_p2="//rstr(xm_p2)//" rel="//rstr(rel))
      call check(error, rel < 0.05_wp, &
                 "picard=2 (unstr) surface momentum deviates from picard=1: "// &
                 rstr(rel))
   end subroutine test_si_keps_picard2

   subroutine run_picard_channel(picard, nz, tau, dt, n_steps, xmom_mid, &
                                 pos_ok, ok)
      !! One SI+keps wind-channel run at the given Picard count; returns the
      !! mid-domain surface-layer x-momentum + a positivity flag.  Closure
      !! knobs are set BEFORE state_unstr_init (the allocation-gating
      !! contract).
      integer, intent(in) :: picard, nz, n_steps
      real(wp), intent(in) :: tau, dt
      real(wp), intent(out) :: xmom_mid
      logical, intent(out) :: pos_ok, ok

      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, n, nco, c_mid
      real(wp) :: xcc, best

      call build_rect_mesh(mesh)
      ! Closure knobs BEFORE init: they gate the Picard save-slot allocation.
      state%multilayer%vmix_use_keps = .true.
      state%multilayer%keps_si_picard = picard
      call state_unstr_init(state, mesh, nz_ml=nz)

      state%multilayer%keps_stability = KEPS_STAB_CONSTANT
      state%multilayer%vmix_nu_max = 0.01_wp
      state%multilayer%vmix_nu_bg = 1.0e-4_wp
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp
      state%barotropic%manning_n = 0.01_wp
      state%barotropic%wind_tau_x = tau
      state%barotropic%wind_tau_y = 0.0_wp

      state%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state, mesh, SI_PRECOND_NONE, 1.0e-10_wp, dt, ok)
      if (.not. ok) then
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         xmom_mid = 0.0_wp
         pos_ok = .false.
         return
      end if

      do c = 1, mesh%ncells_total
         state%barotropic%bed_cv(c) = 0.0_wp
         state%barotropic%stage_cv(c) = H0
         state%barotropic%height_cv(c) = H0
         state%barotropic%xmom_cv(c) = 0.0_wp
         state%barotropic%ymom_cv(c) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_init_layers_unstr(state)
      call ml_seed_standard_tracers_unstr(state)

      do n = 1, n_steps
         call solver_unstr_step(state, mesh, dt)
      end do

      !$omp target update from(state%multilayer%hu_layer, &
      !$omp&            state%multilayer%keps_tke, state%multilayer%keps_eps)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      nco = mesh%ncells_owned
      ! Mid-domain cell: circumcentre closest to the basin centre.
      c_mid = 1
      best = huge(1.0_wp)
      do c = 1, nco
         xcc = abs(mesh%circumcentre_x(c) - 0.5_wp*LX)
         if (xcc < best) then
            best = xcc
            c_mid = c
         end if
      end do
      xmom_mid = state%multilayer%hu_layer(nz, c_mid)

      pos_ok = .true.
      do c = 1, nco
         if (any(state%multilayer%keps_tke(:, c) < K_MIN)) pos_ok = .false.
         if (any(state%multilayer%keps_eps(:, c) < EPS_MIN)) pos_ok = .false.
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine run_picard_channel

   ! ---------------------------------------------------------------------------
   ! Test 4: SI + KPP stratified basin: T within seeded envelope + S conservation.
   ! ---------------------------------------------------------------------------
   subroutine test_si_kpp_stratified(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 120.0_wp  ! SI step
      integer, parameter :: N_STEPS = 40    ! 80 minutes
      real(wp), parameter :: T_SEED_BOT = 10.0_wp, T_SEED_TOP = 20.0_wp
      real(wp), parameter :: S_SEED = 35.0_wp
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      integer :: c, k, nco, nz_ml
      logical :: ok
      real(wp) :: t_mean, s_tot, s_tot0, rel_s_drift
      real(wp) :: t_prof_k

      call build_rect_mesh(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      state%multilayer%vmix_use_kpp = .true.
      ! Scalar KPP params — defaults on the state are fine (kpp_ri_crit=0.3,
      ! cs_nonlocal=6.3, c_vt2=1.8).
      state%multilayer%vmix_nu_max = 0.01_wp
      state%multilayer%vmix_nu_bg = 1.0e-4_wp
      state%multilayer%kappa_S_bg = 1.0e-5_wp
      state%multilayer%kappa_T_bg = 1.0e-5_wp
      ! Linear EOS so T stratification creates density structure.
      state%multilayer%alpha_T = 0.17_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp
      state%barotropic%manning_n = 0.0_wp
      ! Light wind => small BL for stability verification.
      state%barotropic%wind_tau_x = 0.02_wp
      state%barotropic%wind_tau_y = 0.0_wp
      ! No surface heat/salt flux (closed basin, zero evap => S conserved).
      ! state%forcing%heat_flux allocated by state_unstr_init, default 0.
      ! evap_minus_precip default 0 — no need to set.

      state%time_scheme = TIME_SEMI_IMPLICIT
      call state_unstr_enable_si(state, mesh, SI_PRECOND_NONE, 1.0e-10_wp, DT, ok)
      if (.not. ok) then
         call check(error, .false., "si_kpp_stratified_unstr: quality gate failed")
         call state_unstr_destroy(state)
         call mesh_unstr_destroy(mesh)
         return
      end if

      ! Flat bed + rest IC.
      do c = 1, mesh%ncells_total
         state%barotropic%bed_cv(c) = 0.0_wp
         state%barotropic%stage_cv(c) = H0
         state%barotropic%height_cv(c) = H0
         state%barotropic%xmom_cv(c) = 0.0_wp
         state%barotropic%ymom_cv(c) = 0.0_wp
      end do

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_init_layers_unstr(state)
      call ml_seed_standard_tracers_unstr(state)

      ! Overwrite tracers with a linear T profile (hot-surface, cold-bed) and
      ! uniform S.  h_layer was filled ON DEVICE by ml_init_layers_unstr, so
      ! pull it to the host FIRST (on GPU builds the host copy is stale —
      ! seeding from it wrote hTr = 0, a GPU-only failure), then seed on the
      ! host and push the tracers back.  All acc directives use the
      ! associate-leaf pattern (deep `state%multilayer%tracers(i)%hTr`
      ! chains in update clauses are the known NVHPC/ifx hazard).
      nz_ml = state%multilayer%nz_ml
      nco = mesh%ncells_owned
      associate (hl_leaf => state%multilayer%h_layer)
         !$omp target update from(hl_leaf)
      end associate
      associate (ht_leaf => state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                 hs_leaf => state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                 hl => state%multilayer%h_layer)
         do c = 1, mesh%ncells_total
            do k = 1, nz_ml
               ! k=1 bed, k=nz surface (bottom-up convention).
               t_prof_k = T_SEED_BOT + (T_SEED_TOP - T_SEED_BOT) &
                          *(real(k, wp) - 0.5_wp)/real(nz_ml, wp)
               ht_leaf(k, c) = t_prof_k*hl(k, c)
               hs_leaf(k, c) = S_SEED*hl(k, c)
            end do
         end do
         !$omp target update to(ht_leaf, hs_leaf)
      end associate

      ! Record initial total salt (owned cells).  Host hTr/h_layer are
      ! current — just seeded/pulled above; no further update needed.
      s_tot0 = 0.0_wp
      do c = 1, nco
         do k = 1, nz_ml
            s_tot0 = s_tot0 + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(k, c) &
                     *mesh%area(c)
         end do
      end do

      block
         integer :: n
         do n = 1, N_STEPS
            call solver_unstr_step(state, mesh, DT)
         end do
      end block

      associate (ht_leaf => state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                 hs_leaf => state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                 hl_leaf => state%multilayer%h_layer)
         !$omp target update from(ht_leaf, hs_leaf, hl_leaf)
      end associate

      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)

      ! T must remain within [T_SEED_BOT - 1, T_SEED_TOP + 1] at any cell.
      t_mean = 0.0_wp
      do c = 1, nco
         do k = 1, nz_ml
            if (state%multilayer%h_layer(k, c) > 1.0e-4_wp) then
               t_mean = state%multilayer%tracers(state%multilayer%idx_temperature)%hTr(k, c) &
                        /state%multilayer%h_layer(k, c)
               if (t_mean < T_SEED_BOT - 1.0_wp .or. t_mean > T_SEED_TOP + 1.0_wp) then
                  call check(error, .false., &
                             "SI KPP (unstr): T out of seeded envelope: T=" &
                             //rstr(t_mean)//" at k="//rint(k)//" c="//rint(c))
                  call state_unstr_destroy(state)
                  call mesh_unstr_destroy(mesh)
                  return
               end if
            end if
         end do
      end do

      ! Column salt conservation (zero surface flux => S total conserved).
      s_tot = 0.0_wp
      do c = 1, nco
         do k = 1, nz_ml
            s_tot = s_tot + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(k, c) &
                    *mesh%area(c)
         end do
      end do
      rel_s_drift = abs(s_tot - s_tot0)/max(abs(s_tot0), 1.0_wp)
      call global_logger%info("[si_kpp_stratified_unstr] rel_S_drift="//rstr(rel_s_drift))
      call check(error, rel_s_drift < 1.0e-10_wp, &
                 "SI KPP (unstr): column salt drift "//rstr(rel_s_drift)// &
                 " exceeds 1e-10 (closed basin, zero evap)")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_si_kpp_stratified

   ! ---------------------------------------------------------------------------
   ! Test 5: the implicit KPP tracer kernel conserves column S/T content.
   ! Kernel-level on a flat (nz, ncells) column set — no mesh.
   ! ---------------------------------------------------------------------------
   subroutine test_si_kpp_nonlocal(error)
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NZ = 8
      real(wp), parameter :: DT = 480.0_wp
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N)
      real(wp) :: hS(NZ, N), hT(NZ, N), rho(NZ, N)
      real(wp) :: heat(N), evap(N)
      real(wp) :: ustar(N), bflux(N), bld(N)
      real(wp) :: s_col0(N), t_col0(N)
      real(wp) :: s_col1, t_col1, s_err, t_err
      integer :: c, k

      do k = 1, NZ
         h(k, :) = 2.0_wp
         hu(k, :) = h(k, :)*(0.05_wp*real(k, wp))
         hS(k, :) = (34.0_wp + 0.2_wp*real(k, wp))*h(k, :)
         hT(k, :) = (12.0_wp + 0.8_wp*real(k, wp))*h(k, :)
         rho(k, :) = 1026.0_wp - 0.1_wp*real(k, wp)
      end do
      hv = 0.0_wp
      heat = -300.0_wp                  ! strong cooling: nonlocal active
      evap = 2.0e-8_wp                  ! nonzero F_S
      ustar = 0.02_wp
      bflux = -5.0e-7_wp                ! destabilizing
      bld = 10.0_wp                     ! BL covers the upper 5 layers

      do c = 1, N
         s_col0(c) = sum(hS(:, c))
         t_col0(c) = sum(hT(:, c))
      end do

      !$omp target data map(to: h, hu, hv, rho, heat, evap, ustar, bflux, bld) &
      !$omp&     map(tofrom: hS, hT)
      call si_unstr_kpp_diffusion_tracers(NZ, NZ + 1, N, N, h, hu, hv, hS, hT, &
                                          rho, 1025.0_wp, 0.01_wp, 1.0e-5_wp, &
                                          1.0e-5_wp, 5.0_wp, 6.3_wp, heat, evap, &
                                          CP_WATER, ustar, bflux, bld, DT)
      !$omp end target data

      s_err = 0.0_wp
      t_err = 0.0_wp
      do c = 1, N
         s_col1 = sum(hS(:, c))
         t_col1 = sum(hT(:, c))
         s_err = max(s_err, abs(s_col1 - s_col0(c))/s_col0(c))
         t_err = max(t_err, abs(t_col1 - t_col0(c))/t_col0(c))
      end do
      call check(error, s_err < 1.0e-12_wp .and. t_err < 1.0e-12_wp, &
                 "SI KPP tracer kernel (unstr, diffusion+nonlocal) leaks column "// &
                 "content: S "//rstr(s_err)//", T "//rstr(t_err))
   end subroutine test_si_kpp_nonlocal

   ! ---------------------------------------------------------------------------
   ! Formatting helpers.
   ! ---------------------------------------------------------------------------
   function rstr(x) result(str)
      real(wp), intent(in) :: x
      character(len=32) :: str
      write (str, "(es12.4)") x
   end function rstr

   function rint(i) result(str)
      integer, intent(in) :: i
      character(len=16) :: str
      write (str, "(i0)") i
   end function rint

end module test_si_closures_unstr
