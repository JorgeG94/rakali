!! Carrier-Greenspan periodic runup on a plane beach.
module test_carrier_greenspan
   !! Thacker's parabolic bowl is currently the ONLY moving-shoreline
   !! analytical test in the suite, and its planar solution is asserted at
   !! 8 %. This adds the canonical periodic-runup case on a uniform slope.
   !!
   !! ## Setup — a self-forcing standing eigenmode
   !!
   !! Plane beach of slope `s`, offshore wall at the first zero of `J1` so the
   !! wall is an exact velocity node (see `carrier_greenspan`). The J0 mode is
   !! then an eigenmode of the closed basin: seed `eta = A0 J0(xi)`, `u = 0`,
   !! and it oscillates at `w` with the shoreline running up and down the
   !! beach, with **no boundary forcing at all**.
   !!
   !! That is deliberate. Driving the offshore end with a time-dependent
   !! tidal boundary would entangle a runup measurement with the BC
   !! implementation, so a failure would not localise. A closed basin also
   !! makes total mass exactly invariant, which turns the wet/dry cycle into a
   !! conservation assertion for free — and this is the only fixture in the
   !! suite where water repeatedly wets and dries the same cells.
   !!
   !! ## Outcome — what is and is not validated
   !!
   !! **Validated.** The basin oscillates at the Carrier-Greenspan
   !! frequency to within 5 % (measured 1.5 %), total mass is conserved to
   !! 1e-8 through repeated wet/dry cycling, and the shoreline excursion at
   !! detection depths above the retention film (10 cm and 30 cm) is within
   !! 12 % of the exact 20 m and reproducible: 21.5-22.0 m and 18.6-18.75 m
   !! across gfortran -O3, nvfortran -O3 -fast and nvfortran -O0, and across
   !! initial-depth perturbations from 1e-15 to 1e-8 relative.
   !!
   !! **Not validated: the 1-cm detector.**  Its excursion is chaotic at
   !! round-off level: a 1e-15 relative change of the initial depth moves it
   !! from 18.9 m to 17.3 m, 1e-14 gives 16.5 m, 1e-10 gives 20.8 m (gfortran).
   !! The same run measures 17.0 m (nvfortran -O3 -fast), 18.9 m (gfortran)
   !! and 20.9 m (nvfortran -O0).  The detector sits on the ~slope*dx = 12.5 mm
   !! retention film (see below), so a single run is a draw from a spread of
   !! roughly +-2 m and says nothing about the scheme.  Earlier versions
   !! gated on it ("within 10 %", measured 5.6 % on gcc) and failed on
   !! nvfortran only because of that spread, not because of a solver defect;
   !! it is now held to a coarse sanity band that still catches a frozen or
   !! blown-up shoreline.  The exact A0/slope is the h -> 0 contour, so the
   !! finite-depth contours above bracket it (+7.5..+10 % and -6.3..-7 %)
   !! rather than hit it: this is a consistency check, not the 5 % target.
   !!
   !! **Finding (residual) — a millimetre-scale film still pins the 1-mm
   !! detector.**  Pre-fix, water below ~1e-2 m froze outright (momentum
   !! clamps).  Post-fix the film drains down to the hydrostatic-
   !! reconstruction retention scale: with Audusse-type well-balancing a
   !! plane beach retains up to slope*dx = 12.5 mm per cell (measured ~5 mm
   !! mean on the static drain probe), and the freed runup jet overwets the
   !! upper beach (a >1 mm film reaches x ~ 128.4, beyond the exact 120 m
   !! maximum runup), so the 1-mm detector reads an excursion of exactly
   !! zero, pinned at the overshoot high-water mark.  This is a
   !! DISCRETIZATION retention (grid-dependent, ~slope*dx), no longer a
   !! momentum-clamp freeze; it still matters for intertidal stranding and
   !! is asserted two-sidedly below.
   !!
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, PI
   use rki_state, only: state_t
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_init, only: compute_total_mass
   use rki_boundary_types, only: bc_type_from_string
   use pic_logger, only: global_logger
   use carrier_greenspan, only: cg_eta, cg_omega_for_wall, cg_breaking_number, &
                                CG_J1_ZERO1
   implicit none
   private

   public :: collect_carrier_greenspan_tests

   ! ---- Beach ----
   real(wp), parameter :: SLOPE = 0.05_wp          !! 1:20
   real(wp), parameter :: X_SHORE = 100.0_wp
      !! Still-water shoreline, measured from the offshore wall at x = 0.
      !! Depth at the wall is SLOPE*X_SHORE = 5 m.
   real(wp), parameter :: LX = 130.0_wp
      !! Domain length. Max runup reaches x = X_SHORE + A0/SLOPE = 120 m.
      !! Post-fix, a millimetre-scale splash film overshoots to x ~ 128.4,
      !! so the beach is film-wetted almost to the east wall — but the
      !! DYNAMIC shoreline (>= 1 cm) stays below ~121.6 m, well clear.
   real(wp), parameter :: A0 = 1.0_wp
      !! Vertical runup amplitude => 20 m horizontal excursion. Chosen large
      !! enough that the P4 thin-front deficit (~0.03 m depth, i.e. ~0.6 m
      !! horizontal on this slope) would be a few per cent rather than
      !! dominant — and small enough to stay well clear of breaking.

   integer, parameter :: NX = 520                  !! dx = 0.25 m
   integer, parameter :: NY = 4
   integer, parameter :: NGHOST = 2

   integer, parameter :: N_PERIODS = 3
   integer, parameter :: SAMPLES_PER_PERIOD = 400

   real(wp), parameter :: H_SHORE = 1.0e-2_wp
      !! Depth defining the shoreline. The threshold-sensitivity test sweeps
      !! this; the amplitude assertions use it.

   ! ---- Tolerances ----
   real(wp), parameter :: RUNUP_DEEP_A = 1.0e-1_wp
   real(wp), parameter :: RUNUP_DEEP_B = 3.0e-1_wp
      !! The two detection depths the runup amplitude is gated on.  Both sit
      !! well above the ~slope*dx = 12.5 mm retention film, where the
      !! excursion is reproducible: 21.5-22.0 m (10 cm) and 18.6-18.75 m
      !! (30 cm) across gfortran -O3, nvfortran -O3 -fast and nvfortran -O0,
      !! and across initial-depth perturbations from 1e-15 to 1e-8.
   real(wp), parameter :: RUNUP_GUARD = 0.12_wp
      !! Relative tolerance on the deep-detector excursions against the exact
      !! shoreline excursion A0/slope = 20 m.  The exact value is the h->0
      !! contour; a contour at finite depth h* is displaced by O(h*/slope), so
      !! this is a consistency band (measured +7.5..+10 % at 10 cm, -6.3..-7 %
      !! at 30 cm), NOT the 5 % target.  See the module docstring.
   real(wp), parameter :: SHORE_SANE_LO = 0.65_wp
   real(wp), parameter :: SHORE_SANE_HI = 1.35_wp
      !! Coarse sanity band, as a fraction of A0/slope, for the 1-cm detector.
      !! Not a validation: that detector sits on the retention film and its
      !! excursion is chaotic (16.5-20.9 m from 1e-15 relative noise alone),
      !! so only a frozen or blown-up shoreline is detectable there.
   real(wp), parameter :: PERIOD_TOL = 0.05_wp     !! 5 % on oscillation period
   real(wp), parameter :: MASS_TOL = 1.0e-8_wp

contains

   subroutine collect_carrier_greenspan_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cg_runup_amplitude", test_runup_amplitude), &
                  new_unittest("cg_shoreline_period", test_period), &
                  new_unittest("cg_mass_conserved_through_wetdry", test_mass), &
                  new_unittest("cg_thin_film_persists_on_beach", test_threshold), &
                  new_unittest("cg_fixture_is_nonbreaking_and_moves", test_fixture) &
                  ]
   end subroutine collect_carrier_greenspan_tests

   pure function omega() result(w)
      real(wp) :: w
      w = cg_omega_for_wall(X_SHORE, SLOPE, CG_J1_ZERO1)
   end function omega

   pure function bed_at(x) result(b)
      !! Bed rises shoreward; zero at the still-water shoreline.
      real(wp), intent(in) :: x
      real(wp) :: b
      b = SLOPE*(x - X_SHORE)
   end function bed_at

   subroutine setup(state)
      type(state_t), intent(inout) :: state

      integer :: i, j
      real(wp) :: dx, x, b, eta

      dx = LX/real(NX, wp)
      call state%init(NX, NY, NGHOST, dx, dx, 0.45_wp)
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            x = (real(i - NGHOST, wp) - 0.5_wp)*dx
            b = bed_at(x)
            ! Seed the eigenmode at t = 0, where cos(w t) = 1 and u = 0.
            eta = cg_eta(A0, X_SHORE - x, 0.0_wp, omega(), SLOPE)
            state%barotropic%b(i, j) = b
            state%barotropic%h(i, j) = max(eta - b, 0.0_wp)
            state%barotropic%hu(i, j) = 0.0_wp
            state%barotropic%hv(i, j) = 0.0_wp
         end do
      end do
   end subroutine setup

   subroutine shoreline_at(state, h_thresh, x_s)
      !! Shoreline = the most shoreward interior cell still wetter than
      !! `h_thresh`.
      type(state_t), intent(in) :: state
      real(wp), intent(in) :: h_thresh
      real(wp), intent(out) :: x_s

      integer :: i, j_mid
      real(wp) :: dx

      dx = LX/real(NX, wp)
      j_mid = NGHOST + NY/2
      x_s = 0.0_wp
      do i = NGHOST + 1, state%grid%nx_total - NGHOST
         if (state%barotropic%h(i, j_mid) > h_thresh) then
            x_s = (real(i - NGHOST, wp) - 0.5_wp)*dx
         end if
      end do
   end subroutine shoreline_at

   subroutine run_and_sample_multi(h_thresh, x_track, eta_track, t_track, n_track, m_rel)
      !! Run N_PERIODS and record, on a fixed sample cadence, the shoreline
      !! position for EACH detection depth in `h_thresh` (one simulation, many
      !! detectors), plus the relative mass change over the whole run.
      real(wp), intent(in) :: h_thresh(:)
      real(wp), intent(out) :: x_track(:, :)
         !! `x_track(k, d)`: shoreline of detector `d` at sample `k`.
      real(wp), intent(out) :: t_track(:)
      real(wp), intent(out) :: eta_track(:)
         !! Free-surface elevation at the first interior cell — the offshore
         !! end, in ~5 m of water, which never wets or dries.
      integer, intent(out) :: n_track
      real(wp), intent(out) :: m_rel

      type(state_t) :: state
      real(wp) :: period, t_end, dt_sample, t, m0, m1
      integer :: k, d, n_samples

      call setup(state)
      m0 = compute_total_mass(state)

      period = 2.0_wp*PI/omega()
      t_end = real(N_PERIODS, wp)*period
      n_samples = N_PERIODS*SAMPLES_PER_PERIOD
      dt_sample = t_end/real(n_samples, wp)

      call solver_enter_data(state)
      t = 0.0_wp
      n_track = 0
      do k = 1, n_samples
         do while (t < real(k, wp)*dt_sample)
            call solver_step(state, dt_max=real(k, wp)*dt_sample - t, time=t)
            t = t + state%dt
         end do
         !$acc update self(state%barotropic%h)
         n_track = n_track + 1
         t_track(n_track) = t
         do d = 1, size(h_thresh)
            call shoreline_at(state, h_thresh(d), x_track(n_track, d))
         end do
         eta_track(n_track) = state%barotropic%h(NGHOST + 1, NGHOST + NY/2) &
                              + state%barotropic%b(NGHOST + 1, NGHOST + NY/2)
      end do
      call solver_exit_data(state)

      m1 = compute_total_mass(state)
      m_rel = abs(m1 - m0)/max(m0, tiny(1.0_wp))

      call state%destroy()
   end subroutine run_and_sample_multi

   subroutine run_and_sample(h_thresh, x_track, eta_track, t_track, n_track, m_rel)
      !! Single-detector convenience wrapper around `run_and_sample_multi`.
      real(wp), intent(in) :: h_thresh
      real(wp), intent(out) :: x_track(:), t_track(:)
      real(wp), intent(out) :: eta_track(:)
      integer, intent(out) :: n_track
      real(wp), intent(out) :: m_rel

      real(wp) :: x2(size(x_track), 1)

      call run_and_sample_multi([h_thresh], x2, eta_track, t_track, n_track, m_rel)
      x_track = x2(:, 1)
   end subroutine run_and_sample

   subroutine excursion_over_last_period(x_track, n_track, x_min, x_max)
      !! Min/max shoreline over the FINAL period only, so any startup
      !! transient from seeding a discrete grid with a continuous eigenmode
      !! has been left behind.
      real(wp), intent(in) :: x_track(:)
      integer, intent(in) :: n_track
      real(wp), intent(out) :: x_min, x_max

      integer :: k, k0

      k0 = n_track - SAMPLES_PER_PERIOD + 1
      x_min = huge(1.0_wp)
      x_max = -huge(1.0_wp)
      do k = max(k0, 1), n_track
         x_min = min(x_min, x_track(k))
         x_max = max(x_max, x_track(k))
      end do
   end subroutine excursion_over_last_period

   subroutine test_runup_amplitude(error)
      !! Runup excursion vs the exact Carrier-Greenspan value A0/slope = 20 m,
      !! measured with detectors above the retention film (10 cm and 30 cm),
      !! each within RUNUP_GUARD.  The 1-cm detector is only held to a coarse
      !! sanity band — it is chaotic at round-off level (module docstring).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x_track(N_PERIODS*SAMPLES_PER_PERIOD, 3)
      real(wp) :: eta_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: t_track(N_PERIODS*SAMPLES_PER_PERIOD)
      integer :: n_track, d
      real(wp) :: m_rel, x_min, x_max, amp(3), amp_exact
      real(wp), parameter :: DEPTHS(3) = [RUNUP_DEEP_A, RUNUP_DEEP_B, H_SHORE]
      character(len=96) :: msg

      call run_and_sample_multi(DEPTHS, x_track, eta_track, t_track, n_track, m_rel)
      amp_exact = A0/SLOPE
      do d = 1, 3
         call excursion_over_last_period(x_track(:, d), n_track, x_min, x_max)
         amp(d) = 0.5_wp*(x_max - x_min)
         write (msg, '(A,ES9.2,A,F8.3,A,F8.3,A,F8.3,A)') "[cg runup] detector ", DEPTHS(d), &
            " m: x_min=", x_min, " x_max=", x_max, " amp=", amp(d), " m"
         call global_logger%info(trim(msg))
      end do

      call check(error, abs(amp(1) - amp_exact)/amp_exact < RUNUP_GUARD, &
                 "10-cm runup excursion outside 12% of the exact Carrier-Greenspan value")
      if (allocated(error)) return
      call check(error, abs(amp(2) - amp_exact)/amp_exact < RUNUP_GUARD, &
                 "30-cm runup excursion outside 12% of the exact Carrier-Greenspan value")
      if (allocated(error)) return
      call check(error, amp(3) > SHORE_SANE_LO*amp_exact .and. amp(3) < SHORE_SANE_HI*amp_exact, &
                 "1-cm shoreline excursion outside the coarse sanity band (frozen or blown-up shoreline)")
   end subroutine test_runup_amplitude

   subroutine test_period(error)
      !! Oscillation period from the shoreline track, via mean spacing of
      !! upward zero-crossings about the excursion midpoint.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: eta_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: t_track(N_PERIODS*SAMPLES_PER_PERIOD)
      integer :: n_track, k, n_cross
      real(wp) :: m_rel, t_first, t_last
      real(wp) :: period_num, period_exact, rel, eta_mean

      call run_and_sample(H_SHORE, x_track, eta_track, t_track, n_track, m_rel)

      ! Measure the frequency from the OFFSHORE ELEVATION, not the shoreline.
      !
      ! The shoreline track is a cell-index staircase whose rundown extent is
      ! itself affected by the undrained-film defect, so any crossing count
      ! built on it is fragile: it gave 3 crossings on gcc and 2 on nvfortran
      ! — a factor-of-two period difference from a measurement artefact, not a
      ! physics difference. The offshore cell sits in ~5 m of water, never
      ! wets or dries, and carries a clean sinusoid, so simple upward
      ! zero-crossings about the mean suffice and both toolchains agree.
      !
      ! This still tests the physics the phase is after: `omega` is fixed by
      ! the Carrier-Greenspan dispersion relation for the beach, so the basin
      ! oscillating at that frequency is the eigenmode being verified.
      eta_mean = sum(eta_track(1:n_track))/real(n_track, wp)
      n_cross = 0
      t_first = 0.0_wp
      t_last = 0.0_wp
      do k = 2, n_track
         if (eta_track(k - 1) <= eta_mean .and. eta_track(k) > eta_mean) then
            n_cross = n_cross + 1
            if (n_cross == 1) t_first = t_track(k)
            t_last = t_track(k)
         end if
      end do

      call check(error, n_cross >= 2, &
                 "fewer than two shoreline crossings — cannot measure a period")
      if (allocated(error)) return

      period_num = (t_last - t_first)/real(n_cross - 1, wp)
      period_exact = 2.0_wp*PI/omega()
      rel = abs(period_num - period_exact)/period_exact

      call check(error, rel < PERIOD_TOL, &
                 "Carrier-Greenspan shoreline period outside 5%")
   end subroutine test_period

   subroutine test_mass(error)
      !! The basin is closed, so total water mass is invariant through every
      !! wet/dry cycle. This is the only fixture in the suite that repeatedly
      !! wets and dries the same cells.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: eta_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: t_track(N_PERIODS*SAMPLES_PER_PERIOD)
      integer :: n_track
      real(wp) :: m_rel

      call run_and_sample(H_SHORE, x_track, eta_track, t_track, n_track, m_rel)

      call check(error, m_rel < MASS_TOL, &
                 "mass not conserved through the runup wet/dry cycle")
   end subroutine test_mass

   subroutine test_threshold(error)
      !! CHARACTERISATION — a millimetre-scale film is left on the beach.
      !!
      !! Locating the shoreline at 1e-3 m still gives an excursion of
      !! EXACTLY ZERO after the thin-cell fix — but the mechanism changed.
      !! Pre-fix, momentum clamps froze all water below ~1e-2 m at the
      !! high-water mark.  Post-fix the beach drains down to the
      !! hydrostatic-reconstruction retention scale (~slope*dx = 12.5 mm
      !! per cell here) and the freed runup jet leaves an overshoot film,
      !! so the 1-mm detector is pinned at the overshoot mark (~128.4 m)
      !! rather than the high-water mark.  A moving shoreline appears from
      !! ~1e-2 m and deeper (excursion 16.5-20.9 m at 1 cm, chaotic).
      !!
      !! This residual is grid-dependent discretization retention, not a
      !! clamp: it still matters for intertidal work (stranded water each
      !! cycle, consistent with the intertidal salt leak recorded in
      !! `docs/CAPABILITIES_AND_LIMITATIONS.md`).
      !!
      !! Asserted two-sidedly, like the P4 tip guard: the film must still
      !! be there.  If someone adds a thin-film drainage treatment (e.g.
      !! bed-slope-aware reconstruction below the retention scale) this
      !! fails, and that failure is the prompt to retire the guard and
      !! tighten `cg_runup_amplitude` onto the real 5 % target.
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: eta_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: t_track(N_PERIODS*SAMPLES_PER_PERIOD)
      integer :: n_track, m
      real(wp) :: amp(3)
      real(wp) :: m_rel, x_min, x_max, thr, spread
      real(wp), parameter :: THR_LIST(3) = [3.0e-1_wp, 3.0e-2_wp, 1.0e-3_wp]

      do m = 1, 3
         thr = THR_LIST(m)
         call run_and_sample(thr, x_track, eta_track, t_track, n_track, m_rel)
         call excursion_over_last_period(x_track, n_track, x_min, x_max)
         amp(m) = 0.5_wp*(x_max - x_min)
      end do

      spread = (maxval(amp) - minval(amp))/max(maxval(amp), tiny(1.0_wp))

      ! amp(3) is measured at 1e-3 m, where the shoreline is pinned by the
      ! undrained film and the excursion collapses to zero.
      call check(error, amp(3) < 1.0_wp, &
                 "thin film now drains off the beach — the wet/dry deficit "// &
                 "appears fixed; retire this guard and tighten the runup "// &
                 "amplitude assertion onto its 5% target")
      if (allocated(error)) return
      call check(error, amp(1) > 1.0_wp, &
                 "no moving shoreline even at a deep detection threshold")
   end subroutine test_threshold

   subroutine test_fixture(error)
      !! Guard the premises: the mode must be non-breaking (or the
      !! Carrier-Greenspan oracle does not apply), and the shoreline must
      !! actually move (or every assertion above is vacuous).
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: x_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: eta_track(N_PERIODS*SAMPLES_PER_PERIOD)
      real(wp) :: t_track(N_PERIODS*SAMPLES_PER_PERIOD)
      integer :: n_track
      real(wp) :: m_rel, x_min, x_max, bn

      bn = cg_breaking_number(A0, omega(), SLOPE)
      call check(error, bn < 1.0_wp, &
                 "fixture: mode is breaking, the Carrier-Greenspan oracle does not apply")
      if (allocated(error)) return

      call run_and_sample(H_SHORE, x_track, eta_track, t_track, n_track, m_rel)
      call excursion_over_last_period(x_track, n_track, x_min, x_max)

      call check(error, (x_max - x_min) > 1.0_wp, &
                 "fixture: shoreline barely moved, runup assertions are vacuous")
   end subroutine test_fixture

end module test_carrier_greenspan
