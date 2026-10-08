!! WENO tracer face-reconstruction ladder for the structured multilayer solver.
!!
!! Provides a drop-in replacement for PPM tracer advection, offering a
!! rung-adaptive reconstruction that automatically degrades to lower-order
!! schemes near boundaries where the full stencil is unavailable.
!!
!! Schemes (all with Z-weights, Borges et al. 2008):
!!   WENO5  — Jiang & Shu (1996),  3 quadratic candidates, nghost >= 3
!!   WENO7  — Balsara & Shu (2000), 4 cubic     candidates, nghost >= 4
!!   WENO9  — Balsara & Shu (2000), 5 quartic   candidates, nghost >= 5
!!
!! Caveat on WENO9 order: this rung uses simplified 3-term smoothness
!! indicators rather than the full quartic betas, so its effective
!! convergence is ~order 6 (not the nominal 9), and at a sharp front it
!! overshoots slightly MORE than WENO7 — WENO7 gives the best front L1.
!! Do not assume weno9 > weno7. See the function docstring and
!! docs/REFERENCE.md for the measured comparison.
!!
!! At interior cells far from boundaries the requested scheme is used.
!! As a face approaches the boundary the stencil is automatically degraded:
!!   WENO9 → WENO7 → WENO5 → PLM → donor-cell
!!
!! The swept-average form is used throughout: for upwind cell c with
!! CFL sigma, the face value integrates the reconstruction over
!! [c + 1/2 - sigma, c + 1/2], giving a time-averaged flux consistent
!! with the subcycled PPM approach in the rest of the solver.
module rki_recon_weno
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, DRY_TOLERANCE
#else
   use rki_constants, only: NZ_STACK_MAX, wp, DRY_TOLERANCE
#endif
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH, VANISHING_LAYER_TOL
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   ! ---- Public integer parameters (rung codes) ----
   integer, parameter, public :: TRACER_RECON_PPM = 0
      !! Use PPM (the default path — does not call this module).
   integer, parameter, public :: TRACER_RECON_WENO5 = 1
      !! WENO5-Z (Jiang-Shu 1996 + Borges 2008 Z-weights), nghost >= 3.
   integer, parameter, public :: TRACER_RECON_WENO7 = 2
      !! WENO7-Z (Balsara-Shu 2000), nghost >= 4.
   integer, parameter, public :: TRACER_RECON_WENO9 = 3
      !! WENO9-Z (Balsara-Shu 2000), nghost >= 5.

   ! ---- Support-status codes (tracer_recon_support_status) ----
   integer, parameter, public :: TRACER_RECON_OK = 0
      !! Config honours the requested tracer_recon code.
   integer, parameter, public :: TRACER_RECON_REJECT_PQM = 1
      !! "pqm" requested — deferred, not implemented on any config.
   integer, parameter, public :: TRACER_RECON_REJECT_UNKNOWN = 2
      !! Unrecognised tracer_recon string.
   integer, parameter, public :: TRACER_RECON_REJECT_OCEAN = 3
      !! weno* requested on sim_type='ocean' (coastal-only in wave 1).
   integer, parameter, public :: TRACER_RECON_REJECT_NO_ML = 4
      !! weno* requested on a coastal run without multilayer tracers.

   ! ---- Public procedures ----
   public :: parse_tracer_recon
   public :: tracer_recon_required_nghost
   public :: tracer_recon_support_status
   public :: recon_rung_for_face
   public :: ml_advect_tracer_weno
   ! Face helpers exported for the polynomial-exactness unit tests
   ! (test_tracer_recon); production code calls them only from the kernel.
   public :: plm_face_swept
   public :: weno5_face_swept
   public :: weno7_face_swept
   public :: weno9_face_swept

   ! The per-call concentration snapshot is caller-owned scratch
   ! (`tracer_scratch_t%weno_tr_snap`), passed in as an `intent(out)` dummy
   ! so `ml_advect_tracer_weno` stays `pure` — no module state, no per-step
   ! allocation (mirrors ml_advect_tracer_ppm; see the #315 refactor).

contains

   ! ==========================================================================
   !  Utility: parse / query
   ! ==========================================================================

   pure function parse_tracer_recon(str) result(code)
      !! Map a namelist string onto a TRACER_RECON_* code.
      !!
      !! Returns:
      !!   TRACER_RECON_PPM  (0) for "ppm"
      !!   TRACER_RECON_WENO5(1) for "weno5"
      !!   TRACER_RECON_WENO7(2) for "weno7"
      !!   TRACER_RECON_WENO9(3) for "weno9"
      !!   -2                    for "pqm" (deferred)
      !!   -1                    for anything else (unrecognised)
      character(len=*), intent(in) :: str
      integer :: code
      select case (trim(str))
      case ("ppm")
         code = TRACER_RECON_PPM
      case ("weno5")
         code = TRACER_RECON_WENO5
      case ("weno7")
         code = TRACER_RECON_WENO7
      case ("weno9")
         code = TRACER_RECON_WENO9
      case ("pqm")
         code = -2
      case default
         code = -1
      end select
   end function parse_tracer_recon

   pure function tracer_recon_required_nghost(code) result(ng)
      !! Minimum nghost required for a given TRACER_RECON_* code.
      integer, intent(in) :: code
         !! TRACER_RECON_* constant.
      integer :: ng
      select case (code)
      case (TRACER_RECON_WENO5)
         ng = 3
      case (TRACER_RECON_WENO7)
         ng = 4
      case (TRACER_RECON_WENO9)
         ng = 5
      case default
         ng = 2  ! PPM or unrecognised — PPM requires nghost >= 2
      end select
   end function tracer_recon_required_nghost

   pure function tracer_recon_support_status(recon_code, is_ocean, &
                                             has_multilayer_tracers) result(status)
      !! Decide whether a parsed tracer_recon code is honoured on this
      !! config, independent of the multilayer nghost/use_ppm_tracer
      !! quality guards (those apply only once this returns OK for a weno
      !! rung).  Factored out as a pure predicate so the fail-loud
      !! decision is unit-testable without triggering error stop.
      !!
      !! `recon_code` is the parse_tracer_recon result (may be < 0).
      !! The coastal multilayer solver is the ONLY consumer, so a weno*
      !! rung is rejected fail-loud on every other config:
      !!   * "pqm"        → REJECT_PQM      (deferred on all configs)
      !!   * unrecognised → REJECT_UNKNOWN
      !!   * sim_type='ocean' → REJECT_OCEAN  (coastal-only in wave 1)
      !!   * coastal, not multilayer → REJECT_NO_ML (2D barotropic / NH)
      !! "ppm" is always OK (the default path; bit-identical).
      integer, intent(in) :: recon_code
         !! parse_tracer_recon(cfg%tracer_recon) result.
      logical, intent(in) :: is_ocean
         !! .true. when sim_type == 'ocean'.
      logical, intent(in) :: has_multilayer_tracers
         !! cfg%use_multilayer — the coastal multilayer tracer registry.
      integer :: status
      if (recon_code == -2) then
         status = TRACER_RECON_REJECT_PQM
      else if (recon_code < 0) then
         status = TRACER_RECON_REJECT_UNKNOWN
      else if (recon_code == TRACER_RECON_PPM) then
         status = TRACER_RECON_OK
      else if (is_ocean) then
         status = TRACER_RECON_REJECT_OCEAN
      else if (.not. has_multilayer_tracers) then
         status = TRACER_RECON_REJECT_NO_ML
      else
         status = TRACER_RECON_OK
      end if
   end function tracer_recon_support_status

   ! ==========================================================================
   !  Rung-degradation helper
   ! ==========================================================================

   pure function recon_rung_for_face(avail_up, avail_down, rung_max) result(rung)
      !! Return the highest feasible reconstruction rung for one face.
      !!
      !! Internal rung scale (higher = more accurate):
      !!   4 = WENO9 (needs avail_up >= 5, avail_down >= 4)
      !!   3 = WENO7 (needs avail_up >= 4, avail_down >= 3)
      !!   2 = WENO5 (needs avail_up >= 3, avail_down >= 2)
      !!   1 = PLM   (needs avail_up >= 2, avail_down >= 1)
      !!   0 = donor (needs avail_up >= 1)
      !!  -1 = no data (should never happen at interior cells)
      !!
      !! rung_max is the PUBLIC TRACER_RECON_* code plus one:
      !!   WENO5(1) → rung_max_internal = 2
      !!   WENO7(2) → rung_max_internal = 3
      !!   WENO9(3) → rung_max_internal = 4
      !! The kernel maps the public code to rung_max with `recon + 1`.
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      integer, intent(in) :: avail_up
         !! Interior cells from domain-west to face (including upwind cell c).
      integer, intent(in) :: avail_down
         !! Interior cells from face downwind side to domain-east.
      integer, intent(in) :: rung_max
         !! Maximum internal rung allowed (4 = WENO9, 3 = WENO7, 2 = WENO5).
      integer :: rung
      if (rung_max >= 4 .and. avail_up >= 5 .and. avail_down >= 4) then
         rung = 4
      else if (rung_max >= 3 .and. avail_up >= 4 .and. avail_down >= 3) then
         rung = 3
      else if (rung_max >= 2 .and. avail_up >= 3 .and. avail_down >= 2) then
         rung = 2
      else if (avail_up >= 2 .and. avail_down >= 1) then
         rung = 1
      else if (avail_up >= 1) then
         rung = 0
      else
         rung = -1
      end if
   end function recon_rung_for_face

   ! ==========================================================================
   !  Low-level face helpers (pure, !$acc routine seq)
   ! ==========================================================================

   pure function plm_face_swept(qm1, q0, qp1, sigma) result(face)
      !! PLM (piecewise-linear) swept-average face value.
      !!
      !! For u > 0 the upwind cell is q0; stencil: qm1 = i-1, q0 = i, qp1 = i+1.
      !! slope = minmod(q0 - qm1, qp1 - q0)
      !! face  = q0 + (1 - sigma)/2 * slope
      !!
      !! Reference: van Leer (1977) slope limiter.
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      real(wp), intent(in) :: qm1, q0, qp1, sigma
      real(wp) :: face
      real(wp) :: d0, d1, slope
      d0 = q0 - qm1
      d1 = qp1 - q0
      if (d0*d1 > 0.0_wp) then
         if (abs(d0) < abs(d1)) then
            slope = d0
         else
            slope = d1
         end if
      else
         slope = 0.0_wp
      end if
      face = q0 + 0.5_wp*(1.0_wp - sigma)*slope
   end function plm_face_swept

   pure function weno5_face_swept(qm2, qm1, q0, qp1, qp2, sigma) result(face)
      !! WENO5-Z swept-average face value (u > 0, upwind cell = q0).
      !!
      !! Three quadratic candidates reconstructing the right edge of cell q0.
      !! Smoothness indicators: Jiang & Shu (1996), eq. 3.1.
      !! Z-weights: Borges et al. (2008).
      !!
      !! Stencil offsets relative to upwind cell i=q0:
      !!   r=0: {i-2, i-1, i}     d=1/10
      !!   r=1: {i-1, i,   i+1}   d=6/10
      !!   r=2: {i,   i+1, i+2}   d=3/10
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      real(wp), intent(in) :: qm2, qm1, q0, qp1, qp2
         !! Cell averages at i-2, i-1, i, i+1, i+2.
      real(wp), intent(in) :: sigma
         !! CFL of upwind cell: |u_eff|*dt/dx.
      real(wp) :: face

      real(wp), parameter :: d0 = 0.1_wp, d1 = 0.6_wp, d2 = 0.3_wp
      real(wp), parameter :: eps = 1.0e-36_wp
      real(wp), parameter :: c13_12 = 13.0_wp/12.0_wp
      real(wp), parameter :: c1_4 = 0.25_wp

      ! Curvature/slope/intercept for each candidate (quadratic: A + B*xi + C*xi^2)
      ! Cell-centre xi=0, right edge xi=1/2.
      real(wp) :: C0, B0, A0
      real(wp) :: C1, B1, A1
      real(wp) :: C2, B2, A2

      ! Smoothness indicators
      real(wp) :: beta0, beta1, beta2, tau5

      ! Z-weights
      real(wp) :: alpha0, alpha1, alpha2, alpha_sum
      real(wp) :: omega0, omega1, omega2

      ! Swept-average face values per candidate
      real(wp) :: face0, face1, face2

      ! Candidate r=0: stencil {i-2, i-1, i}, offsets (-2,-1,0)
      C0 = 0.5_wp*(qm2 - 2.0_wp*qm1 + q0)
      B0 = C0 - (qm1 - q0)
      A0 = q0 - C0/12.0_wp

      ! Candidate r=1: stencil {i-1, i, i+1}, offsets (-1,0,+1)
      C1 = 0.5_wp*(qm1 - 2.0_wp*q0 + qp1)
      B1 = 0.5_wp*(qp1 - qm1)
      A1 = q0 - C1/12.0_wp

      ! Candidate r=2: stencil {i, i+1, i+2}, offsets (0,+1,+2)
      C2 = 0.5_wp*(q0 - 2.0_wp*qp1 + qp2)
      B2 = qp1 - q0 - C2
      A2 = q0 - C2/12.0_wp

      ! JS96 smoothness indicators
      beta0 = c13_12*(qm2 - 2.0_wp*qm1 + q0)**2 + c1_4*(qm2 - 4.0_wp*qm1 + 3.0_wp*q0)**2
      beta1 = c13_12*(qm1 - 2.0_wp*q0 + qp1)**2 + c1_4*(qm1 - qp1)**2
      beta2 = c13_12*(q0 - 2.0_wp*qp1 + qp2)**2 + c1_4*(3.0_wp*q0 - 4.0_wp*qp1 + qp2)**2

      ! Z-weights (Borges 2008): tau5 = |beta0 - beta2|
      tau5 = abs(beta0 - beta2)
      alpha0 = d0*(1.0_wp + (tau5/(eps + beta0))**2)
      alpha1 = d1*(1.0_wp + (tau5/(eps + beta1))**2)
      alpha2 = d2*(1.0_wp + (tau5/(eps + beta2))**2)
      alpha_sum = alpha0 + alpha1 + alpha2
      omega0 = alpha0/alpha_sum
      omega1 = alpha1/alpha_sum
      omega2 = alpha2/alpha_sum

      ! Swept-average of quadratic A + B*xi + C*xi^2 over [1/2-sigma, 1/2]:
      !   = (A + B/2 + C/4) - sigma/2*(B + C) + sigma^2/6*(2*C)
      face0 = (A0 + 0.5_wp*B0 + 0.25_wp*C0) - 0.5_wp*sigma*(B0 + C0) &
              + sigma**2/3.0_wp*C0
      face1 = (A1 + 0.5_wp*B1 + 0.25_wp*C1) - 0.5_wp*sigma*(B1 + C1) &
              + sigma**2/3.0_wp*C1
      face2 = (A2 + 0.5_wp*B2 + 0.25_wp*C2) - 0.5_wp*sigma*(B2 + C2) &
              + sigma**2/3.0_wp*C2

      face = omega0*face0 + omega1*face1 + omega2*face2
   end function weno5_face_swept

   pure function weno7_face_swept(qm3, qm2, qm1, q0, qp1, qp2, qp3, sigma) result(face)
      !! WENO7-Z swept-average face value (u > 0, upwind cell = q0).
      !!
      !! Four cubic candidates with optimal weights d=(4,18,12,1)/35.
      !! Z-weights: tau7 = |beta0 - beta3|.
      !!
      !! Coefficient matrices from Balsara & Shu (2000), Table 1.
      !! Smoothness indicators from Balsara & Shu (2000), eq. 2.17.
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      real(wp), intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3
         !! Cell averages at i-3 through i+3.
      real(wp), intent(in) :: sigma
         !! CFL of upwind cell: |u_eff|*dt/dx.
      real(wp) :: face

      ! Optimal weights for WENO7 (BS00)
      real(wp), parameter :: d0_w = 4.0_wp/35.0_wp
      real(wp), parameter :: d1_w = 18.0_wp/35.0_wp
      real(wp), parameter :: d2_w = 12.0_wp/35.0_wp
      real(wp), parameter :: d3_w = 1.0_wp/35.0_wp
      real(wp), parameter :: eps = 1.0e-36_wp

      ! Cubic polynomial coefficients A + B*xi + C*xi^2 + D*xi^3
      ! for each of the 4 stencils (via exact rational matrices from BS00)
      real(wp) :: A0, B0, C0, D0
      real(wp) :: A1, B1, C1, D1
      real(wp) :: A2, B2, C2, D2
      real(wp) :: A3, B3, C3, D3

      ! Smoothness indicators (BS00 eq 2.17 exact)
      real(wp) :: beta0, beta1, beta2, beta3, tau7

      ! Z-weights
      real(wp) :: alpha0, alpha1, alpha2, alpha3, alpha_sum
      real(wp) :: omega0, omega1, omega2, omega3

      ! Swept-average face values per candidate
      real(wp) :: face0, face1, face2, face3

      ! Cubic swept-average helpers: p(1/2), p'(1/2), p''(1/2), p'''(1/2)
      real(wp) :: p_half, pp_half, ppp_half, pppp_half

      ! Stencil 0: offsets (-3,-2,-1,0) — q values: qm3,qm2,qm1,q0
      A0 = (1.0_wp/24.0_wp)*qm3 + (-1.0_wp/6.0_wp)*qm2 &
           + (5.0_wp/24.0_wp)*qm1 + (11.0_wp/12.0_wp)*q0
      B0 = (-7.0_wp/24.0_wp)*qm3 + (11.0_wp/8.0_wp)*qm2 &
           + (-23.0_wp/8.0_wp)*qm1 + (43.0_wp/24.0_wp)*q0
      C0 = (-1.0_wp/2.0_wp)*qm3 + 2.0_wp*qm2 &
           + (-5.0_wp/2.0_wp)*qm1 + 1.0_wp*q0
      D0 = (-1.0_wp/6.0_wp)*qm3 + (1.0_wp/2.0_wp)*qm2 &
           + (-1.0_wp/2.0_wp)*qm1 + (1.0_wp/6.0_wp)*q0

      ! Stencil 1: offsets (-2,-1,0,+1) — q values: qm2,qm1,q0,qp1
      A1 = 0.0_wp*qm2 + (-1.0_wp/24.0_wp)*qm1 &
           + (13.0_wp/12.0_wp)*q0 + (-1.0_wp/24.0_wp)*qp1
      B1 = (5.0_wp/24.0_wp)*qm2 + (-9.0_wp/8.0_wp)*qm1 &
           + (5.0_wp/8.0_wp)*q0 + (7.0_wp/24.0_wp)*qp1
      C1 = 0.0_wp*qm2 + (1.0_wp/2.0_wp)*qm1 &
           + (-1.0_wp)*q0 + (1.0_wp/2.0_wp)*qp1
      D1 = (-1.0_wp/6.0_wp)*qm2 + (1.0_wp/2.0_wp)*qm1 &
           + (-1.0_wp/2.0_wp)*q0 + (1.0_wp/6.0_wp)*qp1

      ! Stencil 2: offsets (-1,0,+1,+2) — q values: qm1,q0,qp1,qp2
      A2 = (-1.0_wp/24.0_wp)*qm1 + (13.0_wp/12.0_wp)*q0 &
           + (-1.0_wp/24.0_wp)*qp1 + 0.0_wp*qp2
      B2 = (-7.0_wp/24.0_wp)*qm1 + (-5.0_wp/8.0_wp)*q0 &
           + (9.0_wp/8.0_wp)*qp1 + (-5.0_wp/24.0_wp)*qp2
      C2 = (1.0_wp/2.0_wp)*qm1 + (-1.0_wp)*q0 &
           + (1.0_wp/2.0_wp)*qp1 + 0.0_wp*qp2
      D2 = (-1.0_wp/6.0_wp)*qm1 + (1.0_wp/2.0_wp)*q0 &
           + (-1.0_wp/2.0_wp)*qp1 + (1.0_wp/6.0_wp)*qp2

      ! Stencil 3: offsets (0,+1,+2,+3) — q values: q0,qp1,qp2,qp3
      A3 = (11.0_wp/12.0_wp)*q0 + (5.0_wp/24.0_wp)*qp1 &
           + (-1.0_wp/6.0_wp)*qp2 + (1.0_wp/24.0_wp)*qp3
      B3 = (-43.0_wp/24.0_wp)*q0 + (23.0_wp/8.0_wp)*qp1 &
           + (-11.0_wp/8.0_wp)*qp2 + (7.0_wp/24.0_wp)*qp3
      C3 = 1.0_wp*q0 + (-5.0_wp/2.0_wp)*qp1 &
           + 2.0_wp*qp2 + (-1.0_wp/2.0_wp)*qp3
      D3 = (-1.0_wp/6.0_wp)*q0 + (1.0_wp/2.0_wp)*qp1 &
           + (-1.0_wp/2.0_wp)*qp2 + (1.0_wp/6.0_wp)*qp3

      ! Smoothness indicators: the published integer-coefficient quadratic
      ! forms (Balsara & Shu 2000, eq. 2.17), one per candidate stencil.
      ! These are the forms validated in the numpy prototype
      ! (local_archive/prototypes/tracer_weno/) — do not substitute re-derivations.
      beta0 = qm3*(547.0_wp*qm3 - 3882.0_wp*qm2 + 4642.0_wp*qm1 - 1854.0_wp*q0) &
              + qm2*(7043.0_wp*qm2 - 17246.0_wp*qm1 + 7042.0_wp*q0) &
              + qm1*(11003.0_wp*qm1 - 9402.0_wp*q0) &
              + 2107.0_wp*q0**2
      beta1 = qm2*(267.0_wp*qm2 - 1642.0_wp*qm1 + 1602.0_wp*q0 - 494.0_wp*qp1) &
              + qm1*(2843.0_wp*qm1 - 5966.0_wp*q0 + 1922.0_wp*qp1) &
              + q0*(3443.0_wp*q0 - 2522.0_wp*qp1) &
              + 547.0_wp*qp1**2
      beta2 = qm1*(547.0_wp*qm1 - 2522.0_wp*q0 + 1922.0_wp*qp1 - 494.0_wp*qp2) &
              + q0*(3443.0_wp*q0 - 5966.0_wp*qp1 + 1602.0_wp*qp2) &
              + qp1*(2843.0_wp*qp1 - 1642.0_wp*qp2) &
              + 267.0_wp*qp2**2
      beta3 = q0*(2107.0_wp*q0 - 9402.0_wp*qp1 + 7042.0_wp*qp2 - 1854.0_wp*qp3) &
              + qp1*(11003.0_wp*qp1 - 17246.0_wp*qp2 + 4642.0_wp*qp3) &
              + qp2*(7043.0_wp*qp2 - 3882.0_wp*qp3) &
              + 547.0_wp*qp3**2

      tau7 = abs(beta0 - beta3)
      alpha0 = d0_w*(1.0_wp + (tau7/(eps + beta0))**2)
      alpha1 = d1_w*(1.0_wp + (tau7/(eps + beta1))**2)
      alpha2 = d2_w*(1.0_wp + (tau7/(eps + beta2))**2)
      alpha3 = d3_w*(1.0_wp + (tau7/(eps + beta3))**2)
      alpha_sum = alpha0 + alpha1 + alpha2 + alpha3
      omega0 = alpha0/alpha_sum
      omega1 = alpha1/alpha_sum
      omega2 = alpha2/alpha_sum
      omega3 = alpha3/alpha_sum

      ! Cubic swept-average over [1/2 - sigma, 1/2]:
      !   q_face = p(1/2) - sigma/2 * p'(1/2) + sigma^2/6 * p''(1/2) - sigma^3/24 * p'''(1/2)
      ! where:
      !   p(1/2)   = A + B/2 + C/4 + D/8
      !   p'(1/2)  = B + C + 3D/4
      !   p''(1/2) = 2C + 3D
      !   p'''(1/2)= 6D

      ! Stencil 0
      p_half = A0 + 0.5_wp*B0 + 0.25_wp*C0 + 0.125_wp*D0
      pp_half = B0 + C0 + 0.75_wp*D0
      ppp_half = 2.0_wp*C0 + 3.0_wp*D0
      pppp_half = 6.0_wp*D0
      face0 = p_half - 0.5_wp*sigma*pp_half &
              + sigma**2/6.0_wp*ppp_half - sigma**3/24.0_wp*pppp_half

      ! Stencil 1
      p_half = A1 + 0.5_wp*B1 + 0.25_wp*C1 + 0.125_wp*D1
      pp_half = B1 + C1 + 0.75_wp*D1
      ppp_half = 2.0_wp*C1 + 3.0_wp*D1
      pppp_half = 6.0_wp*D1
      face1 = p_half - 0.5_wp*sigma*pp_half &
              + sigma**2/6.0_wp*ppp_half - sigma**3/24.0_wp*pppp_half

      ! Stencil 2
      p_half = A2 + 0.5_wp*B2 + 0.25_wp*C2 + 0.125_wp*D2
      pp_half = B2 + C2 + 0.75_wp*D2
      ppp_half = 2.0_wp*C2 + 3.0_wp*D2
      pppp_half = 6.0_wp*D2
      face2 = p_half - 0.5_wp*sigma*pp_half &
              + sigma**2/6.0_wp*ppp_half - sigma**3/24.0_wp*pppp_half

      ! Stencil 3
      p_half = A3 + 0.5_wp*B3 + 0.25_wp*C3 + 0.125_wp*D3
      pp_half = B3 + C3 + 0.75_wp*D3
      ppp_half = 2.0_wp*C3 + 3.0_wp*D3
      pppp_half = 6.0_wp*D3
      face3 = p_half - 0.5_wp*sigma*pp_half &
              + sigma**2/6.0_wp*ppp_half - sigma**3/24.0_wp*pppp_half

      face = omega0*face0 + omega1*face1 + omega2*face2 + omega3*face3
   end function weno7_face_swept

   pure function weno9_face_swept(qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, sigma) &
      result(face)
      !! WENO9-Z swept-average face value (u > 0, upwind cell = q0).
      !!
      !! Five quartic candidates with optimal weights
      !! d = (1/126, 10/63, 10/21, 20/63, 5/126).
      !! Z-weights: tau9 = |beta0 - beta4|.
      !! Smoothness indicators use the simplified 3-term form.
      !!
      !! Reference: Balsara & Shu (2000).
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      real(wp), intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4
         !! Cell averages at i-4 through i+4.
      real(wp), intent(in) :: sigma
         !! CFL of upwind cell: |u_eff|*dt/dx.
      real(wp) :: face

      ! Optimal weights (BS00): d = (1/126, 10/63, 10/21, 20/63, 5/126)
      real(wp), parameter :: d0_w = 1.0_wp/126.0_wp
      real(wp), parameter :: d1_w = 10.0_wp/63.0_wp
      real(wp), parameter :: d2_w = 10.0_wp/21.0_wp
      real(wp), parameter :: d3_w = 20.0_wp/63.0_wp
      real(wp), parameter :: d4_w = 5.0_wp/126.0_wp
      real(wp), parameter :: eps = 1.0e-36_wp
      real(wp), parameter :: c13_12 = 13.0_wp/12.0_wp
      real(wp), parameter :: c1_4 = 0.25_wp
      real(wp), parameter :: c1_80 = 1.0_wp/80.0_wp

      ! Quartic polynomial A + B*xi + C*xi^2 + D*xi^3 + E*xi^4 per stencil
      real(wp) :: A0, B0, C0, D0, E0
      real(wp) :: A1, B1, C1, D1, E1
      real(wp) :: A2, B2, C2, D2, E2
      real(wp) :: A3, B3, C3, D3, E3
      real(wp) :: A4, B4, C4, D4, E4

      ! Simplified 3-term beta (centred on stencil centre c):
      !   beta = (13/12)*(q_{c-1}-2*q_c+q_{c+1})^2 + (1/4)*(q_{c-1}-q_{c+1})^2
      !          + (1/80)*(q_{c-2}-4*q_{c-1}+6*q_c-4*q_{c+1}+q_{c+2})^2
      real(wp) :: beta0, beta1, beta2, beta3, beta4, tau9

      ! Z-weights
      real(wp) :: alpha0, alpha1, alpha2, alpha3, alpha4, alpha_sum
      real(wp) :: omega0, omega1, omega2, omega3, omega4

      ! Swept-average face values per candidate
      real(wp) :: face0, face1, face2, face3, face4

      ! Quartic swept-average helpers
      real(wp) :: p_half, pp_half, ppp_half, pppp_half, ppppp_half

      ! -----------------------------------------------------------------------
      ! Quartic coefficients (exact rational from BS00)
      ! -----------------------------------------------------------------------

      ! Stencil 0: offsets (-4,-3,-2,-1,0) → qm4,qm3,qm2,qm1,q0
      A0 = (-71.0_wp/1920.0_wp)*qm4 + (91.0_wp/480.0_wp)*qm3 &
           + (-373.0_wp/960.0_wp)*qm2 + (57.0_wp/160.0_wp)*qm1 &
           + (563.0_wp/640.0_wp)*q0
      B0 = (3.0_wp/16.0_wp)*qm4 + (-25.0_wp/24.0_wp)*qm3 &
           + (5.0_wp/2.0_wp)*qm2 + (-29.0_wp/8.0_wp)*qm1 &
           + (95.0_wp/48.0_wp)*q0
      C0 = (7.0_wp/16.0_wp)*qm4 + (-9.0_wp/4.0_wp)*qm3 &
           + (37.0_wp/8.0_wp)*qm2 + (-17.0_wp/4.0_wp)*qm1 &
           + (23.0_wp/16.0_wp)*q0
      D0 = (1.0_wp/4.0_wp)*qm4 + (-7.0_wp/6.0_wp)*qm3 &
           + 2.0_wp*qm2 + (-3.0_wp/2.0_wp)*qm1 &
           + (5.0_wp/12.0_wp)*q0
      E0 = (1.0_wp/24.0_wp)*qm4 + (-1.0_wp/6.0_wp)*qm3 &
           + (1.0_wp/4.0_wp)*qm2 + (-1.0_wp/6.0_wp)*qm1 &
           + (1.0_wp/24.0_wp)*q0

      ! Stencil 1: offsets (-3,-2,-1,0,+1) → qm3,qm2,qm1,q0,qp1
      A1 = (3.0_wp/640.0_wp)*qm3 + (-3.0_wp/160.0_wp)*qm2 &
           + (-13.0_wp/960.0_wp)*qm1 + (511.0_wp/480.0_wp)*q0 &
           + (-71.0_wp/1920.0_wp)*qp1
      B1 = (-5.0_wp/48.0_wp)*qm3 + (5.0_wp/8.0_wp)*qm2 &
           + (-7.0_wp/4.0_wp)*qm1 + (25.0_wp/24.0_wp)*q0 &
           + (3.0_wp/16.0_wp)*qp1
      C1 = (-1.0_wp/16.0_wp)*qm3 + (1.0_wp/4.0_wp)*qm2 &
           + (1.0_wp/8.0_wp)*qm1 + (-3.0_wp/4.0_wp)*q0 &
           + (7.0_wp/16.0_wp)*qp1
      D1 = (1.0_wp/12.0_wp)*qm3 + (-1.0_wp/2.0_wp)*qm2 &
           + 1.0_wp*qm1 + (-5.0_wp/6.0_wp)*q0 &
           + (1.0_wp/4.0_wp)*qp1
      E1 = (1.0_wp/24.0_wp)*qm3 + (-1.0_wp/6.0_wp)*qm2 &
           + (1.0_wp/4.0_wp)*qm1 + (-1.0_wp/6.0_wp)*q0 &
           + (1.0_wp/24.0_wp)*qp1

      ! Stencil 2: offsets (-2,-1,0,+1,+2) → qm2,qm1,q0,qp1,qp2
      A2 = (3.0_wp/640.0_wp)*qm2 + (-29.0_wp/480.0_wp)*qm1 &
           + (1067.0_wp/960.0_wp)*q0 + (-29.0_wp/480.0_wp)*qp1 &
           + (3.0_wp/640.0_wp)*qp2
      B2 = (5.0_wp/48.0_wp)*qm2 + (-17.0_wp/24.0_wp)*qm1 &
           + 0.0_wp*q0 + (17.0_wp/24.0_wp)*qp1 &
           + (-5.0_wp/48.0_wp)*qp2
      C2 = (-1.0_wp/16.0_wp)*qm2 + (3.0_wp/4.0_wp)*qm1 &
           + (-11.0_wp/8.0_wp)*q0 + (3.0_wp/4.0_wp)*qp1 &
           + (-1.0_wp/16.0_wp)*qp2
      D2 = (-1.0_wp/12.0_wp)*qm2 + (1.0_wp/6.0_wp)*qm1 &
           + 0.0_wp*q0 + (-1.0_wp/6.0_wp)*qp1 &
           + (1.0_wp/12.0_wp)*qp2
      E2 = (1.0_wp/24.0_wp)*qm2 + (-1.0_wp/6.0_wp)*qm1 &
           + (1.0_wp/4.0_wp)*q0 + (-1.0_wp/6.0_wp)*qp1 &
           + (1.0_wp/24.0_wp)*qp2

      ! Stencil 3: offsets (-1,0,+1,+2,+3) → qm1,q0,qp1,qp2,qp3
      A3 = (-71.0_wp/1920.0_wp)*qm1 + (511.0_wp/480.0_wp)*q0 &
           + (-13.0_wp/960.0_wp)*qp1 + (-3.0_wp/160.0_wp)*qp2 &
           + (3.0_wp/640.0_wp)*qp3
      B3 = (-3.0_wp/16.0_wp)*qm1 + (-25.0_wp/24.0_wp)*q0 &
           + (7.0_wp/4.0_wp)*qp1 + (-5.0_wp/8.0_wp)*qp2 &
           + (5.0_wp/48.0_wp)*qp3
      C3 = (7.0_wp/16.0_wp)*qm1 + (-3.0_wp/4.0_wp)*q0 &
           + (1.0_wp/8.0_wp)*qp1 + (1.0_wp/4.0_wp)*qp2 &
           + (-1.0_wp/16.0_wp)*qp3
      D3 = (-1.0_wp/4.0_wp)*qm1 + (5.0_wp/6.0_wp)*q0 &
           + (-1.0_wp)*qp1 + (1.0_wp/2.0_wp)*qp2 &
           + (-1.0_wp/12.0_wp)*qp3
      E3 = (1.0_wp/24.0_wp)*qm1 + (-1.0_wp/6.0_wp)*q0 &
           + (1.0_wp/4.0_wp)*qp1 + (-1.0_wp/6.0_wp)*qp2 &
           + (1.0_wp/24.0_wp)*qp3

      ! Stencil 4: offsets (0,+1,+2,+3,+4) → q0,qp1,qp2,qp3,qp4
      A4 = (563.0_wp/640.0_wp)*q0 + (57.0_wp/160.0_wp)*qp1 &
           + (-373.0_wp/960.0_wp)*qp2 + (91.0_wp/480.0_wp)*qp3 &
           + (-71.0_wp/1920.0_wp)*qp4
      B4 = (-95.0_wp/48.0_wp)*q0 + (29.0_wp/8.0_wp)*qp1 &
           + (-5.0_wp/2.0_wp)*qp2 + (25.0_wp/24.0_wp)*qp3 &
           + (-3.0_wp/16.0_wp)*qp4
      C4 = (23.0_wp/16.0_wp)*q0 + (-17.0_wp/4.0_wp)*qp1 &
           + (37.0_wp/8.0_wp)*qp2 + (-9.0_wp/4.0_wp)*qp3 &
           + (7.0_wp/16.0_wp)*qp4
      D4 = (-5.0_wp/12.0_wp)*q0 + (3.0_wp/2.0_wp)*qp1 &
           + (-2.0_wp)*qp2 + (7.0_wp/6.0_wp)*qp3 &
           + (-1.0_wp/4.0_wp)*qp4
      E4 = (1.0_wp/24.0_wp)*q0 + (-1.0_wp/6.0_wp)*qp1 &
           + (1.0_wp/4.0_wp)*qp2 + (-1.0_wp/6.0_wp)*qp3 &
           + (1.0_wp/24.0_wp)*qp4

      ! -----------------------------------------------------------------------
      ! Simplified 3-term smoothness indicators (centred on each stencil)
      ! Stencil r centre c:
      !   beta_r = (13/12)*(q_{c-1}-2*q_c+q_{c+1})^2
      !          + (1/4)  *(q_{c-1}-q_{c+1})^2
      !          + (1/80) *(q_{c-2}-4*q_{c-1}+6*q_c-4*q_{c+1}+q_{c+2})^2
      ! -----------------------------------------------------------------------
      ! Stencil 0: c = qm2 (index i-2), c±1 = qm3/qm1, c±2 = qm4/q0
      beta0 = c13_12*(qm3 - 2.0_wp*qm2 + qm1)**2 &
              + c1_4*(qm3 - qm1)**2 &
              + c1_80*(qm4 - 4.0_wp*qm3 + 6.0_wp*qm2 - 4.0_wp*qm1 + q0)**2

      ! Stencil 1: c = qm1 (index i-1), c±1 = qm2/q0, c±2 = qm3/qp1
      beta1 = c13_12*(qm2 - 2.0_wp*qm1 + q0)**2 &
              + c1_4*(qm2 - q0)**2 &
              + c1_80*(qm3 - 4.0_wp*qm2 + 6.0_wp*qm1 - 4.0_wp*q0 + qp1)**2

      ! Stencil 2: c = q0 (index i), c±1 = qm1/qp1, c±2 = qm2/qp2
      beta2 = c13_12*(qm1 - 2.0_wp*q0 + qp1)**2 &
              + c1_4*(qm1 - qp1)**2 &
              + c1_80*(qm2 - 4.0_wp*qm1 + 6.0_wp*q0 - 4.0_wp*qp1 + qp2)**2

      ! Stencil 3: c = qp1 (index i+1), c±1 = q0/qp2, c±2 = qm1/qp3
      beta3 = c13_12*(q0 - 2.0_wp*qp1 + qp2)**2 &
              + c1_4*(q0 - qp2)**2 &
              + c1_80*(qm1 - 4.0_wp*q0 + 6.0_wp*qp1 - 4.0_wp*qp2 + qp3)**2

      ! Stencil 4: c = qp2 (index i+2), c±1 = qp1/qp3, c±2 = q0/qp4
      beta4 = c13_12*(qp1 - 2.0_wp*qp2 + qp3)**2 &
              + c1_4*(qp1 - qp3)**2 &
              + c1_80*(q0 - 4.0_wp*qp1 + 6.0_wp*qp2 - 4.0_wp*qp3 + qp4)**2

      tau9 = abs(beta0 - beta4)
      alpha0 = d0_w*(1.0_wp + (tau9/(eps + beta0))**2)
      alpha1 = d1_w*(1.0_wp + (tau9/(eps + beta1))**2)
      alpha2 = d2_w*(1.0_wp + (tau9/(eps + beta2))**2)
      alpha3 = d3_w*(1.0_wp + (tau9/(eps + beta3))**2)
      alpha4 = d4_w*(1.0_wp + (tau9/(eps + beta4))**2)
      alpha_sum = alpha0 + alpha1 + alpha2 + alpha3 + alpha4
      omega0 = alpha0/alpha_sum
      omega1 = alpha1/alpha_sum
      omega2 = alpha2/alpha_sum
      omega3 = alpha3/alpha_sum
      omega4 = alpha4/alpha_sum

      ! Quartic swept-average over [1/2 - sigma, 1/2]:
      !   q_face = p(1/2) - sigma/2 * p'(1/2) + sigma^2/6 * p''(1/2)
      !            - sigma^3/24 * p'''(1/2) + sigma^4/120 * p''''(1/2)
      ! where for p = A + B*xi + C*xi^2 + D*xi^3 + E*xi^4:
      !   p(1/2)    = A + B/2 + C/4 + D/8 + E/16
      !   p'(1/2)   = B + C   + 3D/4 + E/2
      !   p''(1/2)  = 2C + 3D + 3E
      !   p'''(1/2) = 6D + 12E
      !   p''''(1/2)= 24E

      ! Stencil 0
      p_half = A0 + 0.5_wp*B0 + 0.25_wp*C0 + 0.125_wp*D0 + 0.0625_wp*E0
      pp_half = B0 + C0 + 0.75_wp*D0 + 0.5_wp*E0
      ppp_half = 2.0_wp*C0 + 3.0_wp*D0 + 3.0_wp*E0
      pppp_half = 6.0_wp*D0 + 12.0_wp*E0
      ppppp_half = 24.0_wp*E0
      face0 = p_half - 0.5_wp*sigma*pp_half + sigma**2/6.0_wp*ppp_half &
              - sigma**3/24.0_wp*pppp_half + sigma**4/120.0_wp*ppppp_half

      ! Stencil 1
      p_half = A1 + 0.5_wp*B1 + 0.25_wp*C1 + 0.125_wp*D1 + 0.0625_wp*E1
      pp_half = B1 + C1 + 0.75_wp*D1 + 0.5_wp*E1
      ppp_half = 2.0_wp*C1 + 3.0_wp*D1 + 3.0_wp*E1
      pppp_half = 6.0_wp*D1 + 12.0_wp*E1
      ppppp_half = 24.0_wp*E1
      face1 = p_half - 0.5_wp*sigma*pp_half + sigma**2/6.0_wp*ppp_half &
              - sigma**3/24.0_wp*pppp_half + sigma**4/120.0_wp*ppppp_half

      ! Stencil 2
      p_half = A2 + 0.5_wp*B2 + 0.25_wp*C2 + 0.125_wp*D2 + 0.0625_wp*E2
      pp_half = B2 + C2 + 0.75_wp*D2 + 0.5_wp*E2
      ppp_half = 2.0_wp*C2 + 3.0_wp*D2 + 3.0_wp*E2
      pppp_half = 6.0_wp*D2 + 12.0_wp*E2
      ppppp_half = 24.0_wp*E2
      face2 = p_half - 0.5_wp*sigma*pp_half + sigma**2/6.0_wp*ppp_half &
              - sigma**3/24.0_wp*pppp_half + sigma**4/120.0_wp*ppppp_half

      ! Stencil 3
      p_half = A3 + 0.5_wp*B3 + 0.25_wp*C3 + 0.125_wp*D3 + 0.0625_wp*E3
      pp_half = B3 + C3 + 0.75_wp*D3 + 0.5_wp*E3
      ppp_half = 2.0_wp*C3 + 3.0_wp*D3 + 3.0_wp*E3
      pppp_half = 6.0_wp*D3 + 12.0_wp*E3
      ppppp_half = 24.0_wp*E3
      face3 = p_half - 0.5_wp*sigma*pp_half + sigma**2/6.0_wp*ppp_half &
              - sigma**3/24.0_wp*pppp_half + sigma**4/120.0_wp*ppppp_half

      ! Stencil 4
      p_half = A4 + 0.5_wp*B4 + 0.25_wp*C4 + 0.125_wp*D4 + 0.0625_wp*E4
      pp_half = B4 + C4 + 0.75_wp*D4 + 0.5_wp*E4
      ppp_half = 2.0_wp*C4 + 3.0_wp*D4 + 3.0_wp*E4
      pppp_half = 6.0_wp*D4 + 12.0_wp*E4
      ppppp_half = 24.0_wp*E4
      face4 = p_half - 0.5_wp*sigma*pp_half + sigma**2/6.0_wp*ppp_half &
              - sigma**3/24.0_wp*pppp_half + sigma**4/120.0_wp*ppppp_half

      face = omega0*face0 + omega1*face1 + omega2*face2 + omega3*face3 + omega4*face4
   end function weno9_face_swept

   ! ==========================================================================
   !  Main advection kernel
   ! ==========================================================================

   subroutine ml_advect_tracer_weno(h_layer, hu_layer, hv_layer, &
                                         mass_flux_x, mass_flux_y, &
                                         hTr_layer, weno_tr_snap, &
                                         grid, nz, dt, recon)
      !! Per-layer WENO horizontal advection of a scalar tracer hTr = h*Tr.
      !!
      !! Rung-adaptive: near boundaries the stencil automatically degrades
      !! from the requested scheme down to PLM or donor-cell so that every
      !! interior cell is always updated.
      !!
      !! Mass-flux split uses the same vanishing-layer-gated centred CWC
      !! as ml_advect_tracer_ppm — see that kernel's docstring for the
      !! algebra.
      !!
      !! Only called when recon /= TRACER_RECON_PPM (1=WENO5, 2=WENO7,
      !! 3=WENO9).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
         !! Number of vertical layers.
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
         !! Layer depths (m), bottom-up: k=1 bed, k=nz surface.
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
         !! Layer x-momentum flux (m^2/s).
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
         !! Layer y-momentum flux (m^2/s).
      real(wp), intent(in) :: mass_flux_x(grid%nx_total, grid%ny_total)
         !! Barotropic x mass flux (m^2/s).
      real(wp), intent(in) :: mass_flux_y(grid%nx_total, grid%ny_total)
         !! Barotropic y mass flux (m^2/s).
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
         !! Tracer h*Tr (m*[units]), updated in-place.
      real(wp), intent(out) :: weno_tr_snap(grid%nx_total, grid%ny_total, nz)
         !! Caller-owned concentration-snapshot scratch (Tr = hTr/h); filled
         !! then read here.  Passed in (not module state) to keep this `pure`.
      real(wp), intent(in) :: dt
         !! Timestep (s).
      integer, intent(in) :: recon
         !! TRACER_RECON_WENO5/7/9 code.

      integer :: nx, ny, nghost, k, kk, i, j
      real(wp) :: dx, dy, inv_dx, inv_dy
      real(wp) :: h_val, h_col, h_face, f_layer, u_eff, sigma
      real(wp) :: flux_left, flux_right, flux_bot, flux_top
      real(wp) :: tr_face_val
      real(wp) :: corr_xL, corr_xR, corr_yB, corr_yT
      real(wp) :: H_active_xL, H_active_xR, H_active_yB, H_active_yT
      real(wp) :: inv_H_xL, inv_H_xR, inv_H_yB, inv_H_yT
      real(wp) :: h_iL, h_iR, h_jB, h_jT
      real(wp) :: v_xL(NZ_STACK_MAX), v_xR(NZ_STACK_MAX)
      real(wp) :: v_yB(NZ_STACK_MAX), v_yT(NZ_STACK_MAX)
      real(wp) :: hf_xL(NZ_STACK_MAX), hf_xR(NZ_STACK_MAX)
      real(wp) :: hf_yB(NZ_STACK_MAX), hf_yT(NZ_STACK_MAX)
      integer :: rung_max_int, rung
      integer :: avail_up, avail_down

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy
      inv_dx = 1.0_wp/dx
      inv_dy = 1.0_wp/dy

      ! Map public TRACER_RECON_* to internal rung_max:
      !   WENO5(1) -> 2, WENO7(2) -> 3, WENO9(3) -> 4
      rung_max_int = recon + 1

      associate (tr_snap => weno_tr_snap)

         ! ---- Phase 0: concentration snapshot ----
         !$omp parallel do collapse(2) private(k, h_val)
         do j = 1, ny
         do i = 1, nx
         do k = 1, nz
            h_val = h_layer(i, j, k)
            if (h_val > DRY_TOLERANCE) then
               tr_snap(i, j, k) = hTr_layer(i, j, k)/h_val
            else
               tr_snap(i, j, k) = 0.0_wp
            end if
         end do
         end do
         end do
         !$omp end parallel do

         ! ---- Phase 1+2: combined flux + update ----
         !$omp parallel do collapse(2) private(k, kk, h_val, h_col, h_face, f_layer, u_eff, sigma, flux_left, flux_right, flux_bot, flux_top, tr_face_val, corr_xL, corr_xR, corr_yB, corr_yT, H_active_xL, H_active_xR, H_active_yB, H_active_yT, inv_H_xL, inv_H_xR, inv_H_yB, inv_H_yT, h_iL, h_iR, h_jB, h_jT, rung, avail_up, avail_down, v_xL, v_xR, v_yB, v_yT, hf_xL, hf_xR, hf_yB, hf_yT)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost

            h_col = 0.0_wp
            do kk = 1, nz
               h_col = h_col + h_layer(i, j, kk)
            end do

            if (h_col > ML_MIN_DEPTH) then

               ! Vanishing-layer-gated centred CWC (same as ml_advect_tracer_ppm).
               corr_xL = mass_flux_x(i - 1, j)
               corr_xR = mass_flux_x(i, j)
               corr_yB = mass_flux_y(i, j - 1)
               corr_yT = mass_flux_y(i, j)
               H_active_xL = 0.0_wp
               H_active_xR = 0.0_wp
               H_active_yB = 0.0_wp
               H_active_yT = 0.0_wp

               do k = 1, nz
                  h_val = h_layer(i, j, k)
                  h_iL = h_layer(i - 1, j, k)
                  h_iR = h_layer(i + 1, j, k)
                  h_jB = h_layer(i, j - 1, k)
                  h_jT = h_layer(i, j + 1, k)

                  if (h_val > VANISHING_LAYER_TOL .and. h_iL > VANISHING_LAYER_TOL) then
                     v_xL(k) = 0.5_wp*(hu_layer(i - 1, j, k) + hu_layer(i, j, k))
                     hf_xL(k) = 0.5_wp*(h_iL + h_val)
                  else
                     v_xL(k) = 0.0_wp
                     hf_xL(k) = 0.0_wp
                  end if
                  if (h_val > VANISHING_LAYER_TOL .and. h_iR > VANISHING_LAYER_TOL) then
                     v_xR(k) = 0.5_wp*(hu_layer(i, j, k) + hu_layer(i + 1, j, k))
                     hf_xR(k) = 0.5_wp*(h_val + h_iR)
                  else
                     v_xR(k) = 0.0_wp
                     hf_xR(k) = 0.0_wp
                  end if
                  if (h_val > VANISHING_LAYER_TOL .and. h_jB > VANISHING_LAYER_TOL) then
                     v_yB(k) = 0.5_wp*(hv_layer(i, j - 1, k) + hv_layer(i, j, k))
                     hf_yB(k) = 0.5_wp*(h_jB + h_val)
                  else
                     v_yB(k) = 0.0_wp
                     hf_yB(k) = 0.0_wp
                  end if
                  if (h_val > VANISHING_LAYER_TOL .and. h_jT > VANISHING_LAYER_TOL) then
                     v_yT(k) = 0.5_wp*(hv_layer(i, j, k) + hv_layer(i, j + 1, k))
                     hf_yT(k) = 0.5_wp*(h_val + h_jT)
                  else
                     v_yT(k) = 0.0_wp
                     hf_yT(k) = 0.0_wp
                  end if

                  corr_xL = corr_xL - v_xL(k)
                  corr_xR = corr_xR - v_xR(k)
                  corr_yB = corr_yB - v_yB(k)
                  corr_yT = corr_yT - v_yT(k)
                  H_active_xL = H_active_xL + hf_xL(k)
                  H_active_xR = H_active_xR + hf_xR(k)
                  H_active_yB = H_active_yB + hf_yB(k)
                  H_active_yT = H_active_yT + hf_yT(k)
               end do

               if (H_active_xL > 0.0_wp) then
                  inv_H_xL = 1.0_wp/H_active_xL
               else
                  inv_H_xL = 0.0_wp
               end if
               if (H_active_xR > 0.0_wp) then
                  inv_H_xR = 1.0_wp/H_active_xR
               else
                  inv_H_xR = 0.0_wp
               end if
               if (H_active_yB > 0.0_wp) then
                  inv_H_yB = 1.0_wp/H_active_yB
               else
                  inv_H_yB = 0.0_wp
               end if
               if (H_active_yT > 0.0_wp) then
                  inv_H_yT = 1.0_wp/H_active_yT
               else
                  inv_H_yT = 0.0_wp
               end if

               ! Availability counts for the rung ladder.  Stencils may use at
               ! most ONE ghost row beyond the physical interior: ghost rows are
               ! zero-gradient copies of the first interior cell, so deeper rows
               ! carry no information and would bias the smoothness indicators
               ! (this also matches PPM, whose edge interpolation at the first
               ! interior cell reads one ghost row).  For the face between cells
               ! c and c+1 along x (interior cells span [nghost+1, nx-nghost]):
               !   avail_west(c) = c - nghost + 1         (cells nghost .. c)
               !   avail_east(c) = nx - nghost - c + 1    (cells c+1 .. nx-nghost+1)
               ! upwind/downwind assignment follows the flux sign; the y-faces
               ! use the same formulas with (j, ny).  Bounds safety: rung R
               ! engages only when avail_up/down cover its reach, so the deepest
               ! read is >= ghost row nghost (one row in) and <= nx-nghost+1 —
               ! in-array for any nghost >= the configure-enforced minimum.

               do k = 1, nz

                  ! ---- x-direction: LEFT face (between i-1 and i) ----
                  f_layer = v_xL(k) + hf_xL(k)*inv_H_xL*corr_xL
                  h_face = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dx, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     ! Upwind cell c = i-1 (flowing right)
                     avail_up = i - nghost                 ! avail_west(i-1)
                     avail_down = nx - nghost - i + 2      ! avail_east(i-1)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i - 5, j, k), tr_snap(i - 4, j, k), &
                                      tr_snap(i - 3, j, k), tr_snap(i - 2, j, k), &
                                      tr_snap(i - 1, j, k), &
                                      tr_snap(i, j, k), tr_snap(i + 1, j, k), &
                                      tr_snap(i + 2, j, k), tr_snap(i + 3, j, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i - 4, j, k), tr_snap(i - 3, j, k), &
                                      tr_snap(i - 2, j, k), tr_snap(i - 1, j, k), &
                                      tr_snap(i, j, k), tr_snap(i + 1, j, k), &
                                      tr_snap(i + 2, j, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i - 3, j, k), tr_snap(i - 2, j, k), &
                                      tr_snap(i - 1, j, k), &
                                      tr_snap(i, j, k), tr_snap(i + 1, j, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i - 2, j, k), tr_snap(i - 1, j, k), &
                                      tr_snap(i, j, k), sigma)
                     case default
                        tr_face_val = tr_snap(i - 1, j, k)  ! donor
                     end select
                     flux_left = f_layer*tr_face_val
                  else if (f_layer < 0.0_wp) then
                     ! Upwind cell c = i (flowing left — flip stencil)
                     avail_up = nx - nghost - i + 2        ! avail_east(i-1)
                     avail_down = i - nghost               ! avail_west(i-1)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i + 4, j, k), tr_snap(i + 3, j, k), &
                                      tr_snap(i + 2, j, k), tr_snap(i + 1, j, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i - 1, j, k), tr_snap(i - 2, j, k), &
                                      tr_snap(i - 3, j, k), tr_snap(i - 4, j, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i + 3, j, k), tr_snap(i + 2, j, k), &
                                      tr_snap(i + 1, j, k), tr_snap(i, j, k), &
                                      tr_snap(i - 1, j, k), tr_snap(i - 2, j, k), &
                                      tr_snap(i - 3, j, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i + 2, j, k), tr_snap(i + 1, j, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i - 1, j, k), tr_snap(i - 2, j, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i + 1, j, k), tr_snap(i, j, k), &
                                      tr_snap(i - 1, j, k), sigma)
                     case default
                        tr_face_val = tr_snap(i, j, k)  ! donor
                     end select
                     flux_left = f_layer*tr_face_val
                  else
                     flux_left = 0.0_wp
                  end if

                  ! ---- x-direction: RIGHT face (between i and i+1) ----
                  f_layer = v_xR(k) + hf_xR(k)*inv_H_xR*corr_xR
                  h_face = 0.5_wp*(h_layer(i, j, k) + h_layer(i + 1, j, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dx, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     ! Upwind cell c = i (flowing right)
                     avail_up = i - nghost + 1             ! avail_west(i)
                     avail_down = nx - nghost - i + 1      ! avail_east(i)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i - 4, j, k), tr_snap(i - 3, j, k), &
                                      tr_snap(i - 2, j, k), tr_snap(i - 1, j, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i + 1, j, k), tr_snap(i + 2, j, k), &
                                      tr_snap(i + 3, j, k), tr_snap(i + 4, j, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i - 3, j, k), tr_snap(i - 2, j, k), &
                                      tr_snap(i - 1, j, k), tr_snap(i, j, k), &
                                      tr_snap(i + 1, j, k), tr_snap(i + 2, j, k), &
                                      tr_snap(i + 3, j, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i - 2, j, k), tr_snap(i - 1, j, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i + 1, j, k), tr_snap(i + 2, j, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i - 1, j, k), tr_snap(i, j, k), &
                                      tr_snap(i + 1, j, k), sigma)
                     case default
                        tr_face_val = tr_snap(i, j, k)  ! donor
                     end select
                     flux_right = f_layer*tr_face_val
                  else if (f_layer < 0.0_wp) then
                     ! Upwind cell c = i+1 (flowing left — flip stencil)
                     avail_up = nx - nghost - i + 1        ! avail_east(i)
                     avail_down = i - nghost + 1           ! avail_west(i)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i + 5, j, k), tr_snap(i + 4, j, k), &
                                      tr_snap(i + 3, j, k), tr_snap(i + 2, j, k), &
                                      tr_snap(i + 1, j, k), &
                                      tr_snap(i, j, k), tr_snap(i - 1, j, k), &
                                      tr_snap(i - 2, j, k), tr_snap(i - 3, j, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i + 4, j, k), tr_snap(i + 3, j, k), &
                                      tr_snap(i + 2, j, k), tr_snap(i + 1, j, k), &
                                      tr_snap(i, j, k), tr_snap(i - 1, j, k), &
                                      tr_snap(i - 2, j, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i + 3, j, k), tr_snap(i + 2, j, k), &
                                      tr_snap(i + 1, j, k), &
                                      tr_snap(i, j, k), tr_snap(i - 1, j, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i + 2, j, k), tr_snap(i + 1, j, k), &
                                      tr_snap(i, j, k), sigma)
                     case default
                        tr_face_val = tr_snap(i + 1, j, k)  ! donor
                     end select
                     flux_right = f_layer*tr_face_val
                  else
                     flux_right = 0.0_wp
                  end if

                  ! ---- y-direction: BOTTOM face (between j-1 and j) ----
                  f_layer = v_yB(k) + hf_yB(k)*inv_H_yB*corr_yB
                  h_face = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dy, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     ! Upwind cell c = j-1 (flowing north)
                     avail_up = j - nghost                 ! avail_south(j-1)
                     avail_down = ny - nghost - j + 2      ! avail_north(j-1)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i, j - 5, k), tr_snap(i, j - 4, k), &
                                      tr_snap(i, j - 3, k), tr_snap(i, j - 2, k), &
                                      tr_snap(i, j - 1, k), &
                                      tr_snap(i, j, k), tr_snap(i, j + 1, k), &
                                      tr_snap(i, j + 2, k), tr_snap(i, j + 3, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i, j - 4, k), tr_snap(i, j - 3, k), &
                                      tr_snap(i, j - 2, k), tr_snap(i, j - 1, k), &
                                      tr_snap(i, j, k), tr_snap(i, j + 1, k), &
                                      tr_snap(i, j + 2, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i, j - 3, k), tr_snap(i, j - 2, k), &
                                      tr_snap(i, j - 1, k), &
                                      tr_snap(i, j, k), tr_snap(i, j + 1, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i, j - 2, k), tr_snap(i, j - 1, k), &
                                      tr_snap(i, j, k), sigma)
                     case default
                        tr_face_val = tr_snap(i, j - 1, k)  ! donor
                     end select
                     flux_bot = f_layer*tr_face_val
                  else if (f_layer < 0.0_wp) then
                     ! Upwind cell c = j (flowing south — flip stencil)
                     avail_up = ny - nghost - j + 2        ! avail_north(j-1)
                     avail_down = j - nghost               ! avail_south(j-1)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i, j + 4, k), tr_snap(i, j + 3, k), &
                                      tr_snap(i, j + 2, k), tr_snap(i, j + 1, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i, j - 1, k), tr_snap(i, j - 2, k), &
                                      tr_snap(i, j - 3, k), tr_snap(i, j - 4, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i, j + 3, k), tr_snap(i, j + 2, k), &
                                      tr_snap(i, j + 1, k), tr_snap(i, j, k), &
                                      tr_snap(i, j - 1, k), tr_snap(i, j - 2, k), &
                                      tr_snap(i, j - 3, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i, j + 2, k), tr_snap(i, j + 1, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i, j - 1, k), tr_snap(i, j - 2, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i, j + 1, k), tr_snap(i, j, k), &
                                      tr_snap(i, j - 1, k), sigma)
                     case default
                        tr_face_val = tr_snap(i, j, k)  ! donor
                     end select
                     flux_bot = f_layer*tr_face_val
                  else
                     flux_bot = 0.0_wp
                  end if

                  ! ---- y-direction: TOP face (between j and j+1) ----
                  f_layer = v_yT(k) + hf_yT(k)*inv_H_yT*corr_yT
                  h_face = 0.5_wp*(h_layer(i, j, k) + h_layer(i, j + 1, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dy, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     ! Upwind cell c = j (flowing north)
                     avail_up = j - nghost + 1             ! avail_south(j)
                     avail_down = ny - nghost - j + 1      ! avail_north(j)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i, j - 4, k), tr_snap(i, j - 3, k), &
                                      tr_snap(i, j - 2, k), tr_snap(i, j - 1, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i, j + 1, k), tr_snap(i, j + 2, k), &
                                      tr_snap(i, j + 3, k), tr_snap(i, j + 4, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i, j - 3, k), tr_snap(i, j - 2, k), &
                                      tr_snap(i, j - 1, k), tr_snap(i, j, k), &
                                      tr_snap(i, j + 1, k), tr_snap(i, j + 2, k), &
                                      tr_snap(i, j + 3, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i, j - 2, k), tr_snap(i, j - 1, k), &
                                      tr_snap(i, j, k), &
                                      tr_snap(i, j + 1, k), tr_snap(i, j + 2, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i, j - 1, k), tr_snap(i, j, k), &
                                      tr_snap(i, j + 1, k), sigma)
                     case default
                        tr_face_val = tr_snap(i, j, k)  ! donor
                     end select
                     flux_top = f_layer*tr_face_val
                  else if (f_layer < 0.0_wp) then
                     ! Upwind cell c = j+1 (flowing south — flip stencil)
                     avail_up = ny - nghost - j + 1        ! avail_north(j)
                     avail_down = j - nghost + 1           ! avail_south(j)
                     rung = recon_rung_for_face(avail_up, avail_down, rung_max_int)
                     select case (rung)
                     case (4)
                        tr_face_val = weno9_face_swept( &
                                      tr_snap(i, j + 5, k), tr_snap(i, j + 4, k), &
                                      tr_snap(i, j + 3, k), tr_snap(i, j + 2, k), &
                                      tr_snap(i, j + 1, k), &
                                      tr_snap(i, j, k), tr_snap(i, j - 1, k), &
                                      tr_snap(i, j - 2, k), tr_snap(i, j - 3, k), sigma)
                     case (3)
                        tr_face_val = weno7_face_swept( &
                                      tr_snap(i, j + 4, k), tr_snap(i, j + 3, k), &
                                      tr_snap(i, j + 2, k), tr_snap(i, j + 1, k), &
                                      tr_snap(i, j, k), tr_snap(i, j - 1, k), &
                                      tr_snap(i, j - 2, k), sigma)
                     case (2)
                        tr_face_val = weno5_face_swept( &
                                      tr_snap(i, j + 3, k), tr_snap(i, j + 2, k), &
                                      tr_snap(i, j + 1, k), &
                                      tr_snap(i, j, k), tr_snap(i, j - 1, k), sigma)
                     case (1)
                        tr_face_val = plm_face_swept( &
                                      tr_snap(i, j + 2, k), tr_snap(i, j + 1, k), &
                                      tr_snap(i, j, k), sigma)
                     case default
                        tr_face_val = tr_snap(i, j + 1, k)  ! donor
                     end select
                     flux_top = f_layer*tr_face_val
                  else
                     flux_top = 0.0_wp
                  end if

                  ! Conservative update
                  hTr_layer(i, j, k) = hTr_layer(i, j, k) &
                                       - dt*((flux_right - flux_left)*inv_dx &
                                             + (flux_top - flux_bot)*inv_dy)

               end do  ! k loop
            end if  ! h_col > ML_MIN_DEPTH
         end do
         end do
         !$omp end parallel do

      end associate

   end subroutine ml_advect_tracer_weno

end module rki_recon_weno
