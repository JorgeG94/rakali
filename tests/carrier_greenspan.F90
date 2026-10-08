!! Carrier-Greenspan periodic runup on a plane beach.
module carrier_greenspan
   !! Closed-form oracle for wave runup on a uniform slope.
   !!
   !! ## The solution
   !!
   !! On a plane beach of slope `s`, with `X` the distance OFFSHORE from the
   !! still-water shoreline, the still depth is `h = s X`. Time-harmonic
   !! linear shallow water reduces to
   !!
   !!     d/dX ( X d(eta)/dX ) = (w^2 / (g s)) eta
   !!
   !! whose solution regular at the shoreline is a Bessel function of order
   !! zero in the stretched coordinate
   !!
   !!     xi(X) = 2 w sqrt( X / (g s) ),     eta(X, t) = A0 J0(xi) cos(w t)
   !!
   !! `A0` is therefore the **vertical runup amplitude** directly: `J0(0) = 1`,
   !! so `A0` is the shoreline elevation swing. The horizontal shoreline
   !! excursion is `A0 / s`.
   !!
   !! ## Why the linear solution is the right oracle here
   !!
   !! Carrier & Greenspan's (1958) central result is that the hodograph
   !! transform maps the NONLINEAR shallow-water runup problem onto this same
   !! linear equation — so the linear formula predicts runup correctly for
   !! finite-amplitude waves too, right up to the onset of breaking. Using it
   !! is not a small-amplitude approximation being quietly substituted for the
   !! real thing; it is the standard result. What the linear form does NOT
   !! give is the nonlinear distortion of the waveform between the toe and the
   !! shore, which is why this module is used for runup amplitude and phase
   !! rather than for a pointwise profile comparison.
   !!
   !! The non-breaking condition is `A0 w^2 / (g s) < 1`; `cg_breaking_number`
   !! returns that ratio so a fixture can assert it stays in range.
   !!
   !! ## The J1-zero trick
   !!
   !! Driving this with a time-dependent offshore boundary condition would
   !! entangle the test with the BC implementation. Instead, note that the
   !! velocity is proportional to `d(eta)/dX ~ -J1(xi)`, so an offshore WALL
   !! placed exactly at a zero of `J1` is a natural velocity node. The J0 mode
   !! is then an exact standing eigenmode of a closed basin — offshore wall at
   !! one end, moving shoreline at the other — needing no forcing at all: seed
   !! it and it oscillates. `CG_J1_ZERO1` is the first non-trivial zero.
   !!
   !! Bessel functions come from the Fortran 2008 intrinsics `bessel_j0` /
   !! `bessel_j1`, which gfortran 15.1 and nvfortran 26.5 evaluate identically
   !! (verified: `J0(3.8317059702) = -0.402759396` on both), so no series
   !! implementation is needed and there is no cross-toolchain drift.
   use rki_constants, only: wp, GRAVITY
   implicit none
   private

   public :: cg_xi
   public :: cg_eta
   public :: cg_omega_for_wall
   public :: cg_shoreline_elevation
   public :: cg_breaking_number
   public :: CG_J1_ZERO1

   real(wp), parameter :: CG_J1_ZERO1 = 3.8317059702075123_wp
      !! First non-trivial zero of J1. An offshore wall placed where
      !! `xi = CG_J1_ZERO1` is an exact velocity node.

contains

   pure function cg_xi(x_star, omega, slope) result(xi)
      !! Stretched coordinate `xi = 2 w sqrt(X / (g s))`.
      real(wp), intent(in) :: x_star
         !! Distance offshore from the still-water shoreline (m), >= 0
      real(wp), intent(in) :: omega
         !! Angular frequency (rad/s)
      real(wp), intent(in) :: slope
         !! Beach slope (rise/run), > 0
      real(wp) :: xi

      xi = 2.0_wp*omega*sqrt(max(x_star, 0.0_wp)/(GRAVITY*slope))
   end function cg_xi

   pure function cg_omega_for_wall(x_wall, slope, xi_wall) result(omega)
      !! Frequency that puts the offshore wall at a prescribed `xi` — invert
      !! `xi = 2 w sqrt(X/(g s))` for `w`. Pass `CG_J1_ZERO1` to make the wall
      !! a velocity node and the J0 mode an exact eigenmode.
      real(wp), intent(in) :: x_wall
         !! Distance of the wall offshore from the still-water shoreline (m)
      real(wp), intent(in) :: slope
      real(wp), intent(in) :: xi_wall
      real(wp) :: omega

      omega = 0.5_wp*xi_wall*sqrt(GRAVITY*slope/x_wall)
   end function cg_omega_for_wall

   pure function cg_eta(a0, x_star, t, omega, slope) result(eta)
      !! Free-surface elevation above the still level (m).
      real(wp), intent(in) :: a0
         !! Vertical runup amplitude (m) — the shoreline elevation swing
      real(wp), intent(in) :: x_star, t, omega, slope
      real(wp) :: eta

      eta = a0*bessel_j0(cg_xi(x_star, omega, slope))*cos(omega*t)
   end function cg_eta

   pure function cg_shoreline_elevation(a0, t, omega) result(eta_s)
      !! Shoreline elevation, `A0 cos(w t)` — the `X -> 0` limit where
      !! `J0(0) = 1`. Horizontal excursion is this divided by the slope.
      real(wp), intent(in) :: a0, t, omega
      real(wp) :: eta_s

      eta_s = a0*cos(omega*t)
   end function cg_shoreline_elevation

   pure function cg_breaking_number(a0, omega, slope) result(bn)
      !! `A0 w^2 / (g s)`. The Carrier-Greenspan solution stays single-valued
      !! (non-breaking) while this is below 1; at 1 the shoreline acceleration
      !! reaches gravity and the hodograph map folds.
      real(wp), intent(in) :: a0, omega, slope
      real(wp) :: bn

      bn = a0*omega*omega/(GRAVITY*slope)
   end function cg_breaking_number

end module carrier_greenspan
