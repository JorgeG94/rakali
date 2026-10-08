!! Thacker parabolic bowl exact solution for wetting/drying validation
module thacker_bowl
   !! Provides the exact analytical solution for the Thacker (1981)
   !! parabolic bowl problem: periodic sloshing in a parabolic basin
   !! with a moving shoreline.
   !!
   !! Bathymetry: B(r) = h0 * (r^2 / a^2 - 1)
   !! where h0 is the depth at centre and a is the bowl radius.
   !!
   !! The exact solution has the water surface as a tilted plane that
   !! oscillates with period T = 2*pi / omega, where
   !! omega = sqrt(2*g*h0) / a.
   use rki_constants, only: wp, GRAVITY
   implicit none
   private

   public :: thacker_omega
   public :: thacker_eta_exact
   public :: thacker_u_exact

contains

   pure function thacker_omega(h0, a) result(omega)
      !! Natural oscillation frequency of the parabolic bowl
      real(wp), intent(in) :: h0
         !! Depth at centre of bowl (m)
      real(wp), intent(in) :: a
         !! Bowl radius (m)
      real(wp) :: omega

      omega = sqrt(2.0_wp*GRAVITY*h0)/a

   end function thacker_omega

   pure function thacker_eta_exact(x, y, t, h0, a, eta_amp) result(eta)
      !! Exact free surface elevation at position (x,y) and time t
      !!
      !! The solution is a planar surface oscillating in x:
      !! eta(x,t) = h0 * [sqrt(1 - A^2) / (1 - A*cos(omega*t)) - 1]
      !!          + h0 * A * cos(omega*t) * x / a / (1 - A*cos(omega*t))
      !! where A = eta_amp and omega = sqrt(2*g*h0)/a.
      !! Simplified from Thacker (1981).
      real(wp), intent(in) :: x
         !! x-position relative to bowl centre (m)
      real(wp), intent(in) :: y
         !! y-position relative to bowl centre (m)
      real(wp), intent(in) :: t
         !! Time (s)
      real(wp), intent(in) :: h0
         !! Central depth (m)
      real(wp), intent(in) :: a
         !! Bowl radius (m)
      real(wp), intent(in) :: eta_amp
         !! Non-dimensional amplitude parameter (0 < A < 1)
      real(wp) :: eta

      real(wp) :: omega, cot, r_sq, denom

      omega = thacker_omega(h0, a)
      cot = cos(omega*t)
      denom = 1.0_wp - eta_amp*cot
      r_sq = x*x + y*y

      ! Free surface: eta = B + h, with B = h0*(r^2/a^2 - 1)
      ! The exact water depth is:
      ! h = h0 * [sqrt(1-A^2)/(1-A*cos(wt)) - 1 - r^2/a^2
      !          + A*cos(wt)*(2*x)/a / (1-A*cos(wt))^2 * ...]
      ! Using the simpler 1D sloshing form (x-direction only):
      eta = h0*(sqrt(1.0_wp - eta_amp*eta_amp)/denom - 1.0_wp) &
            - h0*r_sq/(a*a) &
            + h0*eta_amp*cot*2.0_wp*x/(a*denom)

   end function thacker_eta_exact

   pure function thacker_u_exact(t, h0, a, eta_amp) result(u_val)
      !! Exact x-velocity (spatially uniform for the planar solution)
      real(wp), intent(in) :: t
         !! Time (s)
      real(wp), intent(in) :: h0
         !! Central depth (m)
      real(wp), intent(in) :: a
         !! Bowl radius (m)
      real(wp), intent(in) :: eta_amp
         !! Non-dimensional amplitude parameter
      real(wp) :: u_val

      real(wp) :: omega, sot, denom

      omega = thacker_omega(h0, a)
      sot = sin(omega*t)
      denom = 1.0_wp - eta_amp*cos(omega*t)

      u_val = eta_amp*omega*a*sot/denom

   end function thacker_u_exact

end module thacker_bowl
