!! Shared k-epsilon stability functions (c_mu, c_mu') — backend-agnostic.
module rki_keps_stability
   !! Stratification-aware k-epsilon stability functions (S^2/N^2 -> c_mu/c_mu'),
   !! shared by structured (`rki_ml_keps`) and unstructured (`rki_ml_keps_unstr`).
   !! CONSTANT: c_mu = 0.09. GALPERIN: Galperin et al. (1988) / Mellor-Yamada
   !! level-2.5 algebraic form. CANUTO: Canuto et al. (2001) Model A in the
   !! rational alpha_M/alpha_N form of Umlauf & Burchard (2003).
   use rki_constants, only: wp, KEPS_STAB_GALPERIN, KEPS_STAB_CANUTO
   implicit none
   private

   public :: keps_stability_funcs

   ! Neutral c_mu + turbulent Prandtl number (the CONSTANT scheme + the guards).
   ! Self-contained here so the routine carries no host dependency on a kernel
   ! module; the kernels keep their own copies for their wall/clamp constants.
   real(wp), parameter :: CMU = 0.09_wp
   real(wp), parameter :: PR_T = 1.0_wp
   real(wp), parameter :: EPS_MIN = 1.0e-12_wp

contains

   pure subroutine keps_stability_funcs(scheme, s2, n2, k, e, cmu_o, cmup_o)
      !$acc routine seq
      !! Stability functions (c_mu, c_mu') for the chosen scheme, given shear S^2
      !! and buoyancy N^2 at the interface.
      !! Result dummies are cmu_o/cmup_o, NOT cmu/cmup: a `cmu` dummy shadows the
      !! module parameter CMU (case-insensitive) -> self-assignment of undefined
      !! intent(out) -> nu_t floored on GPU/gfortran.
      integer, intent(in) :: scheme
      real(wp), intent(in) :: s2, n2, k, e
      real(wp), intent(out) :: cmu_o, cmup_o
      ! Galperin (1988) / Mellor-Yamada level-2.5 constants.
      real(wp), parameter :: A1 = 0.92_wp, A2 = 0.74_wp, B1 = 16.6_wp
      real(wp), parameter :: B2 = 10.1_wp, C1 = 0.08_wp
      real(wp), parameter :: GH_MIN = -0.28_wp, GH_MAX = 0.0233_wp
      ! Canuto (2001) Model A structural constants (L*) + derived numerator
      ! (CA_S*) and denominator (CA_D*) coefficients of the rational form.
      real(wp), parameter :: L1 = 0.107_wp, L2 = 0.0032_wp, L3 = 0.0864_wp
      real(wp), parameter :: L4 = 0.12_wp, L5 = 11.9_wp, L6 = 0.4_wp
      real(wp), parameter :: L7 = 0.0_wp, L8 = 0.48_wp
      real(wp), parameter :: CA_S0 = 1.5_wp*L1*L5**2
      real(wp), parameter :: CA_S1 = -L4*(L6 + L7) + 2.0_wp*L4*L5*(L1 - L2/3.0_wp - L3) &
                             + 1.5_wp*L1*L5*L8
      real(wp), parameter :: CA_S2 = -3.0_wp/8.0_wp*L1*(L6**2 - L7**2)
      real(wp), parameter :: CA_S4 = 2.0_wp*L5
      real(wp), parameter :: CA_S5 = 2.0_wp*L4
      real(wp), parameter :: CA_S6 = 2.0_wp/3.0_wp*L5*(3.0_wp*L3**2 - L2**2) &
                             - 0.5_wp*L5*L1*(3.0_wp*L3 - L2) + 0.75_wp*L1*(L6 - L7)
      real(wp), parameter :: CA_D0 = 3.0_wp*L5**2
      real(wp), parameter :: CA_D1 = L5*(7.0_wp*L4 + 3.0_wp*L8)
      real(wp), parameter :: CA_D2 = L5**2*(3.0_wp*L3**2 - L2**2) - 0.75_wp*(L6**2 - L7**2)
      real(wp), parameter :: CA_D3 = L4*(4.0_wp*L4 + 3.0_wp*L8)
      real(wp), parameter :: CA_D4 = L4*(L2*L6 - 3.0_wp*L3*L7 - L5*(L2**2 - L3**2)) &
                             + L5*L8*(3.0_wp*L3**2 - L2**2)
      real(wp), parameter :: CA_D5 = 0.25_wp*(L2**2 - 3.0_wp*L3**2)*(L6**2 - L7**2)
      real(wp), parameter :: AN_MIN = -1.5_wp, AN_MAX = 19.0_wp, AM_MAX = 200.0_wp
      real(wp) :: alpha_n, gh, sh, sm, alpha_m, tn, ts, denom, eok2

      select case (scheme)
      case (KEPS_STAB_GALPERIN)
         alpha_n = (k/max(e, EPS_MIN))**2*n2
         gh = min(max(-(4.0_wp/(B1*B1))*alpha_n, GH_MIN), GH_MAX)
         sh = A2*(1.0_wp - 6.0_wp*A1/B1) &
              /(1.0_wp - (18.0_wp*A1*A2 + 3.0_wp*A2*B2)*gh)
         sm = (A1*(1.0_wp - 3.0_wp*C1 - 6.0_wp*A1/B1) &
               + sh*(18.0_wp*A1*A1 + 9.0_wp*A1*A2)*gh) &
              /(1.0_wp - 9.0_wp*A1*A2*gh)
         cmu_o = 4.0_wp*sm/B1
         cmup_o = 4.0_wp*sh/B1
      case (KEPS_STAB_CANUTO)
         eok2 = (k/max(e, EPS_MIN))**2
         alpha_n = min(max(eok2*n2, AN_MIN), AN_MAX)
         alpha_m = min(max(eok2*s2, 0.0_wp), AM_MAX)
         tn = 4.0_wp*alpha_n
         ts = 4.0_wp*alpha_m
         denom = CA_D0 + CA_D1*tn + CA_D2*ts + CA_D3*tn*tn + CA_D4*tn*ts + CA_D5*ts*ts
         cmu_o = 2.0_wp*(CA_S0 + CA_S1*tn + CA_S2*ts)/denom
         cmup_o = 2.0_wp*(CA_S4 + CA_S5*tn + CA_S6*ts)/denom
      case default
         cmu_o = CMU
         cmup_o = CMU/PR_T
      end select
   end subroutine keps_stability_funcs

end module rki_keps_stability
