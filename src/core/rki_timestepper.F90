!! Explicit SSP Runge-Kutta time-integrator selection for the coastal solver.
module rki_timestepper
   !! Holds the Shu-Osher coefficient table for the selectable explicit
   !! time integrator (`time_scheme`).  Every scheme is built from the same
   !! forward-Euler stage primitive: stage `i` advances the running iterate by
   !! `beta(i)*dt`, then blends with the start-of-step state `U^n`:
   !!
   !!     E       = U^(i-1) + beta(i)*dt * L(U^(i-1))     ! Euler step
   !!     U^(i)   = alpha(i)*U^n + (1 - alpha(i))*E       ! convex blend
   !!
   !! All four schemes are strong-stability-preserving (SSP) so the per-stage
   !! convex structure keeps the wetting/drying positivity of forward Euler.
   !! `ssp_coeff` is the factor by which the stable timestep exceeds the
   !! forward-Euler CFL limit (1 for euler/rk2/rk3, 2 for rk43).
   !!
   !! Data, not polymorphism: a concrete container + an enum, mirroring the
   !! `vcoord_type` / EOS-variant / `use_hllc` dispatch idiom.  The type is
   !! host-side orchestration only (its weights reach kernels as scalars by
   !! value) so it is never device-mapped.
   use rki_constants, only: wp
   implicit none
   private

   public :: timestepper_t, make_timestepper, parse_time_scheme
   public :: TIME_EULER, TIME_RK2, TIME_RK3, TIME_RK43, TIME_SEMI_IMPLICIT, MAX_STAGES

   integer, parameter :: TIME_EULER = 1
      !! Forward Euler, SSPRK(1,1): 1st order, 1 stage.
   integer, parameter :: TIME_RK2 = 2
      !! Heun / SSPRK(2,2): 2nd order, 2 stages (the default).
   integer, parameter :: TIME_RK3 = 3
      !! Shu-Osher SSPRK(3,3): 3rd order, 3 stages.
   integer, parameter :: TIME_RK43 = 4
      !! SSPRK(4,3): 3rd order, 4 stages, SSP coefficient 2 (2x timestep).
   integer, parameter :: TIME_SEMI_IMPLICIT = 5
      !! Semi-implicit free-surface (Casulli): implicit gravity wave via a CG Helmholtz
      !! eta-solve. Not an SSP stage table — dispatched separately in solver_step.

   integer, parameter :: MAX_STAGES = 4
      !! Largest stage count across the supported schemes.

   type :: timestepper_t
      !! Concrete Shu-Osher coefficient table for one scheme.  Host-side only.
      integer  :: scheme = TIME_RK2
         !! Active scheme tag (`TIME_*`).
      integer  :: n_stages = 2
         !! Number of forward-Euler stages.
      real(wp) :: ssp_coeff = 1.0_wp
         !! Stable-dt multiplier vs the forward-Euler CFL limit.
      real(wp) :: alpha(MAX_STAGES) = 0.0_wp
         !! Per-stage blend weight on `U^n` (combine = alpha*U^n + (1-alpha)*E).
      real(wp) :: beta(MAX_STAGES) = 0.0_wp
         !! Per-stage Euler-step fraction of `dt`.
      real(wp) :: c(MAX_STAGES) = 0.0_wp
         !! Per-stage time abscissa (fraction of `dt`) for time-dependent BCs.
   end type timestepper_t

contains

   pure function make_timestepper(scheme) result(ts)
      !! Build the coefficient table for a scheme tag.
      integer, intent(in) :: scheme
      type(timestepper_t) :: ts

      ts%scheme = scheme
      select case (scheme)
      case (TIME_EULER)
         ts%n_stages = 1
         ts%ssp_coeff = 1.0_wp
         ts%alpha(1) = 0.0_wp
         ts%beta(1) = 1.0_wp
         ts%c(1) = 0.0_wp
      case (TIME_RK3)
         ts%n_stages = 3
         ts%ssp_coeff = 1.0_wp
         ts%alpha(1:3) = [0.0_wp, 0.75_wp, 1.0_wp/3.0_wp]
         ts%beta(1:3) = [1.0_wp, 1.0_wp, 1.0_wp]
         ts%c(1:3) = [0.0_wp, 1.0_wp, 0.5_wp]
      case (TIME_RK43)
         ts%n_stages = 4
         ts%ssp_coeff = 2.0_wp
         ts%alpha(1:4) = [0.0_wp, 0.0_wp, 2.0_wp/3.0_wp, 0.0_wp]
         ts%beta(1:4) = [0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp]
         ts%c(1:4) = [0.0_wp, 0.5_wp, 1.0_wp, 0.5_wp]
      case (TIME_SEMI_IMPLICIT)
         ! Semi-implicit has NO SSP stage table — solver_step dispatches it
         ! to si_step before any coefficient table is consulted.  Reaching
         ! here means a new call site asked for stage coefficients that do
         ! not exist; fail loud rather than silently handing back RK2.
         error stop "make_timestepper: time_scheme='semi_implicit' is not "// &
            "stage-table-driven; it dispatches in solver_step (si_step)."
      case default  ! TIME_RK2 (Heun)
         ts%scheme = TIME_RK2
         ts%n_stages = 2
         ts%ssp_coeff = 1.0_wp
         ts%alpha(1:2) = [0.0_wp, 0.5_wp]
         ts%beta(1:2) = [1.0_wp, 1.0_wp]
         ts%c(1:2) = [0.0_wp, 1.0_wp]
      end select
   end function make_timestepper

   pure function parse_time_scheme(name) result(scheme)
      !! Map a namelist string to a `TIME_*` tag.  Unrecognised -> TIME_RK2.
      character(len=*), intent(in) :: name
      integer :: scheme

      select case (trim(adjustl(name)))
      case ("euler")
         scheme = TIME_EULER
      case ("rk3")
         scheme = TIME_RK3
      case ("rk43")
         scheme = TIME_RK43
      case ("semi_implicit")
         scheme = TIME_SEMI_IMPLICIT
      case default  ! "rk2" and anything unrecognised
         scheme = TIME_RK2
      end select
   end function parse_time_scheme

end module rki_timestepper
