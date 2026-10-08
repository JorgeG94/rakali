!! Semi-implicit settings shared by the explicit config path and the
!! unstructured semi-implicit solver.
module rki_si_common
   !! Sub-bottom fill modes for the Casulli wet/dry subgrid volume
   !! (`&si_nml wetdry_bsub_fill`), their parser, and the conservation
   !! conflict predicate used by the state-enable guards.
   use pic_logger, only: global_logger
   implicit none
   private

   public :: parse_si_bsub_fill
   public :: si_bsub_ml_conflict
   public :: SI_BSUB_FILL_FLAT, SI_BSUB_FILL_LINEAR

   integer, parameter :: SI_BSUB_FILL_FLAT = 0
      !! Sub-bottom fill: all sub-bottoms = the plain bed `b` (bit-identical
      !! to n_sub=1 regardless of n_sub).  The default.
   integer, parameter :: SI_BSUB_FILL_LINEAR = 1
      !! Sub-bottom fill: an analytic per-cell linear sub-cell bed ramp
      !! across `[b-half, b+half]`, `half = 0.5*max` neighbour bed diff.

contains

   function parse_si_bsub_fill(tag) result(mode)
      !! Map the `&si_nml wetdry_bsub_fill` string to its `SI_BSUB_FILL_*`
      !! enum.  Fail-loud allowlist (shared by both backends).
      character(len=*), intent(in) :: tag
      integer :: mode

      select case (trim(adjustl(tag)))
      case ("flat", "FLAT", "")
         mode = SI_BSUB_FILL_FLAT
      case ("linear", "LINEAR")
         mode = SI_BSUB_FILL_LINEAR
      case default
         call global_logger%error("si_wetdry_bsub_fill: unknown fill mode '"// &
                                  trim(adjustl(tag))//"' (allowed: flat, linear)")
         error stop "si_wetdry_bsub_fill: unknown fill mode"
      end select
   end function parse_si_bsub_fill

   pure function si_bsub_ml_conflict(bsub_fill, n_sub, use_multilayer, nz) &
      result(conflict)
      !! The CONSERVATION conflict predicate (shared by both backends' enable
      !! guards): the `linear` subgrid b_sub fill (n_sub >= 2) is incompatible
      !! with the multilayer (nz > 1) wet/dry step.  The free-surface solve
      !! conserves the subgrid wetted volume v_cell(eta) = mean_s
      !! max(eta - b_sub_s, 0), but the ML thickness update + recovery
      !! Dirichlet ghost consume max(eta - b, 0) on the PLAIN bed; under
      !! `linear` these diverge across the partial-flooding band => mass leaks
      !! at the wet/dry front.  `flat + ML` (b_sub == b, conservative) and
      !! `linear + 2D` (barotropic uses v_cell throughout) both return
      !! .false.  A v_cell-consistent ML update is a documented deferred
      !! limitation; until then callers must reject a .true. combination
      !! fail-loud.  Pure so tests can assert it directly (the enable path
      !! error-stops, which cannot be caught in-process).
      integer, intent(in) :: bsub_fill
         !! `SI_BSUB_FILL_*` fill mode.
      integer, intent(in) :: n_sub
         !! Subgrid sub-bottoms per cell.
      logical, intent(in) :: use_multilayer
         !! Multilayer (3D) mode on.
      integer, intent(in) :: nz
         !! Layer count (only matters together with use_multilayer).
      logical :: conflict

      conflict = (bsub_fill == SI_BSUB_FILL_LINEAR) .and. (n_sub >= 2) .and. &
                 use_multilayer .and. (nz > 1)
   end function si_bsub_ml_conflict

end module rki_si_common
