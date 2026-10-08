!! Shared state for the structured multilayer kernels.
module rki_ml_common
   !! State and helpers shared by every multilayer sub-module (dynamics,
   !! tracers, pressure, EOS, momentum).
   use rki_constants, only: wp
   implicit none

   ! Minimum total water depth for vertical layer operations. Columns shallower
   ! than this skip vertical exchange, mixing, tracer advection, and BPG.
   real(wp), protected :: ML_MIN_DEPTH = 0.5_wp
   ! [acc->omp-cpu] dropped device directive (host no-op): declare

   ! Vanishing-layer threshold for VCOORD_ZSTAR_FULL. A layer with
   ! h_layer < VANISHING_LAYER_TOL is treated as a collapsed structural
   ! placeholder by the multilayer operators (BPG/EOS/PP81/friction/
   ! truncation/tracers). 1.5× margin over zstar_h_min (1e-4 m) for round-off.
   ! Other coord types never go below this, so the guards are transparent.
   real(wp), parameter :: VANISHING_LAYER_TOL = 1.5e-4_wp

   public :: ML_MIN_DEPTH
   public :: VANISHING_LAYER_TOL
   public :: ml_set_min_depth
   private

contains

   subroutine ml_set_min_depth(depth)
      !! Set the minimum total water depth for vertical layer operations.
      !! Call once at startup from config.
      real(wp), intent(in) :: depth
      ML_MIN_DEPTH = depth
      ! [acc->omp-cpu] dropped device directive (host no-op): update
   end subroutine ml_set_min_depth

end module rki_ml_common
