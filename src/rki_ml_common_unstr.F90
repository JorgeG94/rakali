!! Shared constants for the unstructured multilayer kernels.
module rki_ml_common_unstr
   !! Parameters shared across the unstructured multilayer kernels (dynamics,
   !! tracers, pressure, EOS, viscosity). Mirrors structured `rki_ml_common`.
   use rki_constants, only: wp
   implicit none
   private

   public :: ML_MIN_DEPTH_UNSTR
   public :: VANISHING_LAYER_TOL_UNSTR
   public :: BPG_NZ_MAX
   public :: bpg_nz_required, bpg_nz_is_sufficient

   ! Minimum total column depth for vertical operations (skip thin columns)
   real(wp), parameter :: ML_MIN_DEPTH_UNSTR = 0.5_wp

   ! Vanishing-layer floor for VCOORD_ZSTAR_FULL; see structured
   ! VANISHING_LAYER_TOL in rki_ml_common.F90 for semantics.
   real(wp), parameter :: VANISHING_LAYER_TOL_UNSTR = 1.5e-4_wp

#ifndef RKI_BPG_NZ_MAX
#define RKI_BPG_NZ_MAX 20
#endif
   integer, parameter :: BPG_NZ_MAX = RKI_BPG_NZ_MAX
      !! Per-column stack-workspace bound for the unstructured multilayer
      !! kernels (BPG, vertical exchange, vertical viscosity, KPP, k-eps,
      !! tracer advection).  Fixed size so the compiler can emit fixed-size
      !! thread-local storage on GPU.  DISTINCT from the structured/ocean
      !! `NZ_STACK_MAX` in `rki_constants`: the unstructured path is
      !! thread-per-cell with many such locals per kernel frame, so this
      !! bound stays deliberately small — raising it multiplies the
      !! per-thread local memory of every unstructured column kernel.
      !! Build-configurable via `-DRAKALI_BPG_NZ_MAX=<n>` (CMake cache
      !! variable, plumbed as the `RKI_BPG_NZ_MAX` preprocessor define;
      !! default 20).
      !!
      !! **The requirement is `BPG_NZ_MAX >= nz + 1`** — see
      !! `bpg_nz_required` below for the derivation and the survey of
      !! every consumer.  Enforced fail-loud by `validate_config`.

contains

   pure function bpg_nz_required(nz) result(req)
      !! Smallest `BPG_NZ_MAX` that safely covers an unstructured
      !! multilayer run of `nz` layers.
      !!
      !! `nz + 1`.  Derived from a survey (2026-08) of every
      !! `BPG_NZ_MAX`-dimensioned local in the unstructured kernels; the
      !! binding consumer is the one that indexes a plainly-`BPG_NZ_MAX`
      !! array to `nz + 1` rather than the `BPG_NZ_MAX + 1` an interface
      !! array would use:
      !!
      !!   * `rki_ml_keps_unstr` `ml_keps_update_unstr_impl`
      !!     `S2`/`N2`/`dz_if`/`a`/`b`/`c`/`rhs`/`cp` — written and read
      !!     at `ni = nz + 1` (interface-indexed; `keps_tridiag_solve_u`
      !!     touches rows `1..ni` of all five solver arrays).
      !!
      !! Every other consumer needs only `nz`:
      !!
      !!   * `rki_ml_pressure_unstr` `p_int`/`p_int_nbr` — max index `nz`
      !!     (the interface arrays `z_if`/`p_if`/`zn_if_e`/`pn_if_e` are
      !!     declared `BPG_NZ_MAX + 1` and carry their own slack).
      !!   * `rki_ml_tracers_unstr` `v_layer`/`hf_layer`/`frac` — max
      !!     `nz`; `flux_iface` — max `nz - 1`.
      !!   * `rki_ml_dynamics_unstr`, `rki_ml_vertical_viscosity_unstr`,
      !!     `rki_ml_kpp_unstr` `flux_*_iface` — max `nz - 1`.
      !!   * `rki_ml_keps_unstr` momentum/tracer applies
      !!     `a`/`b`/`c`/`ru`/`rv`/`r`/`cp` — max `nz`.
      integer, intent(in) :: nz
      integer :: req

      req = nz + 1
   end function bpg_nz_required

   pure function bpg_nz_is_sufficient(nz) result(ok)
      !! `.true.` when the COMPILED `BPG_NZ_MAX` covers `nz` layers on
      !! the unstructured multilayer path.  Callers that get `.false.`
      !! must refuse the run: exceeding the bound is NOT a crash — the
      !! unstructured BPG kernels silently return (no baroclinic
      !! pressure gradient at all) and the unguarded column kernels
      !! silently overrun thread-local storage.
      integer, intent(in) :: nz
      logical :: ok

      ok = bpg_nz_required(nz) <= BPG_NZ_MAX
   end function bpg_nz_is_sufficient

end module rki_ml_common_unstr
