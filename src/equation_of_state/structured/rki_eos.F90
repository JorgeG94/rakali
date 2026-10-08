!! Ocean equation-of-state state + kernel.
module rki_eos
   !! Holds the EOS variant tag + scalar EOS coefficients (linear
   !! Boussinesq) and the kernels that convert (T, S) -> ρ on the
   !! multilayer C-grid.  Coastal uses a sibling linear EOS in
   !! `rki_ml_eos`; the ocean path needs a real nonlinear EOS along
   !! the FV pressure-gradient integration path.  Default is Wright
   !! (1997) — the same EOS MOM6 uses by default — chosen for cost
   !! (one rational expression, no LUT) and accuracy in the open
   !! ocean.  TEOS-10 is a future option.
   !!
   !! Phase 5c status: linear branch live; Wright + LUT scratch
   !! still pending.  The kernel follows the outer-shim + flat-impl
   !! pattern: an outer shim pulls the registered S, T tracer arrays
   !! off the owning state and forwards them as bare 3D arrays to the
   !! inner `_impl` routine.  This avoids the NVHPC stdpar deep-deref
   !! issue with array-of-derived-types tracer registries (see
   !! CLAUDE.md feedback_nvhpc_impl_inlining).  This module is
   !! regime-agnostic — it carries NO ocean-state dependency; the
   !! ocean C-grid shim `ocean_eos_compute` lives in
   !! `rki_ocean_eos_compute` (core/ocean/state/), and the coastal
   !! shim in `rki_ml_eos`.
   use rki_constants, only: wp, H_VANISHED
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: eos_t
   public :: eos_compute_arrays
   public :: eos_validate
   public :: eos_wright_pgf_column_sweep_impl
   public :: eos_specvol_derivs
   public :: eos_density_point
   public :: eos_freezing_point
   public :: parse_eos_variant

   integer, parameter, public :: EOS_VARIANT_LINEAR = 1
      !! Linear T/S (debug / lock-exchange / Eady).
   integer, parameter, public :: EOS_VARIANT_WRIGHT_97 = 2
      !! Wright (1997) "An Equation of State for Use in Ocean
      !! Models: Eckart's Formula Revisited" (JAOT 14:735-740).
      !! Rational form: ρ = (P + p_0) / (λ + α_0*(P + p_0)) with
      !! cubic / linear polynomials in (T, S).  Matches the MOM6
      !! default EOS.  Phase Tier-1 evaluates at P = `eos%p_ref`
      !! (0 by default → surface ρ); Phase 5d adds the FV-PGF
      !! integration path with in-situ pressure.
   integer, parameter, public :: EOS_VARIANT_TEOS10 = 3
      !! TEOS-10 (future).
   integer, parameter, public :: EOS_VARIANT_ROQUET_SPV = 4
      !! Roquet et al. (2015) "Accurate polynomial expressions for the
      !! density and specific volume of seawater using the TEOS-10
      !! standard." Ocean Modelling 90:29-43 — the specific-volume
      !! (SpV) polynomial variant (NEMO / MOM6 `Roquet_SpV`).  ~75-term
      !! polynomial in (CT, SA, p) giving TEOS-10-class accuracy as
      !! self-contained `parameter` arithmetic (no GSW LUT, two sqrt,
      !! no iteration).  Device-clean.  The coefficient assembly + the
      !! PT->CT conversion poly were transcribed from the validated
      !! Python prototype (`local_archive/prototypes/roquet_spv_eos.py`,
      !! 128/128 identical to the published MOM6 `MOM_EOS_Roquet_SpV`
      !! transcription); MOM6 is credited as the published-coefficient
      !! reference.  CONVENTION (pinned 2026-06-16): consumers work in
      !! model (PT, SP); this branch converts SR = SP·(35.16504/35) for
      !! Reference Salinity and CT = ct_from_pt(SR, PT) locally, so it
      !! stays device-callable and returns derivatives w.r.t. the model
      !! (PT, SP) via the analytic chain rule.

   integer, parameter, public :: TS_POT_PRAC = 1
      !! Tracer T/S convention: potential temperature + practical
      !! salinity (the model-prognostic pair).  Identity for the
      !! linear and Wright (1997) branches — they consume (PT, SP)
      !! directly.  Default.
   integer, parameter, public :: TS_CONS_ABS = 2
      !! Reserved: conservative temperature + absolute salinity
      !! (CT, SA) — for a future TEOS-10/Roquet-CT branch that does
      !! the CT/SA->PT/SP conversion locally inside its `acc routine
      !! seq` body.  Not yet device-callable; consumers always pass
      !! the model-prognostic (T, S) regardless of convention.

   ! Seawater freezing-point (liquidus) coefficients — SIS2/MOM6 linear
   ! form `T_Freeze` (sea-ice PR 1, PLAN_SEA_ICE.md "Prerequisites").
   real(wp), parameter, public :: TFR_S_COEFF = -0.054_wp
      !! Liquidus slope dT_f/dS (degC per g/kg) — SIS2's `T_Freeze`
      !! μ = −0.054 °C/(g/kg); T_f(S=35) = −1.89 °C.
   real(wp), parameter, public :: TFR_P_COEFF = -7.53e-8_wp
      !! Pressure depression dT_f/dp (degC/Pa) — MOM6 `DTFREEZE_DP`
      !! reference value (−7.53e-8 °C/Pa ≈ −0.75 °C per 1000 dbar).

   ! Wright (1997) coefficients from Table A1 of the paper.  Units: SI
   ! throughout (T in degC, S in PSU, P in Pa, ρ in kg/m^3).
   real(wp), parameter :: WRIGHT_A0 = 7.057924e-4_wp
   real(wp), parameter :: WRIGHT_A1 = 3.480336e-7_wp
   real(wp), parameter :: WRIGHT_A2 = -1.112733e-7_wp
   real(wp), parameter :: WRIGHT_B0 = 5.790749e8_wp
   real(wp), parameter :: WRIGHT_B1 = 3.516535e6_wp
   real(wp), parameter :: WRIGHT_B2 = -4.002714e4_wp
   real(wp), parameter :: WRIGHT_B3 = 2.084372e2_wp
   real(wp), parameter :: WRIGHT_B4 = 5.944068e5_wp
   real(wp), parameter :: WRIGHT_B5 = -9.643486e3_wp
   real(wp), parameter :: WRIGHT_C0 = 1.704853e5_wp
   real(wp), parameter :: WRIGHT_C1 = 7.904722e2_wp
   real(wp), parameter :: WRIGHT_C2 = -7.984422e0_wp
   real(wp), parameter :: WRIGHT_C3 = 5.140652e-2_wp
   real(wp), parameter :: WRIGHT_C4 = -2.302158e2_wp
   real(wp), parameter :: WRIGHT_C5 = -3.079464e0_wp

   ! ======================================================================
   ! Roquet et al. (2015) specific-volume (SpV) polynomial coefficients.
   ! ----------------------------------------------------------------------
   ! Citation of record:
   !   Roquet, F., Madec, G., McDougall, T. J., Barker, P. M. (2015):
   !   "Accurate polynomial expressions for the density and specific volume
   !    of seawater using the TEOS-10 standard." Ocean Modelling 90:29-43.
   ! Published-coefficient reference (transcription source): the Roquet
   !   et al. (2015) TEOS-10 standard SpV coefficient set.
   !
   ! Every value below is transcribed verbatim from the VERIFIED Python
   ! prototype local_archive/prototypes/roquet_spv_eos.py (128/128 identical
   ! to MOM6).  The I_Ts / Pa2kb scalings are FOLDED INTO the coefficients
   ! here (exactly as the prototype does), so the point math feeds degC and
   ! Pa directly with no extra normalisation.
   !
   ! Normalisation constants (g/kg, degC, Pa).
   real(wp), parameter :: ROQ_PA2KB = 1.0e-8_wp
      !! Pa -> kbar.
   real(wp), parameter :: ROQ_RDELTAS = 24.0_wp
      !! Salinity offset before the sqrt (g/kg).
   real(wp), parameter :: ROQ_R1_S0 = 0.875_wp/35.16504_wp
      !! Inverse plausible salinity range (kg/g).
   real(wp), parameter :: ROQ_I_TS = 0.025_wp
      !! Inverse plausible temperature range (1/degC).
   real(wp), parameter :: ROQ_SR_FACTOR = 35.16504_wp/35.0_wp
      !! SP -> SR (Reference Salinity) conversion factor.

   ! Reference-profile (SV00p) pressure coefficients, in Pa-powers.
   real(wp), parameter :: ROQ_V00 = -4.4015007269e-05_wp*ROQ_PA2KB
   real(wp), parameter :: ROQ_V01 = 6.9232335784e-06_wp*ROQ_PA2KB**2
   real(wp), parameter :: ROQ_V02 = -7.5004675975e-07_wp*ROQ_PA2KB**3
   real(wp), parameter :: ROQ_V03 = 1.7009109288e-08_wp*ROQ_PA2KB**4
   real(wp), parameter :: ROQ_V04 = -1.6884162004e-08_wp*ROQ_PA2KB**5
   real(wp), parameter :: ROQ_V05 = 1.9613503930e-09_wp*ROQ_PA2KB**6

   ! SV(zs,zt,zp) term coefficients  SPV_abc * zs**a * zt**b * zp**c.
   real(wp), parameter :: SPV000 = 1.0772899069e-03_wp
   real(wp), parameter :: SPV100 = -3.1263658781e-04_wp
   real(wp), parameter :: SPV200 = 6.7615860683e-04_wp
   real(wp), parameter :: SPV300 = -8.6127884515e-04_wp
   real(wp), parameter :: SPV400 = 5.9010812596e-04_wp
   real(wp), parameter :: SPV500 = -2.1503943538e-04_wp
   real(wp), parameter :: SPV600 = 3.2678954455e-05_wp
   real(wp), parameter :: SPV010 = -1.4949652640e-05_wp*ROQ_I_TS
   real(wp), parameter :: SPV110 = 3.1866349188e-05_wp*ROQ_I_TS
   real(wp), parameter :: SPV210 = -3.8070687610e-05_wp*ROQ_I_TS
   real(wp), parameter :: SPV310 = 2.9818473563e-05_wp*ROQ_I_TS
   real(wp), parameter :: SPV410 = -1.0011321965e-05_wp*ROQ_I_TS
   real(wp), parameter :: SPV510 = 1.0751931163e-06_wp*ROQ_I_TS
   real(wp), parameter :: SPV020 = 2.7546851539e-05_wp*ROQ_I_TS**2
   real(wp), parameter :: SPV120 = -3.6597334199e-05_wp*ROQ_I_TS**2
   real(wp), parameter :: SPV220 = 3.4489154625e-05_wp*ROQ_I_TS**2
   real(wp), parameter :: SPV320 = -1.7663254122e-05_wp*ROQ_I_TS**2
   real(wp), parameter :: SPV420 = 3.5965131935e-06_wp*ROQ_I_TS**2
   real(wp), parameter :: SPV030 = -1.6506828994e-05_wp*ROQ_I_TS**3
   real(wp), parameter :: SPV130 = 2.4412359055e-05_wp*ROQ_I_TS**3
   real(wp), parameter :: SPV230 = -1.4606740723e-05_wp*ROQ_I_TS**3
   real(wp), parameter :: SPV330 = 2.3293406656e-06_wp*ROQ_I_TS**3
   real(wp), parameter :: SPV040 = 6.7896174634e-06_wp*ROQ_I_TS**4
   real(wp), parameter :: SPV140 = -8.7951832993e-06_wp*ROQ_I_TS**4
   real(wp), parameter :: SPV240 = 4.4249040774e-06_wp*ROQ_I_TS**4
   real(wp), parameter :: SPV050 = -7.2535743349e-07_wp*ROQ_I_TS**5
   real(wp), parameter :: SPV150 = -3.4680559205e-07_wp*ROQ_I_TS**5
   real(wp), parameter :: SPV060 = 1.9041365570e-07_wp*ROQ_I_TS**6
   real(wp), parameter :: SPV001 = -1.6889436589e-05_wp*ROQ_PA2KB
   real(wp), parameter :: SPV101 = 2.1106556158e-05_wp*ROQ_PA2KB
   real(wp), parameter :: SPV201 = -2.1322804368e-05_wp*ROQ_PA2KB
   real(wp), parameter :: SPV301 = 1.7347655458e-05_wp*ROQ_PA2KB
   real(wp), parameter :: SPV401 = -4.3209400767e-06_wp*ROQ_PA2KB
   real(wp), parameter :: SPV011 = 1.5355844621e-05_wp*(ROQ_I_TS*ROQ_PA2KB)
   real(wp), parameter :: SPV111 = 2.0914122241e-06_wp*(ROQ_I_TS*ROQ_PA2KB)
   real(wp), parameter :: SPV211 = -5.7751479725e-06_wp*(ROQ_I_TS*ROQ_PA2KB)
   real(wp), parameter :: SPV311 = 1.0767234341e-06_wp*(ROQ_I_TS*ROQ_PA2KB)
   real(wp), parameter :: SPV021 = -9.6659393016e-06_wp*(ROQ_I_TS**2*ROQ_PA2KB)
   real(wp), parameter :: SPV121 = -7.0686982208e-07_wp*(ROQ_I_TS**2*ROQ_PA2KB)
   real(wp), parameter :: SPV221 = 1.4488066593e-06_wp*(ROQ_I_TS**2*ROQ_PA2KB)
   real(wp), parameter :: SPV031 = 3.1134283336e-06_wp*(ROQ_I_TS**3*ROQ_PA2KB)
   real(wp), parameter :: SPV131 = 7.9562529879e-08_wp*(ROQ_I_TS**3*ROQ_PA2KB)
   real(wp), parameter :: SPV041 = -5.6590253863e-07_wp*(ROQ_I_TS**4*ROQ_PA2KB)
   real(wp), parameter :: SPV002 = 1.0500241168e-06_wp*ROQ_PA2KB**2
   real(wp), parameter :: SPV102 = 1.9600661704e-06_wp*ROQ_PA2KB**2
   real(wp), parameter :: SPV202 = -2.1666693382e-06_wp*ROQ_PA2KB**2
   real(wp), parameter :: SPV012 = -3.8541359685e-06_wp*(ROQ_I_TS*ROQ_PA2KB**2)
   real(wp), parameter :: SPV112 = 1.0157632247e-06_wp*(ROQ_I_TS*ROQ_PA2KB**2)
   real(wp), parameter :: SPV022 = 1.7178343158e-06_wp*(ROQ_I_TS**2*ROQ_PA2KB**2)
   real(wp), parameter :: SPV003 = -4.1503454190e-07_wp*ROQ_PA2KB**3
   real(wp), parameter :: SPV103 = 3.5627020989e-07_wp*ROQ_PA2KB**3
   real(wp), parameter :: SPV013 = -1.1293871415e-07_wp*(ROQ_I_TS*ROQ_PA2KB**3)

   ! dSV/dCT coefficient table (ALP = d/dzt of the SV table; b-power -> b*coef).
   real(wp), parameter :: ALP000 = SPV010, ALP100 = SPV110, ALP200 = SPV210
   real(wp), parameter :: ALP300 = SPV310, ALP400 = SPV410, ALP500 = SPV510
   real(wp), parameter :: ALP010 = 2.0_wp*SPV020, ALP110 = 2.0_wp*SPV120
   real(wp), parameter :: ALP210 = 2.0_wp*SPV220, ALP310 = 2.0_wp*SPV320
   real(wp), parameter :: ALP410 = 2.0_wp*SPV420
   real(wp), parameter :: ALP020 = 3.0_wp*SPV030, ALP120 = 3.0_wp*SPV130
   real(wp), parameter :: ALP220 = 3.0_wp*SPV230, ALP320 = 3.0_wp*SPV330
   real(wp), parameter :: ALP030 = 4.0_wp*SPV040, ALP130 = 4.0_wp*SPV140
   real(wp), parameter :: ALP230 = 4.0_wp*SPV240
   real(wp), parameter :: ALP040 = 5.0_wp*SPV050, ALP140 = 5.0_wp*SPV150
   real(wp), parameter :: ALP050 = 6.0_wp*SPV060
   real(wp), parameter :: ALP001 = SPV011, ALP101 = SPV111, ALP201 = SPV211
   real(wp), parameter :: ALP301 = SPV311
   real(wp), parameter :: ALP011 = 2.0_wp*SPV021, ALP111 = 2.0_wp*SPV121
   real(wp), parameter :: ALP211 = 2.0_wp*SPV221
   real(wp), parameter :: ALP021 = 3.0_wp*SPV031, ALP121 = 3.0_wp*SPV131
   real(wp), parameter :: ALP031 = 4.0_wp*SPV041
   real(wp), parameter :: ALP002 = SPV012, ALP102 = SPV112
   real(wp), parameter :: ALP012 = 2.0_wp*SPV022
   real(wp), parameter :: ALP003 = SPV013

   ! dSV/dSA coefficient table (BET; folds the 0.5*r1_S0 zs chain factor
   ! per-coef as the prototype does; the residual 1/zs is applied at the
   ! call site).
   real(wp), parameter :: BET000 = 0.5_wp*SPV100*ROQ_R1_S0, BET100 = SPV200*ROQ_R1_S0
   real(wp), parameter :: BET200 = 1.5_wp*SPV300*ROQ_R1_S0, BET300 = 2.0_wp*SPV400*ROQ_R1_S0
   real(wp), parameter :: BET400 = 2.5_wp*SPV500*ROQ_R1_S0, BET500 = 3.0_wp*SPV600*ROQ_R1_S0
   real(wp), parameter :: BET010 = 0.5_wp*SPV110*ROQ_R1_S0, BET110 = SPV210*ROQ_R1_S0
   real(wp), parameter :: BET210 = 1.5_wp*SPV310*ROQ_R1_S0, BET310 = 2.0_wp*SPV410*ROQ_R1_S0
   real(wp), parameter :: BET410 = 2.5_wp*SPV510*ROQ_R1_S0
   real(wp), parameter :: BET020 = 0.5_wp*SPV120*ROQ_R1_S0, BET120 = SPV220*ROQ_R1_S0
   real(wp), parameter :: BET220 = 1.5_wp*SPV320*ROQ_R1_S0, BET320 = 2.0_wp*SPV420*ROQ_R1_S0
   real(wp), parameter :: BET030 = 0.5_wp*SPV130*ROQ_R1_S0, BET130 = SPV230*ROQ_R1_S0
   real(wp), parameter :: BET230 = 1.5_wp*SPV330*ROQ_R1_S0
   real(wp), parameter :: BET040 = 0.5_wp*SPV140*ROQ_R1_S0, BET140 = SPV240*ROQ_R1_S0
   real(wp), parameter :: BET050 = 0.5_wp*SPV150*ROQ_R1_S0
   real(wp), parameter :: BET001 = 0.5_wp*SPV101*ROQ_R1_S0, BET101 = SPV201*ROQ_R1_S0
   real(wp), parameter :: BET201 = 1.5_wp*SPV301*ROQ_R1_S0, BET301 = 2.0_wp*SPV401*ROQ_R1_S0
   real(wp), parameter :: BET011 = 0.5_wp*SPV111*ROQ_R1_S0, BET111 = SPV211*ROQ_R1_S0
   real(wp), parameter :: BET211 = 1.5_wp*SPV311*ROQ_R1_S0
   real(wp), parameter :: BET021 = 0.5_wp*SPV121*ROQ_R1_S0, BET121 = SPV221*ROQ_R1_S0
   real(wp), parameter :: BET031 = 0.5_wp*SPV131*ROQ_R1_S0
   real(wp), parameter :: BET002 = 0.5_wp*SPV102*ROQ_R1_S0, BET102 = SPV202*ROQ_R1_S0
   real(wp), parameter :: BET012 = 0.5_wp*SPV112*ROQ_R1_S0
   real(wp), parameter :: BET003 = 0.5_wp*SPV103*ROQ_R1_S0

   ! PT->CT conversion (7-term gsw_CT_from_pt surface polynomial, in SR, PT).
   real(wp), parameter :: ROQ_CP0 = 3991.86795711963_wp
   real(wp), parameter :: ROQ_CT_SFAC = 0.0248826675584615_wp
      !! (35.16504/35)/40  [(g/kg)^-1] — normalises SR for the poly.

   type :: eos_t
      logical :: is_init = .false.
         !! True between `init` and `destroy`.  Prefer this to
         !! `allocated(...)` — tracks GPU device attachment too.
      integer  :: variant = EOS_VARIANT_LINEAR
         !! Active EOS variant.  Defaults to linear two-tracer
         !! (the simplest implemented branch).  Wright (1997)
         !! rational EOS is also shipped — opt in by setting
         !! `variant = EOS_VARIANT_WRIGHT_97`.  TEOS-10 reserved.
      real(wp) :: rho0 = 1035.0_wp
         !! Boussinesq reference density (kg/m^3).
      real(wp) :: T_ref = 10.0_wp
         !! Reference temperature for linear EOS (degC).
      real(wp) :: S_ref = 35.0_wp
         !! Reference salinity for linear EOS (PSU).
      real(wp) :: alpha_T = 1.7e-4_wp
         !! Thermal expansion coeff (linear EOS), kg/m^3 per degC.
         !! NOTE: this default is ~1000× smaller than the standard
         !! seawater value (~0.17 kg/m³/K).  The small value
         !! suppresses baroclinic feedback from PPM round-off in
         !! tests that don't care about realistic density gradients.
         !! Density-driven tests (e.g. `lock_exchange_2layer`) MUST
         !! override to the realistic value.
      real(wp) :: beta_S = 7.6e-4_wp
         !! Haline contraction coeff (linear EOS), kg/m^3 per PSU.
         !! See `alpha_T` for the rationale on this small default.
      real(wp) :: p_ref = 0.0_wp
         !! Reference pressure for the Wright EOS evaluation (Pa).
         !! Phase Tier-1 evaluates ρ at p_ref = 0 (surface) — gives
         !! the nonlinear-in-(T, S) ρ that cabbeling and densification
         !! tests need.  The FV-PGF integration path uses the in-situ
         !! hydrostatic pressure at each layer centre (its own column
         !! scratch on the PGF state, not an EOS field).
      integer :: ts_convention = TS_POT_PRAC
         !! Tracer T/S convention this EOS expects (TS_POT_PRAC /
         !! TS_CONS_ABS).  Identity for linear + Wright.  This type is
         !! a flat POD — NO allocatable component — so it copies into
         !! registers for free when passed by value into the device
         !! point routines (`!$acc routine seq`).
   contains
      procedure, non_overridable :: init => eos_init
      procedure, non_overridable :: destroy => eos_destroy
   end type eos_t

contains

   subroutine eos_init(this, grid)
      class(eos_t), intent(inout) :: this
      type(hgrid_t), intent(in) :: grid
      if (.false.) this%rho0 = real(grid%nx_total, wp)
      ! NOTE: the device-callable-variant gate is NOT run here anymore.
      ! `eos%variant` is still at its default at `init` time; the config
      ! knob (`&ocean_eos_nml eos=...`) is applied later, in
      ! `configure_ocean_pgf`, which calls `eos_validate` POST-config
      ! (so it sees the real variant) plus the FV_WRIGHT+Roquet fail-loud.
      this%is_init = .true.
   end subroutine eos_init

   subroutine eos_validate(eos)
      !! Host-side fail-loud gate over the device-supported variant
      !! set.  Device point routines (`!$acc routine seq`) cannot
      !! `error stop`, so membership in the device-callable set is
      !! guaranteed HERE at configure time; the device `else` branch
      !! is then unreachable-by-contract.
      use pic_logger, only: global_logger
      type(eos_t), intent(in) :: eos
      select case (eos%variant)
      case (EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, EOS_VARIANT_ROQUET_SPV)
         ! Device-callable — OK.
      case (EOS_VARIANT_TEOS10)
         call global_logger%error("ocean_eos: TEOS-10 variant is not "// &
                                  "yet device-callable")
         error stop "ocean_eos: TEOS-10 not yet device-callable"
      case default
         call global_logger%error("ocean_eos: unknown eos%variant")
         error stop "ocean_eos: unknown eos%variant"
      end select
   end subroutine eos_validate

   subroutine eos_destroy(this)
      class(eos_t), intent(inout) :: this
      this%is_init = .false.
   end subroutine eos_destroy

   pure function parse_eos_variant(name) result(code)
      !! Translate a `&ocean_eos_nml eos=...` string into an
      !! `EOS_VARIANT_*` code.  Unrecognised values fall back to
      !! `EOS_VARIANT_LINEAR` (the default = bit-identical to runs that
      !! omit the knob).  Membership in the device-callable set is
      !! enforced separately by `eos_validate` at configure time.
      character(len=*), intent(in) :: name
      integer :: code
      select case (trim(adjustl(name)))
      case ("linear", "LINEAR")
         code = EOS_VARIANT_LINEAR
      case ("wright", "WRIGHT", "wright_97", "wright97")
         code = EOS_VARIANT_WRIGHT_97
      case ("roquet_spv", "ROQUET_SPV", "roquet", "ROQUET")
         code = EOS_VARIANT_ROQUET_SPV
      case ("teos10", "TEOS10", "teos-10")
         code = EOS_VARIANT_TEOS10
      case default
         code = EOS_VARIANT_LINEAR
      end select
   end function parse_eos_variant

   subroutine eos_compute_arrays(eos, h_layer, hS_layer, hT_layer, &
                                 rho_layer, nx, ny, nz)
      !! Variant dispatch on bare 3D arrays — the regime-agnostic body of
      !! `ocean_eos_compute` (the ocean C-grid shim, hosted in
      !! `rki_ocean_eos_compute`).  That shim forwards the C-grid
      !! multilayer state's registry arrays here; the coastal opt-in path
      !! (`rki_ml_eos::ml_equation_of_state_shared`, `&tracer_nml
      !! eos_form`) forwards the A-grid multilayer arrays — the two
      !! regimes share one EOS through this single dispatch.  All impls
      !! evaluate at the fixed reference pressure `eos%p_ref` (default 0
      !! ⇒ surface/potential density) with the `H_VANISHED` vanishing-
      !! layer fallback to `eos%rho0`.
      type(eos_t), intent(in) :: eos
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: hS_layer(nx, ny, nz)
      real(wp), intent(in) :: hT_layer(nx, ny, nz)
      real(wp), intent(out) :: rho_layer(nx, ny, nz)

      select case (eos%variant)
      case (EOS_VARIANT_LINEAR)
         call eos_linear_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                              eos%rho0, eos%beta_S, eos%S_ref, &
                              eos%alpha_T, eos%T_ref, &
                              nx, ny, nz)
      case (EOS_VARIANT_WRIGHT_97)
         call eos_wright_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                              eos%rho0, eos%p_ref, &
                              nx, ny, nz)
      case (EOS_VARIANT_ROQUET_SPV)
         call eos_roquet_spv_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                                  eos%rho0, eos%p_ref, &
                                  nx, ny, nz)
      case default
         error stop "eos_compute_arrays: unknown eos%variant"
      end select
   end subroutine eos_compute_arrays

   subroutine eos_linear_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                                   rho_0, beta_S, S_ref, &
                                   alpha_T, T_ref, &
                                   nx, ny, nz)
      !! Linear two-tracer EOS, flat-impl form:
      !!
      !!   rho_k = rho_0 + beta_S * (S_k - S_ref) - alpha_T * (T_k - T_ref)
      !!
      !! where S_k = hS_k / h_k, T_k = hT_k / h_k.  `beta_S` and
      !! `alpha_T` are pre-multiplied sensitivities in kg/m³ per unit
      !! S / T — standard seawater values are 0.78 and 0.17.
      !!
      !! For vanishing layers (`h_layer <= H_VANISHED`) the cell falls
      !! back to rho_0 — same defensive branch as the coastal kernel; keeps
      !! the EOS finite under ZSTAR_FULL when bed-side layers can pinch out.
      !! The gate is `> H_VANISHED` (not `> 0`): during an active drain the
      !! PPM positivity limiter guarantees `h >= 0` but NOT `h >= H_VANISHED`,
      !! so a layer at e.g. `h = 1e-8` with `hS ≈ 35·1e-8` would pass a
      !! `> 0` gate and give `S = hS/h ≈ 5e2 PSU` ⇒ corrupted ρ ⇒ garbage
      !! PGF.  `> H_VANISHED` (the D4 vanished-layer role) returns rho_0 for
      !! any layer in `(0, H_VANISHED]`.  Bit-identical for any config whose
      !! layers all exceed H_VANISHED.  CAVEAT: ZSTAR_FULL floors vanishing
      !! bed layers to `zstar_h_min` (default 1.5e-4 == H_VANISHED exactly),
      !! so such a bed layer now takes the rho_0 fallback instead of the
      !! computed density — a dynamically negligible change on a 0.15 mm
      !! layer (PGF contribution ~1e-4 of a normal layer); the shipped
      !! anchors (ocean_analytical 8/8, dyn_split, baroclinic_longrun,
      !! double-gyre helpers) pass unchanged on both toolchains.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: hS_layer(nx, ny, nz)
      real(wp), intent(in) :: hT_layer(nx, ny, nz)
      real(wp), intent(out) :: rho_layer(nx, ny, nz)
      real(wp), intent(in) :: rho_0, beta_S, S_ref, alpha_T, T_ref

      integer :: i, j, k
      real(wp) :: inv_h, S_k, T_k

      !$omp parallel do collapse(3) private(inv_h, S_k, T_k)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (h_layer(i, j, k) > H_VANISHED) then
            inv_h = 1.0_wp/h_layer(i, j, k)
            S_k = hS_layer(i, j, k)*inv_h
            T_k = hT_layer(i, j, k)*inv_h
         else
            S_k = S_ref
            T_k = T_ref
         end if
         rho_layer(i, j, k) = rho_0 + beta_S*(S_k - S_ref) - alpha_T*(T_k - T_ref)
      end do
      end do
      end do
      !$omp end parallel do
   end subroutine eos_linear_impl

   subroutine eos_wright_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                                   rho_0, p_ref, nx, ny, nz)
      !! Wright (1997) rational EOS evaluated at a single reference
      !! pressure `p_ref`.  Same outer-shim signature as
      !! `eos_linear_impl` — bare 3D arrays, vanishing-layer
      !! fallback to `rho_0`.
      !!
      !!   α_0(T, S) = a0 + a1*T + a2*S
      !!   p_0(T, S) = b0 + b1*T + b2*T^2 + b3*T^3 + b4*S + b5*S*T
      !!   λ  (T, S) = c0 + c1*T + c2*T^2 + c3*T^3 + c4*S + c5*S*T
      !!
      !!   ρ = (P + p_0) / (λ + α_0 * (P + p_0))
      !!
      !! Reference behaviour (T=10, S=35, P=0): ρ ≈ 1027.3 kg/m^3,
      !! within 0.01 kg/m^3 of the surface seawater density used
      !! across the MOM6 test suite.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: hS_layer(nx, ny, nz)
      real(wp), intent(in) :: hT_layer(nx, ny, nz)
      real(wp), intent(out) :: rho_layer(nx, ny, nz)
      real(wp), intent(in) :: rho_0, p_ref

      integer :: i, j, k
      real(wp) :: inv_h, S_k, T_k, T_sq, T_cu
      real(wp) :: alpha_0, p_0, lambda, p_plus_p0, denom

      !$omp parallel do collapse(3) private(inv_h, S_k, T_k, T_sq, T_cu, alpha_0, p_0, lambda, p_plus_p0, denom)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (h_layer(i, j, k) > H_VANISHED) then
            inv_h = 1.0_wp/h_layer(i, j, k)
            S_k = hS_layer(i, j, k)*inv_h
            T_k = hT_layer(i, j, k)*inv_h
            T_sq = T_k*T_k
            T_cu = T_sq*T_k

            alpha_0 = WRIGHT_A0 + WRIGHT_A1*T_k + WRIGHT_A2*S_k
            p_0 = WRIGHT_B0 + WRIGHT_B1*T_k + WRIGHT_B2*T_sq + WRIGHT_B3*T_cu + &
                  WRIGHT_B4*S_k + WRIGHT_B5*S_k*T_k
            lambda = WRIGHT_C0 + WRIGHT_C1*T_k + WRIGHT_C2*T_sq + WRIGHT_C3*T_cu + &
                     WRIGHT_C4*S_k + WRIGHT_C5*S_k*T_k

            p_plus_p0 = p_ref + p_0
            denom = lambda + alpha_0*p_plus_p0
            rho_layer(i, j, k) = p_plus_p0/denom
         else
            rho_layer(i, j, k) = rho_0
         end if
      end do
      end do
      end do
      !$omp end parallel do
   end subroutine eos_wright_impl

   subroutine eos_roquet_spv_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                                       rho_0, p_ref, nx, ny, nz)
      !! Roquet et al. (2015) SpV EOS evaluated at a single reference
      !! pressure `p_ref`.  Same outer-shim signature as
      !! `eos_linear_impl`/`eos_wright_impl` — bare 3D arrays, model
      !! (PT, SP) tracers, vanishing-layer fallback to `rho_0`.
      !!
      !! Density = 1 / SV(CT(SR,PT), SR, p_ref) with SR = SP·(35.16504/35)
      !! and CT = ct_from_pt(SR, PT) — the conversions live inside the
      !! fused `roquet_spv_point` point routine (the only sqrt path).
      !!
      !! FV_LITE / FV_MOM6 PGF read `rho_layer` generically, so they pick
      !! up Roquet automatically.  (FV_WRIGHT re-evaluates Wright in its
      !! Picard sweep → unsupported with Roquet; gated fail-loud at
      !! configure in `configure_ocean_pgf`.)
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: hS_layer(nx, ny, nz)
      real(wp), intent(in) :: hT_layer(nx, ny, nz)
      real(wp), intent(out) :: rho_layer(nx, ny, nz)
      real(wp), intent(in) :: rho_0, p_ref

      integer :: i, j, k
      real(wp) :: inv_h, S_k, T_k, sv_k, d_dum1, d_dum2

      !$omp parallel do collapse(3) private(inv_h, S_k, T_k, sv_k, d_dum1, d_dum2)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (h_layer(i, j, k) > H_VANISHED) then
            inv_h = 1.0_wp/h_layer(i, j, k)
            S_k = hS_layer(i, j, k)*inv_h
            T_k = hT_layer(i, j, k)*inv_h
            call roquet_spv_point(T_k, S_k, p_ref, sv_k, d_dum1, d_dum2)
            rho_layer(i, j, k) = 1.0_wp/sv_k
         else
            rho_layer(i, j, k) = rho_0
         end if
      end do
      end do
      end do
      !$omp end parallel do
   end subroutine eos_roquet_spv_impl

   subroutine eos_wright_pgf_column_sweep_impl(h_layer, hS_layer, hT_layer, &
                                                    rho_layer_seed, p_edge_out, &
                                                    rho_insitu_out, gravity, rho_0, &
                                                    nx, ny, nz)
      !! FV-Wright PGF column sweep.  Top-down per-column traversal
      !! that simultaneously produces the hydrostatic pressure stack
      !! `p_edge_out` and the in-situ density `rho_insitu_out` at each
      !! layer centre.
      !!
      !! For each layer k from the surface (k=nz) down to the bed
      !! (k=1) we use a single Picard step:
      !!
      !!   1. Seed the half-layer pressure with the surface-evaluated
      !!      density `rho_layer_seed` (= ms%rho_layer, the EOS output
      !!      at p_ref=0):
      !!         p_centre_seed = p_above + 0.5 * g * rho_seed * h
      !!   2. Re-evaluate Wright at the seed pressure:
      !!         rho_insitu(k) = ρ(T_k, S_k, p_centre_seed)
      !!   3. Accumulate the bottom edge of layer k using the in-situ
      !!      density:
      !!         p_edge(k) = p_above + g * rho_insitu(k) * h
      !!
      !! This is one Picard iteration of the implicit
      !!   p_centre(k) = p_above + 0.5*g*ρ(T, S, p_centre(k))*h.
      !! For ocean conditions (Δρ along path << ρ_0) one iteration is
      !! within ~1e-5 of the converged value.
      !!
      !! `rho_layer_seed` must be `ms%rho_layer` from
      !! `ocean_eos_compute` with `EOS_VARIANT_WRIGHT_97` —
      !! the seed is the *full* nonlinear ρ at p_ref=0, not a Boussinesq
      !! constant; this avoids a second Picard iteration in 99% of
      !! cases.
      !!
      !! Vanishing-layer fallback: if `h_layer(k) <= 0` the Wright eval
      !! is skipped and `rho_insitu(k) = rho_0` — matches the existing
      !! `eos_wright_impl` defensive branch.
      !!
      !! Loop order: outer `do concurrent (j, i)` for GPU parallelism;
      !! inner serial k loop for the column recurrence (same shape as
      !! the existing PGF Pass 1 + vdiff column kernels).
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(in) :: h_layer(nx, ny, nz)
      real(wp), intent(in) :: hS_layer(nx, ny, nz)
      real(wp), intent(in) :: hT_layer(nx, ny, nz)
      real(wp), intent(in) :: rho_layer_seed(nx, ny, nz)
      real(wp), intent(out) :: p_edge_out(nx, ny, nz + 1)
      real(wp), intent(out) :: rho_insitu_out(nx, ny, nz)
      real(wp), intent(in) :: gravity, rho_0

      integer :: i, j, k
      real(wp) :: p_above, p_centre_seed, inv_h, S_k, T_k, T_sq, T_cu
      real(wp) :: alpha_0, p_0, lambda, p_plus_p0, denom, rho_k

      !$omp parallel do collapse(2) private(k, p_above, p_centre_seed, inv_h, S_k, T_k, T_sq, T_cu, alpha_0, p_0, lambda, p_plus_p0, denom, rho_k)
      do j = 1, ny
      do i = 1, nx
         p_edge_out(i, j, nz + 1) = 0.0_wp
         p_above = 0.0_wp
         do k = nz, 1, -1
            if (h_layer(i, j, k) > H_VANISHED) then
               inv_h = 1.0_wp/h_layer(i, j, k)
               S_k = hS_layer(i, j, k)*inv_h
               T_k = hT_layer(i, j, k)*inv_h
               T_sq = T_k*T_k
               T_cu = T_sq*T_k

               p_centre_seed = p_above + 0.5_wp*gravity*rho_layer_seed(i, j, k)*h_layer(i, j, k)

               alpha_0 = WRIGHT_A0 + WRIGHT_A1*T_k + WRIGHT_A2*S_k
               p_0 = WRIGHT_B0 + WRIGHT_B1*T_k + WRIGHT_B2*T_sq + WRIGHT_B3*T_cu + &
                     WRIGHT_B4*S_k + WRIGHT_B5*S_k*T_k
               lambda = WRIGHT_C0 + WRIGHT_C1*T_k + WRIGHT_C2*T_sq + WRIGHT_C3*T_cu + &
                        WRIGHT_C4*S_k + WRIGHT_C5*S_k*T_k

               p_plus_p0 = p_centre_seed + p_0
               denom = lambda + alpha_0*p_plus_p0
               rho_k = p_plus_p0/denom
            else
               rho_k = rho_0
            end if
            rho_insitu_out(i, j, k) = rho_k
            p_above = p_above + gravity*rho_k*h_layer(i, j, k)
            p_edge_out(i, j, k) = p_above
         end do
      end do
      end do
      !$omp end parallel do
   end subroutine eos_wright_pgf_column_sweep_impl

   pure subroutine roquet_spv_point(T_pt, S_sp, p, sv, dsv_dt_model, dsv_ds_model)
      !! Fused Roquet et al. (2015) SpV evaluation at a point, in MODEL
      !! variables (potential temperature `T_pt` degC, practical salinity
      !! `S_sp` PSU, pressure `p` Pa).  Returns specific volume `sv`
      !! (m^3/kg) and the analytic sensitivities w.r.t. the MODEL
      !! variables (`dsv_dt_model` = dSV/dPT, `dsv_ds_model` = dSV/dSP)
      !! so the variant-agnostic consumers (which work in PT, SP) get the
      !! correct chain-ruled derivatives.
      !!
      !! Convention (pinned 2026-06-16 — see EOS_VARIANT_ROQUET_SPV):
      !!   SR = S_sp * (35.16504/35)         [Reference Salinity, g/kg]
      !!   CT = ct_from_pt(SR, T_pt)         [Conservative Temperature]
      !! The polynomial is fit in (CT, SA, p); we feed SR for SA (bounded
      !! anomaly deviation) and CT via the local conversion poly.
      !!
      !! Chain rule back to the model variables.  CT = ct_from_pt(SR, PT)
      !! depends on BOTH PT and SR (= SP·factor), so the SP derivative
      !! carries TWO routes into SV — the direct SA route AND the
      !! CT-via-SR route:
      !!   dSV/dPT = (dSV/dCT) · (dCT/dPT)
      !!   dSV/dSP = [ (dSV/dSA) + (dSV/dCT)·(dCT/dSR) ] · (35.16504/35)
      !! The CT-via-SR term is ~0.5 % of dSV/dSP (dCT/dSR ≈ −0.02 degC per
      !! g/kg); dropping it (the simplified spec formula) leaves a real
      !! ~5e-3 relative error against the true total derivative, so we
      !! keep the full chain — this is what the FD regression locks.
      !!
      !! dCT/dPT and dCT/dSR are analytic derivatives of the ct_from_pt
      !! poly: with yy = 0.025·PT, x2 = sfac·SR, xx = sqrt(x2), each
      !! yy-coefficient c_k(xx) = A + B·xx² + C·xx³ + D·xx⁴ + E·xx⁵, so
      !!   dCT/dPT = 0.025·(dhh/dyy)/cp0
      !!   dCT/dSR = (dhh/dxx)·(sfac/(2·xx))/cp0.
      !!
      !! Transcribed from the verified prototype roquet_spv_eos.py.  One
      !! sqrt for zs (shared by SV and both derivatives) + one sqrt for
      !! the ct_from_pt poly normalisation.
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      real(wp), intent(in) :: T_pt, S_sp, p
      real(wp), intent(out) :: sv, dsv_dt_model, dsv_ds_model

      real(wp) :: SR, CT, dct_dpt, dct_dsr
      real(wp) :: zt, zs, zp
      real(wp) :: sv_ts0, sv_ts1, sv_ts2, sv_ts3, sv_0s0, sv_00p
      real(wp) :: dvdzt0, dvdzt1, dvdzt2, dvdzt3, dsv_dct
      real(wp) :: dvdzs0, dvdzs1, dvdzs2, dvdzs3, dsv_dsa
      real(wp) :: x2, xx, yy, hh, dh_dy, dh_dx
      real(wp) :: c0, c1, c2, c3, c4, c5, c6, c7
      real(wp) :: d0, d1, d2, d3, d4, d5, d6, d7

      SR = S_sp*ROQ_SR_FACTOR

      ! CT = ct_from_pt(SR, PT): the 7-term gsw surface poly is degree 7 in
      ! yy = 0.025*PT.  Each yy-coefficient c_k is a polynomial in xx =
      ! sqrt(sfac*SR): c_k = A + B*xx^2 + C*xx^3 + D*xx^4 + E*xx^5 (x2 = xx^2).
      ! Building c_k and dc_k/dxx (= d_k) gives exact analytic dCT/dPT and
      ! dCT/dSR by Horner — no risk of mis-differentiating the published
      ! nesting.
      ! Floor x2 to a tiny positive so the `dct_dsr = dh_dx/(2*xx)` below
      ! never hits 0/0 at SP=0 (fresh water): dh_dx is itself proportional
      ! to xx (lowest term 2B*xx), so dh_dx/xx has a finite limit — the
      ! floor reproduces it (xx tiny-but-nonzero) instead of NaN.  The
      ! floor only bites at SP < ~1e-9 PSU, so it is bit-identical for any
      ! real-ocean salinity.
      x2 = max(ROQ_CT_SFAC*SR, 1.0e-20_wp)
      xx = sqrt(x2)
      yy = T_pt*0.025_wp
      c0 = 61.01362420681071_wp &
           + x2*(268.5520265845071_wp &
                 + xx*(937.2099110620707_wp &
                       + xx*(-1687.914374187449_wp + xx*246.9598888781377_wp)))
      c1 = 168776.46138048015_wp &
           + x2*(-12019.028203559312_wp &
                 + xx*(588.1802812170108_wp &
                       + xx*(936.3206544460336_wp + xx*123.59576582457964_wp)))
      c2 = -2735.2785605119625_wp &
           + x2*(3734.858026725145_wp &
                 + xx*(248.39476522971285_wp &
                       + xx*(-942.7827304544439_wp + xx*(-48.5891069025409_wp))))
      c3 = 2574.2164453821433_wp &
           + x2*(-2046.7671145057618_wp &
                 + xx*(-3.871557904936333_wp + xx*369.4389437509002_wp))
      c4 = -1536.6644434977543_wp &
           + x2*(465.28655623126450_wp + xx*(-2.6268019854268356_wp + xx*(-33.83664947895248_wp)))
      c5 = 545.7340497931629_wp &
           + x2*(-0.6370820302831379_wp + xx*(-9.987880382780322_wp))
      c6 = -50.91091728474331_wp + x2*(-10.650848542359153_wp)
      c7 = -18.30489878927802_wp
      ! dc_k/dxx = 2B*xx + 3C*xx^2 + 4D*xx^3 + 5E*xx^4 (A and the xx^0/xx^1
      ! terms vanish; B,C,D,E read off the c_k expansions above).
      d0 = xx*(2.0_wp*268.5520265845071_wp &
               + xx*(3.0_wp*937.2099110620707_wp &
                     + xx*(4.0_wp*(-1687.914374187449_wp) + xx*5.0_wp*246.9598888781377_wp)))
      d1 = xx*(2.0_wp*(-12019.028203559312_wp) &
               + xx*(3.0_wp*588.1802812170108_wp &
                     + xx*(4.0_wp*936.3206544460336_wp + xx*5.0_wp*123.59576582457964_wp)))
      d2 = xx*(2.0_wp*3734.858026725145_wp &
               + xx*(3.0_wp*248.39476522971285_wp &
                     + xx*(4.0_wp*(-942.7827304544439_wp) + xx*5.0_wp*(-48.5891069025409_wp))))
      d3 = xx*(2.0_wp*(-2046.7671145057618_wp) &
               + xx*(3.0_wp*(-3.871557904936333_wp) + xx*4.0_wp*369.4389437509002_wp))
      d4 = xx*(2.0_wp*465.28655623126450_wp &
               + xx*(3.0_wp*(-2.6268019854268356_wp) + xx*4.0_wp*(-33.83664947895248_wp)))
      d5 = xx*(2.0_wp*(-0.6370820302831379_wp) + xx*3.0_wp*(-9.987880382780322_wp))
      d6 = xx*2.0_wp*(-10.650848542359153_wp)
      d7 = 0.0_wp
      hh = c0 + yy*(c1 + yy*(c2 + yy*(c3 + yy*(c4 + yy*(c5 + yy*(c6 + yy*c7))))))
      dh_dy = c1 + yy*(2.0_wp*c2 + yy*(3.0_wp*c3 + yy*(4.0_wp*c4 &
                                                       + yy*(5.0_wp*c5 + yy*(6.0_wp*c6 + yy*7.0_wp*c7)))))
      dh_dx = d0 + yy*(d1 + yy*(d2 + yy*(d3 + yy*(d4 + yy*(d5 + yy*(d6 + yy*d7))))))
      CT = hh/ROQ_CP0
      dct_dpt = 0.025_wp*dh_dy/ROQ_CP0
      dct_dsr = dh_dx*(ROQ_CT_SFAC/(2.0_wp*xx))/ROQ_CP0

      zt = CT
      zs = sqrt(abs(S_sp*ROQ_SR_FACTOR + ROQ_RDELTAS)*ROQ_R1_S0)
      zp = p

      ! --- specific volume SV(zs, zt, zp) ---
      sv_ts3 = SPV003 + (zs*SPV103 + zt*SPV013)
      sv_ts2 = SPV002 + (zs*(SPV102 + zs*SPV202) &
                         + zt*(SPV012 + (zs*SPV112 + zt*SPV022)))
      sv_ts1 = SPV001 + (zs*(SPV101 + zs*(SPV201 + zs*(SPV301 + zs*SPV401))) &
                         + zt*(SPV011 + (zs*(SPV111 + zs*(SPV211 + zs*SPV311)) &
                                         + zt*(SPV021 + (zs*(SPV121 + zs*SPV221) &
                                                         + zt*(SPV031 + (zs*SPV131 + zt*SPV041)))))))
      sv_ts0 = zt*(SPV010 &
                   + (zs*(SPV110 + zs*(SPV210 + zs*(SPV310 + zs*(SPV410 + zs*SPV510)))) &
                      + zt*(SPV020 + (zs*(SPV120 + zs*(SPV220 + zs*(SPV320 + zs*SPV420))) &
                                      + zt*(SPV030 + (zs*(SPV130 + zs*(SPV230 + zs*SPV330)) &
                                                      + zt*(SPV040 + (zs*(SPV140 + zs*SPV240) &
                                                                      + zt*(SPV050 + (zs*SPV150 + zt*SPV060))))))))))
      sv_0s0 = SPV000 + zs*(SPV100 + zs*(SPV200 + zs*(SPV300 + zs*(SPV400 &
                                                                   + zs*(SPV500 + zs*SPV600)))))
      sv_00p = zp*(ROQ_V00 + zp*(ROQ_V01 + zp*(ROQ_V02 + zp*(ROQ_V03 &
                                                             + zp*(ROQ_V04 + zp*ROQ_V05)))))
      sv = ((sv_ts0 + sv_0s0) + zp*(sv_ts1 + zp*(sv_ts2 + zp*sv_ts3))) + sv_00p

      ! --- dSV/dCT ---
      dvdzt3 = ALP003
      dvdzt2 = ALP002 + (zs*ALP102 + zt*ALP012)
      dvdzt1 = ALP001 + (zs*(ALP101 + zs*(ALP201 + zs*ALP301)) &
                         + zt*(ALP011 + (zs*(ALP111 + zs*ALP211) &
                                         + zt*(ALP021 + (zs*ALP121 + zt*ALP031)))))
      dvdzt0 = ALP000 + (zs*(ALP100 + zs*(ALP200 + zs*(ALP300 + zs*(ALP400 + zs*ALP500)))) &
                         + zt*(ALP010 + (zs*(ALP110 + zs*(ALP210 + zs*(ALP310 + zs*ALP410))) &
                                         + zt*(ALP020 + (zs*(ALP120 + zs*(ALP220 + zs*ALP320)) &
                                                         + zt*(ALP030 + (zt*(ALP040 + (zs*ALP140 + zt*ALP050)) &
                                                                         + zs*(ALP130 + zs*ALP230))))))))
      dsv_dct = dvdzt0 + zp*(dvdzt1 + zp*(dvdzt2 + zp*dvdzt3))

      ! --- dSV/dSA (per-coef 0.5*r1_S0 folded into BET; residual /zs here) ---
      dvdzs3 = BET003
      dvdzs2 = BET002 + (zs*BET102 + zt*BET012)
      dvdzs1 = BET001 + (zs*(BET101 + zs*(BET201 + zs*BET301)) &
                         + zt*(BET011 + (zs*(BET111 + zs*BET211) &
                                         + zt*(BET021 + (zs*BET121 + zt*BET031)))))
      dvdzs0 = BET000 + (zs*(BET100 + zs*(BET200 + zs*(BET300 + zs*(BET400 + zs*BET500)))) &
                         + zt*(BET010 + (zs*(BET110 + zs*(BET210 + zs*(BET310 + zs*BET410))) &
                                         + zt*(BET020 + (zs*(BET120 + zs*(BET220 + zs*BET320)) &
                                                         + zt*(BET030 + (zt*(BET040 + (zs*BET140 + zt*BET050)) &
                                                                         + zs*(BET130 + zs*BET230))))))))
      dsv_dsa = (dvdzs0 + zp*(dvdzs1 + zp*(dvdzs2 + zp*dvdzs3)))/zs

      ! Chain rule to the model variables (PT, SP).  The SP route includes
      ! the CT-via-SR coupling (CT depends on SR = SP·factor).
      dsv_dt_model = dsv_dct*dct_dpt
      dsv_ds_model = (dsv_dsa + dsv_dct*dct_dsr)*ROQ_SR_FACTOR
   end subroutine roquet_spv_point

   pure subroutine eos_specvol_derivs(eos, T, S, p, dsv_dt, dsv_ds)
      !! Analytic specific-volume sensitivities dSV/dT and dSV/dS
      !! (SV = 1/rho) at a point.  Needed by the EPBL energy
      !! bookkeeping (pressure-weighted PE-per-unit-tracer-change
      !! weights) and kappa-shear buoyancy — dSV/dX = -(1/rho^2)
      !! d(rho)/dX.
      !!
      !! Takes the shared `eos_t` handle BY VALUE (flat POD,
      !! no allocatable) so the device copy is register-resident.
      !! The `select case (eos%variant)` body is warp-uniform (one
      !! variant per run) — ~free.
      !!
      !! Linear variant: rho = rho0 + beta_S (S - S_ref) - alpha_T
      !! (T - T_ref) gives constant dSV/dT = +alpha_T/rho0^2 and
      !! dSV/dS = -beta_S/rho0^2 (evaluated at the reference density,
      !! consistent with the Boussinesq weights that consume them).
      !!
      !! Wright (1997) variant: SV = alpha_0(T,S) + lambda(T,S)/P with
      !! P = p + p_0(T,S), differentiable in closed form from the
      !! Table A1 polynomials:
      !!   dSV/dX = d(alpha_0)/dX + d(lambda)/dX / P
      !!            - lambda * d(p_0)/dX / P^2.
      !! The `else` is unreachable-by-contract: `eos_validate`
      !! guarantees `eos%variant` is in the device-callable set at
      !! configure time (device code cannot `error stop`).
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle (variant + scalar coeffs), by value.
      real(wp), intent(in) :: T, S
         !! In-situ temperature (degC) and salinity (PSU).
      real(wp), intent(in) :: p
         !! Pressure (Pa), hydrostatic surface-relative.
      real(wp), intent(out) :: dsv_dt
         !! dSV/dT (m^3/kg/degC); > 0 for warm-expands water.
      real(wp), intent(out) :: dsv_ds
         !! dSV/dS (m^3/kg/PSU); < 0 (salt contracts).

      real(wp) :: T_sq, p_0, lambda, big_p, inv_p, inv_p2
      real(wp) :: dp0_dt, dlam_dt, dp0_ds, dlam_ds
      real(wp) :: sv_roq

      if (eos%variant == EOS_VARIANT_ROQUET_SPV) then
         call roquet_spv_point(T, S, p, sv_roq, dsv_dt, dsv_ds)
      else if (eos%variant == EOS_VARIANT_WRIGHT_97) then
         T_sq = T*T
         p_0 = WRIGHT_B0 + WRIGHT_B1*T + WRIGHT_B2*T_sq + WRIGHT_B3*T_sq*T + &
               WRIGHT_B4*S + WRIGHT_B5*S*T
         lambda = WRIGHT_C0 + WRIGHT_C1*T + WRIGHT_C2*T_sq + WRIGHT_C3*T_sq*T + &
                  WRIGHT_C4*S + WRIGHT_C5*S*T
         dp0_dt = WRIGHT_B1 + 2.0_wp*WRIGHT_B2*T + 3.0_wp*WRIGHT_B3*T_sq + &
                  WRIGHT_B5*S
         dlam_dt = WRIGHT_C1 + 2.0_wp*WRIGHT_C2*T + 3.0_wp*WRIGHT_C3*T_sq + &
                   WRIGHT_C5*S
         dp0_ds = WRIGHT_B4 + WRIGHT_B5*T
         dlam_ds = WRIGHT_C4 + WRIGHT_C5*T
         big_p = p + p_0
         inv_p = 1.0_wp/big_p
         inv_p2 = inv_p*inv_p
         dsv_dt = WRIGHT_A1 + dlam_dt*inv_p - lambda*dp0_dt*inv_p2
         dsv_ds = WRIGHT_A2 + dlam_ds*inv_p - lambda*dp0_ds*inv_p2
      else
         dsv_dt = eos%alpha_T/(eos%rho0*eos%rho0)
         dsv_ds = -eos%beta_S/(eos%rho0*eos%rho0)
      end if
   end subroutine eos_specvol_derivs

   pure function eos_density_point(eos, T, S, p) result(rho)
      !! Scalar density evaluation at a point — the same formulas the
      !! 3D `eos_*_impl` kernels apply, exposed for finite-difference
      !! verification of `eos_specvol_derivs` and for host-side
      !! diagnostics.  Takes the shared `eos_t` handle by value.
      !! The `else` is unreachable-by-contract (see
      !! `eos_validate`).
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      type(eos_t), intent(in) :: eos
      real(wp), intent(in) :: T, S, p
      real(wp) :: rho

      real(wp) :: T_sq, alpha_0, p_0, lambda, p_plus_p0
      real(wp) :: sv_roq, d_dum1, d_dum2

      if (eos%variant == EOS_VARIANT_ROQUET_SPV) then
         call roquet_spv_point(T, S, p, sv_roq, d_dum1, d_dum2)
         rho = 1.0_wp/sv_roq
      else if (eos%variant == EOS_VARIANT_WRIGHT_97) then
         T_sq = T*T
         alpha_0 = WRIGHT_A0 + WRIGHT_A1*T + WRIGHT_A2*S
         p_0 = WRIGHT_B0 + WRIGHT_B1*T + WRIGHT_B2*T_sq + WRIGHT_B3*T_sq*T + &
               WRIGHT_B4*S + WRIGHT_B5*S*T
         lambda = WRIGHT_C0 + WRIGHT_C1*T + WRIGHT_C2*T_sq + WRIGHT_C3*T_sq*T + &
                  WRIGHT_C4*S + WRIGHT_C5*S*T
         p_plus_p0 = p + p_0
         rho = p_plus_p0/(lambda + alpha_0*p_plus_p0)
      else
         rho = eos%rho0 + eos%beta_S*(S - eos%S_ref) - eos%alpha_T*(T - eos%T_ref)
      end if
   end function eos_density_point

   pure elemental function eos_freezing_point(eos, S, p) result(T_f)
      !! Seawater freezing point T_f (degC) at salinity `S` and pressure
      !! `p` — the ocean-side prerequisite for the sea-ice port
      !! (PLAN_SEA_ICE.md, PR 1).  Same point-function style as
      !! `eos_density_point` (flat-POD `eos_t` by value, device-
      !! callable), plus `elemental` so callers can evaluate whole
      !! salinity arrays in one reference.
      !!
      !! v1 is the SIS2/MOM6 LINEAR liquidus for EVERY `eos%variant`:
      !!
      !!   T_f = TFR_S_COEFF·S + TFR_P_COEFF·p
      !!
      !! (μ = −0.054 °C per g/kg, MOM6 `TFREEZE_FORM = "LINEAR"` — MOM6
      !! keeps the linear form as its default even under the Wright /
      !! TEOS-10 density branches, so applying it across variants is
      !! parity, not a shortcut.)  A TEOS-10 `t_freezing(SA, p)`
      !! polynomial slots in later as an `eos%variant ==
      !! EOS_VARIANT_ROQUET_SPV / EOS_VARIANT_TEOS10` branch, mirroring
      !! the dispatch in `eos_density_point`.
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle (variant tag), by value.  Read for the
         !! future per-variant dispatch; every variant takes the linear
         !! branch in v1.
      real(wp), intent(in) :: S
         !! Salinity (PSU / g/kg).
      real(wp), intent(in) :: p
         !! Pressure (Pa), hydrostatic surface-relative (0 at the
         !! surface — the frazil kernel's use case).
      real(wp) :: T_f

      ! Variant dispatch seam: v1 = the same linear liquidus for all
      ! variants (see docstring), so `eos%variant` is read only by the
      ! future TEOS-10 branch — which slots in as
      !   if (eos%variant == EOS_VARIANT_ROQUET_SPV) then ... else
      ! without touching this default expression.
      T_f = TFR_S_COEFF*S + TFR_P_COEFF*p
   end function eos_freezing_point

end module rki_eos
