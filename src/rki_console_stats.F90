!! Regime-agnostic MOM6-style console statistics formatter.
module rki_console_stats
   !! The shared print layout for the periodic conservation + stability
   !! console line, used by BOTH the ocean and coastal drivers so the two
   !! regimes report identically.  Each regime computes its own totals
   !! (area-weighted device reductions over its own state), then hands the
   !! scalars here; this module owns the reference-snapshot latching, the
   !! relative-drift arithmetic, the exact line format, and the NaN / CFL
   !! panic guards.
   !!
   !! Layout (mirrors MOM6's "MOM Day N:" status line):
   !!
   !!   [stats] Day D  step N  En E  MaxCFL C  Mass M  [Salt S  Temp T]
   !!       Mass : <total>  Error <drift-from-t0>
   !!       Salt : <total>  Error <drift>     (thermodynamics on)
   !!       Heat : <total>  Error <drift>     (thermodynamics on)
   !!       Age  : <mean>  days               (ideal-age tracer on)
   !!       En   : <total KE>  Growth <ratio vs t0 KE>x
   !!
   !! Initial values are captured on the first call (`is_initialised`);
   !! every later call reports drift relative to them.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use rki_constants, only: wp
   use rki_efp, only: efp_t, efp_real_diff
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: console_stats_t
   public :: conservation_budget_t
   public :: console_stats_report

   real(wp), parameter :: CFL_PANIC_THRESHOLD = 0.9_wp
      !! MaxCFL above this triggers an early-stop panic. Just below the
      !! stability limit (1.0) so we abort before NaN on doomed runs.

   type :: console_stats_t
      !! Holds initial-snapshot values for relative-drift reporting.
      !! `is_initialised` flips on the first `report` call.
      logical :: is_initialised = .false.
      real(wp) :: mass0 = 0.0_wp
      real(wp) :: salt0 = 0.0_wp
      real(wp) :: heat0 = 0.0_wp
      real(wp) :: ke0 = 0.0_wp
      ! Cumulative budget terms at the reference (latched with the totals on
      ! the first report), so the residual is measured over the SAME window
      ! as the drift — the accumulators run from t=0 but the reference report
      ! may fire after a warmup step.
      real(wp) :: mass_out0 = 0.0_wp, salt_out0 = 0.0_wp, heat_out0 = 0.0_wp
      real(wp) :: mass_src0 = 0.0_wp, salt_src0 = 0.0_wp, heat_src0 = 0.0_wp
      type(efp_t) :: mass0_efp, salt0_efp, heat0_efp
         !! PR-32: extended-fixed-point (EFP) reference snapshots, latched
         !! alongside `mass0`/`salt0`/`heat0` on the first report when the
         !! caller supplies `mass_efp`/`salt_efp`/`heat_efp` (the ocean
         !! `reproducing_sums = .true.` path).  Their sole purpose is
         !! `efp_real_diff(current_efp, mass0_efp)` -- a difference formed
         !! in FIXED POINT, not `current - mass0` in double.  This is the
         !! SS2.2 fix: latching the reference as a `real(wp)` alone (as
         !! `mass0` does) already quantises it to ~1 ulp of a ~1e21 total
         !! (~1.3e5 kg) before any subtraction happens, so an exact SUM
         !! feeding an unchanged double-difference would not move the
         !! drift floor at all.
      logical :: exact_sums = .false.
         !! `.true.` once `mass_efp`/`salt_efp`/`heat_efp` have been
         !! latched (the ocean EFP path is active); gates whether
         !! `emit_drift_line`'s caller may rely on `mass0_efp` etc. being
         !! meaningful.  Coastal callers never set this (stays `.false.`
         !! for the lifetime of a coastal run's `console_stats_t`).
      logical :: report_thermodynamics = .true.
         !! `.false.` suppresses the Salt / Temp columns + detail blocks
         !! (a 2D-barotropic / adiabatic run has no meaningful S/T, so it
         !! doesn't print misleading lines).
      logical :: panic_on_nan = .true.
         !! Early-stop abort when any stat is NaN (a NaN is always fatal).
      logical :: panic_on_cfl = .true.
         !! Early-stop abort when MaxCFL exceeds the panic threshold.  Right
         !! for the ocean (advective CFL IS its stability criterion), but
         !! the coastal path sets this `.false.`: coastal stability is
         !! gravity-wave-limited, so a high *advective* CFL is normal (fast
         !! shallow flows) and must not trip a spurious abort.
   end type console_stats_t

   type :: conservation_budget_t
      !! Cumulative (time-integrated, since t=0) budget terms that CLOSE the
      !! domain integral of each extensive quantity, so the reported `Error`
      !! is a true numerical-leak residual even with open boundaries + surface
      !! forcing.  For quantity Q the residual is
      !!
      !!   leak_Q = (Q_total − Q0) + out_Q − src_Q
      !!
      !! (change in the domain integral, plus what left through the open
      !! boundaries, minus what the surface added), which should stay
      !! ~round-off for a conservative scheme.  Same units as the matching
      !! total.  A quantity whose `*_active` flag is `.false.` falls back to
      !! raw drift `(Q_total − Q0)/Q0` — byte-identical to the pre-budget
      !! behaviour (used for the not-yet-instrumented quantities/paths).
      real(wp) :: mass_out = 0.0_wp
      real(wp) :: salt_out = 0.0_wp
      real(wp) :: heat_out = 0.0_wp
         !! Net OUTFLUX through open boundaries (positive = left the domain).
      real(wp) :: mass_src = 0.0_wp
      real(wp) :: salt_src = 0.0_wp
      real(wp) :: heat_src = 0.0_wp
         !! Net SURFACE source (positive = added to the domain: precip,
         !! surface heat/salt flux).
      logical :: mass_active = .false.
      logical :: salt_active = .false.
      logical :: heat_active = .false.
         !! Per-quantity: `.true.` ⇒ that quantity's `Error` is the closed
         !! budget residual + shows `out`/`src`; `.false.` ⇒ plain drift.
         !! Per-quantity so a phased rollout (mass first) doesn't print
         !! misleading `out 0` for a quantity whose flux isn't tracked yet.
   end type conservation_budget_t

contains

   subroutine console_stats_report(this, t, step, total_mass, total_ke, &
                                   mean_S, mean_T, max_cfl, total_salt, total_heat, &
                                   has_salt, has_temp, mean_age, has_age, budget, &
                                   mass_efp, salt_efp, heat_efp)
      !! Emit one MOM6-style console block from pre-computed totals.  The
      !! caller has already area-weighted + reduced each scalar over its
      !! own state; this routine only latches the t=0 reference (first
      !! call), formats the lines, and runs the panic guards.  When a
      !! `budget` with `active=.true.` is passed, the `Error` column becomes
      !! the boundary-flux + surface-source-corrected residual (see
      !! `conservation_budget_t`) and the outflux / source are shown.
      type(console_stats_t), intent(inout) :: this
      real(wp), intent(in) :: t
         !! Simulation time (s); the compact line reports `Day = t/86400`.
      integer, intent(in) :: step
      real(wp), intent(in) :: total_mass
      real(wp), intent(in) :: total_ke
         !! Absolute kinetic energy (J); the compact "En" column reports
         !! KE per unit mass (m²/s²), the detail block the absolute value.
      real(wp), intent(in) :: mean_S, mean_T
         !! Volume-mean salinity / temperature (compact-line columns).
      real(wp), intent(in) :: max_cfl
      real(wp), intent(in) :: total_salt, total_heat
         !! Column-integrated salt / heat (drift-block totals).
      logical, intent(in) :: has_salt, has_temp
         !! Whether a salinity / temperature tracer is registered.
      real(wp), intent(in) :: mean_age
         !! Volume-mean ideal age (s); printed in days when `has_age`.
      logical, intent(in) :: has_age
      type(conservation_budget_t), intent(in), optional :: budget
         !! Cumulative boundary outflux + surface source that close the
         !! extensive budgets; absent / inactive ⇒ raw drift.
      type(efp_t), intent(in), optional :: mass_efp, salt_efp, heat_efp
         !! PR-32: order-invariant EFP totals for the SAME quantities as
         !! `total_mass`/`total_salt`/`total_heat` (already globally
         !! combined by the caller via `halo_allreduce_efp_list`).  Absent
         !! (default — every existing caller) ⇒ `exact_sums` stays
         !! `.false.` and the console block is byte-identical to pre-PR-32.
         !! When present, latches `mass0_efp`/etc. on the first call and
         !! makes `emit_drift_line` form its residual via `efp_real_diff`
         !! against the EFP reference rather than `(total - ref)` in
         !! double — the SS2.2 fix.

      real(wp) :: day, en_per_mass
      real(wp) :: m_out, s_out, h_out, m_src, s_src, h_src
      logical :: mass_on, salt_on, heat_on
      character(len=256) :: line
      character(len=32) :: cfl_tok

      day = t/86400.0_wp

      mass_on = .false.
      salt_on = .false.
      heat_on = .false.
      m_out = 0.0_wp
      s_out = 0.0_wp
      h_out = 0.0_wp
      m_src = 0.0_wp
      s_src = 0.0_wp
      h_src = 0.0_wp
      if (present(budget)) then
         mass_on = budget%mass_active
         salt_on = budget%salt_active
         heat_on = budget%heat_active
         m_out = budget%mass_out
         s_out = budget%salt_out
         h_out = budget%heat_out
         m_src = budget%mass_src
         s_src = budget%salt_src
         h_src = budget%heat_src
      end if

      ! MOM6-comparable "En" — kinetic energy per unit mass (m²/s²).
      if (total_mass > 0.0_wp) then
         en_per_mass = total_ke/total_mass
      else
         en_per_mass = 0.0_wp
      end if

      if (.not. this%is_initialised) then
         this%mass0 = total_mass
         this%salt0 = total_salt
         this%heat0 = total_heat
         this%ke0 = total_ke
         this%mass_out0 = m_out
         this%salt_out0 = s_out
         this%heat_out0 = h_out
         this%mass_src0 = m_src
         this%salt_src0 = s_src
         this%heat_src0 = h_src
         ! PR-32: latch the EFP reference alongside the double one, iff the
         ! caller supplied it (the ocean `reproducing_sums = .true.` path).
         ! `exact_sums` then gates whether the drift lines below use
         ! `efp_real_diff` (fixed-point difference) instead of the double
         ! subtraction `(total - ref)` — the SS2.2 fix.
         this%exact_sums = present(mass_efp)
         if (present(mass_efp)) this%mass0_efp = mass_efp
         if (present(salt_efp)) this%salt0_efp = salt_efp
         if (present(heat_efp)) this%heat0_efp = heat_efp
         this%is_initialised = .true.
      end if

      ! Measure the budget terms over the same window as the drift: subtract
      ! the reference latched with the totals above.
      m_out = m_out - this%mass_out0
      s_out = s_out - this%salt_out0
      h_out = h_out - this%heat_out0
      m_src = m_src - this%mass_src0
      s_src = s_src - this%salt_src0
      h_src = h_src - this%heat_src0

      ! Compact single-line summary (mirrors MOM6's "MOM Day N:" line).  The
      ! CFL column is a self-contained token (label + value, own leading
      ! spacing).
      write (cfl_tok, "('  MaxCFL ',F8.5)") max_cfl
      if (this%report_thermodynamics) then
         write (line, &
                "('[stats] Day ',F8.3,'  step ',I8,&
                 &'  En ',ES10.3,A,&
                 &'  Mass ',ES12.5,'  Salt ',F8.3,'  Temp ',F8.3)") &
            day, step, en_per_mass, trim(cfl_tok), total_mass, mean_S, mean_T
      else
         write (line, &
                "('[stats] Day ',F8.3,'  step ',I8,&
                 &'  En ',ES10.3,A,&
                 &'  Mass ',ES12.5)") &
            day, step, en_per_mass, trim(cfl_tok), total_mass
      end if
      call logger%info(trim(line))

      ! Conservation detail block — the `Error` is the relative residual vs
      ! the IC snapshot: raw drift `(Q − Q0)/Q0` normally, or the closed
      ! budget `((Q − Q0) + out − src)/Q0` when a budget is active (so open
      ! boundaries + surface forcing don't masquerade as a leak).  Should
      ! stay ~1e-15; growth beyond means real non-conservation.
      if (this%exact_sums .and. present(mass_efp)) then
         call emit_drift_line("Mass", total_mass, this%mass0, m_out, m_src, mass_on, &
                              residual_exact=efp_real_diff(mass_efp, this%mass0_efp))
      else
         call emit_drift_line("Mass", total_mass, this%mass0, m_out, m_src, mass_on)
      end if
      if (this%report_thermodynamics) then
         if (has_salt) then
            if (this%exact_sums .and. present(salt_efp)) then
               call emit_drift_line("Salt", total_salt, this%salt0, s_out, s_src, salt_on, &
                                    residual_exact=efp_real_diff(salt_efp, this%salt0_efp))
            else
               call emit_drift_line("Salt", total_salt, this%salt0, s_out, s_src, salt_on)
            end if
         end if
         if (has_temp) then
            if (this%exact_sums .and. present(heat_efp)) then
               call emit_drift_line("Heat", total_heat, this%heat0, h_out, h_src, heat_on, &
                                    residual_exact=efp_real_diff(heat_efp, this%heat0_efp))
            else
               call emit_drift_line("Heat", total_heat, this%heat0, h_out, h_src, heat_on)
            end if
         end if
      end if

      ! Ideal-age line — volume-mean age in days (passive tracer).
      if (has_age) then
         write (line, "('    Age  : ',F12.4,'  days (volume-mean ideal age)')") &
            mean_age/86400.0_wp
         call logger%info(trim(line))
      end if

      ! Energy block — absolute Joules + growth ratio vs initial KE.
      ! NOT a conservation invariant (wind in, visc/drag out); "Growth"
      ! tracks spin-up. ke0 = 0 for at-rest starts ⇒ ratio falls back to "—".
      if (this%ke0 > tiny(0.0_wp)) then
         write (line, "('    En   : ',ES16.9,'  Growth ',ES10.3,'x')") &
            total_ke, total_ke/this%ke0
      else
         write (line, "('    En   : ',ES16.9,'  Growth         —')") total_ke
      end if
      call logger%info(trim(line))

      ! Early-stop guards — printed after the line above so the last good
      ! (or first-bad) numbers are visible before the abort.
      ! Budget-scalar NaN guard: a poisoned cumulative accumulator (e.g. from
      ! a NaN in the kernel that fills mass_out / salt_src / etc.) would print
      ! a silently-wrong Error line and then recurse into the panic on the next
      ! report.  Test the active budget scalars now so a poisoned accumulator
      ! triggers the panic immediately — and the pre-abort line above shows
      ! which term is NaN.
      if (this%panic_on_nan .and. &
          (any_nan(total_mass, total_ke, total_salt, total_heat, max_cfl) .or. &
           (mass_on .and. (ieee_is_nan(m_out) .or. ieee_is_nan(m_src))) .or. &
           (salt_on .and. (ieee_is_nan(s_out) .or. ieee_is_nan(s_src))) .or. &
           (heat_on .and. (ieee_is_nan(h_out) .or. ieee_is_nan(h_src))))) then
         call logger%error("============================================")
         call logger%error("[panic] NaN detected in console stats — aborting")
         call logger%error("        last good t = "//to_string(t)// &
                           " s (day "//to_string(day)//")")
         call logger%error("        step "//to_string(step))
         call logger%error("============================================")
         error stop "console stats: NaN detected"
      end if
      if (this%panic_on_cfl .and. max_cfl > CFL_PANIC_THRESHOLD) then
         call logger%error("============================================")
         call logger%error("[panic] MaxCFL "//to_string(max_cfl)// &
                           " exceeded threshold "//to_string(CFL_PANIC_THRESHOLD))
         call logger%error("        Model on the cliff — aborting before NaN.")
         call logger%error("        day "//to_string(day)//"  step "//to_string(step))
         call logger%error("============================================")
         error stop "console stats: CFL > panic threshold"
      end if
   end subroutine console_stats_report

   pure function any_nan(a, b, c, d, e) result(res)
      !! `ieee_is_nan` OR-reduced across the five stats scalars.
      real(wp), intent(in) :: a, b, c, d, e
      logical :: res
      res = ieee_is_nan(a) .or. ieee_is_nan(b) .or. ieee_is_nan(c) &
            .or. ieee_is_nan(d) .or. ieee_is_nan(e)
   end function any_nan

   subroutine emit_drift_line(label, total, ref, q_out, q_src, budget_on, residual_exact)
      !! One `<Label> : <total>  Error <residual>` conservation line.  With
      !! `budget_on` the residual closes the budget — `(total − ref) + out −
      !! src` — and the cumulative outflux / surface source are appended;
      !! without it the residual is raw drift `total − ref` and the line is
      !! byte-identical to the pre-budget format.  `label` is a 4-char tag
      !! (Mass / Salt / Heat) so the colons align.
      !!
      !! Robustness note: when `budget_on = .false.` the `out`/`src` arguments
      !! are zeroed locally before the residual is formed, so a caller that
      !! passes non-zero accumulators for an inactive quantity cannot silently
      !! corrupt the raw-drift line with the budget correction.
      !!
      !! PR-32 `residual_exact`: when present, the caller has already formed
      !! `(total − ref)` in FIXED POINT (`efp_real_diff`) rather than double
      !! subtraction — used in place of `(total - ref)` here.  The
      !! `budget_on` out/src correction, the format strings, and
      !! `relative_drift`'s guard are otherwise UNCHANGED (SS5.4: "the
      !! format strings, the budget_on zeroing logic and the relative_drift
      !! guard are unchanged").  Absent ⇒ byte-identical to pre-PR-32.
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: total, ref, q_out, q_src
      logical, intent(in) :: budget_on
      real(wp), intent(in), optional :: residual_exact
      character(len=256) :: line
      real(wp) :: residual, out_eff, src_eff

      ! Zero the budget terms when the budget is inactive so the raw-drift
      ! residual is exact even if the caller passes non-zero accumulators.
      if (budget_on) then
         out_eff = q_out
         src_eff = q_src
      else
         out_eff = 0.0_wp
         src_eff = 0.0_wp
      end if
      if (present(residual_exact)) then
         residual = residual_exact + out_eff - src_eff
      else
         residual = (total - ref) + out_eff - src_eff
      end if
      if (budget_on) then
         write (line, "('    ',A,' : ',ES16.9,'  Error ',ES10.3,&
                       &'  out ',ES10.3,'  src ',ES10.3)") &
            label, total, relative_drift(residual, ref), out_eff, src_eff
      else
         write (line, "('    ',A,' : ',ES16.9,'  Error ',ES10.3)") &
            label, total, relative_drift(residual, ref)
      end if
      call logger%info(trim(line))
   end subroutine emit_drift_line

   pure function relative_drift(change, init) result(rel)
      !! `change / init` with a zero guard for a never-set baseline.
      real(wp), intent(in) :: change, init
      real(wp) :: rel
      real(wp), parameter :: TINY_INIT = 1.0e-30_wp
      if (abs(init) > TINY_INIT) then
         rel = change/init
      else
         rel = 0.0_wp
      end if
   end function relative_drift

end module rki_console_stats
