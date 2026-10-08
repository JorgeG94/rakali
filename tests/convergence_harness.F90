!! Grid-refinement error norms and observed convergence order.
module convergence_harness
   !! Test-only infrastructure for verifying that a scheme achieves its design
   !! order on a smooth solution.  Two pieces:
   !!
   !!   * **Error norms** — `error_norms_uniform` (structured, constant cell
   !!     measure) and `error_norms_weighted` (unstructured, per-cell areas).
   !!     Both return L1, L2 and Linf as INTEGRALS over the domain, i.e.
   !!     `L1 = sum |e| w`, `L2 = sqrt(sum e^2 w)`, `Linf = max |e|`.
   !!   * **Observed order** — `convergence_rates`, giving the pairwise rates
   !!     between consecutive refinement levels and the overall coarse-to-fine
   !!     rate, both from `p = log(E_c/E_f) / log(h_c/h_f)`.
   !!
   !! ## Why both pairwise and overall
   !!
   !! They fail differently, and reporting only one hides real problems.
   !!
   !! A single **pairwise** rate is noisy: it is formed from two numbers, so
   !! one contaminated level (a limiter kicking in, a boundary effect, an
   !! under-resolved coarsest grid) can drag it far from the design order in
   !! either direction — including *above* it, which looks like good news and
   !! is not. The **overall** coarse-to-fine rate averages that noise out, but
   !! it will happily report a healthy number for a series that plateaus at
   !! the fine end, which is the signature of hitting a different error floor
   !! (round-off, an unconverged solver tolerance, a fixed time-step error).
   !! Asserting the overall rate and *inspecting* the pairwise series catches
   !! both; that is the pattern the campaign's phases are expected to follow.
   !!
   !! ## Normalisation does not affect the rate
   !!
   !! The norms are integrals rather than means. Any constant factor cancels
   !! in `E_c/E_f`, so a consumer that prefers means may divide by the domain
   !! measure without changing any reported order.
   !!
   !! ## Refinement ratio
   !!
   !! `h` is the grid SPACING and must be strictly decreasing (coarse first).
   !! Consumers that hold a list of cell counts should pass `L/n` — or use
   !! `spacing_from_counts`, which does it and keeps the ordering convention
   !! visible at the call site. The ratio need not be 2, and need not be
   !! uniform across the series.
   !!
   !! Assumed-shape dummies are used throughout: this is host-only test
   !! infrastructure, never a `do concurrent` device kernel, so the
   !! explicit-shape rule in CLAUDE.md does not apply.
   use rki_constants, only: wp
   implicit none
   private

   public :: error_norms_uniform
   public :: error_norms_weighted
   public :: convergence_rates
   public :: spacing_from_counts

contains

   pure subroutine error_norms_uniform(numeric, exact, dcell, l1, l2, linf)
      !! Error norms on a uniform grid: every cell carries the same measure
      !! `dcell` (dx in 1-D, dx*dy in 2-D, dx*dy*dz in 3-D).
      real(wp), intent(in) :: numeric(:)   ! assumed-shape-ok: host test infra
      real(wp), intent(in) :: exact(:)     ! assumed-shape-ok: host test infra
      real(wp), intent(in) :: dcell
         !! Measure of one cell (m, m^2 or m^3).
      real(wp), intent(out) :: l1, l2, linf

      real(wp) :: e
      integer :: i

      l1 = 0.0_wp
      l2 = 0.0_wp
      linf = 0.0_wp
      do i = 1, min(size(numeric), size(exact))
         e = abs(numeric(i) - exact(i))
         l1 = l1 + e*dcell
         l2 = l2 + e*e*dcell
         linf = max(linf, e)
      end do
      l2 = sqrt(l2)
   end subroutine error_norms_uniform

   pure subroutine error_norms_weighted(numeric, exact, weight, l1, l2, linf)
      !! Error norms with a per-cell measure — the unstructured case, where
      !! triangle areas vary.
      real(wp), intent(in) :: numeric(:)   ! assumed-shape-ok: host test infra
      real(wp), intent(in) :: exact(:)     ! assumed-shape-ok: host test infra
      real(wp), intent(in) :: weight(:)    ! assumed-shape-ok: host test infra
      real(wp), intent(out) :: l1, l2, linf

      real(wp) :: e
      integer :: i

      l1 = 0.0_wp
      l2 = 0.0_wp
      linf = 0.0_wp
      do i = 1, min(size(numeric), size(exact), size(weight))
         e = abs(numeric(i) - exact(i))
         l1 = l1 + e*weight(i)
         l2 = l2 + e*e*weight(i)
         linf = max(linf, e)
      end do
      l2 = sqrt(l2)
   end subroutine error_norms_weighted

   pure subroutine spacing_from_counts(domain_length, counts, h)
      !! `h(i) = domain_length / counts(i)`.  A one-liner, but having it named
      !! keeps the coarse-first ordering convention explicit at call sites and
      !! stops each consumer re-deriving whether to use counts or spacings —
      !! the two existing hand-rolled rate computations in this suite disagree
      !! on exactly that.
      real(wp), intent(in) :: domain_length
      integer, intent(in) :: counts(:)     ! assumed-shape-ok: host test infra
      real(wp), intent(out) :: h(:)        ! assumed-shape-ok: host test infra

      integer :: i

      do i = 1, min(size(counts), size(h))
         if (counts(i) > 0) then
            h(i) = domain_length/real(counts(i), wp)
         else
            h(i) = 0.0_wp
         end if
      end do
   end subroutine spacing_from_counts

   pure subroutine convergence_rates(h, err, p_pair, p_overall, ok)
      !! Observed order from a refinement series.
      !!
      !!   `p = log(E_coarse / E_fine) / log(h_coarse / h_fine)`
      !!
      !! `p_pair(i)` is the rate between levels `i` and `i+1`; `p_overall` is
      !! the rate from the coarsest to the finest level.
      !!
      !! `ok` is `.false.`, and the outputs are zeroed, if the input cannot
      !! support a rate: fewer than two levels, a non-positive or non-finite
      !! error (an exactly-zero error means the scheme reproduced the oracle
      !! and the order is undefined, not infinite), or an `h` series that is
      !! not strictly decreasing. Callers should assert `ok` before the rate —
      !! silently returning a plausible number from degenerate input is how a
      !! convergence test starts certifying nothing.
      real(wp), intent(in) :: h(:)         ! assumed-shape-ok: host test infra
         !! Grid spacing per level, strictly DECREASING (coarsest first).
      real(wp), intent(in) :: err(:)       ! assumed-shape-ok: host test infra
         !! Error norm per level, same ordering and length as `h`.
      real(wp), intent(out) :: p_pair(:)   ! assumed-shape-ok: host test infra
         !! Pairwise rates; needs at least `size(h) - 1` elements.
      real(wp), intent(out) :: p_overall
      logical, intent(out) :: ok

      integer :: i, n

      p_pair = 0.0_wp
      p_overall = 0.0_wp
      ok = .false.

      n = min(size(h), size(err))
      if (n < 2) return
      if (size(p_pair) < n - 1) return

      do i = 1, n
         if (.not. (err(i) > 0.0_wp)) return   ! also rejects NaN
         if (.not. (h(i) > 0.0_wp)) return
      end do
      do i = 1, n - 1
         if (.not. (h(i) > h(i + 1))) return   ! strictly decreasing
      end do

      do i = 1, n - 1
         p_pair(i) = log(err(i)/err(i + 1))/log(h(i)/h(i + 1))
      end do
      p_overall = log(err(1)/err(n))/log(h(1)/h(n))
      ok = .true.
   end subroutine convergence_rates

end module convergence_harness
