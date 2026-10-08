!! Unit tests for the PQM (piecewise-quartic) ALE vertical-remap reconstruction
module test_remap_pqm
   !! Analytical tests for `remap_column_pqm` (REMAP_PQM, White & Adcroft 2008):
   !! conservation, high-order accuracy vs PPM, monotonicity, the N<5 PPM
   !! fallback, identity remap, dispatch, and an on-device `do concurrent` run.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, REMAP_PQM, REMAP_PPM
   use rki_remap_column, only: remap_column, remap_column_pqm, remap_column_ppm
   use rki_vcoord, only: parse_remap_method
   implicit none
   private

   public :: collect_remap_pqm_tests

   real(wp), parameter :: PI = 3.14159265358979323846_wp

contains

   subroutine collect_remap_pqm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("pqm_conserves", test_pqm_conserves), &
                  new_unittest("pqm_high_order", test_pqm_high_order), &
                  new_unittest("pqm_monotone", test_pqm_monotone), &
                  new_unittest("pqm_limiter_curvature", test_pqm_limiter_curvature), &
                  new_unittest("pqm_smallN_fallback", test_pqm_smalln_fallback), &
                  new_unittest("pqm_identity", test_pqm_identity), &
                  new_unittest("pqm_dispatch", test_pqm_dispatch), &
                  new_unittest("pqm_on_device", test_pqm_on_device), &
                  new_unittest("pqm_parses_from_namelist_string", test_pqm_parses_from_namelist_string) &
                  ]
   end subroutine collect_remap_pqm_tests

   ! ----------------------------------------------------------------------
   ! T0 — reachability: the "pqm" namelist string reaches REMAP_PQM.
   !      This is the exact string->enum link that was missing for the
   !      life of the feature (PR-9 §2A) — the kernel above had 8 green
   !      tests while nothing in a namelist could select it.
   ! ----------------------------------------------------------------------
   subroutine test_pqm_parses_from_namelist_string(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_remap_method("pqm") == REMAP_PQM, &
                 "parse_remap_method('pqm') must yield REMAP_PQM")
   end subroutine test_pqm_parses_from_namelist_string

   ! ----------------------------------------------------------------------
   ! Helpers
   ! ----------------------------------------------------------------------

   subroutine make_nonuniform_dz(nz, amp, seed, dz)
      !! Deterministic (LCG) non-uniform thicknesses, normalised to sum = 1.
      integer, intent(in) :: nz, seed
      real(wp), intent(in) :: amp
      real(wp), intent(out) :: dz(nz)
      integer :: k
      integer(kind=8) :: lstate
      real(wp) :: r, tot
      lstate = int(seed, 8)
      tot = 0.0_wp
      do k = 1, nz
         lstate = mod(lstate*6364136223846793005_8 + 1442695040888963407_8, &
                      9223372036854775783_8)
         r = real(modulo(lstate, 1000000_8), wp)/1.0e6_wp
         dz(k) = 1.0_wp + amp*r
         tot = tot + dz(k)
      end do
      dz = dz/tot
   end subroutine make_nonuniform_dz

   pure function cell_avg_sin(z_lo, z_hi) result(u)
      !! Cell-average of sin(2*pi*z) over [z_lo, z_hi] (H = 1).
      real(wp), intent(in) :: z_lo, z_hi
      real(wp) :: u
      u = (-cos(2.0_wp*PI*z_hi) + cos(2.0_wp*PI*z_lo))/(2.0_wp*PI*(z_hi - z_lo))
   end function cell_avg_sin

   pure function cell_avg_exp(z_lo, z_hi) result(u)
      !! Cell-average of exp(1.5*z) over [z_lo, z_hi] — smooth, monotone, and
      !! curved.  Monotone (no interior extremum) so the W&A limiter never
      !! flattens the interior to PCM: this isolates the genuine high-order
      !! reconstruction (sin's interior extrema get PCM-flattened by BOTH PQM
      !! and PPM, masking the order with equal max-error there).
      real(wp), intent(in) :: z_lo, z_hi
      real(wp) :: u
      u = (exp(1.5_wp*z_hi) - exp(1.5_wp*z_lo))/(1.5_wp*(z_hi - z_lo))
   end function cell_avg_exp

   ! ----------------------------------------------------------------------
   ! T1 — conservation: column integral preserved to round-off (single remap).
   ! ----------------------------------------------------------------------
   subroutine test_pqm_conserves(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 24
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz), u_new(nz)
      real(wp) :: zc, z_mid, mass_old, mass_new, rel
      integer :: k

      call make_nonuniform_dz(nz, 0.5_wp, 123, dz_old)
      call make_nonuniform_dz(nz, 0.5_wp, 321, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      zc = 0.0_wp
      do k = 1, nz
         z_mid = zc + 0.5_wp*dz_old(k)
         u_old(k) = 1.0_wp + sin(2.0_wp*PI*z_mid)
         zc = zc + dz_old(k)
      end do

      call remap_column_pqm(nz, dz_old, dz_new, u_old, u_new)
      mass_old = sum(u_old*dz_old)
      mass_new = sum(u_new*dz_new)
      rel = abs(mass_new - mass_old)/abs(mass_old)
      call check(error, rel <= 1.0e-13_wp, "PQM must conserve column integral to roundoff")
   end subroutine test_pqm_conserves

   ! ----------------------------------------------------------------------
   ! T2 — high order: PQM is strictly more accurate than PPM on a smooth
   !      profile, the ratio shrinks with N (steeper convergence), and at
   !      N=64 the PQM/PPM error ratio is < 0.1.
   ! ----------------------------------------------------------------------
   subroutine test_pqm_high_order(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nN = 3
      integer :: N_vals(nN), iN, N, k, lo, hi
      real(wp) :: err_pqm(nN), err_ppm(nN), ratio(nN)
      real(wp), allocatable :: dz_s(:), dz_t(:), u_s(:), u_ref(:), u_q(:), u_p(:)
      real(wp), allocatable :: z_s(:), z_t(:)

      N_vals = [16, 32, 64]
      do iN = 1, nN
         N = N_vals(iN)
         allocate (dz_s(N), dz_t(N), u_s(N), u_ref(N), u_q(N), u_p(N))
         allocate (z_s(0:N), z_t(0:N))
         call make_nonuniform_dz(N, 0.4_wp, 11, dz_s)
         call make_nonuniform_dz(N, 0.4_wp, 97, dz_t)
         dz_t(N) = dz_t(N) + (sum(dz_s) - sum(dz_t))   ! exact equal depths
         z_s(0) = 0.0_wp; z_t(0) = 0.0_wp
         do k = 1, N
            z_s(k) = z_s(k - 1) + dz_s(k)
            z_t(k) = z_t(k - 1) + dz_t(k)
         end do
         do k = 1, N
            u_s(k) = cell_avg_exp(z_s(k - 1), z_s(k))
            u_ref(k) = cell_avg_exp(z_t(k - 1), z_t(k))
         end do
         call remap_column_pqm(N, dz_s, dz_t, u_s, u_q)
         call remap_column_ppm(N, dz_s, dz_t, u_s, u_p)
         ! Deep-interior max error (boundary closure decays inward).
         lo = N/4 + 1
         hi = 3*N/4
         err_pqm(iN) = maxval(abs(u_q(lo:hi) - u_ref(lo:hi)))
         err_ppm(iN) = maxval(abs(u_p(lo:hi) - u_ref(lo:hi)))
         ratio(iN) = err_pqm(iN)/err_ppm(iN)
         deallocate (dz_s, dz_t, u_s, u_ref, u_q, u_p, z_s, z_t)
      end do

      ! PQM is strictly more accurate than PPM at every resolution.
      do iN = 1, nN
         call check(error, err_pqm(iN) < err_ppm(iN), &
                    "PQM interior error must be smaller than PPM")
         if (allocated(error)) return
      end do
      ! The accuracy advantage grows with resolution (steeper order).
      call check(error, ratio(nN) < ratio(1), &
                 "PQM/PPM error ratio must shrink as N grows (higher order)")
      if (allocated(error)) return
      ! At N=64, PQM is at least 10x more accurate than PPM.
      call check(error, ratio(nN) < 0.1_wp, &
                 "PQM/PPM error ratio must be < 0.1 at N=64")
   end subroutine test_pqm_high_order

   ! ----------------------------------------------------------------------
   ! T3 — monotonicity: a step profile produces no overshoot beyond the
   !      source range (the W&A limiter holds).
   ! ----------------------------------------------------------------------
   subroutine test_pqm_monotone(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 30
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz), u_new(nz)
      real(wp) :: zc, z_mid, qmn, qmx

      integer :: k
      call make_nonuniform_dz(nz, 0.3_wp, 77, dz_old)
      call make_nonuniform_dz(nz, 0.3_wp, 707, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      zc = 0.0_wp
      do k = 1, nz
         z_mid = zc + 0.5_wp*dz_old(k)
         if (z_mid < 0.5_wp) then
            u_old(k) = 1.0_wp
         else
            u_old(k) = 5.0_wp
         end if
         zc = zc + dz_old(k)
      end do
      qmn = minval(u_old); qmx = maxval(u_old)

      call remap_column_pqm(nz, dz_old, dz_new, u_old, u_new)
      call check(error, minval(u_new) >= qmn - 1.0e-9_wp .and. &
                 maxval(u_new) <= qmx + 1.0e-9_wp, &
                 "PQM must not overshoot beyond the source range on a step profile")
   end subroutine test_pqm_monotone

   subroutine test_pqm_limiter_curvature(error)
      !! Drives the limiter's curvature-collapse + post-collapse-reset paths
      !! (the bug-prone L268-318 block).  A STEEP but MONOTONE sigmoid on a
      !! COARSE non-uniform grid: the extremum-flatten path does NOT fire
      !! (no interior extremum), but the high-curvature quartic in the ramp
      !! cells is non-monotone -> the curvature-root test collapses the
      !! inflexion.  A buggy collapse/reset would inject a NON-monotone dip
      !! (overshoot or a local reversal) -> caught here.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 12
      real(wp) :: dz_old(nz), dz_new(nz), u_old(nz), u_new(nz)
      real(wp) :: zc, z_mid, qmn, qmx, s0, s1
      integer :: k
      logical :: monotone
      call make_nonuniform_dz(nz, 0.5_wp, 131, dz_old)   ! coarse + non-uniform
      call make_nonuniform_dz(nz, 0.5_wp, 919, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      zc = 0.0_wp
      do k = 1, nz
         z_mid = zc + 0.5_wp*dz_old(k)
         ! steep monotone sigmoid (no interior extremum, high mid-curvature)
         u_old(k) = 1.0_wp + 4.0_wp/(1.0_wp + exp(-14.0_wp*(z_mid - 0.5_wp)))
         zc = zc + dz_old(k)
      end do
      qmn = minval(u_old); qmx = maxval(u_old)
      s0 = sum(u_old*dz_old)
      call remap_column_pqm(nz, dz_old, dz_new, u_old, u_new)
      s1 = sum(u_new*dz_new)
      call check(error, abs(s1 - s0) < 1.0e-12_wp, &
                 "PQM limiter-collapse: column integral conserved")
      if (allocated(error)) return
      call check(error, minval(u_new) >= qmn - 1.0e-9_wp .and. &
                 maxval(u_new) <= qmx + 1.0e-9_wp, &
                 "PQM limiter-collapse: no overshoot beyond source range")
      if (allocated(error)) return
      ! monotonicity preserved: a buggy collapse/reset would create a dip
      monotone = .true.
      do k = 2, nz
         if (u_new(k) < u_new(k - 1) - 1.0e-9_wp) monotone = .false.
      end do
      call check(error, monotone, &
                 "PQM limiter-collapse: remap of a monotone profile stays monotone")
   end subroutine test_pqm_limiter_curvature

   ! ----------------------------------------------------------------------
   ! T4 — N<5 fallback: bit-identical to remap_column_ppm.
   ! ----------------------------------------------------------------------
   subroutine test_pqm_smalln_fallback(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: nz, k
      real(wp) :: dz_old(4), dz_new(4), q_old(4), a(4), b(4)

      ! Exercise nz = 2, 3, 4 (all < 5 → PPM fallback).  nz=1 is the trivial
      ! identity in both, also covered.
      do nz = 2, 4
         call make_nonuniform_dz(nz, 0.5_wp, 42 + nz, dz_old(1:nz))
         call make_nonuniform_dz(nz, 0.5_wp, 424 + nz, dz_new(1:nz))
         dz_new(nz) = dz_new(nz) + (sum(dz_old(1:nz)) - sum(dz_new(1:nz)))
         do k = 1, nz
            q_old(k) = sin(real(k, wp))
         end do
         call remap_column_pqm(nz, dz_old(1:nz), dz_new(1:nz), q_old(1:nz), a(1:nz))
         call remap_column_ppm(nz, dz_old(1:nz), dz_new(1:nz), q_old(1:nz), b(1:nz))
         do k = 1, nz
            call check(error, a(k) == b(k), &
                       "PQM with N<5 must be bit-identical to PPM")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_pqm_smalln_fallback

   ! ----------------------------------------------------------------------
   ! T5 — identity remap: same grid recovers the input.
   ! ----------------------------------------------------------------------
   subroutine test_pqm_identity(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 9
      real(wp) :: dz(nz), q_old(nz), q_new(nz)
      integer :: k

      dz = [1.0_wp, 2.0_wp, 1.5_wp, 3.0_wp, 2.0_wp, 1.0_wp, 1.2_wp, 0.8_wp, 1.4_wp]
      q_old = [10.0_wp, 20.0_wp, 30.0_wp, 25.0_wp, 15.0_wp, 18.0_wp, 22.0_wp, 19.0_wp, 21.0_wp]
      call remap_column_pqm(nz, dz, dz, q_old, q_new)
      do k = 1, nz
         call check(error, abs(q_new(k) - q_old(k)) < 1.0e-11_wp, &
                    "PQM identity remap must recover q exactly")
         if (allocated(error)) return
      end do
   end subroutine test_pqm_identity

   ! ----------------------------------------------------------------------
   ! T6 — dispatch: REMAP_PQM routes to remap_column_pqm bit-identically.
   ! ----------------------------------------------------------------------
   subroutine test_pqm_dispatch(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nz = 10
      real(wp) :: dz_old(nz), dz_new(nz), q_old(nz), a(nz), b(nz)
      integer :: k

      call make_nonuniform_dz(nz, 0.5_wp, 42, dz_old)
      call make_nonuniform_dz(nz, 0.5_wp, 424, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      do k = 1, nz
         q_old(k) = sin(real(k, wp))
      end do
      call remap_column(REMAP_PQM, nz, dz_old, dz_new, q_old, a)
      call remap_column_pqm(nz, dz_old, dz_new, q_old, b)
      do k = 1, nz
         call check(error, a(k) == b(k), "dispatch REMAP_PQM must equal direct call")
         if (allocated(error)) return
      end do
   end subroutine test_pqm_dispatch

   ! ----------------------------------------------------------------------
   ! T7 — on-device: run a batch of PQM remaps inside `do concurrent` on the
   !      NVHPC build; output must be finite and conservative per column.
   ! ----------------------------------------------------------------------
   subroutine test_pqm_on_device(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ncol = 8, nz = 12
      real(wp) :: dz_old(nz), dz_new(nz)
      real(wp) :: q_in(nz, ncol), q_out(nz, ncol)
      real(wp) :: mass_in(ncol), mass_out(ncol)
      real(wp) :: zc, z_mid, rel
      integer :: i, k
      logical :: ok

      call make_nonuniform_dz(nz, 0.5_wp, 555, dz_old)
      call make_nonuniform_dz(nz, 0.5_wp, 999, dz_new)
      dz_new(nz) = dz_new(nz) + (sum(dz_old) - sum(dz_new))
      do i = 1, ncol
         zc = 0.0_wp
         do k = 1, nz
            z_mid = zc + 0.5_wp*dz_old(k)
            q_in(k, i) = real(i, wp) + sin(2.0_wp*PI*z_mid + 0.3_wp*real(i, wp))
            zc = zc + dz_old(k)
         end do
      end do

      !$omp target data map(to: dz_old, dz_new, q_in) map(from: q_out)
      !$omp target teams distribute parallel do
      !$omp target teams distribute parallel do
      do i = 1, ncol
         call remap_column_pqm(nz, dz_old, dz_new, q_in(:, i), q_out(:, i))
      end do
      !$omp end target teams distribute parallel do
      !$omp end target data

      ok = .true.
      do i = 1, ncol
         mass_in(i) = sum(q_in(:, i)*dz_old)
         mass_out(i) = sum(q_out(:, i)*dz_new)
         do k = 1, nz
            if (q_out(k, i) /= q_out(k, i)) ok = .false.   ! NaN check
         end do
         rel = abs(mass_out(i) - mass_in(i))/abs(mass_in(i))
         if (rel > 1.0e-12_wp) ok = .false.
      end do
      call check(error, ok, "on-device PQM remap must be finite and conservative per column")
   end subroutine test_pqm_on_device

end module test_remap_pqm
