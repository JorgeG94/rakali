!! Conservative vertical remapping for ALE coordinates
module rki_remap_column
   !! Pure per-column conservative vertical remap, shared by both backends.
   !! Reconstruction methods:
   !!   PCM    — piecewise constant (donor cell)
   !!   PLM    — piecewise linear, minmod limiter
   !!   PPM    — piecewise parabolic (Colella & Woodward 1984)
   !!   PPM_H4 — PPM with non-uniform 4th-order edge values (White & Adcroft 2008)
   !!   PQM    — piecewise quartic (White & Adcroft 2008)
   !!
   !! All routines are pure with NZ_STACK_MAX stack workspace so they run
   !! inside `do concurrent` (one GPU thread per column, sequential O(nz)
   !! vertical sweep). Conservation: sum(q_new*dz_new) = sum(q_old*dz_old)
   !! to machine precision when sum(dz_old) = sum(dz_new).
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, REMAP_PCM, REMAP_PLM, REMAP_PPM, &
                            REMAP_PPM_H4, REMAP_PQM
#else
   use rki_constants, only: NZ_STACK_MAX, wp, REMAP_PCM, REMAP_PLM, REMAP_PPM, &
                            REMAP_PPM_H4, REMAP_PQM
#endif
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: remap_column
   public :: remap_column_pcm
   public :: remap_column_plm
   public :: remap_column_ppm
   public :: remap_column_ppm_h4
   public :: remap_column_pqm

   real(wp), parameter :: H_NEGLECT = 1.0e-30_wp
      !! Tiny thickness floor for the PPM_H4 edge singularity guard; only fires
      !! on degenerate (vanishing-layer) columns. White & Adcroft 2008.
   real(wp), parameter :: H_MIN_FRAC = 1.0e-5_wp
      !! Relative thickness floor in the PPM_H4 edge stencil when a consecutive
      !! thickness pair sums to ~0.
   real(wp), parameter :: H_RATIO_FLOOR = 1.0e-12_wp
      !! Relative thickness-ratio floor in the PQM implicit-h4 edge-value
      !! stencil (White & Adcroft 2008, roundoff-safe).
   real(wp), parameter :: PQM_MIN_FRAC = 1.0e-6_wp
      !! Boundary-closure guard for the PQM 4th-order one-sided end fit
      !! (White & Adcroft 2008).

contains

   pure subroutine remap_column(method, nz, dz_old, dz_new, q_old, q_new)
      !$acc routine seq
      !! Dispatch to the requested remapping method.
      integer, intent(in) :: method
         !! REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, or REMAP_PQM
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz_old(nz)
      real(wp), intent(in) :: dz_new(nz)
      real(wp), intent(in) :: q_old(nz)
      real(wp), intent(out) :: q_new(nz)

      select case (method)
      case (REMAP_PCM)
         call remap_column_pcm(nz, dz_old, dz_new, q_old, q_new)
      case (REMAP_PLM)
         call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)
      case (REMAP_PPM)
         call remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)
      case (REMAP_PPM_H4)
         call remap_column_ppm_h4(nz, dz_old, dz_new, q_old, q_new)
      case (REMAP_PQM)
         call remap_column_pqm(nz, dz_old, dz_new, q_old, q_new)
      case default
         call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)
      end select
   end subroutine remap_column

   pure subroutine remap_column_pcm(nz, dz_old, dz_new, q_old, q_new)
      !$acc routine seq
      !! Piecewise-constant (donor cell) remap. Diffusive, guaranteed monotone.
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz_old(nz)
         !! Old layer thicknesses (must sum to same total as dz_new)
      real(wp), intent(in) :: dz_new(nz)
         !! New layer thicknesses
      real(wp), intent(in) :: q_old(nz)
         !! Old cell-average scalar values
      real(wp), intent(out) :: q_new(nz)
         !! New cell-average scalar values (conservative)

      real(wp) :: z_old(0:NZ_STACK_MAX), z_new(0:NZ_STACK_MAX)
      real(wp) :: z_lo, z_hi, overlap, integral
      integer :: k, ko, ko_start

      ! Build interface positions (cumulative sum from bottom)
      z_old(0) = 0.0_wp
      z_new(0) = 0.0_wp
      do k = 1, nz
         z_old(k) = z_old(k - 1) + dz_old(k)
         z_new(k) = z_new(k - 1) + dz_new(k)
      end do

      ! Sweep: for each new layer, integrate PCM from old layers
      ko_start = 1
      do k = 1, nz
         if (dz_new(k) <= 0.0_wp) then
            q_new(k) = 0.0_wp
            cycle
         end if

         integral = 0.0_wp
         do ko = ko_start, nz
            z_lo = max(z_new(k - 1), z_old(ko - 1))
            z_hi = min(z_new(k), z_old(ko))
            overlap = z_hi - z_lo

            if (overlap <= 0.0_wp) then
               if (z_old(ko) > z_new(k)) exit
               cycle
            end if

            integral = integral + q_old(ko)*overlap

            ! Advance scan: if old layer fully consumed, next new layer
            ! can start from the next old layer
            if (z_old(ko) <= z_new(k)) ko_start = ko
         end do

         q_new(k) = integral/dz_new(k)
      end do
   end subroutine remap_column_pcm

   pure subroutine remap_column_plm(nz, dz_old, dz_new, q_old, q_new)
      !$acc routine seq
      !! Piecewise-linear (minmod-limited) remap. Monotone (no new extrema).
      !! Per old layer k: q_hat(xi) = q(k) + slope(k)*(2*xi - 1), xi in [0,1],
      !! slope(k) = 0.5*minmod(q(k+1)-q(k), q(k)-q(k-1)).
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz_old(nz)
         !! Old layer thicknesses
      real(wp), intent(in) :: dz_new(nz)
         !! New layer thicknesses
      real(wp), intent(in) :: q_old(nz)
         !! Old cell-average scalar values
      real(wp), intent(out) :: q_new(nz)
         !! New cell-average scalar values (conservative)

      real(wp) :: z_old(0:NZ_STACK_MAX), z_new(0:NZ_STACK_MAX)
      real(wp) :: slope(NZ_STACK_MAX)
      real(wp) :: z_lo, z_hi, overlap, integral
      real(wp) :: xi_lo, xi_hi, dq_l, dq_r
      integer :: k, ko, ko_start

      ! Single-layer case: identity remap
      if (nz == 1) then
         q_new(1) = q_old(1)
         return
      end if

      ! Build interface positions
      z_old(0) = 0.0_wp
      z_new(0) = 0.0_wp
      do k = 1, nz
         z_old(k) = z_old(k - 1) + dz_old(k)
         z_new(k) = z_new(k - 1) + dz_new(k)
      end do

      ! Compute minmod-limited slopes
      ! slope(k) = half the limited difference across the layer
      slope(1) = 0.0_wp
      do k = 2, nz - 1
         dq_l = q_old(k) - q_old(k - 1)
         dq_r = q_old(k + 1) - q_old(k)
         if (dq_l*dq_r > 0.0_wp) then
            slope(k) = 0.5_wp*sign(min(abs(dq_l), abs(dq_r)), dq_l)
         else
            slope(k) = 0.0_wp
         end if
      end do
      slope(nz) = 0.0_wp

      ! Sweep: for each new layer, integrate PLM from old layers
      ko_start = 1
      do k = 1, nz
         if (dz_new(k) <= 0.0_wp) then
            q_new(k) = 0.0_wp
            cycle
         end if

         integral = 0.0_wp
         do ko = ko_start, nz
            z_lo = max(z_new(k - 1), z_old(ko - 1))
            z_hi = min(z_new(k), z_old(ko))
            overlap = z_hi - z_lo

            if (overlap <= 0.0_wp) then
               if (z_old(ko) > z_new(k)) exit
               cycle
            end if

            if (dz_old(ko) > 0.0_wp) then
               ! Normalised coordinates within old layer ko
               xi_lo = (z_lo - z_old(ko - 1))/dz_old(ko)
               xi_hi = (z_hi - z_old(ko - 1))/dz_old(ko)

               ! Integral of q_hat(xi) = q + slope*(2*xi - 1) over [xi_lo, xi_hi]
               ! = (xi_hi - xi_lo) * (q + slope*(xi_hi + xi_lo - 1))
               ! scaled to physical space: * dz_old(ko)
               integral = integral + overlap* &
                          (q_old(ko) + slope(ko)*(xi_lo + xi_hi - 1.0_wp))
            else
               integral = integral + q_old(ko)*overlap
            end if

            if (z_old(ko) <= z_new(k)) ko_start = ko
         end do

         q_new(k) = integral/dz_new(k)
      end do
   end subroutine remap_column_plm

   pure subroutine remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)
      !$acc routine seq
      !! Piecewise-parabolic (Colella & Woodward 1984) remap.
      !! Per old layer k, xi in [0,1]:
      !!   q_hat(xi) = q_L + xi*(q_R - q_L + q6*(1 - xi)), q6 = 6*q_bar - 3*(q_L+q_R)
      !! Edge values: 4th-order interp + CW monotonicity limiting; boundary
      !! layers fall back to PLM-quality edges.
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz_old(nz)
         !! Old layer thicknesses
      real(wp), intent(in) :: dz_new(nz)
         !! New layer thicknesses
      real(wp), intent(in) :: q_old(nz)
         !! Old cell-average scalar values
      real(wp), intent(out) :: q_new(nz)
         !! New cell-average scalar values (conservative)

      real(wp) :: z_old(0:NZ_STACK_MAX), z_new(0:NZ_STACK_MAX)
      real(wp) :: q_L(NZ_STACK_MAX), q_R(NZ_STACK_MAX), q6(NZ_STACK_MAX)
      real(wp) :: z_lo, z_hi, overlap, integral
      real(wp) :: xi_lo, xi_hi
      real(wp) :: edge, dq, dq_l, dq_r, q_min, q_max
      integer :: k, ko, ko_start

      ! Trivial cases
      if (nz == 1) then
         q_new(1) = q_old(1)
         return
      end if
      if (nz == 2) then
         ! With only 2 layers, PPM reduces to PLM
         call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)
         return
      end if

      ! Build interface positions
      z_old(0) = 0.0_wp
      z_new(0) = 0.0_wp
      do k = 1, nz
         z_old(k) = z_old(k - 1) + dz_old(k)
         z_new(k) = z_new(k - 1) + dz_new(k)
      end do

      ! ---- Step 1: Compute unlimited edge values via 4th-order interp ----
      ! Interior edges (between layers k and k+1) use the 4-cell stencil.
      ! Store in q_R(k) = right edge of layer k = left edge of layer k+1.

      ! Boundary: PCM (layer 1 left, layer nz right)
      q_L(1) = q_old(1)
      q_R(nz) = q_old(nz)

      ! Layer 1 right edge = layer 2 left edge: use 3-cell stencil (one-sided)
      q_R(1) = 0.5_wp*(q_old(1) + q_old(2))
      ! Layer nz left edge = layer nz-1 right edge: use 3-cell stencil
      q_L(nz) = 0.5_wp*(q_old(nz - 1) + q_old(nz))

      ! Interior edges: 4th-order Colella-Woodward interpolation
      ! For uniform layers this gives (7/12)(q_k + q_{k+1}) - (1/12)(q_{k-1} + q_{k+2})
      ! For non-uniform layers, use the simpler weighted average
      do k = 2, nz - 1
         edge = 0.5_wp*(q_old(k) + q_old(k + 1))
         if (k >= 2 .and. k + 1 <= nz) then
            ! Add 4th-order correction when stencil is available
            dq_l = q_old(k) - q_old(k - 1)
            dq_r = q_old(k + 1) - q_old(k)
            if (k - 1 >= 1 .and. k + 2 <= nz) then
               edge = (7.0_wp/12.0_wp)*(q_old(k) + q_old(k + 1)) &
                      - (1.0_wp/12.0_wp)*(q_old(k - 1) + q_old(k + 2))
            end if
         end if
         q_R(k) = edge
         q_L(k + 1) = edge
      end do

      ! Layer 2 left edge (if nz >= 3, was set above; otherwise use average)
      if (nz >= 3) then
         q_L(2) = q_R(1)
      end if

      ! ---- Step 2: Colella-Woodward monotonicity limiting ----
      do k = 1, nz
         q_min = q_old(k)
         q_max = q_old(k)
         if (k > 1) then
            q_min = min(q_min, q_old(k - 1))
            q_max = max(q_max, q_old(k - 1))
         end if
         if (k < nz) then
            q_min = min(q_min, q_old(k + 1))
            q_max = max(q_max, q_old(k + 1))
         end if

         ! Clip edges to local bounds
         q_L(k) = max(q_min, min(q_max, q_L(k)))
         q_R(k) = max(q_min, min(q_max, q_R(k)))

         ! CW monotonicity: if the cell is a local extremum, flatten
         dq = q_R(k) - q_L(k)
         dq_l = q_old(k) - q_L(k)
         dq_r = q_R(k) - q_old(k)
         if (dq_l*dq_r <= 0.0_wp) then
            ! Local extremum: flatten to PCM
            q_L(k) = q_old(k)
            q_R(k) = q_old(k)
         else
            ! Check if parabola overshoots
            ! q6 = 6*q_bar - 3*(q_L + q_R)
            ! The parabola has an extremum inside [0,1] if q6*(q_R - q_L) < 0
            ! and the extremum value exceeds the local bounds.
            q6(k) = 6.0_wp*q_old(k) - 3.0_wp*(q_L(k) + q_R(k))
            if (abs(q6(k)) > abs(dq)) then
               if (q6(k)*dq > 0.0_wp) then
                  ! Overshoot near left edge: adjust q_L
                  q_L(k) = 3.0_wp*q_old(k) - 2.0_wp*q_R(k)
               else
                  ! Overshoot near right edge: adjust q_R
                  q_R(k) = 3.0_wp*q_old(k) - 2.0_wp*q_L(k)
               end if
            end if
         end if

         ! Recompute q6 after limiting
         q6(k) = 6.0_wp*q_old(k) - 3.0_wp*(q_L(k) + q_R(k))
      end do

      ! ---- Step 3: Integrate parabolic reconstruction over new layers ----
      ko_start = 1
      do k = 1, nz
         if (dz_new(k) <= 0.0_wp) then
            q_new(k) = 0.0_wp
            cycle
         end if

         integral = 0.0_wp
         do ko = ko_start, nz
            z_lo = max(z_new(k - 1), z_old(ko - 1))
            z_hi = min(z_new(k), z_old(ko))
            overlap = z_hi - z_lo

            if (overlap <= 0.0_wp) then
               if (z_old(ko) > z_new(k)) exit
               cycle
            end if

            if (dz_old(ko) > 0.0_wp) then
               xi_lo = (z_lo - z_old(ko - 1))/dz_old(ko)
               xi_hi = (z_hi - z_old(ko - 1))/dz_old(ko)

               ! Parabolic integral
               integral = integral + dz_old(ko)*( &
                          (xi_hi - xi_lo)*q_L(ko) &
                          + 0.5_wp*(xi_hi*xi_hi - xi_lo*xi_lo)*(q_R(ko) - q_L(ko) + q6(ko)) &
                          - (xi_hi*xi_hi*xi_hi - xi_lo*xi_lo*xi_lo)*q6(ko)/3.0_wp)
            else
               integral = integral + q_old(ko)*overlap
            end if

            if (z_old(ko) <= z_new(k)) ko_start = ko
         end do

         q_new(k) = integral/dz_new(k)
      end do
   end subroutine remap_column_ppm

   pure subroutine remap_column_ppm_h4(nz, dz_old, dz_new, q_old, q_new)
      !$acc routine seq
      !! PPM with non-uniform 4th-order (H4) edge values (White & Adcroft 2008).
      !! As `remap_column_ppm` but the interior edge estimate is the
      !! thickness-weighted exactly-4th-order stencil (reduces to PPM's
      !! (7/12,-1/12) on uniform layers), cutting spurious diapycnal mixing per
      !! remap (Ilicak et al. 2012). Limiter/reconstruction/integration are
      !! identical to PPM. H4 edge algebra is inlined (helper extraction costs
      !! 4-6% on this hot kernel). Boundary edges: outermost = PCM,
      !! second-from-boundary = non-uniform 3-cell (H3) quadratic; CW limiter
      !! clamps all edges to local monotone bounds.
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz_old(nz)
         !! Old layer thicknesses
      real(wp), intent(in) :: dz_new(nz)
         !! New layer thicknesses
      real(wp), intent(in) :: q_old(nz)
         !! Old cell-average scalar values
      real(wp), intent(out) :: q_new(nz)
         !! New cell-average scalar values (conservative)

      real(wp) :: z_old(0:NZ_STACK_MAX), z_new(0:NZ_STACK_MAX)
      real(wp) :: q_L(NZ_STACK_MAX), q_R(NZ_STACK_MAX), q6(NZ_STACK_MAX)
      real(wp) :: z_lo, z_hi, overlap, integral
      real(wp) :: xi_lo, xi_hi
      real(wp) :: dq, dq_l, dq_r, q_min, q_max
      real(wp) :: h0, h1, h2, h3, hf, h_sum
      real(wp) :: h01, h12, h23, h012, h123, h0123
      real(wp) :: f1, f2, f3, et1, et2, et3
      real(wp) :: m11, m12, m13, m21, m22, m23, m31, m32, m33, det
      real(wp) :: z1, z2, z3, ca, cb, cc, edge_val
      integer :: k, ko, ko_start

      ! Trivial cases (identical to PPM)
      if (nz == 1) then
         q_new(1) = q_old(1)
         return
      end if
      if (nz == 2) then
         ! With only 2 layers, PPM reduces to PLM
         call remap_column_plm(nz, dz_old, dz_new, q_old, q_new)
         return
      end if

      ! Build interface positions
      z_old(0) = 0.0_wp
      z_new(0) = 0.0_wp
      do k = 1, nz
         z_old(k) = z_old(k - 1) + dz_old(k)
         z_new(k) = z_new(k - 1) + dz_new(k)
      end do

      ! ---- Step 1 (PPM_H4): non-uniform 4th-order edge values ----
      ! Store edges in q_R(k) = right edge of layer k = q_L(k+1).

      ! Outermost edges: PCM (layer 1 left, layer nz right)
      q_L(1) = q_old(1)
      q_R(nz) = q_old(nz)

      ! Second-from-boundary edges: non-uniform 3-cell H3 quadratic.
      ! Bed side: edge between layer 1 and 2 — fit a quadratic through layers
      ! 1,2,3 (cell-average constraints) and evaluate at the 1|2 interface.
      ! z-coordinate with the left of layer 1 at 0:
      h1 = dz_old(1)
      h2 = dz_old(2)
      h3 = dz_old(3)
      ! H3 boundary singularity guard (mirror of the interior floor)
      if (h1 <= 0.0_wp .or. h2 <= 0.0_wp .or. h3 <= 0.0_wp) then
         hf = H_MIN_FRAC*max(H_NEGLECT, dz_old(1) + dz_old(2) + dz_old(3))
         h1 = max(h1, hf)
         h2 = max(h2, hf)
         h3 = max(h3, hf)
      end if
      z1 = h1
      z2 = h1 + h2
      z3 = h1 + h2 + h3
      ! Cell moments [int z^2 / h, int z / h, 1] over each layer:
      m11 = (z1*z1*z1)/(3.0_wp*h1)
      m12 = (z1*z1)/(2.0_wp*h1)
      m13 = 1.0_wp
      m21 = (z2*z2*z2 - z1*z1*z1)/(3.0_wp*h2)
      m22 = (z2*z2 - z1*z1)/(2.0_wp*h2)
      m23 = 1.0_wp
      m31 = (z3*z3*z3 - z2*z2*z2)/(3.0_wp*h3)
      m32 = (z3*z3 - z2*z2)/(2.0_wp*h3)
      m33 = 1.0_wp
      det = m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31)
      ! Cramer's rule for the quadratic coefficients ca*z^2 + cb*z + cc:
      ca = (q_old(1)*(m22*m33 - m23*m32) - m12*(q_old(2)*m33 - m23*q_old(3)) &
            + m13*(q_old(2)*m32 - m22*q_old(3)))/det
      cb = (m11*(q_old(2)*m33 - m23*q_old(3)) - q_old(1)*(m21*m33 - m23*m31) &
            + m13*(m21*q_old(3) - q_old(2)*m31))/det
      cc = (m11*(m22*q_old(3) - q_old(2)*m32) - m12*(m21*q_old(3) - q_old(2)*m31) &
            + q_old(1)*(m21*m32 - m22*m31))/det
      edge_val = ca*z1*z1 + cb*z1 + cc   ! evaluate at 1|2 interface
      q_R(1) = edge_val
      q_L(2) = edge_val

      ! Surface side: edge between layers nz-1 and nz — quadratic through
      ! layers nz-2, nz-1, nz, evaluated at the (nz-1)|nz interface.
      h1 = dz_old(nz - 2)
      h2 = dz_old(nz - 1)
      h3 = dz_old(nz)
      ! H3 boundary singularity guard (mirror of the interior floor)
      if (h1 <= 0.0_wp .or. h2 <= 0.0_wp .or. h3 <= 0.0_wp) then
         hf = H_MIN_FRAC*max(H_NEGLECT, dz_old(nz - 2) + dz_old(nz - 1) + dz_old(nz))
         h1 = max(h1, hf)
         h2 = max(h2, hf)
         h3 = max(h3, hf)
      end if
      z1 = h1
      z2 = h1 + h2
      z3 = h1 + h2 + h3
      m11 = (z1*z1*z1)/(3.0_wp*h1)
      m12 = (z1*z1)/(2.0_wp*h1)
      m13 = 1.0_wp
      m21 = (z2*z2*z2 - z1*z1*z1)/(3.0_wp*h2)
      m22 = (z2*z2 - z1*z1)/(2.0_wp*h2)
      m23 = 1.0_wp
      m31 = (z3*z3*z3 - z2*z2*z2)/(3.0_wp*h3)
      m32 = (z3*z3 - z2*z2)/(2.0_wp*h3)
      m33 = 1.0_wp
      det = m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31)
      ca = (q_old(nz - 2)*(m22*m33 - m23*m32) - m12*(q_old(nz - 1)*m33 - m23*q_old(nz)) &
            + m13*(q_old(nz - 1)*m32 - m22*q_old(nz)))/det
      cb = (m11*(q_old(nz - 1)*m33 - m23*q_old(nz)) - q_old(nz - 2)*(m21*m33 - m23*m31) &
            + m13*(m21*q_old(nz) - q_old(nz - 1)*m31))/det
      cc = (m11*(m22*q_old(nz) - q_old(nz - 1)*m32) - m12*(m21*q_old(nz) - q_old(nz - 1)*m31) &
            + q_old(nz - 2)*(m21*m32 - m22*m31))/det
      edge_val = ca*z2*z2 + cb*z2 + cc   ! evaluate at (nz-1)|nz interface
      q_R(nz - 1) = edge_val
      q_L(nz) = edge_val

      ! Interior edges: full non-uniform H4 stencil (4 cells k-1..k+2).
      ! q_R(k) is the edge between layer k and k+1; stencil thicknesses
      ! h0=dz(k-1), h1=dz(k), h2=dz(k+1), h3=dz(k+2).
      do k = 2, nz - 2
         h0 = dz_old(k - 1)
         h1 = dz_old(k)
         h2 = dz_old(k + 1)
         h3 = dz_old(k + 2)

         ! Conditional singularity guard: only floor when a consecutive
         ! thickness pair sums to ~0 (vanishing layers under ZSTAR_FULL).
         h_sum = h0 + h1 + h2 + h3
         if (h0 + h1 <= 0.0_wp .or. h1 + h2 <= 0.0_wp .or. h2 + h3 <= 0.0_wp) then
            hf = H_MIN_FRAC*max(H_NEGLECT, h_sum)
            h0 = max(h0, hf)
            h1 = max(h1, hf)
            h2 = max(h2, hf)
            h3 = max(h3, hf)
         end if

         h01 = h0 + h1
         h12 = h1 + h2
         h23 = h2 + h3
         h012 = h0 + h1 + h2
         h123 = h1 + h2 + h3
         h0123 = h0 + h1 + h2 + h3

         f1 = h01*h23/h12
         f2 = h2*q_old(k) + h1*q_old(k + 1)
         f3 = 1.0_wp/h012 + 1.0_wp/h123
         et1 = f1*f2*f3
         et2 = (h2*h23/(h012*h01))*((h0 + 2.0_wp*h1)*q_old(k) - h1*q_old(k - 1))
         et3 = (h1*h01/(h123*h23))*((2.0_wp*h2 + h3)*q_old(k + 1) - h2*q_old(k + 2))

         q_R(k) = (et1 + et2 + et3)/h0123
         q_L(k + 1) = q_R(k)
      end do

      ! ---- Step 2: Colella-Woodward monotonicity limiting (verbatim) ----
      do k = 1, nz
         q_min = q_old(k)
         q_max = q_old(k)
         if (k > 1) then
            q_min = min(q_min, q_old(k - 1))
            q_max = max(q_max, q_old(k - 1))
         end if
         if (k < nz) then
            q_min = min(q_min, q_old(k + 1))
            q_max = max(q_max, q_old(k + 1))
         end if

         ! Clip edges to local bounds
         q_L(k) = max(q_min, min(q_max, q_L(k)))
         q_R(k) = max(q_min, min(q_max, q_R(k)))

         ! CW monotonicity: if the cell is a local extremum, flatten
         dq = q_R(k) - q_L(k)
         dq_l = q_old(k) - q_L(k)
         dq_r = q_R(k) - q_old(k)
         if (dq_l*dq_r <= 0.0_wp) then
            ! Local extremum: flatten to PCM
            q_L(k) = q_old(k)
            q_R(k) = q_old(k)
         else
            ! Check if parabola overshoots
            q6(k) = 6.0_wp*q_old(k) - 3.0_wp*(q_L(k) + q_R(k))
            if (abs(q6(k)) > abs(dq)) then
               if (q6(k)*dq > 0.0_wp) then
                  ! Overshoot near left edge: adjust q_L
                  q_L(k) = 3.0_wp*q_old(k) - 2.0_wp*q_R(k)
               else
                  ! Overshoot near right edge: adjust q_R
                  q_R(k) = 3.0_wp*q_old(k) - 2.0_wp*q_L(k)
               end if
            end if
         end if

         ! Recompute q6 after limiting
         q6(k) = 6.0_wp*q_old(k) - 3.0_wp*(q_L(k) + q_R(k))
      end do

      ! ---- Step 3: Integrate parabolic reconstruction over new layers (verbatim) ----
      ko_start = 1
      do k = 1, nz
         if (dz_new(k) <= 0.0_wp) then
            q_new(k) = 0.0_wp
            cycle
         end if

         integral = 0.0_wp
         do ko = ko_start, nz
            z_lo = max(z_new(k - 1), z_old(ko - 1))
            z_hi = min(z_new(k), z_old(ko))
            overlap = z_hi - z_lo

            if (overlap <= 0.0_wp) then
               if (z_old(ko) > z_new(k)) exit
               cycle
            end if

            if (dz_old(ko) > 0.0_wp) then
               xi_lo = (z_lo - z_old(ko - 1))/dz_old(ko)
               xi_hi = (z_hi - z_old(ko - 1))/dz_old(ko)

               ! Parabolic integral
               integral = integral + dz_old(ko)*( &
                          (xi_hi - xi_lo)*q_L(ko) &
                          + 0.5_wp*(xi_hi*xi_hi - xi_lo*xi_lo)*(q_R(ko) - q_L(ko) + q6(ko)) &
                          - (xi_hi*xi_hi*xi_hi - xi_lo*xi_lo*xi_lo)*q6(ko)/3.0_wp)
            else
               integral = integral + q_old(ko)*overlap
            end if

            if (z_old(ko) <= z_new(k)) ko_start = ko
         end do

         q_new(k) = integral/dz_new(k)
      end do
   end subroutine remap_column_ppm_h4

   pure subroutine pqm_solve_diag_dominant(n, al, ac, au, r, x)
      !$acc routine seq
      !! Diagonally-dominant tridiagonal solve; central diagonal supplied as the
      !! OFFSET `ac` from `al + au` (full pivot = ac + al + au). Never divides by
      !! zero for positive-definite ac, al, au (White & Adcroft 2008).
      integer, intent(in) :: n
         !! Number of unknowns (= number of edges = nz+1)
      real(wp), intent(in) :: al(n)
         !! Lower diagonal (al(1) unused)
      real(wp), intent(in) :: ac(n)
         !! Central-diagonal OFFSET from al+au (full diagonal = ac+al+au)
      real(wp), intent(in) :: au(n)
         !! Upper diagonal (au(n) unused)
      real(wp), intent(in) :: r(n)
         !! Right-hand side
      real(wp), intent(out) :: x(n)
         !! Solution vector

      real(wp) :: c1(NZ_STACK_MAX + 1)
      real(wp) :: d1, i_pivot, denom_t1
      integer :: k

      i_pivot = 1.0_wp/(ac(1) + au(1))
      d1 = ac(1)*i_pivot
      c1(1) = au(1)*i_pivot
      x(1) = r(1)*i_pivot
      do k = 2, n - 1
         denom_t1 = ac(k) + d1*al(k)
         i_pivot = 1.0_wp/(denom_t1 + au(k))
         d1 = denom_t1*i_pivot
         c1(k) = au(k)*i_pivot
         x(k) = (r(k) - al(k)*x(k - 1))*i_pivot
      end do
      i_pivot = 1.0_wp/(ac(n) + d1*al(n))
      x(n) = (r(n) - al(n)*x(n - 1))*i_pivot
      do k = n - 1, 1, -1
         x(k) = x(k) - c1(k)*x(k + 1)
      end do
   end subroutine pqm_solve_diag_dominant

   pure subroutine pqm_end_value_h4(dz, u, csys)
      !$acc routine seq
      !! One-sided 4th-order polynomial fit of the cell averages `u` to the
      !! four boundary layers `dz` (thicknesses, must be positive), returning
      !! the four coefficients `csys` of the fit (White & Adcroft 2008,
      !! appendix; roundoff-safe closed form).  `csys(1)` is the edge VALUE at
      !! the boundary interface and `csys(2)` is the edge SLOPE there.
      real(wp), intent(in) :: dz(4)
         !! Thicknesses of the 4 boundary layers, starting at the edge
      real(wp), intent(in) :: u(4)
         !! Cell averages of the 4 boundary layers, starting at the edge
      real(wp), intent(out) :: csys(4)
         !! Coefficients of the 4th-order fit polynomial in z

      real(wp) :: wt(3, 4)
      real(wp) :: h1, h2, h3, h4
      real(wp) :: h12, h23, h34, h123, h234, h1234
      real(wp) :: i_h12, i_h23, i_h34, i_h123, i_h234, i_h1234
      real(wp) :: i_denom, i_denb3
      real(wp) :: du1, du2, du3

      h1 = dz(1)
      h2 = dz(2)
      h3 = dz(3)
      h4 = dz(4)
      ! Bound the thickness ratios so property differences at the level of
      ! roundoff are not amplified to order one.
      if ((h2 + h3) < PQM_MIN_FRAC*h1) h3 = PQM_MIN_FRAC*h1 - h2
      if ((h3 + h4) < PQM_MIN_FRAC*h1) h4 = PQM_MIN_FRAC*h1 - h3

      h12 = h1 + h2
      h23 = h2 + h3
      h34 = h3 + h4
      h123 = h12 + h3
      h234 = h2 + h34
      h1234 = h12 + h34
      ! Three reciprocals from a single division each, for efficiency.
      i_denb3 = 1.0_wp/(h123*h12*h23)
      i_h12 = (h123*h23)*i_denb3
      i_h23 = (h12*h123)*i_denb3
      i_h123 = (h12*h23)*i_denb3
      i_denom = 1.0_wp/(h1234*(h234*h34))
      i_h34 = (h1234*h234)*i_denom
      i_h234 = (h1234*h34)*i_denom
      i_h1234 = (h234*h34)*i_denom

      wt(1, 1) = -h1*(i_h1234 + i_h123 + i_h12)
      wt(2, 1) = h1*h12*(i_h234*i_h1234 + i_h23*(i_h234 + i_h123))
      wt(3, 1) = -h1*h12*h123*i_denom

      wt(1, 2) = 2.0_wp*(i_h12*(1.0_wp + (h1 + h12)*(i_h1234 + i_h123)) + h1*i_h1234*i_h123)
      wt(2, 2) = -2.0_wp*((h1*h12*i_h1234)*(i_h23*(i_h234 + i_h123)) + &
                          (h1 + h12)*(i_h1234*i_h234 + i_h23*(i_h234 + i_h123)))
      wt(3, 2) = 2.0_wp*((h1 + h12)*h123 + h1*h12)*i_denom

      wt(1, 3) = -3.0_wp*i_h12*i_h123*(1.0_wp + i_h1234*((h1 + h12) + h123))
      wt(2, 3) = 3.0_wp*i_h23*(i_h123 + i_h1234*((h1 + h12) + h123)*(i_h123 + i_h234))
      wt(3, 3) = -3.0_wp*((h1 + h12) + h123)*i_denom

      wt(1, 4) = 4.0_wp*i_h1234*i_h123*i_h12
      wt(2, 4) = -4.0_wp*i_h1234*(i_h23*(i_h123 + i_h234))
      wt(3, 4) = 4.0_wp*i_denom

      du1 = u(2) - u(1)
      du2 = u(3) - u(2)
      du3 = u(4) - u(3)
      csys(1) = ((u(1) + (wt(1, 1)*du1)) + (wt(2, 1)*du2)) + (wt(3, 1)*du3)
      csys(2) = ((wt(1, 2)*du1) + (wt(2, 2)*du2)) + (wt(3, 2)*du3)
      csys(3) = ((wt(1, 3)*du1) + (wt(2, 3)*du2)) + (wt(3, 3)*du3)
      csys(4) = ((wt(1, 4)*du1) + (wt(2, 4)*du2)) + (wt(3, 4)*du3)
   end subroutine pqm_end_value_h4

   pure subroutine remap_column_pqm(nz, dz_old, dz_new, q_old, q_new)
      !$acc routine seq
      !! Piecewise-quartic (PQM_IH4IH3) conservative remap (White & Adcroft 2008).
      !! Implicit-h4 edge VALUES + implicit-h3 edge SLOPES (each a
      !! diagonally-dominant tridiagonal solve with one-sided 4-cell boundary
      !! closure), per-cell quartic, W&A monotonicity limiter, conservative
      !! quartic overlap integral. Cuts diapycnal mixing per remap vs PPM/PPM_H4
      !! (Ilicak et al. 2012).
      !! Per cell k, xi in [0,1]: q_hat(xi) = a + b*xi + c*xi^2 + d*xi^3 + e*xi^4.
      !! Boundary cells reconstruct as PCM. nz < 5 falls back to REMAP_PPM
      !! (W&A boundary closure needs >= 4 cells).
      integer, intent(in) :: nz
      real(wp), intent(in) :: dz_old(nz)
         !! Old layer thicknesses (must sum to same total as dz_new)
      real(wp), intent(in) :: dz_new(nz)
         !! New layer thicknesses
      real(wp), intent(in) :: q_old(nz)
         !! Old cell-average scalar values
      real(wp), intent(out) :: q_new(nz)
         !! New cell-average scalar values (conservative)

      real(wp) :: z_old(0:NZ_STACK_MAX), z_new(0:NZ_STACK_MAX)
      ! Edge values / slopes, two per cell: index 1 = left, 2 = right.
      real(wp) :: ev_l(NZ_STACK_MAX), ev_r(NZ_STACK_MAX)
      real(wp) :: es_l(NZ_STACK_MAX), es_r(NZ_STACK_MAX)
      ! Per-cell quartic coefficients a..e.
      real(wp) :: pa(NZ_STACK_MAX), pb(NZ_STACK_MAX), pc(NZ_STACK_MAX)
      real(wp) :: pd(NZ_STACK_MAX), pe(NZ_STACK_MAX)
      ! Tridiagonal workspace (N+1 edges).
      real(wp) :: tri_l(NZ_STACK_MAX + 1), tri_c(NZ_STACK_MAX + 1)
      real(wp) :: tri_u(NZ_STACK_MAX + 1), tri_b(NZ_STACK_MAX + 1)
      real(wp) :: tri_x(NZ_STACK_MAX + 1)
      real(wp) :: dz4(4), u4(4), csys(4)
      real(wp) :: h0, h1, i_h2, alpha, beta, abmix, aco, bco
      real(wp) :: i_h, h0h1, h0_2, h1_2, h0_3, h1_3, i_d
      real(wp) :: z_lo, z_hi, overlap, integral, xi_lo, xi_hi
      real(wp) :: h_c, u0_l, u0_r, u1_l, u1_r, u_l, u_c, u_r, h_l, h_r
      real(wp) :: sigma_l, sigma_c, sigma_r, slope, slope_x_h
      real(wp) :: u0_avg
      real(wp) :: a, b, cco, dco, eco, alpha1, alpha2, alpha3
      real(wp) :: rho, sqrt_rho, x1, x2, grad1, grad2
      integer :: k, ko, ko_start, np1, inflexion_l, inflexion_r

      ! Trivial / degenerate cases — fall back to lower-order safe paths.
      if (nz == 1) then
         q_new(1) = q_old(1)
         return
      end if
      if (nz < 5) then
         call remap_column_ppm(nz, dz_old, dz_new, q_old, q_new)
         return
      end if

      np1 = nz + 1

      ! Interface positions.
      z_old(0) = 0.0_wp
      z_new(0) = 0.0_wp
      do k = 1, nz
         z_old(k) = z_old(k - 1) + dz_old(k)
         z_new(k) = z_new(k - 1) + dz_new(k)
      end do

      ! ---- Step 1: implicit-h4 edge VALUES (tridiagonal, N+1 edges) ----
      ! Interior edges i+1 (between cells i and i+1) — roundoff-safe stencil
      ! al*x(i)+x(i+1)+be*x(i+2) with diagonal OFFSET tri_c = 2*abmix.
      do k = 1, nz - 1
         h0 = max(dz_old(k), H_NEGLECT)
         h1 = max(dz_old(k + 1), H_NEGLECT)
         if (abs(h0) < H_RATIO_FLOOR*abs(h1)) h0 = H_RATIO_FLOOR*h1
         if (abs(h1) < H_RATIO_FLOOR*abs(h0)) h1 = H_RATIO_FLOOR*h0
         i_h2 = 1.0_wp/((h0 + h1)**2)
         alpha = (h1*h1)*i_h2
         beta = (h0*h0)*i_h2
         abmix = (h0*h1)*i_h2
         aco = 2.0_wp*alpha*(alpha + 2.0_wp*beta + 3.0_wp*abmix)
         bco = 2.0_wp*beta*(beta + 2.0_wp*alpha + 3.0_wp*abmix)
         tri_l(k + 1) = alpha
         tri_c(k + 1) = 2.0_wp*abmix
         tri_u(k + 1) = beta
         tri_b(k + 1) = aco*q_old(k) + bco*q_old(k + 1)
      end do
      ! Top boundary (edge 1): exact 4-cell one-sided fit, value = csys(1).
      do k = 1, 4
         dz4(k) = max(H_NEGLECT, dz_old(k))
         u4(k) = q_old(k)
      end do
      call pqm_end_value_h4(dz4, u4, csys)
      tri_b(1) = csys(1)
      tri_c(1) = 1.0_wp
      tri_u(1) = 0.0_wp
      ! Bottom boundary (edge N+1): layers REVERSED, value = csys(1).
      do k = 1, 4
         dz4(k) = max(H_NEGLECT, dz_old(nz + 1 - k))
         u4(k) = q_old(nz + 1 - k)
      end do
      call pqm_end_value_h4(dz4, u4, csys)
      tri_b(np1) = csys(1)
      tri_c(np1) = 1.0_wp
      tri_l(np1) = 0.0_wp

      call pqm_solve_diag_dominant(np1, tri_l, tri_c, tri_u, tri_b, tri_x)

      ! Scatter the N+1 edge values into per-cell (left,right) pairs.
      ev_l(1) = tri_x(1)
      do k = 2, nz
         ev_l(k) = tri_x(k)
         ev_r(k - 1) = tri_x(k)
      end do
      ev_r(nz) = tri_x(np1)

      ! ---- Step 2: implicit-h3 edge SLOPES (tridiagonal, N+1 edges) ----
      ! Nondimensionalised stencil; diagonal OFFSET tri_c (NOT 1.0 — the
      ! diagonal-dominant solver adds tri_l+tri_u to form the pivot).
      do k = 1, nz - 1
         h0 = max(dz_old(k), H_NEGLECT)
         h1 = max(dz_old(k + 1), H_NEGLECT)
         i_h = 1.0_wp/(h0 + h1)
         h0 = h0*i_h
         h1 = h1*i_h
         h0h1 = h0*h1
         h0_2 = h0*h0
         h1_2 = h1*h1
         h0_3 = h0_2*h0
         h1_3 = h1_2*h1
         i_d = 1.0_wp/(4.0_wp*h0h1*(h0 + h1) + h1_3 + h0_3)
         tri_l(k + 1) = (h1*((h0_2 + h0h1) - h1_2))*i_d
         tri_c(k + 1) = 2.0_wp*((h0_2 + h1_2)*(h0 + h1))*i_d
         tri_u(k + 1) = (h0*((h1_2 + h0h1) - h0_2))*i_d
         tri_b(k + 1) = 12.0_wp*(h0h1*i_d)*((q_old(k + 1) - q_old(k))*i_h)
      end do
      ! Top boundary slope = csys(2) of the 4-cell fit.
      do k = 1, 4
         dz4(k) = max(H_NEGLECT, dz_old(k))
         u4(k) = q_old(k)
      end do
      call pqm_end_value_h4(dz4, u4, csys)
      tri_b(1) = csys(2)
      tri_c(1) = 1.0_wp
      tri_u(1) = 0.0_wp
      ! Bottom boundary slope = -csys(2) (layers reversed → sign flip).
      do k = 1, 4
         dz4(k) = max(H_NEGLECT, dz_old(nz + 1 - k))
         u4(k) = q_old(nz + 1 - k)
      end do
      call pqm_end_value_h4(dz4, u4, csys)
      tri_b(np1) = -csys(2)
      tri_c(np1) = 1.0_wp
      tri_l(np1) = 0.0_wp

      call pqm_solve_diag_dominant(np1, tri_l, tri_c, tri_u, tri_b, tri_x)

      es_l(1) = tri_x(1)
      do k = 2, nz
         es_l(k) = tri_x(k)
         es_r(k - 1) = tri_x(k)
      end do
      es_r(nz) = tri_x(np1)

      ! ---- Step 3: PQM limiter (White & Adcroft 2008) ----
      ! 3a. bound_edge_values: van-Leer edge limiting + neighbour-mean clamp.
      do k = 1, nz
         u_l = q_old(max(1, k - 1))
         u_c = q_old(k)
         u_r = q_old(min(k + 1, nz))
         h_l = dz_old(max(1, k - 1))
         h_c = dz_old(k)
         h_r = dz_old(min(k + 1, nz))
         slope_x_h = 0.0_wp
         if (((h_l + h_r) + 2.0_wp*h_c) > 0.0_wp) then
            sigma_l = (u_c - u_l)
            sigma_c = (u_r - u_l)*(h_c/((h_l + h_r) + 2.0_wp*h_c))
            sigma_r = (u_r - u_c)
            if ((sigma_l*sigma_r) > 0.0_wp) then
               slope_x_h = sign(min(abs(sigma_l), abs(sigma_c), abs(sigma_r)), sigma_c)
            end if
         end if
         if ((u_l - ev_l(k))*(ev_l(k) - u_c) < 0.0_wp) then
            ev_l(k) = u_c - sign(min(abs(slope_x_h), abs(ev_l(k) - u_c)), slope_x_h)
         end if
         if ((u_r - ev_r(k))*(ev_r(k) - u_c) < 0.0_wp) then
            ev_r(k) = u_c + sign(min(abs(slope_x_h), abs(ev_r(k) - u_c)), slope_x_h)
         end if
         ev_l(k) = max(min(ev_l(k), max(u_l, u_c)), min(u_l, u_c))
         ev_r(k) = max(min(ev_r(k), max(u_r, u_c)), min(u_r, u_c))
      end do

      ! 3b. check_discontinuous_edge_values: average non-monotonic collocated
      ! edges.  Sweep low→high; ev_km1_r holds the (possibly updated) right
      ! edge of cell k so the pair update stays consistent.
      do k = 1, nz - 1
         if ((ev_l(k + 1) - ev_r(k))*(q_old(k + 1) - q_old(k)) < 0.0_wp) then
            u0_avg = 0.5_wp*(ev_r(k) + ev_l(k + 1))
            u0_avg = max(min(u0_avg, max(q_old(k), q_old(k + 1))), &
                         min(q_old(k), q_old(k + 1)))
            ev_r(k) = u0_avg
            ev_l(k + 1) = u0_avg
         end if
      end do

      ! 3c. interior cells: PLM-slope consistency, extremum flatten, quartic
      ! curvature / inflexion test, collapse + post-collapse resets.
      do k = 2, nz - 1
         inflexion_l = 0
         inflexion_r = 0
         u0_l = ev_l(k)
         u0_r = ev_r(k)
         u1_l = es_l(k)
         u1_r = es_r(k)
         h_l = dz_old(k - 1)
         h_c = dz_old(k)
         h_r = dz_old(k + 1)
         u_l = q_old(k - 1)
         u_c = q_old(k)
         u_r = q_old(k + 1)

         sigma_l = 2.0_wp*(u_c - u_l)/(h_c + H_NEGLECT)
         sigma_c = 2.0_wp*(u_r - u_l)/(h_l + 2.0_wp*h_c + h_r + H_NEGLECT)
         sigma_r = 2.0_wp*(u_r - u_c)/(h_c + H_NEGLECT)
         if ((sigma_l*sigma_r) > 0.0_wp) then
            slope = sign(min(abs(sigma_l), abs(sigma_c), abs(sigma_r)), sigma_c)
         else
            slope = 0.0_wp
         end if

         if (u1_l*slope <= 0.0_wp) u1_l = slope
         if (u1_r*slope <= 0.0_wp) u1_r = slope

         if ((u0_r - u_c)*(u_c - u0_l) <= 0.0_wp) then
            u0_l = u_c
            u0_r = u_c
            u1_l = 0.0_wp
            u1_r = 0.0_wp
            inflexion_l = -1
            inflexion_r = -1
         end if

         if ((inflexion_l == 0) .and. (inflexion_r == 0)) then
            a = u0_l
            b = h_c*u1_l
            cco = 30.0_wp*u_c - 12.0_wp*u0_r - 18.0_wp*u0_l + 1.5_wp*h_c*(u1_r - 3.0_wp*u1_l)
            dco = -60.0_wp*u_c + h_c*(6.0_wp*u1_l - 4.0_wp*u1_r) + 28.0_wp*u0_r + 32.0_wp*u0_l
            eco = 30.0_wp*u_c + 2.5_wp*h_c*(u1_r - u1_l) - 15.0_wp*(u0_l + u0_r)

            alpha1 = 6.0_wp*eco
            alpha2 = 3.0_wp*dco
            alpha3 = cco
            rho = alpha2*alpha2 - 4.0_wp*alpha1*alpha3

            if ((alpha1 /= 0.0_wp) .and. (rho >= 0.0_wp)) then
               sqrt_rho = sqrt(rho)
               x1 = 0.5_wp*(-alpha2 - sqrt_rho)/alpha1
               x2 = 0.5_wp*(-alpha2 + sqrt_rho)/alpha1
               if ((x1 >= 0.0_wp) .and. (x1 <= 1.0_wp) .and. &
                   (x2 >= 0.0_wp) .and. (x2 <= 1.0_wp)) then
                  grad1 = 4.0_wp*eco*(x1**3) + 3.0_wp*dco*(x1**2) + 2.0_wp*cco*x1 + b
                  grad2 = 4.0_wp*eco*(x2**3) + 3.0_wp*dco*(x2**2) + 2.0_wp*cco*x2 + b
                  if ((grad1*slope < 0.0_wp) .or. (grad2*slope < 0.0_wp)) then
                     if (abs(sigma_l) < abs(sigma_r)) then
                        inflexion_l = 1
                     else
                        inflexion_r = 1
                     end if
                  end if
               else if ((x1 >= 0.0_wp) .and. (x1 <= 1.0_wp)) then
                  grad1 = 4.0_wp*eco*(x1**3) + 3.0_wp*dco*(x1**2) + 2.0_wp*cco*x1 + b
                  if (grad1*slope < 0.0_wp) then
                     if (abs(sigma_l) < abs(sigma_r)) then
                        inflexion_l = 1
                     else
                        inflexion_r = 1
                     end if
                  end if
               else if ((x2 >= 0.0_wp) .and. (x2 <= 1.0_wp)) then
                  grad2 = 4.0_wp*eco*(x2**3) + 3.0_wp*dco*(x2**2) + 2.0_wp*cco*x2 + b
                  if (grad2*slope < 0.0_wp) then
                     if (abs(sigma_l) < abs(sigma_r)) then
                        inflexion_l = 1
                     else
                        inflexion_r = 1
                     end if
                  end if
               end if
            end if

            if ((alpha1 == 0.0_wp) .and. (alpha2 /= 0.0_wp)) then
               x1 = -alpha3/alpha2
               if ((x1 >= 0.0_wp) .and. (x1 <= 1.0_wp)) then
                  grad1 = 4.0_wp*eco*(x1**3) + 3.0_wp*dco*(x1**2) + 2.0_wp*cco*x1 + b
                  if (grad1*slope < 0.0_wp) then
                     if (abs(sigma_l) < abs(sigma_r)) then
                        inflexion_l = 1
                     else
                        inflexion_r = 1
                     end if
                  end if
               end if
            end if
         end if

         if (inflexion_l == 1) then
            ! Collapse both inflexion points onto the LEFT edge.
            u1_l = (10.0_wp*u_c - 2.0_wp*u0_r - 8.0_wp*u0_l)/(3.0_wp*h_c + H_NEGLECT)
            u1_r = (-10.0_wp*u_c + 6.0_wp*u0_r + 4.0_wp*u0_l)/(h_c + H_NEGLECT)
            if (u1_l*slope < 0.0_wp) then
               u1_l = 0.0_wp
               u0_r = 5.0_wp*u_c - 4.0_wp*u0_l
               u1_r = 20.0_wp*(u_c - u0_l)/(h_c + H_NEGLECT)
            else if (u1_r*slope < 0.0_wp) then
               u1_r = 0.0_wp
               u0_l = (5.0_wp*u_c - 3.0_wp*u0_r)/2.0_wp
               u1_l = 10.0_wp*(-u_c + u0_r)/(3.0_wp*h_c + H_NEGLECT)
            end if
         else if (inflexion_r == 1) then
            ! Collapse both inflexion points onto the RIGHT edge.
            u1_r = (-10.0_wp*u_c + 8.0_wp*u0_r + 2.0_wp*u0_l)/(3.0_wp*h_c + H_NEGLECT)
            u1_l = (10.0_wp*u_c - 4.0_wp*u0_r - 6.0_wp*u0_l)/(h_c + H_NEGLECT)
            if (u1_l*slope < 0.0_wp) then
               u1_l = 0.0_wp
               u0_r = (5.0_wp*u_c - 3.0_wp*u0_l)/2.0_wp
               u1_r = 10.0_wp*(u_c - u0_l)/(3.0_wp*h_c + H_NEGLECT)
            else if (u1_r*slope < 0.0_wp) then
               u1_r = 0.0_wp
               u0_l = 5.0_wp*u_c - 4.0_wp*u0_r
               u1_l = 20.0_wp*(-u_c + u0_r)/(h_c + H_NEGLECT)
            end if
         end if

         ev_l(k) = u0_l
         ev_r(k) = u0_r
         es_l(k) = u1_l
         es_r(k) = u1_r
      end do

      ! Boundary cells: PCM (constant reconstruction).
      ev_l(1) = q_old(1)
      ev_r(1) = q_old(1)
      es_l(1) = 0.0_wp
      es_r(1) = 0.0_wp
      ev_l(nz) = q_old(nz)
      ev_r(nz) = q_old(nz)
      es_l(nz) = 0.0_wp
      es_r(nz) = 0.0_wp

      ! ---- Step 4: per-cell quartic coefficients (xi in [0,1]) ----
      do k = 1, nz
         h_c = dz_old(k)
         u0_l = ev_l(k)
         u0_r = ev_r(k)
         u1_l = es_l(k)
         u1_r = es_r(k)
         u_c = q_old(k)
         pa(k) = u0_l
         pb(k) = h_c*u1_l
         pc(k) = 30.0_wp*u_c - 12.0_wp*u0_r - 18.0_wp*u0_l + 1.5_wp*h_c*(u1_r - 3.0_wp*u1_l)
         pd(k) = -60.0_wp*u_c + h_c*(6.0_wp*u1_l - 4.0_wp*u1_r) + 28.0_wp*u0_r + 32.0_wp*u0_l
         pe(k) = 30.0_wp*u_c + 2.5_wp*h_c*(u1_r - u1_l) - 15.0_wp*(u0_l + u0_r)
      end do

      ! ---- Step 5: conservative quartic overlap integration ----
      ! For old cell ko, the mean of the quartic over [xi_lo,xi_hi] is the
      ! analytic average_value_ppoly form:
      !   a + b*<xi> + c*<xi^2> + d*<xi^3> + e*<xi^4>
      ! with the symmetric power-mean expressions over the sub-interval.
      ! Multiplying the mean by the overlap thickness gives the contribution.
      ko_start = 1
      do k = 1, nz
         if (dz_new(k) <= 0.0_wp) then
            q_new(k) = 0.0_wp
            cycle
         end if
         integral = 0.0_wp
         do ko = ko_start, nz
            z_lo = max(z_new(k - 1), z_old(ko - 1))
            z_hi = min(z_new(k), z_old(ko))
            overlap = z_hi - z_lo
            if (overlap <= 0.0_wp) then
               if (z_old(ko) > z_new(k)) exit
               cycle
            end if
            if (dz_old(ko) > 0.0_wp) then
               xi_lo = (z_lo - z_old(ko - 1))/dz_old(ko)
               xi_hi = (z_hi - z_old(ko - 1))/dz_old(ko)
               integral = integral + overlap*( &
                          pa(ko) &
                          + pb(ko)*0.5_wp*(xi_lo + xi_hi) &
                          + pc(ko)*(1.0_wp/3.0_wp)*(xi_lo*xi_lo + xi_hi*xi_hi + xi_lo*xi_hi) &
                          + pd(ko)*0.25_wp*((xi_lo*xi_lo + xi_hi*xi_hi)*(xi_lo + xi_hi)) &
                          + pe(ko)*0.2_wp*((xi_hi**3 + xi_lo**3)*(xi_lo + xi_hi) &
                                           + xi_lo*xi_lo*xi_hi*xi_hi))
            else
               integral = integral + q_old(ko)*overlap
            end if
            if (z_old(ko) <= z_new(k)) ko_start = ko
         end do
         q_new(k) = integral/dz_new(k)
      end do
   end subroutine remap_column_pqm

end module rki_remap_column
