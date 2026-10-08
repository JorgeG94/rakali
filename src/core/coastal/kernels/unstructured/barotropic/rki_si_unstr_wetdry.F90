!! Casulli-style wetting/drying machinery for the unstructured SI surface.
module rki_si_unstr_wetdry
   !! Casulli (2009) subgrid wetting/drying for the unstructured
   !! semi-implicit free-surface mode — ported per-cell (`(n_sub, ncells)` sub-bottom
   !! layout) exactly as the prototype confirmed mesh-agnostic
   !! (`si_unstr_proto.py` study 2: V/a/Newton + the edge sill rule ran
   !! UNCHANGED on triangles).  Replaces the linear solver's constant
   !! `1.0` cell-mass term by the exact wetted-volume function per unit
   !! cell area (depth units):
   !!
   !!   v_c(eta) = (1/n_sub) sum_s max(eta - b_c^(s), 0)   (mean wet depth)
   !!   a_c(eta) = v'_c(eta) in [0,1]                      (wet fraction)
   !!
   !! driving a PLAIN Newton outer loop (V convex => monotone finite
   !! convergence) around the SAME area-weighted TPFA CG, with the
   !! operator's mass term replaced per cell by `a(eta) + floor` (the
   !! `mass_diag` argument already threaded through
   !! `si_unstr_apply_helmholtz` / `si_unstr_cg_solve`).
   !!
   !! FROZEN edge rule (Casulli, rule R1) — the GATED-MEAN conveyance:
   !!   Hf_e = 0.5*(max(eta_L-b_L,0) + max(eta_R-b_R,0))
   !!          gated to 0 unless  max(eta_L,eta_R) > max(b_L,b_R)
   !! built ONCE from eta^n and held over the step, shared by the
   !! operator, back-substitution and RHS (R2).  Where both cells are wet
   !! this reduces EXACTLY to the linear edge depth 0.5*(h_L+h_R) — the
   !! all-wet self-consistency invariant.  Dead edges (higher surface
   !! below the sill) become walls; their U^n edge transports are zeroed
   !! at step start (R3 — structural on this backend: the RHS consumes
   !! per-edge transports directly, so zeroing u_star/u_n_edge on Hf=0
   !! edges closes the dead-edge flux exactly; no cell-momentum
   !! reconstruction masking is needed, unlike the structured path).
   !! Dry-cell eta is NOT clamped at the bed — Newton lets it settle
   !! below (V = 0 regardless).
   !!
   !! References: Casulli (2009) IJNMF 60:391-408; Casulli & Zanolli
   !! (2010) SISC 32(4):2255-2273.
   !!
   !! All kernels explicit-shape `(n_sub|3, ncells)`, device-resident,
   !! cell-index-last per the unstructured convention.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   implicit none
   private

   public :: si_unstr_wetdry_fill_bsub
   public :: si_unstr_wetdry_volume
   public :: si_unstr_wetdry_wetfrac
   public :: si_unstr_wetdry_mass_diag
   public :: si_unstr_wetdry_face_depths
   public :: si_unstr_wetdry_newton_residual
   public :: si_unstr_wetdry_wetcount

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.
   real(wp), parameter :: WETDRY_FLOOR_TOL = 1.0e-14_wp
      !! Threshold below which a cell's wet fraction AND its three edge
      !! stiffness coefficients count as zero => a fully-decoupled dry row
      !! that receives the +1 identity floor (Casulli rule R5): keeps the
      !! Newton Jacobian SPD with O(1) conditioning where a would be 0.

contains

   pure subroutine si_unstr_wetdry_fill_bsub(mesh, nsub, nct, nco, bed, &
                                             bsub_fill, b_sub, b_min)
      !! Fill the per-cell sub-bottom table from the plain bed.  v1
      !! (n_sub=1): b_sub(1,:) = bed, b_min = bed — always, regardless of
      !! mode.  Subgrid (n_sub>=2):
      !!   FLAT (default): all sub-bottoms = bed => BIT-IDENTICAL to n_sub=1.
      !!   LINEAR: `nsub` equal-area sub-bottoms evenly spaced across the
      !!     per-cell relief band `[bed-half, bed+half]`, `half = 0.5*max`
      !!     edge-neighbour |Δbed| (a linear sub-cell bed ramp; Casulli 2009
      !!     subgrid bathymetry, prototype `subgrid_bathy_from_bed`).
      !!     Sub-bottoms ascending; mean == bed (still-water volume
      !!     unchanged); wet fraction staircases 1/nsub per breakpoint so a
      !!     coarse tidal flat drains/floods SMOOTHLY.  Boundary edges (no
      !!     neighbour) contribute 0 to the relief so a shore cell falls back
      !!     toward FLAT.  On flat bed (half=0) => == FLAT.  b_min = the
      !!     deepest (first, sill) sub-bottom.
      use rki_si_common, only: SI_BSUB_FILL_LINEAR
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nsub
         !! Sub-bottoms per cell (>= 1).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: bed(nct)
      integer, intent(in) :: bsub_fill
         !! `SI_BSUB_FILL_*`.  FLAT => flat fill; LINEAR => the ramp.
      real(wp), intent(out) :: b_sub(nsub, nct)
      real(wp), intent(out) :: b_min(nct)

      integer :: c, s, i, kn
      real(wp) :: half, b_lo, span, frac, dmax

      if (bsub_fill /= SI_BSUB_FILL_LINEAR .or. nsub < 2) then
         do concurrent(c=1:nco) local(s)
            do s = 1, nsub
               b_sub(s, c) = bed(c)
            end do
            b_min(c) = bed(c)
         end do
         return
      end if

      do concurrent(c=1:nco) local(s, i, kn, half, b_lo, span, frac, dmax)
         dmax = 0.0_wp
         do i = 1, EDGES_PER_TRI
            kn = mesh%neighbours(i, c)
            if (kn > 0) dmax = max(dmax, abs(bed(kn) - bed(c)))
         end do
         half = 0.5_wp*dmax
         b_lo = bed(c) - half
         span = 2.0_wp*half
         do s = 1, nsub
            frac = (real(s, wp) - 0.5_wp)/real(nsub, wp)
            b_sub(s, c) = b_lo + span*frac
         end do
         b_min(c) = b_lo
      end do
   end subroutine si_unstr_wetdry_fill_bsub

   pure subroutine si_unstr_wetdry_volume(nsub, nct, nco, eta, b_sub, v_cell)
      !! Per-cell wetted mean depth v_c(eta) = (1/n_sub) sum_s max(eta-b^s,0)
      !! (m; per unit cell area).  n_sub=1: v = max(eta-b,0) = the
      !! still-water total depth.
      integer, intent(in) :: nsub
         !! Sub-bottoms per cell.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta(nct)
      real(wp), intent(in) :: b_sub(nsub, nct)
      real(wp), intent(out) :: v_cell(nct)

      integer :: c, s
      real(wp) :: acc, inv_nsub

      inv_nsub = 1.0_wp/real(nsub, wp)
      do concurrent(c=1:nco) local(s, acc)
         acc = 0.0_wp
         do s = 1, nsub
            acc = acc + max(eta(c) - b_sub(s, c), 0.0_wp)
         end do
         v_cell(c) = acc*inv_nsub
      end do
   end subroutine si_unstr_wetdry_volume

   pure subroutine si_unstr_wetdry_wetfrac(nsub, nct, nco, eta, b_sub, a_cell)
      !! Per-cell wet fraction a_c(eta) = v' in [0,1] — the Newton diagonal.
      !! RIGHT-continuous at breakpoints (`eta >= b_sub` counts wet), the
      !! sub-derivative that keeps Newton-on-convex-V monotone.
      integer, intent(in) :: nsub
         !! Sub-bottoms per cell.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta(nct)
      real(wp), intent(in) :: b_sub(nsub, nct)
      real(wp), intent(out) :: a_cell(nct)

      integer :: c, s, cnt
      real(wp) :: inv_nsub

      inv_nsub = 1.0_wp/real(nsub, wp)
      do concurrent(c=1:nco) local(s, cnt)
         cnt = 0
         do s = 1, nsub
            if (eta(c) >= b_sub(s, c)) cnt = cnt + 1
         end do
         a_cell(c) = real(cnt, wp)*inv_nsub
      end do
   end subroutine si_unstr_wetdry_wetfrac

   pure subroutine si_unstr_wetdry_mass_diag(nct, nco, a_cell, Hf_edge, mass_diag)
      !! Assemble the Newton mass diagonal M = a + floor (Casulli R5): the
      !! wet fraction a plus a +1 identity for a fully-decoupled dry row
      !! (a <= tol AND all three edge stiffness coefficients Hf <= tol).
      !! The `+1` matches the linear operator's constant mass exactly at
      !! wet cells where a == 1 and no floor fires.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: a_cell(nct)
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(out) :: mass_diag(nct)

      integer :: c, i
      real(wp) :: stiff_sum, floor_c

      do concurrent(c=1:nco) local(i, stiff_sum, floor_c)
         stiff_sum = 0.0_wp
         do i = 1, EDGES_PER_TRI
            stiff_sum = stiff_sum + Hf_edge(i, c)
         end do
         if (a_cell(c) <= WETDRY_FLOOR_TOL .and. stiff_sum <= WETDRY_FLOOR_TOL) then
            floor_c = 1.0_wp
         else
            floor_c = 0.0_wp
         end if
         mass_diag(c) = a_cell(c) + floor_c
      end do
   end subroutine si_unstr_wetdry_mass_diag

   pure subroutine si_unstr_wetdry_face_depths(mesh, nct, nco, n_bc, eta, &
                                               b_min, bc_type, bc_eta_ext, &
                                               Hf_edge, eta_ghost_edge)
      !! FROZEN wet/dry edge depths — the GATED-MEAN conveyance rule per
      !! edge (Casulli rule R1, prototype-verified mesh-agnostic):
      !!   wet  <=>  max(eta_L, eta_R) > max(b_L, b_R)     (max-eta/sill)
      !!   Hf   = 0.5*(max(eta_L-b_L,0) + max(eta_R-b_R,0))  where wet
      !! Where BOTH cells are wet this is EXACTLY the linear edge depth
      !! 0.5*(h_L+h_R) — an all-wet basin reproduces the linear operator's
      !! Hf bit-for-bit (NOT a raw max-eta or bare sill depth, which would
      !! change the all-wet answer — the load-bearing structured lesson).
      !! Dead edges (higher surface below the higher bed — a wet cell
      !! cannot "see" through a dry one over the sill) are walls (Hf=0).
      !!
      !! Boundary edges: WALL => 0; TIDAL/CLAMPED (Dirichlet) => the same
      !! gated-mean rule against the seeded ghost (ghost bed mirrors the
      !! cell bed; ghost surface = the prescribed eta_bc, also stashed in
      !! `eta_ghost_edge` for the operator/RHS/back-sub).
      use rki_boundary_unstr, only: BC_WALL, BC_TIDAL, BC_CLAMPED
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      integer, intent(in) :: n_bc
         !! BC group table size.
      real(wp), intent(in) :: eta(nct)
         !! eta^n at circumcentres (the FROZEN surface).
      real(wp), intent(in) :: b_min(nct)
         !! Sill bed elevation (deepest sub-bottom) at cells (m).
      integer, intent(in) :: bc_type(n_bc)
      real(wp), intent(in) :: bc_eta_ext(n_bc)
      real(wp), intent(out) :: Hf_edge(EDGES_PER_TRI, nct)
      real(wp), intent(out) :: eta_ghost_edge(EDGES_PER_TRI, nct)

      integer :: c, i, kn, g, bt
      real(wp) :: eta_bc

      do concurrent(c=1:nco) local(i, kn, g, bt, eta_bc)
         do i = 1, EDGES_PER_TRI
            eta_ghost_edge(i, c) = 0.0_wp
            kn = mesh%neighbours(i, c)
            if (kn > 0) then
               Hf_edge(i, c) = edge_depth_wetdry(eta(c), eta(kn), &
                                                 b_min(c), b_min(kn))
            else
               g = mesh%bc_group(i, c)
               bt = BC_WALL
               if (g >= 1 .and. g <= n_bc) bt = bc_type(g)
               if ((bt == BC_TIDAL .or. bt == BC_CLAMPED) .and. g >= 1) then
                  eta_bc = bc_eta_ext(g)
                  Hf_edge(i, c) = edge_depth_wetdry(eta(c), eta_bc, &
                                                    b_min(c), b_min(c))
                  eta_ghost_edge(i, c) = eta_bc
               else
                  Hf_edge(i, c) = 0.0_wp
               end if
            end if
         end do
      end do
   end subroutine si_unstr_wetdry_face_depths

   pure function edge_depth_wetdry(eta_l, eta_r, b_l, b_r) result(hf)
      !! Frozen wet/dry edge conveyance: the arithmetic-mean water depth of
      !! the two adjacent cells, gated to 0 by the max-eta/sill wet test.
      !! `!$acc routine seq` so device do-concurrent bodies can call it.
      !$omp declare target
      real(wp), intent(in) :: eta_l, eta_r, b_l, b_r
      real(wp) :: hf

      if (max(eta_l, eta_r) > max(b_l, b_r)) then
         hf = 0.5_wp*(max(eta_l - b_l, 0.0_wp) + max(eta_r - b_r, 0.0_wp))
      else
         hf = 0.0_wp
      end if
   end function edge_depth_wetdry

   subroutine si_unstr_wetdry_newton_residual(nct, nco, eta, Jeta, mass_diag, &
                                              v_cell, rhs, res, res2_out)
      !! Form the Newton nonlinear residual and its squared interior norm:
      !!   F_c = v_c(eta) + [T eta]_c - rhs_c,
      !! with the stiffness recovered matmul-free from the full mass-diag
      !! operator apply:  T*eta = Jeta - mass_diag*eta.  Returns sum(F^2)
      !! over the owned cells (the Newton convergence norm — plain, not
      !! area-weighted, matching the prototype/structured gate; it is a
      !! convergence heuristic, not the CG inner product).
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta(nct), Jeta(nct), mass_diag(nct)
      real(wp), intent(in) :: v_cell(nct), rhs(nct)
      real(wp), intent(out) :: res(nct)
      real(wp), intent(out) :: res2_out

      integer :: c
      real(wp) :: res2

      do concurrent(c=1:nco)
         res(c) = v_cell(c) + (Jeta(c) - mass_diag(c)*eta(c)) - rhs(c)
      end do

      res2 = 0.0_wp
      do concurrent(c=1:nco) reduce(+:res2)
         res2 = res2 + res(c)*res(c)
      end do
      res2_out = res2
   end subroutine si_unstr_wetdry_newton_residual

   subroutine si_unstr_wetdry_wetcount(nsub, nct, nco, eta, b_sub, count_out)
      !! Total owned wet-sub-bottom count — the active-set fingerprint the
      !! Newton loop watches (frozen set + small residual = converged).
      !! Integer so equality is exact.
      integer, intent(in) :: nsub
         !! Sub-bottoms per cell.
      integer, intent(in) :: nct
         !! Cell array extent (ncells_total).
      integer, intent(in) :: nco
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta(nct)
      real(wp), intent(in) :: b_sub(nsub, nct)
      integer, intent(out) :: count_out

      integer :: c, s, cnt

      cnt = 0
      do concurrent(c=1:nco) reduce(+:cnt)
         do s = 1, nsub
            if (eta(c) >= b_sub(s, c)) cnt = cnt + 1
         end do
      end do
      count_out = cnt
   end subroutine si_unstr_wetdry_wetcount

end module rki_si_unstr_wetdry
