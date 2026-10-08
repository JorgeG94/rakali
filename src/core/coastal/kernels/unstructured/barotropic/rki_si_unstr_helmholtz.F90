!! Matrix-free TPFA Helmholtz operator + RHS for the unstructured SI surface.
module rki_si_unstr_helmholtz
   !! The elliptic (wave-continuity) machinery for the semi-implicit
   !! free-surface mode on the UNSTRUCTURED (triangular) backend — the
   !! two-point flux (TPFA, Casulli & Walters 2000) scheme.  Substituting the semi-implicit edge momentum into
   !! depth-integrated continuity gives a Helmholtz system for eta^{n+1} at
   !! the circumcentres:
   !!
   !!   [ M - g*theta^2*dt^2 * div_c( Hf grad_e ) ] eta^{n+1} = RHS
   !!
   !! with the CELL-CENTRED FV divergence (÷area_c, depth units) and the
   !! per-edge TPFA transmissibility `K_e = Hf_e * len_e / delta_e`
   !! (`delta_e` = circumcentre-to-circumcentre distance).  `M` is the cell
   !! mass coefficient (1 on the all-wet 2D path; the wet/dry Newton diagonal
   !! in P3.3).
   !!
   !! THE LOAD-BEARING DETAIL: the operator (T eta)_c = -(coef/area_c) *
   !! sum_e sign_e K_e (eta_nbr - eta_c) is self-adjoint in the AREA-WEIGHTED
   !! inner product <x,y> = sum_c area_c x_c y_c, NOT the plain dot — variable
   !! cell areas scale the rows.  Every CG/Jacobi reduction here carries the
   !! `area` weight; uniform area recovers the structured harness.  Get it
   !! wrong and CG stalls (prototype §7.5).
   !!
   !! CELL-CENTRIC DOUBLE-VISIT: every apply computes each cell's own row
   !! from its own 3 edges (reads neighbour eta, writes own Ap) — the house
   !! flux pattern, GPU-race-free.  Interior edges are visited twice (once
   !! from each side); no scatter, no atomics.  Single-rank, owned-cell loops.
   use rki_constants, only: wp, GRAVITY
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_si_unstr_mg, only: agglom_hierarchy_unstr_t, &
                              agglom_unstr_vcycle_precond
   implicit none
   private

   public :: si_unstr_edge_flux
   public :: si_unstr_face_depths
   public :: si_unstr_apply_helmholtz
   public :: si_unstr_compute_diag
   public :: si_unstr_compute_rhs
   public :: si_unstr_cg_solve

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.
   real(wp), parameter :: DELTA_EPS = 1.0e-30_wp
      !! 1/0 armour on the circumcentre-line distance (never fires on a
      !! gated mesh — the quality gate rejects delta<=0 at configure).
   integer, parameter :: SI_MAX_ITER = 4000
      !! CG iteration cap (triangle plain-CG is ~7x the structured count at
      !! matched Courant — prototype study 4; the cap leaves headroom).

contains

   pure function si_unstr_edge_flux(hf, eta_c, eta_nbr, inv_delta) result(flux)
      !! The single discrete edge-normal transport used by the TPFA
      !! operator's `div(Hf grad eta)`, the edge predictor's explicit PGF,
      !! and the edge back-substitution:
      !!   flux = Hf * (eta_nbr - eta_c) / delta   (along the outward normal)
      !! Sharing ONE flux function across the operator/back-sub is what makes
      !! the back-substituted edge transports satisfy discrete continuity
      !! identically (mass consistency).  `!$acc routine seq` for device
      !! do-concurrent bodies.
      !$omp declare target
      real(wp), intent(in) :: hf
         !! Frozen edge effective depth Hf_e.
      real(wp), intent(in) :: eta_c, eta_nbr
         !! eta at this cell / the across-edge neighbour (or Dirichlet ghost).
      real(wp), intent(in) :: inv_delta
         !! 1/delta_e (circumcentre-line distance reciprocal).
      real(wp) :: flux

      flux = hf*(eta_nbr - eta_c)*inv_delta
   end function si_unstr_edge_flux

   pure function si_unstr_face_depth_linear(h_c, h_nbr) result(hf)
      !! All-wet (2D barotropic) edge effective depth: the arithmetic mean of
      !! the two cells' still-water depths, `Hf = 0.5*(h_c + h_nbr)`.  This is
      !! the gated-mean rule (`face_depth_wetdry`) in its
      !! all-wet limit — it reduces EXACTLY to the linear edge depth.  Wall edges never call
      !! this (Hf stays 0); Dirichlet edges pass the ghost `h`.
      !$omp declare target
      real(wp), intent(in) :: h_c, h_nbr
      real(wp) :: hf

      hf = 0.5_wp*(h_c + h_nbr)
   end function si_unstr_face_depth_linear

   subroutine si_unstr_face_depths(mesh, ncell_total, ncell_owned, n_bc, &
                                        height, bed, bc_type, bc_eta_ext, &
                                        Hf_edge, eta_ghost_edge)
      !! Build the per-edge frozen effective depth Hf_e (2D all-wet path) AND
      !! the per-edge Dirichlet ghost surface `eta_ghost_edge`:
      !!   interior/partition edge: Hf = 0.5*(h_c + h_nbr); ghost unused (0).
      !!   WALL edge: Hf = 0 (no coupling — no-flux Neumann in the coeffs).
      !!   TIDAL/CLAMPED (Dirichlet) edge: the ghost surface is the prescribed
      !!     eta_bc; the ghost bed mirrors the cell bed so the ghost depth is
      !!     h_ghost = eta_bc - b_c, and Hf = 0.5*(h_c + h_ghost).  The ghost
      !!     eta is stashed in `eta_ghost_edge(i,k)` for the operator + RHS +
      !!     back-sub (the coefficient-space Dirichlet route).
      !!
      !! Cell-centric double-visit `(3, ncell)`: each cell fills its own 3
      !! edge values from its own `height`/`bed` + the neighbour's `height`
      !! (read-only); interior edges are visited twice with the same result.
      use rki_boundary_unstr, only: BC_WALL, BC_TIDAL, BC_CLAMPED
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: ncell_total
         !! Cell array extent (owned + halo) — the explicit-shape dim.
      integer, intent(in) :: ncell_owned
         !! Owned-cell loop bound.
      integer, intent(in) :: n_bc
         !! BC group table size (`MAX_BC_GROUPS`).
      real(wp), intent(in) :: height(ncell_total)
         !! Still-water depth h at cells (m).
      real(wp), intent(in) :: bed(ncell_total)
         !! Bed elevation b at cells (m).
      integer, intent(in) :: bc_type(n_bc)
         !! Per-group BC type (`state%bc_type`), the boundary lookup.
      real(wp), intent(in) :: bc_eta_ext(n_bc)
         !! Per-group Dirichlet external surface eta (m) (`state%bc_eta_ext`).
      real(wp), intent(out) :: Hf_edge(EDGES_PER_TRI, ncell_total)
         !! Edge effective depth (m).
      real(wp), intent(out) :: eta_ghost_edge(EDGES_PER_TRI, ncell_total)
         !! Dirichlet ghost surface eta (m) on Dirichlet boundary edges;
         !! 0 elsewhere (unread — Hf=0 there).

      integer :: k, i, kn, g, bt
      real(wp) :: h_ghost, eta_bc

      !$omp target teams distribute parallel do private(i, kn, g, bt, h_ghost, eta_bc)
      do k = 1, ncell_owned
         do i = 1, EDGES_PER_TRI
            eta_ghost_edge(i, k) = 0.0_wp
            kn = mesh%neighbours(i, k)
            if (kn > 0) then
               ! Interior / partition edge: arithmetic-mean still depth.
               Hf_edge(i, k) = si_unstr_face_depth_linear(height(k), height(kn))
            else
               g = mesh%bc_group(i, k)
               bt = BC_WALL
               if (g >= 1 .and. g <= n_bc) bt = bc_type(g)
               if ((bt == BC_TIDAL .or. bt == BC_CLAMPED) .and. g >= 1) then
                  eta_bc = bc_eta_ext(g)
                  h_ghost = max(eta_bc - bed(k), 0.0_wp)
                  Hf_edge(i, k) = si_unstr_face_depth_linear(height(k), h_ghost)
                  eta_ghost_edge(i, k) = eta_bc
               else
                  Hf_edge(i, k) = 0.0_wp
               end if
            end if
         end do
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_face_depths

   subroutine si_unstr_apply_helmholtz(mesh, ncell_total, ncell_owned, &
                                            eta, Ap, Hf_edge, g_theta2_dt2, &
                                            mass_diag)
      !! Matrix-free apply `Ap = A*eta` for the SPD TPFA Helmholtz operator:
      !!
      !!   Ap_c = M_c*eta_c - (g*theta^2*dt^2/area_c) *
      !!            sum_{e in c} sign_e * K_e * (eta_nbr - eta_c)
      !!
      !! with sign_e = +1 (flux LEAVES the cell along the outward normal),
      !! K_e = Hf_e*len_e/delta_e, and the divergence in depth units (÷area).
      !! Wall edges carry Hf=0 => no coupling (no-flux Neumann).  Dirichlet
      !! edges read the ghost eta the caller seeded into eta(kn) — the
      !! coefficient-space route (no separate elimination pass).
      !!
      !! `M_c` is the cell mass: absent `mass_diag` => M=1 (the linear
      !! all-wet operator); present => the wet/dry Newton diagonal a+floor.
      !!
      !! Cell-centric double-visit: reads neighbour eta, writes own Ap.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: ncell_total
         !! Cell array extent (owned + halo) — the explicit-shape dim.
      integer, intent(in) :: ncell_owned
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta(ncell_total)
      real(wp), intent(out) :: Ap(ncell_total)
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, ncell_total)
      real(wp), intent(in) :: g_theta2_dt2
      real(wp), intent(in), optional :: mass_diag(ncell_total)
         !! Per-cell Newton mass diagonal (wet/dry); absent => 1.0.

      integer :: k, i, kn
      real(wp) :: div_flux, ke, mass_k
      logical :: have_mass

      have_mass = present(mass_diag)

      !$omp target teams distribute parallel do private(i, kn, div_flux, ke, mass_k)
      do k = 1, ncell_owned
         div_flux = 0.0_wp
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, k) > 0.0_wp) then
               ke = Hf_edge(i, k)*mesh%edge_len(i, k) &
                    /max(mesh%edge_delta(i, k), DELTA_EPS)
               kn = mesh%neighbours(i, k)
               if (kn > 0) then
                  ! Interior edge: read the neighbour eta.
                  div_flux = div_flux + ke*(eta(kn) - eta(k))
               else
                  ! Dirichlet boundary edge (Hf>0 <=> tidal/clamped): the
                  ! FIXED ghost eta is NOT a CG unknown, so it reads as 0 in
                  ! the operator (A*p) — this stiffens the diagonal; the
                  ! ghost's real value is folded into the RHS (coefficient-
                  ! space Dirichlet route).  Contribution: ke*(0 - eta(k)).
                  div_flux = div_flux - ke*eta(k)
               end if
            end if
         end do
         if (have_mass) then
            mass_k = mass_diag(k)
         else
            mass_k = 1.0_wp
         end if
         Ap(k) = mass_k*eta(k) - g_theta2_dt2*div_flux/mesh%area(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_apply_helmholtz

   subroutine si_unstr_compute_diag(mesh, ncell_total, ncell_owned, &
                                         Hf_edge, eta_diag, g_theta2_dt2, &
                                         mass_diag)
      !! Precompute the reciprocal of the TPFA Helmholtz operator diagonal
      !! (Jacobi preconditioner):
      !!   diag_c = M_c + (g*theta^2*dt^2/area_c) * sum_{e in c} K_e
      !! stored as 1/diag.  M_c = 1 (linear) or a+floor (wet/dry, present
      !! `mass_diag`).  Wall edges (Hf=0) drop out => diag -> M there.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: ncell_total
         !! Cell array extent (owned + halo) — the explicit-shape dim.
      integer, intent(in) :: ncell_owned
         !! Owned-cell loop bound.
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, ncell_total)
      real(wp), intent(out) :: eta_diag(ncell_total)
      real(wp), intent(in) :: g_theta2_dt2
      real(wp), intent(in), optional :: mass_diag(ncell_total)

      integer :: k, i
      real(wp) :: ksum, ke, mass_k, diag_k
      logical :: have_mass

      have_mass = present(mass_diag)

      !$omp target teams distribute parallel do private(i, ksum, ke, mass_k, diag_k)
      do k = 1, ncell_owned
         ksum = 0.0_wp
         do i = 1, EDGES_PER_TRI
            if (Hf_edge(i, k) > 0.0_wp) then
               ke = Hf_edge(i, k)*mesh%edge_len(i, k) &
                    /max(mesh%edge_delta(i, k), DELTA_EPS)
               ksum = ksum + ke
            end if
         end do
         if (have_mass) then
            mass_k = mass_diag(k)
         else
            mass_k = 1.0_wp
         end if
         diag_k = mass_k + g_theta2_dt2*ksum/mesh%area(k)
         eta_diag(k) = 1.0_wp/diag_k
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_compute_diag

   subroutine si_unstr_compute_rhs(mesh, ncell_total, ncell_owned, &
                                   eta_n, u_star, u_n, rhs, dt, theta, &
                                   Hf_edge, eta_ghost_edge, g_theta2_dt2)
      !! Assemble the theta-blend RHS from the CELL edge transports:
      !!   RHS_c = eta^n_c - (dt/area_c) * sum_{e in c} sign_e *
      !!             [theta*U*_e + (1-theta)*U^n_e] * len_e
      !!             + (coef/area_c) * sum_{Dirichlet e in c} K_e * eta_ghost_e
      !! U*_e / U^n_e are per-edge transports (double-visit (3,ncell)); the
      !! blend flux telescopes over interior edges (each shared by 2 cells
      !! with opposite outward normals) and every wall edge carries 0
      !! (Hf=0 kills the predictor PGF, and U^n is zeroed on dead edges), so
      !! total mass is conserved to CG tolerance + round-off by construction.
      !! The Dirichlet term is the coefficient-space route: the operator
      !! (`si_unstr_apply_helmholtz`) treats the fixed ghost eta as 0, so its
      !! real value is carried here.
      !!
      !! The theta-blend is a STABILITY condition (not accuracy): bare
      !! div(U*) over-weights the explicit Laplacian by 1/theta and diverges
      !! at large Courant (structured proto A/B: ~1e11x growth) — same trap
      !! edge-wise.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: ncell_total
         !! Cell array extent (owned + halo) — the explicit-shape dim.
      integer, intent(in) :: ncell_owned
         !! Owned-cell loop bound.
      real(wp), intent(in) :: eta_n(ncell_total)
      real(wp), intent(in) :: u_star(EDGES_PER_TRI, ncell_total)
      real(wp), intent(in) :: u_n(EDGES_PER_TRI, ncell_total)
      real(wp), intent(out) :: rhs(ncell_total)
      real(wp), intent(in) :: dt, theta
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, ncell_total)
      real(wp), intent(in) :: eta_ghost_edge(EDGES_PER_TRI, ncell_total)
      real(wp), intent(in) :: g_theta2_dt2

      integer :: k, i, kn
      real(wp) :: one_m_theta, divU, flux_e, ghost_sum, ke

      one_m_theta = 1.0_wp - theta

      !$omp target teams distribute parallel do private(i, kn, divU, flux_e, ghost_sum, ke)
      do k = 1, ncell_owned
         divU = 0.0_wp
         ghost_sum = 0.0_wp
         do i = 1, EDGES_PER_TRI
            ! Blended edge transport (m^2/s); U^n already zeroed on dead
            ! edges by the caller (R3).  +sign: leaves the cell.
            flux_e = theta*u_star(i, k) + one_m_theta*u_n(i, k)
            divU = divU + flux_e*mesh%edge_len(i, k)
            ! Dirichlet ghost coupling (boundary edge with Hf>0).
            kn = mesh%neighbours(i, k)
            if (kn < 0 .and. Hf_edge(i, k) > 0.0_wp) then
               ke = Hf_edge(i, k)*mesh%edge_len(i, k) &
                    /max(mesh%edge_delta(i, k), DELTA_EPS)
               ghost_sum = ghost_sum + ke*eta_ghost_edge(i, k)
            end if
         end do
         rhs(k) = eta_n(k) - dt*divU/mesh%area(k) &
                  + g_theta2_dt2*ghost_sum/mesh%area(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine si_unstr_compute_rhs

   subroutine si_unstr_cg_solve(mesh, ncell_total, ncell_owned, &
                                eta, cg_r, cg_p, cg_Ap, cg_z, eta_diag, &
                                Hf_edge, g_theta2_dt2, precond, rtol, &
                                iter_out, mass_res_out, mass_diag, mg)
      !! Area-weighted (optionally preconditioned) CG for the TPFA Helmholtz
      !! eta-solve, warm-started from eta^n.  On entry `cg_r` holds the RHS and
      !! `eta` the warm start.  EVERY dot product carries the `mesh%area`
      !! weight (the operator is self-adjoint in <x,y>_A only — prototype
      !! §7.5); the convergence gate is the area-weighted relative residual
      !! `||r||_A^2 < rtol^2 ||rhs||_A^2` on the UNPRECONDITIONED residual in
      !! every mode (comparable iteration counts).  Preconditioners:
      !! NONE (plain), JACOBI (`z = eta_diag*r`), MG (aggregation-multigrid
      !! V-cycle via `mg`, `z = Vcycle_Ahat(area*r)` — the same `area*r`
      !! symmetric-space RHS the Jacobi apply implicitly uses).
      use rki_si_unstr_state, only: SI_PRECOND_NONE, SI_PRECOND_JACOBI, &
                                    SI_PRECOND_MG
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: ncell_total
         !! Cell array extent (owned + halo) — the explicit-shape dim.
      integer, intent(in) :: ncell_owned
         !! Owned-cell loop bound.
      real(wp), intent(inout) :: eta(ncell_total)
      real(wp), intent(inout) :: cg_r(ncell_total), cg_p(ncell_total)
      real(wp), intent(inout) :: cg_Ap(ncell_total), cg_z(ncell_total)
      real(wp), intent(inout) :: eta_diag(ncell_total)
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, ncell_total)
      real(wp), intent(in) :: g_theta2_dt2, rtol
      integer, intent(in) :: precond
      integer, intent(out) :: iter_out
      real(wp), intent(out) :: mass_res_out
         !! Area-weighted interior residual sum at solve end (mass diagnostic).
      real(wp), intent(in), optional :: mass_diag(ncell_total)
      type(agglom_hierarchy_unstr_t), intent(inout), optional :: mg
         !! Aggregation-MG hierarchy (required when precond == MG); its
         !! operator must already be refreshed for this step's frozen Hf.

      integer :: k, iter
      real(wp) :: rz, rr, rz_new, pAp, alpha, beta, tol2, rhs2, res_sum
      logical :: use_pc, use_mg, have_mass

      iter_out = 0
      mass_res_out = 0.0_wp
      use_mg = (precond == SI_PRECOND_MG)
      use_pc = (precond == SI_PRECOND_JACOBI) .or. use_mg
      have_mass = present(mass_diag)

      ! Jacobi needs the operator diagonal; MG carries its own per-level diag
      ! (under MG `eta_diag` stays STALE and is threaded through
      ! si_unstr_precond_apply UNREAD — do not consume it on that branch).
      if (precond == SI_PRECOND_JACOBI) then
         if (have_mass) then
            call si_unstr_compute_diag(mesh, ncell_total, ncell_owned, Hf_edge, &
                                       eta_diag, g_theta2_dt2, mass_diag=mass_diag)
         else
            call si_unstr_compute_diag(mesh, ncell_total, ncell_owned, Hf_edge, &
                                       eta_diag, g_theta2_dt2)
         end if
      end if

      ! ||rhs||_A^2 for the relative test (cg_r still holds the RHS).
      rhs2 = 0.0_wp
      !$omp target teams distribute parallel do
      do k = 1, ncell_owned
         rhs2 = rhs2 + mesh%area(k)*cg_r(k)*cg_r(k)
      end do
      !$omp end target teams distribute parallel do
      tol2 = max(rtol*rtol*rhs2, 1.0e-60_wp)

      ! Zero the full search direction (ghost/halo band stays 0 so a
      ! Dirichlet ghost is never a CG unknown).
      !$omp target teams distribute parallel do
      do k = 1, ncell_total
         cg_p(k) = 0.0_wp
      end do
      !$omp end target teams distribute parallel do

      ! r = rhs - A*eta (warm start eta = eta^n).
      if (have_mass) then
         call si_unstr_apply_helmholtz(mesh, ncell_total, ncell_owned, eta, &
                                       cg_Ap, Hf_edge, g_theta2_dt2, &
                                       mass_diag=mass_diag)
      else
         call si_unstr_apply_helmholtz(mesh, ncell_total, ncell_owned, eta, &
                                       cg_Ap, Hf_edge, g_theta2_dt2)
      end if
      !$omp target teams distribute parallel do
      do k = 1, ncell_owned
         cg_r(k) = cg_r(k) - cg_Ap(k)
      end do
      !$omp end target teams distribute parallel do

      rr = 0.0_wp
      !$omp target teams distribute parallel do
      do k = 1, ncell_owned
         rr = rr + mesh%area(k)*cg_r(k)*cg_r(k)
      end do
      !$omp end target teams distribute parallel do

      if (rr < tol2) then
         res_sum = 0.0_wp
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            res_sum = res_sum + mesh%area(k)*cg_r(k)
         end do
         !$omp end target teams distribute parallel do
         mass_res_out = res_sum
         return
      end if

      if (use_pc) then
         call si_unstr_precond_apply(mesh, ncell_total, ncell_owned, cg_r, cg_z, &
                                     eta_diag, use_mg, mg)
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            cg_p(k) = cg_z(k)
         end do
         !$omp end target teams distribute parallel do
         rz = 0.0_wp
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            rz = rz + mesh%area(k)*cg_r(k)*cg_z(k)
         end do
         !$omp end target teams distribute parallel do
      else
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            cg_p(k) = cg_r(k)
         end do
         !$omp end target teams distribute parallel do
         rz = rr
      end if

      do iter = 1, SI_MAX_ITER
         if (have_mass) then
            call si_unstr_apply_helmholtz(mesh, ncell_total, ncell_owned, cg_p, &
                                          cg_Ap, Hf_edge, g_theta2_dt2, &
                                          mass_diag=mass_diag)
         else
            call si_unstr_apply_helmholtz(mesh, ncell_total, ncell_owned, cg_p, &
                                          cg_Ap, Hf_edge, g_theta2_dt2)
         end if

         pAp = 0.0_wp
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            pAp = pAp + mesh%area(k)*cg_p(k)*cg_Ap(k)
         end do
         !$omp end target teams distribute parallel do

         if (abs(pAp) < 1.0e-30_wp) exit
         alpha = rz/pAp

         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            eta(k) = eta(k) + alpha*cg_p(k)
            cg_r(k) = cg_r(k) - alpha*cg_Ap(k)
         end do
         !$omp end target teams distribute parallel do

         rr = 0.0_wp
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            rr = rr + mesh%area(k)*cg_r(k)*cg_r(k)
         end do
         !$omp end target teams distribute parallel do

         iter_out = iter
         if (rr < tol2) exit

         if (use_pc) then
            call si_unstr_precond_apply(mesh, ncell_total, ncell_owned, cg_r, &
                                        cg_z, eta_diag, use_mg, mg)
            rz_new = 0.0_wp
            !$omp target teams distribute parallel do
            do k = 1, ncell_owned
               rz_new = rz_new + mesh%area(k)*cg_r(k)*cg_z(k)
            end do
            !$omp end target teams distribute parallel do
            beta = rz_new/rz
            rz = rz_new
            !$omp target teams distribute parallel do
            do k = 1, ncell_owned
               cg_p(k) = cg_z(k) + beta*cg_p(k)
            end do
            !$omp end target teams distribute parallel do
         else
            beta = rr/rz
            rz = rr
            !$omp target teams distribute parallel do
            do k = 1, ncell_owned
               cg_p(k) = cg_r(k) + beta*cg_p(k)
            end do
            !$omp end target teams distribute parallel do
         end if
      end do

      res_sum = 0.0_wp
      !$omp target teams distribute parallel do
      do k = 1, ncell_owned
         res_sum = res_sum + mesh%area(k)*cg_r(k)
      end do
      !$omp end target teams distribute parallel do
      mass_res_out = res_sum
   end subroutine si_unstr_cg_solve

   subroutine si_unstr_precond_apply(mesh, ncell_total, ncell_owned, r, z, &
                                     eta_diag, use_mg, mg)
      !! Preconditioner apply z = M^{-1} r for the area-weighted PCG.
      !! JACOBI: z = eta_diag*r (diagonal).  MG: one aggregation-multigrid
      !! V-cycle of the symmetric Ahat on the area-weighted RHS
      !! (`agglom_unstr_vcycle_precond` folds the `area*r`), which reduces to
      !! the Jacobi apply when the V-cycle is a single diagonal sweep — the
      !! generalisation is exact (module header).  Owned cells only; the ghost
      !! band of z is left as the caller set it (pre-zeroed).
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: ncell_total, ncell_owned
      real(wp), intent(in) :: r(ncell_total)
      real(wp), intent(inout) :: z(ncell_total)
      real(wp), intent(in) :: eta_diag(ncell_total)
      logical, intent(in) :: use_mg
      type(agglom_hierarchy_unstr_t), intent(inout), optional :: mg

      integer :: k

      if (use_mg) then
         ! Belt for direct callers: the wetdry-Newton / ML call sites omit
         ! `mg` and can never carry use_mg (enable_si fences mg to the 2D
         ! all-wet path) — but an absent optional here would be a silent
         ! deref, not a diagnosis.
         if (.not. present(mg)) then
            error stop "si_unstr_precond_apply: use_mg without the mg "// &
               "hierarchy argument (mg precond is 2D-barotropic-only)"
         end if
         call agglom_unstr_vcycle_precond(mg, r, z, mesh%area, ncell_total, &
                                          ncell_owned)
      else
         !$omp target teams distribute parallel do
         do k = 1, ncell_owned
            z(k) = eta_diag(k)*r(k)
         end do
         !$omp end target teams distribute parallel do
      end if
   end subroutine si_unstr_precond_apply

end module rki_si_unstr_helmholtz
