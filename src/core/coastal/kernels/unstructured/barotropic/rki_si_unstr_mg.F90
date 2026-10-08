!! Agglomeration-multigrid V-cycle preconditioner for the unstructured SI
!! eta-solve.
module rki_si_unstr_mg
   !! Aggregation (agglomeration) multigrid V(2,2)-cycle preconditioner for the
   !! TRIANGLE semi-implicit free-surface Helmholtz system — geometric MG has no
   !! triangle hierarchy, so we coarsen the CELL GRAPH by greedy pairwise
   !! aggregation instead.  The design is
   !! numpy-prototype-validated (`local_archive/prototypes/q12_agglom_mg`,
   !! STRONG GO: plain-CG plateau of thousands -> tens of MG-PCG iters).
   !!
   !! WHAT IT PRECONDITIONS.  The Fortran CG (`si_unstr_cg_solve`) runs the
   !! AREA-WEIGHTED PCG on the operator `A` (`si_unstr_apply_helmholtz`), inner
   !! product `<x,y>_area = sum_c area_c x_c y_c`.  This module works in the
   !! genuinely-SYMMETRIC space `Ahat = diag(area)*A` (SPD M-matrix):
   !!   Ahat[c,c]   = area_c + coef*sum_e K_e      (coef = g*theta^2*dt^2)
   !!   Ahat[c,nbr] = -coef*K_e ,   K_e = Hf_e*len_e/delta_e
   !! The Jacobi preconditioner the code already ships is, in this space,
   !! exactly `z = Ahat_diag^{-1} (area*r)`; the MG apply is the direct
   !! generalisation `z = Vcycle_Ahat(area*r)` — one V(2,2) cycle from a zero
   !! guess (SPD + linear => a valid CG preconditioner).  Passing `area*r`
   !! (not `r`) into the symmetric V-cycle is what makes it consistent with the
   !! area-weighted outer PCG (prototype driver §Ahat).
   !!
   !! HIERARCHY.  Cell graph coarsened ~2x/level by greedy pairwise
   !! aggregation (strongest-transmissibility partner).  Built ONCE at init
   !! (static mesh).  Every level's operator is stored in a COMPACT
   !! variable-valence edge form `nbr(mv,nc)` / `ke(mv,nc)` / `diag(nc)`; a
   !! coarse level's operator is the GALERKIN `P^T Ahat P` with a 0/1
   !! piecewise-constant P — which for aggregation equals a rediscretized
   !! agglomerate-TPFA (coarse mass = sum member areas, coarse transmissibility
   !! = sum of fine transmissibilities across the cut) EXACTLY, so the coarse
   !! STENCIL DOES NOT GROW (prototype: Galerkin == rediscretized to
   !! round-off).  Transfers are P (inject) / R = P^T (sum-gather) — exact
   !! transpose by construction, no scaling factor.
   !!
   !! REFRESH.  Only the transmissibilities `ke` + diagonals depend on the
   !! per-step frozen `Hf` (coef is fixed).  The STRUCTURE (aggregates, coarse
   !! adjacency, member lists) is static; `ke`/`diag`/`dinv` are rebuilt ONCE
   !! PER STEP on the HOST (fine level from the mesh geometry + current Hf,
   !! coarse levels by the Galerkin cascade) and pushed to the device — the
   !! same host-side per-step operator refresh the structured MG uses.  NOT per
   !! V-cycle.
   !!
   !! GPU.  Each level owns explicit-shape device arrays + enter/exit_data.
   !! The level loop runs on the HOST; every device kernel (apply / damped-
   !! Jacobi smooth / residual / restrict / prolong) is called with FLAT
   !! EXPLICIT-SHAPE arguments — the level array is NEVER dereferenced inside a
   !! `do concurrent` (the array-of-derived-types device-indirection gotcha).
   !! Kernels are cell-centric (each cell writes only its own row / entry;
   !! restriction sums a coarse cell's own members) — no neighbour scatter, no
   !! atomics.  Damped Jacobi ping-pongs (read old, write new) — race-free.
   !! Single-rank v1 (no per-cycle halo).
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use pic_logger, only: global_logger
   implicit none
   private

   public :: agglom_level_unstr_t, agglom_hierarchy_unstr_t
   public :: agglom_hierarchy_unstr_build, agglom_hierarchy_unstr_refresh
   public :: agglom_unstr_vcycle_precond
   public :: AGGLOM_MIN_SIZE, AGGLOM_MAX_LEVELS, AGGLOM_OMEGA
   public :: AGGLOM_NU_PRE, AGGLOM_NU_POST, AGGLOM_COARSE_SWEEPS

   integer, parameter :: EDGES_PER_TRI = 3
   integer, parameter :: AGGLOM_MIN_SIZE = 150
      !! Stop coarsening when a level has <= this many cells (prototype-tuned:
      !! 2x pairwise coarsening down to ~150 cells, then a coarse Jacobi solve).
   integer, parameter :: AGGLOM_MAX_LEVELS = 16
      !! Hard cap on the level count (fixed-size level array bound).
   real(wp), parameter :: AGGLOM_OMEGA = 0.8_wp
      !! Damped-Jacobi relaxation factor (matches the structured MG).
   ! NU_PRE == NU_POST is LOAD-BEARING: the V-cycle is self-adjoint (and
   ! therefore a valid plain-PCG preconditioner) only with symmetric
   ! pre/post smoothing — an unequal edit silently breaks SPD.
   integer, parameter :: AGGLOM_NU_PRE = 2
      !! Pre-smoothing Jacobi sweeps per level.
   integer, parameter :: AGGLOM_NU_POST = 2
      !! Post-smoothing Jacobi sweeps per level.
   integer, parameter :: AGGLOM_COARSE_SWEEPS = 40
      !! Coarsest-level Jacobi sweeps (approximate solve from zero).  FIXED (not
      !! convergence-gated) on purpose: a fixed sweep count makes the whole
      !! V-cycle a FIXED LINEAR SPD operator, which the outer STANDARD (non-
      !! flexible) area-weighted PCG requires — a residual-gated coarse solve
      !! would vary the sweep count with the RHS, making M RHS-dependent
      !! (flexible preconditioning) and breaking the constant-M PCG assumption.
      !! EFFICIENCY CAVEAT: for a very stiff Helmholtz (large SI dt) 40 sweeps
      !! can under-resolve the ~<=150-cell coarse system and lift the outer PCG
      !! count on larger/stiffer meshes; the funnel CI case is well inside the
      !! budget.  Robustify (if ever needed) by raising this constant, not by
      !! gating — keep the operator fixed + linear.
   real(wp), parameter :: DELTA_EPS = 1.0e-30_wp
      !! 1/0 armour on the circumcentre-line distance (matches the operator).

   type :: agglom_level_unstr_t
      !! One aggregation-MG level: the compact variable-valence Helmholtz
      !! operator (symmetric Ahat) on `nc` cells + the V-cycle work buffers +
      !! the aggregation map to the next coarser level.  Structure arrays
      !! (nbr / n_nbr / agg_of / member CSR) are STATIC (built once); the
      !! coefficient arrays (ke / diag / dinv) are refreshed per step.  All
      !! arrays explicit-shape; the level loop dereferences on the HOST and
      !! passes flat args to the device kernels.
      integer :: nc = 0
         !! Cells on this level.
      integer :: mv = 0
         !! Max valence (fixed leading dim of nbr/ke; unused slots padded
         !! nbr=1, ke=0 so the apply loop is branchless).

      ! --- STATIC operator STRUCTURE (built once at init) ---
      integer, allocatable :: n_nbr(:)
         !! Valence per cell, shape (nc).
      integer, allocatable :: nbr(:, :)
         !! Neighbour cell id on THIS level, shape (mv, nc); padded slots = 1.
      integer, allocatable :: agg_of(:)
         !! Aggregate id in the NEXT COARSER level for each of this level's
         !! cells, shape (nc).  Drives prolongation (fine cell reads its
         !! parent).  Unallocated on the coarsest level.

      ! Member CSR: which cells of the NEXT FINER level belong to each of THIS
      ! level's cells (this level is the coarse side).  Drives restriction (a
      ! coarse cell sums its own members) + the host Galerkin build.
      integer, allocatable :: mstart(:)
         !! Member CSR offsets, shape (nc+1).  Unallocated on the finest level.
      integer, allocatable :: mlist(:)
         !! Member cell ids on the finer level, shape (n_finer).

      ! --- PER-STEP coefficients (host + device; refreshed each step) ---
      real(wp), allocatable :: ke(:, :)
         !! Transmissibility coef*K_e per edge, shape (mv, nc); off-diag = -ke.
      real(wp), allocatable :: diag(:)
         !! Operator diagonal, shape (nc).
      real(wp), allocatable :: dinv(:)
         !! 1/diag (Jacobi), shape (nc).

      ! --- V-cycle work buffers (device) ---
      real(wp), allocatable :: e(:)
         !! Error iterate, shape (nc).
      real(wp), allocatable :: r(:)
         !! Level RHS / residual, shape (nc).
      real(wp), allocatable :: tmp(:)
         !! Jacobi ping-pong / residual scratch, shape (nc).
   contains
      procedure, non_overridable :: enter_data => agglom_level_enter_data
      procedure, non_overridable :: exit_data => agglom_level_exit_data
      procedure, non_overridable :: destroy => agglom_level_destroy
   end type agglom_level_unstr_t

   type :: agglom_hierarchy_unstr_t
      !! Fine-to-coarse level stack for the aggregation-MG preconditioner.
      !! Allocated + built only when `precond == mg`.
      integer :: n_levels = 0
      type(agglom_level_unstr_t), allocatable :: lev(:)
   contains
      procedure, non_overridable :: enter_data => agglom_hierarchy_enter_data
      procedure, non_overridable :: exit_data => agglom_hierarchy_exit_data
      procedure, non_overridable :: destroy => agglom_hierarchy_destroy
   end type agglom_hierarchy_unstr_t

contains

   ! ======================================================================
   !  BUILD (init, host) — aggregate the cell graph, size + fill the static
   !  structure of every level.  Uses the STILL-DEPTH Hf as the aggregation
   !  strength; the per-step refresh overwrites ke/diag with the live Hf.
   ! ======================================================================

   subroutine agglom_hierarchy_unstr_build(this, mesh, nco)
      !! Build the level stack for `nco` owned cells.  Only the mesh GEOMETRY
      !! is needed: the aggregation strength is the geometric transmissibility
      !! `len/delta` (a uniform-`Hf`, `coef=1` proxy).  Greedy pairwise
      !! aggregation picks the MAX-strength neighbour, which is invariant to a
      !! uniform positive scale, so the live per-step `coef*Hf` (folded in at
      !! refresh) never changes the STRUCTURE — the hierarchy is static.
      !! Host-only (runs before enter_data); the diag/ke are placeholders here
      !! and are overwritten by the first `refresh` before any solve.
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nco

      integer :: l, nc, l_alloc, k, i, kn
      real(wp), allocatable :: str_l(:, :)
      real(wp), allocatable :: hf0(:, :)

      ! ---- Level 1: the fine mesh operator structure (mv = 3) ----
      call build_level1_structure(this, mesh, nco)
      call alloc_level_coeffs(this%lev(1))
      allocate (str_l(this%lev(1)%mv, this%lev(1)%nc))
      ! Uniform per-edge Hf=1 on interior edges => aggregation strength =
      ! len/delta (the geometric transmissibility).  Boundary edges left 0.
      allocate (hf0(EDGES_PER_TRI, nco), source=0.0_wp)
      do k = 1, nco
         do i = 1, EDGES_PER_TRI
            kn = mesh%neighbours(i, k)
            if (kn > 0 .and. kn <= nco) hf0(i, k) = 1.0_wp
         end do
      end do
      call level1_operator(mesh, nco, 1.0_wp, hf0, this%lev(1), str_l)
      deallocate (hf0)

      l = 1
      do
         nc = this%lev(l)%nc
         if (nc <= AGGLOM_MIN_SIZE .or. l >= AGGLOM_MAX_LEVELS) exit
         ! Aggregate level l (strength = |ke| = str_l) -> define level l+1.
         call aggregate_and_build_coarse(this, l, str_l)
         if (this%lev(l + 1)%nc >= nc) then
            ! No coarsening achieved: drop the stub level, stop.
            call this%lev(l + 1)%destroy()
            this%n_levels = l
            exit
         end if
         call alloc_level_coeffs(this%lev(l + 1))
         if (allocated(str_l)) deallocate (str_l)
         allocate (str_l(this%lev(l + 1)%mv, this%lev(l + 1)%nc))
         call galerkin_operator(this%lev(l), this%lev(l + 1), str_l)
         l = l + 1
         this%n_levels = l
      end do
      if (allocated(str_l)) deallocate (str_l)

      ! Shrink the level array to n_levels (it was over-allocated).
      l_alloc = size(this%lev)
      if (this%n_levels < l_alloc) call trim_levels(this)

      call global_logger%info("si_unstr_mg: built "// &
                              itoa(this%n_levels)//" levels (fine nc="// &
                              itoa(this%lev(1)%nc)//", coarse nc="// &
                              itoa(this%lev(this%n_levels)%nc)//")")
   end subroutine agglom_hierarchy_unstr_build

   subroutine build_level1_structure(this, mesh, nco)
      !! Allocate the level array + fill level 1's static interior-neighbour
      !! structure from the mesh (mv = EDGES_PER_TRI = 3).
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nco

      integer :: k, i, kn, ns

      allocate (this%lev(AGGLOM_MAX_LEVELS))
      this%n_levels = 1
      this%lev(1)%nc = nco
      this%lev(1)%mv = EDGES_PER_TRI
      allocate (this%lev(1)%n_nbr(nco), source=0)
      allocate (this%lev(1)%nbr(EDGES_PER_TRI, nco), source=1)

      do k = 1, nco
         ns = 0
         do i = 1, EDGES_PER_TRI
            kn = mesh%neighbours(i, k)
            if (kn > 0 .and. kn <= nco) then
               ns = ns + 1
               this%lev(1)%nbr(ns, k) = kn
            end if
         end do
         this%lev(1)%n_nbr(k) = ns
      end do
   end subroutine build_level1_structure

   subroutine level1_operator(mesh, nco, coef, Hf_edge, lev, strength)
      !! Fill level 1's ke/diag/dinv from the mesh geometry + the per-edge
      !! frozen Hf (still-water synthetic at build, the live `Hf_edge` at
      !! refresh).  Mirrors `si_unstr_compute_diag` + the interior off-diagonals
      !! of `si_unstr_apply_helmholtz`, so the symmetric Ahat row is exact:
      !!   diag(k) = area_k + coef*sum_{Hf>0} K_e   (interior + Dirichlet)
      !!   off-diag ke to interior neighbours only.
      !! A Dirichlet boundary edge carries Hf>0 with no CG-unknown neighbour =>
      !! it lands on the diagonal (ghost = homogeneous 0 in the error
      !! operator).  Wall edges carry Hf=0 => dropped.  `strength` returns |ke|
      !! per interior edge for aggregation.  HOST; caller pushes to the device.
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nco
      real(wp), intent(in) :: coef
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nco)
      type(agglom_level_unstr_t), intent(inout) :: lev
      real(wp), intent(out), optional :: strength(lev%mv, nco)
         !! |ke| per edge for the NEXT aggregation.  Needed at BUILD only;
         !! the per-step refresh omits it (no scratch allocation in the
         !! per-step path).

      integer :: k, i, kn, ns
      real(wp) :: kloc, diagk

      if (present(strength)) strength = 0.0_wp
      lev%ke = 0.0_wp
      do k = 1, nco
         diagk = mesh%area(k)
         ns = 0
         do i = 1, EDGES_PER_TRI
            kn = mesh%neighbours(i, k)
            ! Diagonal: every Hf>0 edge contributes coef*K_e (interior OR
            ! Dirichlet boundary edge — the latter lands on the diagonal with
            ! no CG-unknown neighbour slot).  Hf<=0 (wall / dried edge) => 0.
            if (Hf_edge(i, k) > 0.0_wp) then
               kloc = coef*Hf_edge(i, k)*mesh%edge_len(i, k) &
                      /max(mesh%edge_delta(i, k), DELTA_EPS)
               diagk = diagk + kloc
            else
               kloc = 0.0_wp
            end if
            ! Off-diagonal slot: advance `ns` on the SAME interior-edge
            ! predicate as `build_level1_structure` so the ke slot index stays
            ! LOCKED to the nbr slot index — independent of whether this edge's
            ! Hf is zero.  An interior edge that dries (Hf<=0) writes a ZERO
            ! (wall-like) transmissibility into its own aligned slot rather than
            ! skipping it, keeping the operator symmetric + SPD even with mixed
            ! zero/nonzero interior Hf (the plain 2D-barotropic dry-cell path).
            if (kn > 0 .and. kn <= nco) then
               ns = ns + 1
               lev%ke(ns, k) = kloc
               if (present(strength)) strength(ns, k) = kloc
            end if
         end do
         lev%diag(k) = diagk
         lev%dinv(k) = 1.0_wp/diagk
      end do
   end subroutine level1_operator

   subroutine aggregate_and_build_coarse(this, l, strength)
      !! Greedy pairwise aggregation of level l (strongest-transmissibility
      !! partner), then build level l+1's static structure: nc, member CSR,
      !! coarse adjacency (nbr / n_nbr / mv).  `strength(mv,nc)` = |ke| per
      !! edge on level l.
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      integer, intent(in) :: l
      real(wp), intent(in) :: strength(:, :)

      integer :: nf, k, s, j, best, nc, i
      real(wp) :: bw, w
      integer, allocatable :: aggof(:)

      nf = this%lev(l)%nc
      allocate (aggof(nf), source=0)
      nc = 0
      do k = 1, nf
         if (aggof(k) /= 0) cycle
         best = -1
         bw = 0.0_wp
         do s = 1, this%lev(l)%n_nbr(k)
            j = this%lev(l)%nbr(s, k)
            if (aggof(j) /= 0) cycle
            w = strength(s, k)
            if (w > bw) then
               bw = w
               best = j
            end if
         end do
         nc = nc + 1
         aggof(k) = nc
         if (best > 0) aggof(best) = nc
      end do

      ! Store agg_of on level l (drives prolongation).
      allocate (this%lev(l)%agg_of(nf))
      this%lev(l)%agg_of = aggof

      ! Level l+1 skeleton.
      this%lev(l + 1)%nc = nc

      ! Member CSR (which level-l cells belong to each coarse cell), STANDARD
      ! 1-based: members of coarse cell I are mlist(mstart(I):mstart(I+1)-1).
      allocate (this%lev(l + 1)%mstart(nc + 1), source=0)
      do k = 1, nf
         this%lev(l + 1)%mstart(aggof(k) + 1) = &
            this%lev(l + 1)%mstart(aggof(k) + 1) + 1   ! count(I) into slot I+1
      end do
      this%lev(l + 1)%mstart(1) = 1
      do i = 1, nc
         this%lev(l + 1)%mstart(i + 1) = this%lev(l + 1)%mstart(i) &
                                         + this%lev(l + 1)%mstart(i + 1)
      end do
      allocate (this%lev(l + 1)%mlist(nf))
      block
         integer, allocatable :: cursor(:)
         allocate (cursor(nc))
         cursor = this%lev(l + 1)%mstart(1:nc)   ! cursor(I) = mstart(I)
         do k = 1, nf
            i = aggof(k)
            this%lev(l + 1)%mlist(cursor(i)) = k
            cursor(i) = cursor(i) + 1
         end do
         deallocate (cursor)
      end block

      call build_coarse_adjacency(this%lev(l), this%lev(l + 1), aggof)
      deallocate (aggof)
   end subroutine aggregate_and_build_coarse

   subroutine build_coarse_adjacency(levf, levc, aggof)
      !! Build coarse-level neighbour lists from the fine adjacency: coarse
      !! cells I,J are neighbours iff some fine edge crosses the I-J aggregate
      !! boundary.  Sets levc%mv / n_nbr / nbr (unique neighbours per coarse
      !! cell).  Host, O(fine edges).
      type(agglom_level_unstr_t), intent(in) :: levf
      type(agglom_level_unstr_t), intent(inout) :: levc
      integer, intent(in) :: aggof(:)

      integer :: nc, k, s, j, capI, capJ, mv, i
      integer, allocatable :: cnt(:), tmpnbr(:, :)

      nc = levc%nc
      allocate (cnt(nc), source=0)
      ! First pass: over-estimate max valence = sum of member valences.
      mv = 0
      do capI = 1, nc
         capJ = 0
         do i = levc%mstart(capI), levc%mstart(capI + 1) - 1
            k = levc%mlist(i)
            capJ = capJ + levf%n_nbr(k)
         end do
         mv = max(mv, capJ)
      end do
      mv = max(mv, 1)
      allocate (tmpnbr(mv, nc), source=0)

      do capI = 1, nc
         do i = levc%mstart(capI), levc%mstart(capI + 1) - 1
            k = levc%mlist(i)
            do s = 1, levf%n_nbr(k)
               j = aggof(levf%nbr(s, k))
               if (j == capI) cycle
               call push_unique(tmpnbr(:, capI), cnt(capI), j)
            end do
         end do
      end do

      levc%mv = max(maxval(cnt), 1)
      allocate (levc%n_nbr(nc))
      allocate (levc%nbr(levc%mv, nc), source=1)
      do capI = 1, nc
         levc%n_nbr(capI) = cnt(capI)
         do s = 1, cnt(capI)
            levc%nbr(s, capI) = tmpnbr(s, capI)
         end do
      end do
      deallocate (cnt, tmpnbr)
   end subroutine build_coarse_adjacency

   subroutine push_unique(list, n, val)
      !! Append val to list(1:n) if not already present; bump n.
      integer, intent(inout) :: list(:)
      integer, intent(inout) :: n
      integer, intent(in) :: val
      integer :: s
      do s = 1, n
         if (list(s) == val) return
      end do
      n = n + 1
      list(n) = val
   end subroutine push_unique

   subroutine galerkin_operator(levf, levc, strength)
      !! Build the coarse operator ke/diag/dinv by the Galerkin cascade from
      !! the finer level's ke/diag (== rediscretized agglomerate-TPFA for 0/1
      !! aggregation).  `strength(mvc,ncc)` returns |ke| for the NEXT
      !! aggregation.  HOST; caller pushes to the device.
      type(agglom_level_unstr_t), intent(in) :: levf
      type(agglom_level_unstr_t), intent(inout) :: levc
      real(wp), intent(out), optional :: strength(levc%mv, levc%nc)
         !! |ke| per coarse edge for the NEXT aggregation (BUILD only; the
         !! per-step refresh omits it).

      integer :: capI, i, k, s, j, capJ, t
      real(wp) :: diagI, v

      if (present(strength)) strength = 0.0_wp
      levc%ke = 0.0_wp
      do capI = 1, levc%nc
         diagI = 0.0_wp
         do i = levc%mstart(capI), levc%mstart(capI + 1) - 1
            k = levc%mlist(i)
            diagI = diagI + levf%diag(k)
            do s = 1, levf%n_nbr(k)
               j = levf%nbr(s, k)
               v = levf%ke(s, k)
               capJ = levf%agg_of(j)
               if (capJ == capI) then
                  diagI = diagI - v            ! internal edge folds to diag
               else
                  t = find_slot(levc%nbr(:, capI), levc%n_nbr(capI), capJ)
                  ! ke >= 0 invariant (Hf >= 0 => level-1 ke >= 0 => Galerkin
                  ! sums stay >= 0), so the signed ke doubles as the next
                  ! aggregation strength; a negative ke would silently
                  ! degrade coarsening to singletons.
                  levc%ke(t, capI) = levc%ke(t, capI) + v
                  if (present(strength)) strength(t, capI) = levc%ke(t, capI)
               end if
            end do
         end do
         levc%diag(capI) = diagI
         levc%dinv(capI) = 1.0_wp/diagI
      end do
   end subroutine galerkin_operator

   integer function find_slot(nbr, n, tgt) result(t)
      !! Locate the coarse-edge slot for neighbour `tgt` in nbr(1:n).  The
      !! target is ALWAYS present — the coarse adjacency + the Galerkin fold
      !! are built from the same aggregation map, so every folded fine edge has
      !! a coarse slot.  A not-found here means that invariant broke; fail loud
      !! (host-only build path) rather than silently mislocate a coefficient
      !! into slot 1.
      integer, intent(in) :: nbr(:)
      integer, intent(in) :: n, tgt
      integer :: s
      t = 0
      do s = 1, n
         if (nbr(s) == tgt) then
            t = s
            return
         end if
      end do
      call global_logger%error("si_unstr_mg: find_slot could not locate coarse "// &
                               "neighbour "//itoa(tgt)//" in the coarse adjacency "// &
                               "— Galerkin fold inconsistent with the aggregation map")
      error stop "si_unstr_mg: find_slot: coarse neighbour not found"
   end function find_slot

   subroutine alloc_level_coeffs(lev)
      !! Allocate the per-step coefficient + V-cycle work arrays of a level.
      type(agglom_level_unstr_t), intent(inout) :: lev
      allocate (lev%ke(lev%mv, lev%nc), source=0.0_wp)
      allocate (lev%diag(lev%nc), source=1.0_wp)
      allocate (lev%dinv(lev%nc), source=1.0_wp)
      allocate (lev%e(lev%nc), source=0.0_wp)
      allocate (lev%r(lev%nc), source=0.0_wp)
      allocate (lev%tmp(lev%nc), source=0.0_wp)
   end subroutine alloc_level_coeffs

   subroutine trim_levels(this)
      !! Move-shrink the over-allocated level array to n_levels.
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      type(agglom_level_unstr_t), allocatable :: keep(:)
      integer :: l
      allocate (keep(this%n_levels))
      do l = 1, this%n_levels
         call move_level(this%lev(l), keep(l))
      end do
      call move_alloc(keep, this%lev)
   end subroutine trim_levels

   subroutine move_level(src, dst)
      !! move_alloc every component from src to dst (host, pre-enter_data).
      type(agglom_level_unstr_t), intent(inout) :: src
      type(agglom_level_unstr_t), intent(inout) :: dst
      dst%nc = src%nc
      dst%mv = src%mv
      if (allocated(src%n_nbr)) call move_alloc(src%n_nbr, dst%n_nbr)
      if (allocated(src%nbr)) call move_alloc(src%nbr, dst%nbr)
      if (allocated(src%agg_of)) call move_alloc(src%agg_of, dst%agg_of)
      if (allocated(src%mstart)) call move_alloc(src%mstart, dst%mstart)
      if (allocated(src%mlist)) call move_alloc(src%mlist, dst%mlist)
      if (allocated(src%ke)) call move_alloc(src%ke, dst%ke)
      if (allocated(src%diag)) call move_alloc(src%diag, dst%diag)
      if (allocated(src%dinv)) call move_alloc(src%dinv, dst%dinv)
      if (allocated(src%e)) call move_alloc(src%e, dst%e)
      if (allocated(src%r)) call move_alloc(src%r, dst%r)
      if (allocated(src%tmp)) call move_alloc(src%tmp, dst%tmp)
   end subroutine move_level

   ! ======================================================================
   !  REFRESH (per step, host build + device push)
   ! ======================================================================

   subroutine agglom_hierarchy_unstr_refresh(this, mesh, nco, coef, Hf_edge)
      !! Rebuild every level's ke/diag/dinv from the current frozen per-edge
      !! `Hf_edge` (the same array the operator + RHS use) + push to the
      !! device.  Fine level from the mesh; coarse levels by the Galerkin
      !! cascade.  Structures are static.  Runs ONCE PER STEP (host loop), NOT
      !! per V-cycle.  `Hf_edge` is device-resident; the caller pulls it self
      !! before calling (host build reads it).
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: nco
      real(wp), intent(in) :: coef
      real(wp), intent(in) :: Hf_edge(EDGES_PER_TRI, nco)

      integer :: l
      call level1_operator(mesh, nco, coef, Hf_edge, this%lev(1))
      call level_coeffs_to_device(this%lev(1))

      do l = 2, this%n_levels
         call galerkin_operator(this%lev(l - 1), this%lev(l))
         call level_coeffs_to_device(this%lev(l))
      end do
   end subroutine agglom_hierarchy_unstr_refresh

   subroutine level_coeffs_to_device(lev)
      !! Push the host-refreshed coefficient arrays to the device.
      type(agglom_level_unstr_t), intent(inout) :: lev
      !$acc update device(lev%ke, lev%diag, lev%dinv)
   end subroutine level_coeffs_to_device

   ! ======================================================================
   !  DEVICE KERNELS (flat explicit-shape args; cell-centric, race-free)
   ! ======================================================================

   subroutine agglom_jacobi_one(e_old, e_new, r, nbr, ke, diag, dinv, &
                                mv, nc, omega)
      !! One damped-Jacobi sweep e_new = e_old + omega*(r - Ahat*e_old)/diag —
      !! reads e_old only, writes e_new only (race-free).
      integer, intent(in) :: mv, nc
      real(wp), intent(in) :: e_old(nc)
      real(wp), intent(out) :: e_new(nc)
      real(wp), intent(in) :: r(nc)
      integer, intent(in) :: nbr(mv, nc)
      real(wp), intent(in) :: ke(mv, nc)
      real(wp), intent(in) :: diag(nc), dinv(nc)
      real(wp), intent(in) :: omega

      integer :: k, s
      real(wp) :: acc

      do concurrent(k=1:nc) local(s, acc)
         acc = diag(k)*e_old(k)
         do s = 1, mv
            acc = acc - ke(s, k)*e_old(nbr(s, k))
         end do
         e_new(k) = e_old(k) + omega*(r(k) - acc)*dinv(k)
      end do
   end subroutine agglom_jacobi_one

   subroutine agglom_jacobi_sweeps(e, r, tmp, nbr, ke, diag, dinv, &
                                   mv, nc, omega, nu)
      !! nu damped-Jacobi sweeps, ping-pong (e,tmp); latest iterate in `e`.
      integer, intent(in) :: mv, nc, nu
      real(wp), intent(inout) :: e(nc)
      real(wp), intent(in) :: r(nc)
      real(wp), intent(inout) :: tmp(nc)
      integer, intent(in) :: nbr(mv, nc)
      real(wp), intent(in) :: ke(mv, nc)
      real(wp), intent(in) :: diag(nc), dinv(nc)
      real(wp), intent(in) :: omega

      integer :: s, np

      np = nu/2
      do s = 1, np
         call agglom_jacobi_one(e, tmp, r, nbr, ke, diag, dinv, mv, nc, omega)
         call agglom_jacobi_one(tmp, e, r, nbr, ke, diag, dinv, mv, nc, omega)
      end do
      if (mod(nu, 2) == 1) then
         call agglom_jacobi_one(e, tmp, r, nbr, ke, diag, dinv, mv, nc, omega)
         call agglom_copy(tmp, e, nc)
      end if
   end subroutine agglom_jacobi_sweeps

   subroutine agglom_residual(e, r, res, nbr, ke, diag, mv, nc)
      !! res = r - Ahat*e.
      integer, intent(in) :: mv, nc
      real(wp), intent(in) :: e(nc), r(nc)
      real(wp), intent(out) :: res(nc)
      integer, intent(in) :: nbr(mv, nc)
      real(wp), intent(in) :: ke(mv, nc)
      real(wp), intent(in) :: diag(nc)

      integer :: k, s
      real(wp) :: acc

      do concurrent(k=1:nc) local(s, acc)
         acc = diag(k)*e(k)
         do s = 1, mv
            acc = acc - ke(s, k)*e(nbr(s, k))
         end do
         res(k) = r(k) - acc
      end do
   end subroutine agglom_residual

   subroutine agglom_restrict(rf, rc, mstart, mlist, ncc, ncf)
      !! rc = P^T rf: each coarse cell sums its own fine members (cell-centric,
      !! race-free — a coarse cell writes only its own entry).
      integer, intent(in) :: ncc, ncf
      real(wp), intent(in) :: rf(ncf)
      real(wp), intent(out) :: rc(ncc)
      integer, intent(in) :: mstart(ncc + 1)
      integer, intent(in) :: mlist(ncf)

      integer :: capI, i
      real(wp) :: acc

      do concurrent(capI=1:ncc) local(i, acc)
         acc = 0.0_wp
         do i = mstart(capI), mstart(capI + 1) - 1
            acc = acc + rf(mlist(i))
         end do
         rc(capI) = acc
      end do
   end subroutine agglom_restrict

   subroutine agglom_prolong_add(ec, ef, agg_of, ncf, ncc)
      !! ef += P ec: each fine cell adds its parent aggregate's correction
      !! (piecewise-constant inject; fine cell writes only its own entry).
      integer, intent(in) :: ncf, ncc
      real(wp), intent(in) :: ec(ncc)
      real(wp), intent(inout) :: ef(ncf)
      integer, intent(in) :: agg_of(ncf)

      integer :: k

      do concurrent(k=1:ncf)
         ef(k) = ef(k) + ec(agg_of(k))
      end do
   end subroutine agglom_prolong_add

   subroutine agglom_copy(src, dst, nc)
      integer, intent(in) :: nc
      real(wp), intent(in) :: src(nc)
      real(wp), intent(out) :: dst(nc)
      integer :: k
      do concurrent(k=1:nc)
         dst(k) = src(k)
      end do
   end subroutine agglom_copy

   subroutine agglom_zero(a, nc)
      integer, intent(in) :: nc
      real(wp), intent(out) :: a(nc)
      integer :: k
      do concurrent(k=1:nc)
         a(k) = 0.0_wp
      end do
   end subroutine agglom_zero

   ! ======================================================================
   !  V-CYCLE + PRECONDITIONER APPLY
   ! ======================================================================

   recursive subroutine agglom_vcycle(this, l)
      !! One V(nu_pre,nu_post) cycle solving Ahat*e = r on level l from ZERO.
      !! On entry lev(l)%r holds the RHS; on exit lev(l)%e holds the solution.
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      integer, intent(in) :: l

      integer :: nc, ncc

      nc = this%lev(l)%nc
      call agglom_zero(this%lev(l)%e, nc)

      if (l == this%n_levels) then
         call agglom_jacobi_sweeps(this%lev(l)%e, this%lev(l)%r, this%lev(l)%tmp, &
                                   this%lev(l)%nbr, this%lev(l)%ke, &
                                   this%lev(l)%diag, this%lev(l)%dinv, &
                                   this%lev(l)%mv, nc, &
                                   AGGLOM_OMEGA, AGGLOM_COARSE_SWEEPS)
         return
      end if

      call agglom_jacobi_sweeps(this%lev(l)%e, this%lev(l)%r, this%lev(l)%tmp, &
                                this%lev(l)%nbr, this%lev(l)%ke, &
                                this%lev(l)%diag, this%lev(l)%dinv, &
                                this%lev(l)%mv, nc, &
                                AGGLOM_OMEGA, AGGLOM_NU_PRE)

      call agglom_residual(this%lev(l)%e, this%lev(l)%r, this%lev(l)%tmp, &
                           this%lev(l)%nbr, this%lev(l)%ke, this%lev(l)%diag, &
                           this%lev(l)%mv, nc)

      ncc = this%lev(l + 1)%nc
      call agglom_restrict(this%lev(l)%tmp, this%lev(l + 1)%r, &
                           this%lev(l + 1)%mstart, this%lev(l + 1)%mlist, ncc, nc)

      call agglom_vcycle(this, l + 1)

      call agglom_prolong_add(this%lev(l + 1)%e, this%lev(l)%e, &
                              this%lev(l)%agg_of, nc, ncc)

      call agglom_jacobi_sweeps(this%lev(l)%e, this%lev(l)%r, this%lev(l)%tmp, &
                                this%lev(l)%nbr, this%lev(l)%ke, &
                                this%lev(l)%diag, this%lev(l)%dinv, &
                                this%lev(l)%mv, nc, &
                                AGGLOM_OMEGA, AGGLOM_NU_POST)
   end subroutine agglom_vcycle

   subroutine agglom_unstr_vcycle_precond(this, r, z, area, nct, nco)
      !! z = M^{-1} r for the AREA-WEIGHTED outer PCG: form rhs_hat = area*r,
      !! run ONE V(2,2) cycle of the symmetric Ahat from zero, copy the
      !! correction into z (owned cells).  `r`/`z`/`area` are the fine
      !! (ncell_total-shaped) device arrays; only the first `nco` owned entries
      !! are touched (ghost band left as the caller set it — pre-zeroed).
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      integer, intent(in) :: nct, nco
      real(wp), intent(in) :: r(nct)
      real(wp), intent(inout) :: z(nct)
      real(wp), intent(in) :: area(nct)

      ! Direct-API belt: enable_si always builds; a caller that sets
      ! precond=mg without build_mg would deref unallocated levels.
      if (this%n_levels < 1) then
         error stop "agglom_unstr_vcycle_precond: hierarchy not built "// &
            "(state_unstr_enable_si builds it; direct callers must too)"
      end if

      ! rhs_hat = area * r  ->  lev(1)%r  (symmetric-space RHS)
      call agglom_scale_into_shim(r, area, this%lev(1)%r, nct, nco)
      call agglom_vcycle(this, 1)
      call agglom_copy_out(this%lev(1)%e, z, nct, nco)
   end subroutine agglom_unstr_vcycle_precond

   subroutine agglom_scale_into_shim(src, fac, dst, nct, nco)
      !! dst(1:nco) = fac*src (explicit-shape ncell_total state arrays -> the
      !! level-1 RHS buffer).
      integer, intent(in) :: nct, nco
      real(wp), intent(in) :: src(nct), fac(nct)
      real(wp), intent(out) :: dst(nco)
      integer :: k
      do concurrent(k=1:nco)
         dst(k) = fac(k)*src(k)
      end do
   end subroutine agglom_scale_into_shim

   subroutine agglom_copy_out(src, dst, nct, nco)
      !! dst(1:nco) = src (level-1 buffer -> state z, owned cells).
      integer, intent(in) :: nct, nco
      real(wp), intent(in) :: src(nco)
      real(wp), intent(inout) :: dst(nct)
      integer :: k
      do concurrent(k=1:nco)
         dst(k) = src(k)
      end do
   end subroutine agglom_copy_out

   ! ======================================================================
   !  enter/exit/destroy
   ! ======================================================================

   subroutine agglom_level_enter_data(this)
      class(agglom_level_unstr_t), intent(inout) :: this
      select type (this)
      type is (agglom_level_unstr_t)
         call agglom_level_enter_data_impl(this)
      end select
   end subroutine agglom_level_enter_data

   subroutine agglom_level_enter_data_impl(lev)
      type(agglom_level_unstr_t), intent(inout) :: lev
      !$acc enter data copyin(lev%nbr, lev%n_nbr, lev%ke, lev%diag, lev%dinv)
      !$acc enter data copyin(lev%e, lev%r, lev%tmp)
      if (allocated(lev%agg_of)) then
         !$acc enter data copyin(lev%agg_of)
      end if
      if (allocated(lev%mstart)) then
         !$acc enter data copyin(lev%mstart, lev%mlist)
      end if
   end subroutine agglom_level_enter_data_impl

   subroutine agglom_level_exit_data(this)
      class(agglom_level_unstr_t), intent(inout) :: this
      select type (this)
      type is (agglom_level_unstr_t)
         call agglom_level_exit_data_impl(this)
      end select
   end subroutine agglom_level_exit_data

   subroutine agglom_level_exit_data_impl(lev)
      type(agglom_level_unstr_t), intent(inout) :: lev
      if (allocated(lev%mstart)) then
         !$acc exit data delete(lev%mstart, lev%mlist)
      end if
      if (allocated(lev%agg_of)) then
         !$acc exit data delete(lev%agg_of)
      end if
      !$acc exit data delete(lev%e, lev%r, lev%tmp)
      !$acc exit data delete(lev%nbr, lev%n_nbr, lev%ke, lev%diag, lev%dinv)
   end subroutine agglom_level_exit_data_impl

   subroutine agglom_level_destroy(this)
      class(agglom_level_unstr_t), intent(inout) :: this
      if (allocated(this%n_nbr)) deallocate (this%n_nbr)
      if (allocated(this%nbr)) deallocate (this%nbr)
      if (allocated(this%agg_of)) deallocate (this%agg_of)
      if (allocated(this%mstart)) deallocate (this%mstart)
      if (allocated(this%mlist)) deallocate (this%mlist)
      if (allocated(this%ke)) deallocate (this%ke)
      if (allocated(this%diag)) deallocate (this%diag)
      if (allocated(this%dinv)) deallocate (this%dinv)
      if (allocated(this%e)) deallocate (this%e)
      if (allocated(this%r)) deallocate (this%r)
      if (allocated(this%tmp)) deallocate (this%tmp)
      this%nc = 0
      this%mv = 0
   end subroutine agglom_level_destroy

   subroutine agglom_hierarchy_enter_data(this)
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      integer :: l
      do l = 1, this%n_levels
         call this%lev(l)%enter_data()
      end do
   end subroutine agglom_hierarchy_enter_data

   subroutine agglom_hierarchy_exit_data(this)
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      integer :: l
      do l = 1, this%n_levels
         call this%lev(l)%exit_data()
      end do
   end subroutine agglom_hierarchy_exit_data

   subroutine agglom_hierarchy_destroy(this)
      class(agglom_hierarchy_unstr_t), intent(inout) :: this
      integer :: l
      if (allocated(this%lev)) then
         do l = 1, size(this%lev)
            call this%lev(l)%destroy()
         end do
         deallocate (this%lev)
      end if
      this%n_levels = 0
   end subroutine agglom_hierarchy_destroy

   pure function itoa(i) result(s)
      !! Minimal integer->string for the build log line.
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, "(I0)") i
      s = trim(buf)
   end function itoa

end module rki_si_unstr_mg
