!! Semi-implicit free-surface solver state (unstructured triangular 2D barotropic).
module rki_si_unstr_state
   !! Persistent workspace for the semi-implicit (Casulli) free-surface
   !! mode on the UNSTRUCTURED (triangular) backend
   !! (`time_scheme="semi_implicit"`, `grid_type="unstructured"`).  The
   !! semi-implicit workspace for the unstructured backend; allocated only when the mode is
   !! on and living at `state%si` on `unstr_state_t`.
   !!
   !! GEOMETRY OF RECORD (cell-centred FV on a triangulation, prototype
   !! `si_unstr_proto.py`): eta lives at the triangle CIRCUMCENTRE (mesh
   !! `circumcentre_x/y`, NOT the centroid the KNP path uses); the SI face
   !! transports live on EDGES.  A triangle's three edges are stored in the
   !! `(3, ncell)` per-edge-of-cell layout (edge i opposite vertex i) — the
   !! same double-visit layout as the KNP flux path: each cell owns its 3
   !! edge values, interior edges are computed twice (no scatter, no
   !! atomics, GPU-race-free).
   !!
   !! THE LOAD-BEARING CG DETAIL: the cell-centred FV Helmholtz operator is
   !! self-adjoint in the AREA-WEIGHTED inner product `<x,y> = sum_c area_c
   !! x_c y_c`, NOT the plain dot (variable cell areas scale the rows).  The
   !! CG/Jacobi dot products carry the `mesh%area` weight; the residual gate
   !! uses the area-weighted norm.  Get it wrong and the operator is
   !! non-symmetric and CG stalls.  Uniform area recovers the structured
   !! harness.  (See `rki_si_unstr_helmholtz`.)
   !!
   !! All CG vectors are `(ncell)`; the per-edge workspaces are `(3, ncell)`.
   !! The `enter_data`/`exit_data` methods attach the arrays to the device;
   !! the parent `unstr_state_t` struct is mapped by `state_unstr_enter_data`
   !! BEFORE this runs (same discipline as the structured `si` slot).
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_si_unstr_mg, only: agglom_hierarchy_unstr_t, &
                              agglom_hierarchy_unstr_build
   use pic_logger, only: global_logger
   implicit none
   private

   public :: semi_implicit_state_unstr_t
   public :: parse_si_precond_unstr
   public :: SI_PRECOND_NONE, SI_PRECOND_JACOBI, SI_PRECOND_MG

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.

   integer, parameter :: SI_PRECOND_NONE = 0
      !! Unpreconditioned CG (area-weighted plain CG).
   integer, parameter :: SI_PRECOND_JACOBI = 1
      !! Diagonal (Jacobi) PCG: z = eta_diag*r.  Matters MORE on triangles
      !! than on the structured grid — variable cell areas put real
      !! disparity on the Helmholtz diagonal even at uniform depth
      !! (measured diag spread O(100-1000)x on a jittered basin;
      !! ~2.5x iteration win — prototype study 4).
   integer, parameter :: SI_PRECOND_MG = 2
      !! Aggregation (agglomeration) multigrid V(2,2)-cycle PCG — the triangle
      !! plateau-breaker (`rki_si_unstr_mg`).  Greedy pairwise cell-graph
      !! coarsening + Galerkin coarse operators + damped-Jacobi smoothing;
      !! numpy-prototype-validated to collapse the plain-CG plateau of
      !! thousands to tens of iters (Q12).  Opt-in; the hierarchy `mg` slot is
      !! built at `enable_si` and refreshed once per step.

   type :: semi_implicit_state_unstr_t
      !! CG workspace + edge face-depths for the 2D free-surface Helmholtz
      !! solve on the triangulation.

      integer :: ncells_total = 0
         !! Cell count cached at init (all CG vectors are this long).
      integer :: nz_si = 0
         !! Layer count for the 3D SI path (0 or 1 = 2D-only; > 1 = the
         !! P3.2 multilayer column slots below are allocated).

      ! Free-surface elevation solved implicitly, at the CIRCUMCENTRE.
      real(wp), allocatable :: eta(:)
         !! Surface elevation eta^{n+1} (m), size ncells_total.  eta = h + b
         !! at circumcentres.  Warm-started from eta^n each step; the CG
         !! iterate.  Boundary/halo entries hold the seeded Dirichlet ghost.

      ! CG vectors (matrix-free), all size ncells_total.
      real(wp), allocatable :: eta_r(:)
         !! Residual r = rhs - A*eta.
      real(wp), allocatable :: eta_p(:)
         !! Search direction p.
      real(wp), allocatable :: eta_Ap(:)
         !! Operator apply A*p.
      real(wp), allocatable :: eta_z(:)
         !! Preconditioned residual z = M^{-1} r (Jacobi); unread on plain CG.
      real(wp), allocatable :: eta_diag(:)
         !! 1/diag of the Helmholtz operator (Jacobi preconditioner).

      ! Per-edge-of-cell prognostic face transport U_e (m^2/s), the
      ! staggered (Casulli) quantity.  Double-visit `(3, ncell)`: transport
      ! along the OUTWARD edge-normal of edge i of cell k.  Holds the
      ! predictor U*_e between the predictor and the back-substitution, then
      ! the final U^{n+1}_e.  The A-grid cell momenta xmom_cv/ymom_cv are a
      ! step-end least-squares reconstruction for the KNP seam.
      real(wp), allocatable :: u_edge(:, :)
         !! Edge-normal transport (m^2/s), shape (3, ncells_total); + along
         !! the stored outward edge-normal.
      real(wp), allocatable :: u_star(:, :)
         !! Predictor U*_e (m^2/s), shape (3, ncells_total).

      ! Per-edge frozen effective depth Hf weighting the div(Hf grad eta)
      ! TPFA stencil (wall edges 0; open/Dirichlet edges real).  Shape
      ! (3, ncell); shared by the operator, back-sub and recovery (frozen
      ! over the step, rule R2).
      real(wp), allocatable :: Hf_edge(:, :)
         !! Edge effective depth (m), shape (3, ncells_total).
      real(wp), allocatable :: eta_ghost_edge(:, :)
         !! Dirichlet ghost surface eta (m) on tidal/clamped boundary edges,
         !! shape (3, ncells_total); 0 elsewhere (unread — Hf=0 there).
         !! Filled alongside Hf_edge each step (frozen, rule R2).
      real(wp), allocatable :: u_n_edge(:, :)
         !! Saved U^n edge transport (m^2/s), shape (3, ncells_total) — the
         !! (1-theta) side of the theta-blend RHS; dead edges zeroed (R3).

      ! CG controls (`&si_nml`; host scalars, set at init from cfg).
      integer :: precond = SI_PRECOND_NONE
         !! Preconditioner (`SI_PRECOND_*`): NONE = plain CG (default),
         !! JACOBI = diagonal PCG, MG = aggregation-multigrid PCG.
      real(wp) :: rtol = 1.0e-10_wp
         !! CG relative residual tolerance (area-weighted
         !! ||r||^2 < rtol^2*||rhs||^2).

      ! Aggregation-multigrid preconditioner hierarchy — built at enable_si
      ! (host) only when precond == MG, refreshed once per step, attached to
      ! the device by this state's enter_data.
      logical :: use_mg = .false.
         !! MG preconditioner active (host scalar; set by enable_si).
      type(agglom_hierarchy_unstr_t) :: mg
         !! Fine-to-coarse aggregation-MG level stack (`rki_si_unstr_mg`).

      ! Diagnostics (host scalars).
      integer :: cg_iters = 0
         !! CG iterations on the most recent step (go/no-go cost number).
      real(wp) :: mass_res = 0.0_wp
         !! Area-weighted interior residual sum at solve end (m^3) — the
         !! per-step mass-consistency diagnostic.
      integer :: newton_iters = 0
         !! Wet/dry Newton outer iterations on the most recent step (host
         !! scalar; read by tests / logging).

      ! --- Casulli wetting/drying (P3.3) — allocated only when `wetdry`
      !     is on.  The V(eta)/a(eta)/Newton machinery; per-cell,
      !     mesh-agnostic (prototype-confirmed).  All size (ncells) with
      !     the sub-bottom table (n_sub, ncells), cell index LAST. ---
      logical :: wetdry = .false.
         !! Wet/dry mode on (host scalar; set by `state_unstr_enable_si`).
      integer :: n_sub = 1
         !! Subgrid sub-bottoms per cell (1 = v1 flat bottom).
      integer :: bsub_fill = 0
         !! Sub-bottom fill mode (`SI_BSUB_FILL_*` from rki_si_common; 0=FLAT
         !! default ⇒ bit-identical, 1=LINEAR ⇒ analytic subgrid ramp).  Set
         !! by `state_unstr_enable_si`.
      real(wp), allocatable :: b_sub(:, :)
         !! Per-cell SORTED sub-bottom bed elevations (n_sub, ncells) (m).
      real(wp), allocatable :: b_min(:)
         !! Deepest sub-bottom per cell (m) — the sill bed at edges.
      real(wp), allocatable :: a_wet(:)
         !! Wet fraction a(eta) in [0,1] at the current Newton iterate.
      real(wp), allocatable :: mass_diag(:)
         !! Newton mass diagonal a+floor — the operator's `M` term.
      real(wp), allocatable :: v_cell(:)
         !! Wetted mean depth v(eta) at the current iterate (m).
      real(wp), allocatable :: rhs_d(:)
         !! Depth-units Newton RHS v(eta^n) - dt*div[blend] — held across
         !! the Newton loop.
      real(wp), allocatable :: newton_dz(:)
         !! Newton increment: the inner CG solves J*dz = -F from zero.
      real(wp), allocatable :: newton_F(:)
         !! Nonlinear residual F = v + T*eta - rhs.
      real(wp), allocatable :: newton_Jeta(:)
         !! Scratch for the full mass-diag operator apply J*eta (the
         !! stiffness is recovered as J*eta - mass_diag*eta).

      ! --- P3.2 (3D multilayer) column slots — allocated only when
      !     nz_si > 1.  Layer convention: (nz, ncells) axis order — k=1
      !     bed, k=nz surface, cell index LAST (unstructured convention,
      !     OPPOSITE of the structured (nx,ny,nz)). ---
      real(wp), allocatable :: g_col_u(:, :)
         !! Column x-velocity response G_k^u = solve(u^n + tendencies)
         !! (nz, ncells_total) [m/s].
      real(wp), allocatable :: g_col_v(:, :)
         !! Column y-velocity response G_k^v (nz, ncells_total) [m/s].
      real(wp), allocatable :: r_col(:, :)
         !! Column response R_k = solve(-1) (nz, ncells_total) [m/s]
         !! (<= 0); axis-independent — u and v share the tridiagonal.
      real(wp), allocatable :: nu_iface(:, :)
         !! PP81 interface viscosity (nz+1, ncells_total) [m^2/s]; index k
         !! = interface BELOW layer k (bed row 1 and surface row nz+1 = 0).
      real(wp), allocatable :: Hhat_cell(:)
         !! Effective column depth Hhat = -sum_k h_k R_k (ncells_total) [m].
      real(wp), allocatable :: Gbar_cell_x(:)
         !! Column-integrated predictor transport sum_k h_k G_k^u
         !! (ncells_total) [m^2/s].
      real(wp), allocatable :: Gbar_cell_y(:)
         !! Column-integrated y predictor transport (ncells_total) [m^2/s].
      real(wp), allocatable :: u_edge_layer(:, :, :)
         !! Recovered per-layer edge-normal transports
         !! (nz, EDGES_PER_TRI, ncells_total) [m^2/s] — double-visit, same
         !! outward-normal convention as `u_edge`.  Column sums telescope to
         !! the 2D back-substituted edge transport (mass consistency).
      real(wp), allocatable :: mass_flux_edge(:, :)
         !! Blended column edge mass flux theta*sum_k u_edge_layer +
         !! (1-theta)*U^n_e (EDGES_PER_TRI, ncells_total) [m^2/s] — the
         !! EXACT flux the eta update used; drives the dsig tracer split.

   contains
      procedure, non_overridable :: init => si_unstr_state_init
      procedure, non_overridable :: destroy => si_unstr_state_destroy
      procedure, non_overridable :: enter_data => si_unstr_state_enter_data
      procedure, non_overridable :: exit_data => si_unstr_state_exit_data
      procedure, non_overridable :: build_mg => si_unstr_state_build_mg
   end type semi_implicit_state_unstr_t

contains

   subroutine si_unstr_state_init(this, mesh, nz, wetdry, n_sub)
      !! Allocate workspace arrays at the mesh size, zero-filled.  When `nz`
      !! is present and > 1 the P3.2 (3D multilayer) column slots are also
      !! allocated; when `wetdry` is present and .true. the P3.3 Casulli
      !! wet/dry slots are (with `n_sub` sub-bottoms, default 1); otherwise
      !! only the 2D barotropic slots exist and the 2D path is unchanged.
      !! Plain host allocation (runs before `enter_data`, so a device loop
      !! would round-trip).
      class(semi_implicit_state_unstr_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in), optional :: nz
         !! Number of vertical layers.  Absent or <= 1 selects the 2D path.
      logical, intent(in), optional :: wetdry
         !! Casulli wetting/drying on.  Absent/.false. => slots unallocated.
      integer, intent(in), optional :: n_sub
         !! Subgrid sub-bottoms per cell (>= 1).  Absent => 1 (flat bottom).

      integer :: n, nz_loc, ns
      logical :: wd

      n = mesh%ncells_total
      this%ncells_total = n
      nz_loc = 0
      if (present(nz)) nz_loc = nz
      wd = .false.
      if (present(wetdry)) wd = wetdry
      ns = 1
      if (present(n_sub)) ns = max(1, n_sub)

      allocate (this%eta(n), source=0.0_wp)
      allocate (this%eta_r(n), source=0.0_wp)
      allocate (this%eta_p(n), source=0.0_wp)
      allocate (this%eta_Ap(n), source=0.0_wp)
      allocate (this%eta_z(n), source=0.0_wp)
      allocate (this%eta_diag(n), source=0.0_wp)
      allocate (this%u_edge(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%u_star(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%Hf_edge(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%eta_ghost_edge(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%u_n_edge(EDGES_PER_TRI, n), source=0.0_wp)

      ! --- P3.2 (3D multilayer) column slots (only when nz > 1) ---
      if (nz_loc > 1) then
         this%nz_si = nz_loc
         allocate (this%g_col_u(nz_loc, n), source=0.0_wp)
         allocate (this%g_col_v(nz_loc, n), source=0.0_wp)
         allocate (this%r_col(nz_loc, n), source=0.0_wp)
         allocate (this%nu_iface(nz_loc + 1, n), source=0.0_wp)
         allocate (this%Hhat_cell(n), source=0.0_wp)
         allocate (this%Gbar_cell_x(n), source=0.0_wp)
         allocate (this%Gbar_cell_y(n), source=0.0_wp)
         allocate (this%u_edge_layer(nz_loc, EDGES_PER_TRI, n), source=0.0_wp)
         allocate (this%mass_flux_edge(EDGES_PER_TRI, n), source=0.0_wp)
      end if

      ! --- P3.3 Casulli wet/dry slots (only when wetdry on) ---
      if (wd) then
         this%wetdry = .true.
         this%n_sub = ns
         allocate (this%b_sub(ns, n), source=0.0_wp)
         allocate (this%b_min(n), source=0.0_wp)
         allocate (this%a_wet(n), source=1.0_wp)
         allocate (this%mass_diag(n), source=1.0_wp)
         allocate (this%v_cell(n), source=0.0_wp)
         allocate (this%rhs_d(n), source=0.0_wp)
         allocate (this%newton_dz(n), source=0.0_wp)
         allocate (this%newton_F(n), source=0.0_wp)
         allocate (this%newton_Jeta(n), source=0.0_wp)
      end if
   end subroutine si_unstr_state_init

   subroutine si_unstr_state_build_mg(this, mesh)
      !! Build the aggregation-MG preconditioner hierarchy for this state.
      !! Host-only; call from `enable_si` (after `init`, before `enter_data`)
      !! when `precond == SI_PRECOND_MG`.  Structure only (geometric
      !! aggregation); the operator coefficients are filled by the first
      !! per-step `refresh`.
      class(semi_implicit_state_unstr_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in) :: mesh

      call agglom_hierarchy_unstr_build(this%mg, mesh, mesh%ncells_owned)
      this%use_mg = .true.
   end subroutine si_unstr_state_build_mg

   subroutine si_unstr_state_destroy(this)
      !! Deallocate all workspace arrays.
      class(semi_implicit_state_unstr_t), intent(inout) :: this

      if (this%use_mg) then
         call this%mg%destroy()
         this%use_mg = .false.
      end if

      if (allocated(this%eta)) deallocate (this%eta)
      if (allocated(this%eta_r)) deallocate (this%eta_r)
      if (allocated(this%eta_p)) deallocate (this%eta_p)
      if (allocated(this%eta_Ap)) deallocate (this%eta_Ap)
      if (allocated(this%eta_z)) deallocate (this%eta_z)
      if (allocated(this%eta_diag)) deallocate (this%eta_diag)
      if (allocated(this%u_edge)) deallocate (this%u_edge)
      if (allocated(this%u_star)) deallocate (this%u_star)
      if (allocated(this%Hf_edge)) deallocate (this%Hf_edge)
      if (allocated(this%eta_ghost_edge)) deallocate (this%eta_ghost_edge)
      if (allocated(this%u_n_edge)) deallocate (this%u_n_edge)

      ! P3.2 (3D multilayer) column slots.
      if (allocated(this%g_col_u)) deallocate (this%g_col_u)
      if (allocated(this%g_col_v)) deallocate (this%g_col_v)
      if (allocated(this%r_col)) deallocate (this%r_col)
      if (allocated(this%nu_iface)) deallocate (this%nu_iface)
      if (allocated(this%Hhat_cell)) deallocate (this%Hhat_cell)
      if (allocated(this%Gbar_cell_x)) deallocate (this%Gbar_cell_x)
      if (allocated(this%Gbar_cell_y)) deallocate (this%Gbar_cell_y)
      if (allocated(this%u_edge_layer)) deallocate (this%u_edge_layer)
      if (allocated(this%mass_flux_edge)) deallocate (this%mass_flux_edge)

      ! P3.3 Casulli wet/dry slots.
      if (allocated(this%b_sub)) deallocate (this%b_sub)
      if (allocated(this%b_min)) deallocate (this%b_min)
      if (allocated(this%a_wet)) deallocate (this%a_wet)
      if (allocated(this%mass_diag)) deallocate (this%mass_diag)
      if (allocated(this%v_cell)) deallocate (this%v_cell)
      if (allocated(this%rhs_d)) deallocate (this%rhs_d)
      if (allocated(this%newton_dz)) deallocate (this%newton_dz)
      if (allocated(this%newton_F)) deallocate (this%newton_F)
      if (allocated(this%newton_Jeta)) deallocate (this%newton_Jeta)
      this%wetdry = .false.
      this%n_sub = 1
      this%bsub_fill = 0

      ! Reset so a destroy -> re-init-as-2D cycle doesn't leave nz_si > 1
      ! pointing enter_data at unallocated slots.
      this%nz_si = 0
      this%ncells_total = 0
   end subroutine si_unstr_state_destroy

   subroutine si_unstr_state_enter_data(this)
      !! Attach the workspace arrays to the device.  Arrays-only — the parent
      !! `unstr_state_t` is mapped by `state_unstr_enter_data` before this
      !! runs.  TBP wraps a non-polymorphic `_impl` (the `class(...)` box
      !! maps to stack under AMD libomptarget; delegate to `type(...)`).
      class(semi_implicit_state_unstr_t), intent(inout) :: this
      select type (this)
      type is (semi_implicit_state_unstr_t)
         call si_unstr_state_enter_data_impl(this)
      end select
   end subroutine si_unstr_state_enter_data

   subroutine si_unstr_state_enter_data_impl(this)
      type(semi_implicit_state_unstr_t), intent(inout) :: this

      !$omp target enter data map(to: this%eta, this%eta_r, this%eta_p, this%eta_Ap, this%eta_z)
      !$omp target enter data map(to: this%eta_diag, this%u_edge, this%u_star, this%Hf_edge)
      !$omp target enter data map(to: this%eta_ghost_edge, this%u_n_edge)

      ! P3.2 (3D multilayer) column slots — only when allocated.
      if (this%nz_si > 1) then
         !$omp target enter data map(to: this%g_col_u, this%g_col_v, this%r_col)
         !$omp target enter data map(to: this%nu_iface, this%Hhat_cell)
         !$omp target enter data map(to: this%Gbar_cell_x, this%Gbar_cell_y)
         !$omp target enter data map(to: this%u_edge_layer, this%mass_flux_edge)
      end if

      ! P3.3 Casulli wet/dry slots — only when allocated.
      if (this%wetdry) then
         !$omp target enter data map(to: this%b_sub, this%b_min, this%a_wet)
         !$omp target enter data map(to: this%mass_diag, this%v_cell, this%rhs_d)
         !$omp target enter data map(to: this%newton_dz, this%newton_F, this%newton_Jeta)
      end if

      ! Aggregation-MG hierarchy — only when built (precond == MG).
      if (this%use_mg) call this%mg%enter_data()
   end subroutine si_unstr_state_enter_data_impl

   subroutine si_unstr_state_exit_data(this)
      !! Detach the workspace arrays from the device.
      class(semi_implicit_state_unstr_t), intent(inout) :: this
      select type (this)
      type is (semi_implicit_state_unstr_t)
         call si_unstr_state_exit_data_impl(this)
      end select
   end subroutine si_unstr_state_exit_data

   subroutine si_unstr_state_exit_data_impl(this)
      type(semi_implicit_state_unstr_t), intent(inout) :: this

      ! Aggregation-MG hierarchy — only when built (precond == MG).
      if (this%use_mg) call this%mg%exit_data()

      ! P3.3 Casulli wet/dry slots — reverse order of enter.
      if (this%wetdry) then
         !$omp target exit data map(delete: this%newton_dz, this%newton_F, this%newton_Jeta)
         !$omp target exit data map(delete: this%mass_diag, this%v_cell, this%rhs_d)
         !$omp target exit data map(delete: this%b_sub, this%b_min, this%a_wet)
      end if

      ! P3.2 (3D multilayer) column slots — reverse order of enter.
      if (this%nz_si > 1) then
         !$omp target exit data map(delete: this%u_edge_layer, this%mass_flux_edge)
         !$omp target exit data map(delete: this%Gbar_cell_x, this%Gbar_cell_y)
         !$omp target exit data map(delete: this%nu_iface, this%Hhat_cell)
         !$omp target exit data map(delete: this%g_col_u, this%g_col_v, this%r_col)
      end if

      !$omp target exit data map(delete: this%eta_ghost_edge, this%u_n_edge)
      !$omp target exit data map(delete: this%eta_diag, this%u_edge, this%u_star, this%Hf_edge)
      !$omp target exit data map(delete: this%eta, this%eta_r, this%eta_p, this%eta_Ap, this%eta_z)
   end subroutine si_unstr_state_exit_data_impl

   function parse_si_precond_unstr(tag) result(precond)
      !! Map the `&si_nml precond` string to its `SI_PRECOND_*` enum on the
      !! unstructured path.  Fail-loud allowlist: an unknown tag aborts.  MG is
      !! aggregation (agglomeration) multigrid here (`rki_si_unstr_mg`) — the
      !! triangle plateau-breaker — NOT the structured geometric MG.
      character(len=*), intent(in) :: tag
      integer :: precond

      select case (trim(adjustl(tag)))
      case ("none", "NONE", "")
         precond = SI_PRECOND_NONE
      case ("jacobi", "JACOBI")
         precond = SI_PRECOND_JACOBI
      case ("mg", "MG")
         precond = SI_PRECOND_MG
      case default
         call global_logger%error("si_precond (unstructured): unknown "// &
                                  "preconditioner '"//trim(adjustl(tag))// &
                                  "' (allowed: none, jacobi, mg)")
         error stop "si_precond: unknown preconditioner"
      end select
   end function parse_si_precond_unstr

end module rki_si_unstr_state
