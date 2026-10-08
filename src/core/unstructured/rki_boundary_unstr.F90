!! Boundary condition type constants for the unstructured solver
module rki_boundary_unstr
   !! Integer constants identifying the boundary condition type for
   !! each boundary edge group, plus a fixed maximum on the number of
   !! groups a single mesh can carry. Mirrors the role of
   !! `rki_boundary_types` on the structured side. Seven BC types are
   !! fully dispatched inside `compute_boundary_state_impl`
   !! (inlined into `rki_kernel_flux_unstr`): `BC_WALL`, `BC_OPEN`, `BC_INFLOW`,
   !! `BC_DISCHARGE`, `BC_OUTFLOW`, `BC_TIDAL`, `BC_CLAMPED`. `BC_CHAPMAN`
   !! is dispatched one level higher in the flux kernel
   !! (`rki_kernel_flux_unstr`) because it needs per-edge mesh access to
   !! compute the radiating ghost target before falling through into the
   !! Flather/`BC_TIDAL` algebra. The multilayer tracer advection kernel
   !! (`ml_advect_tracer_unstr`) honours `bc_S_ext` / `bc_T_ext` at
   !! inflow edges, so freshwater / clamped-tracer river inflows work
   !! end-to-end. Types that live on the structured side but are not yet
   !! ported: `BC_NESTED`, `BC_SPONGE`.
   !!
   !! `BC_TIDAL` is a Flather radiating boundary that clamps the ghost
   !! stage to `bc_eta_ext(g)` and reconstructs the ghost outward-normal
   !! velocity from the outgoing Riemann invariant. `BC_CLAMPED` does
   !! the same stage clamp but copies momentum transmissively (no
   !! radiation term) -- pick it when you want a pure Dirichlet stage
   !! without the wave-radiation coupling.
   !!
   !! How a mesh wires up its boundary edges:
   !!
   !!   * Every boundary edge of a triangle has `mesh%neighbours(i, k) == -1`
   !!     and `mesh%bc_group(i, k) >= 0`. Group 0 means "no group set
   !!     -- treat as default wall".
   !!   * `state%bc_type(g)` for `g = 1..MAX_BC_GROUPS` says what BC
   !!     type group `g` represents. Defaults to `BC_WALL` so a fresh
   !!     state with no explicit BC setup behaves like an all-walls
   !!     domain (which is what we want for lake-at-rest tests).
   !!   * The flux kernel reads the bc_type via this lookup and
   !!     dispatches inside `compute_boundary_state_impl` (inlined in
   !!     `rki_kernel_flux_unstr`) to build the right-side ghost state.
   !!   * Tracer inflow values live in `state%bc_S_ext(g)` /
   !!     `state%bc_T_ext(g)` and are consumed by
   !!     `ml_advect_tracer_unstr` on edges whose mass flux is into
   !!     the domain.
   !!
   !! Integer values match the structured side
   !! (`rki_boundary_types`): WALL=1, OPEN=2, TIDAL=3, NESTED=4,
   !! INFLOW=5, DISCHARGE=6, CLAMPED=7, SPONGE=8, CHAPMAN=9. Only the
   !! types actually dispatched here are re-declared; the rest of that
   !! range is reserved. `BC_OUTFLOW = 10` is unstructured-only and
   !! sits past the structured range to avoid future collisions.
   implicit none
   private

   public :: BC_INTERIOR
   public :: BC_WALL
   public :: BC_OPEN
   public :: BC_NESTED
   public :: BC_INFLOW
   public :: BC_DISCHARGE
   public :: BC_OUTFLOW
   public :: BC_TIDAL
   public :: BC_CLAMPED
   public :: BC_SPONGE
   public :: BC_CHAPMAN
   public :: MAX_BC_GROUPS

   integer, parameter :: BC_INTERIOR = 0
      !! Sentinel for interior edges. Never appears in `state%bc_type(:)`
      !! itself -- it's the value of `bc_group(i, k)` on edges that
      !! were not assigned to a group at mesh-build time.
   integer, parameter :: BC_WALL = 1
      !! Reflective wall. The right-side ghost state mirrors the left
      !! state across the edge normal: stage / depth / bed are the
      !! same, momentum has its normal component negated. Implemented
      !! in `compute_boundary_state`.
   integer, parameter :: BC_OPEN = 2
      !! Transmissive (zero-gradient) radiation. The right-side ghost
      !! state copies the left state, so the central-upwind solver
      !! emits the physical flux of the interior with no reflection.
      !! Stateless: uses no per-group parameters.
   integer, parameter :: BC_NESTED = 4
      !! Nested boundary: Flather radiation driven by a parent-model
      !! prescribed elevation per BC group.  At smoke level this is
      !! identical to `BC_TIDAL` — per-group `bc_eta_ext(g)` supplies
      !! the external stage and the outward-normal ghost velocity is
      !! reconstructed from the Riemann invariant.  Full parent-model
      !! coupling (live eta from a coarser grid) is deferred. Matches
      !! structured `BC_NESTED = 4`.
   integer, parameter :: BC_TIDAL = 3
      !! Prescribed free-surface elevation with Flather radiation.
      !! The per-group parameter `state%bc_eta_ext(g)` is the absolute
      !! stage (m, in the same datum as the bed) the kernel clamps the
      !! ghost stage to. Ghost outward-normal velocity is rebuilt as
      !!   u_n_ghost = u_n_int + sqrt(g/h_int)*(eta_int - eta_ext)
      !! so interior-going waves enter driven by the stage anomaly and
      !! outgoing waves leave cleanly. Update the value between steps
      !! to drive a tide. Matches structured `BC_TIDAL` semantics.
   integer, parameter :: BC_INFLOW = 5
      !! Prescribed normal velocity. The per-group parameter
      !! `state%bc_un_ext(g)` is the velocity component along the
      !! outward normal -- so a value of -2.0 means flow into the
      !! domain at 2 m/s on a boundary whose outward normal points
      !! out of the domain. The Python wrapper exposes
      !! `set_bc_inflow_speed(g, v)` which negates internally so
      !! positive `v` always means inflow. Depth at the ghost state
      !! is the extrapolated interior depth.
   integer, parameter :: BC_DISCHARGE = 6
      !! Prescribed momentum flux per unit width (m^2/s) along the
      !! outward normal. The per-group parameter `state%bc_Qn_ext(g)`
      !! carries the specific discharge; positive means *out of* the
      !! domain (along the outward normal). Ghost momentum is rebuilt
      !! from `(-bc_Qn_ext * n1, -bc_Qn_ext * n2)` so a negative Qn
      !! pumps mass in without a separate sign flip at the call site.
      !! Depth at the ghost follows the interior extrapolation. Useful
      !! when the driving quantity is a measured river flow per unit
      !! width rather than a velocity. Matches structured BC_DISCHARGE.
   integer, parameter :: BC_CLAMPED = 7
      !! Prescribed free-surface elevation with transmissive momentum
      !! (no Flather radiation term). The per-group parameter
      !! `state%bc_eta_ext(g)` is the absolute stage (m); the kernel
      !! clamps the ghost stage to it and copies interior momentum
      !! across unchanged. Pick this when you need a Dirichlet stage
      !! without wave-radiation coupling. Matches structured
      !! `BC_CLAMPED=7`.
   integer, parameter :: BC_CHAPMAN = 9
      !! Chapman / Orlanski radiating boundary. No per-group external
      !! parameter; instead the kernel keeps a per-boundary-edge
      !! persistent state `state%chapman_eta_old(i, k)` (the previous
      !! step's eta at edge `i` of cell `k`) and reconstructs the ghost
      !! stage each call from a phase-speed estimate
      !!   c = -(deta/dt) / (deta/dn_inward)
      !! clamped to `[0, sqrt(g*h_int)]`. Spatial gradient is
      !! approximated from the cell's centroid stage and the extrapolated
      !! edge-midpoint stage; temporal gradient uses `state%dt`. The
      !! computed ghost target is then dispatched through the
      !! Flather/`BC_TIDAL` algebra in `compute_boundary_state` so the
      !! ghost momentum reconstruction matches `BC_TIDAL`. Matches
      !! structured `BC_CHAPMAN=9` semantics.
   integer, parameter :: BC_SPONGE = 8
      !! Sponge layer: explicit relaxation of h and momentum toward a
      !! prescribed reference state.  Per-group parameters
      !! `state%bc_sponge_strength(g)` [1/s] and `state%bc_sponge_h_ref(g)`
      !! [m] are consumed by `apply_sponge_unstr` (separate kernel, not
      !! the flux kernel): the sponge falls through to `BC_WALL` ghost
      !! state in the flux pass, and relaxation is applied after each
      !! RK2 stage update.  Matches structured `BC_SPONGE = 8`.
   integer, parameter :: BC_OUTFLOW = 10
      !! Prescribed depth at the boundary. The per-group parameter
      !! `state%bc_h_ext(g)` is the depth (m) the kernel clamps to;
      !! momentum is transmissive. Useful as a downstream water level
      !! condition for steady-state channel runs. Unstructured-only;
      !! no structured analog, numbered past the structured range.

   integer, parameter :: MAX_BC_GROUPS = 32
      !! Maximum number of boundary groups per mesh. Plenty for any
      !! realistic coastal case (most have 2-5: open ocean, river
      !! mouths, walls). Fixed-size so the table lives inside the
      !! state struct without needing a separate `map(to:)`.

end module rki_boundary_unstr
