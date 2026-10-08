!! Tracer machinery for the unstructured multilayer solver.
module rki_ml_tracers_unstr
   !! Generic tracer pipeline (advection, distribute, RK2 save/average,
   !! clamp, BCs) plus tracer-side mixing (vertical exchange, vertical
   !! diffusion, horizontal diffusion).  All entry points take a tracer
   !! field as an actual argument (`hTr_layer`) and update it in place,
   !! so the same kernel serves salinity, temperature, and any future
   !! passive tracer.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, CP_WATER, &
                            RHO_WATER, LATENT_HEAT_VAPORIZATION
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   !use rki_boundary_unstr, only: MAX_BC_GROUPS
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR, BPG_NZ_MAX, &
                                  VANISHING_LAYER_TOL_UNSTR
   implicit none
   private

   integer, parameter :: MAX_BC_GROUPS = 32

   public :: ml_init_tracer_unstr
   public :: ml_init_tracer_stratified_unstr
   public :: ml_seed_tracer_unstr
   public :: ml_seed_standard_tracers_unstr
   public :: ml_save_tracer_unstr
   public :: ml_rk2_average_tracer_unstr
   public :: ml_distribute_tracer_unstr
   public :: ml_clamp_tracer_unstr
   public :: ml_advect_tracer_unstr
   public :: ml_advect_tracer_muscl_unstr
   public :: ml_apply_tracer_bcs_unstr
   public :: ml_apply_surface_tracer_fluxes_unstr
   public :: ml_vertical_exchange_tracer_unstr
   public :: ml_horizontal_diffusion_tracer_unstr
   public :: ml_vertical_diffusion_tracer_unstr

contains

   subroutine ml_init_tracer_unstr(state, hTr_layer, tr_init)
      !! Initialise a tracer field from a uniform scalar value:
      !!   hTr_k(cell) = h_layer_k(cell) * tr_init
      !! Called once per tracer after `ml_init_layers_unstr` has
      !! populated `h_layer`. Future phases will accept per-column
      !! profiles via a separate entry point.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — ncells varies per mesh partition;
      ! called once at init (see CLAUDE.md outer-shim + flat-impl pattern).
      real(wp), intent(inout) :: hTr_layer(:, :)
      real(wp), intent(in) :: tr_init

      integer :: cell, k, nz, n

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      !$omp parallel do collapse(2)
      do cell = 1, n
      do k = 1, nz
         hTr_layer(k, cell) = state%multilayer%h_layer(k, cell)*tr_init
      end do
      end do
      !$omp end parallel do
   end subroutine ml_init_tracer_unstr

   subroutine ml_init_tracer_stratified_unstr(state, hTr_layer, tr_bottom, tr_surface)
      !! Seed a tracer with a linear vertical profile: `tr_bottom` at k=1
      !! (bed) to `tr_surface` at k=nz (surface), per the bottom-up layer
      !! convention.  The per-column entry point anticipated by
      !! `ml_init_tracer_unstr`.  The profile is computed inline from the
      !! two scalars (no host-side per-layer array) so the seed runs
      !! correctly on the device after `enter_data`.  Caller guarantees
      !! nz > 1.
      !!   hTr_k(cell) = h_layer_k(cell) * [tr_bottom + dtr*(k - 1)]
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout) :: hTr_layer(:, :)  ! assumed-shape-ok: unstructured outer-shim; ncells varies per mesh partition
      real(wp), intent(in) :: tr_bottom, tr_surface

      integer :: cell, k, nz, n
      real(wp) :: dtr

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total
      dtr = (tr_surface - tr_bottom)/real(nz - 1, wp)

      !$omp parallel do collapse(2)
      do cell = 1, n
      do k = 1, nz
         hTr_layer(k, cell) = state%multilayer%h_layer(k, cell) &
                              *(tr_bottom + dtr*real(k - 1, wp))
      end do
      end do
      !$omp end parallel do
   end subroutine ml_init_tracer_stratified_unstr

   subroutine ml_seed_tracer_unstr(state, hTr_layer, tr_uniform, &
                                        tr_bottom, tr_surface)
      !! Seed one tracer's `hTr` field, dispatching between the uniform and
      !! the linear-stratified profile from the tracer's IC scalars.  When
      !! both `tr_bottom` and `tr_surface` are set (non-zero) and nz > 1 a
      !! linear profile is laid down from `tr_bottom` at the bed (k=1) to
      !! `tr_surface` at the surface (k=nz); otherwise every layer is filled
      !! with `tr_uniform`.  The single stratified-vs-uniform gate shared by
      !! the unstructured IC seed site (solver_unstr_init_multilayer).
      !! Mirrors the structured `ml_seed_tracer`.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout) :: hTr_layer(:, :)
      real(wp), intent(in) :: tr_uniform, tr_bottom, tr_surface

      if (tr_bottom /= 0.0_wp .and. tr_surface /= 0.0_wp .and. &
          state%multilayer%nz_ml > 1) then
         call ml_init_tracer_stratified_unstr(state, hTr_layer, tr_bottom, tr_surface)
      else
         call ml_init_tracer_unstr(state, hTr_layer, tr_uniform)
      end if

   end subroutine ml_seed_tracer_unstr

   subroutine ml_seed_standard_tracers_unstr(state)
      !! Seed the two standard tracers (salinity, temperature) from the IC
      !! scalars carried on `state%multilayer`, each dispatching uniform-vs-
      !! stratified via `ml_seed_tracer_unstr`.  Call once after
      !! `ml_init_layers_unstr` and before EOS.  The unstructured twin of the
      !! structured `multilayer_state_t%seed_tracers`; a free helper (not a
      !! type-bound method) because the unstructured tracer kernels take the
      !! top-level `unstr_state_t`, so binding on the state would close a
      !! `state <-> tracers_unstr` module cycle.
      type(unstr_state_t), intent(inout) :: state

      call ml_seed_tracer_unstr(state, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%initial_salinity, &
                                state%multilayer%S_init_bottom, state%multilayer%S_init_surface)
      call ml_seed_tracer_unstr(state, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%initial_temperature, &
                                state%multilayer%T_init_bottom, state%multilayer%T_init_surface)

   end subroutine ml_seed_standard_tracers_unstr

   subroutine ml_save_tracer_unstr(state, hTr_layer, hTr_layer0)
      !! Snapshot a tracer field at step start for SSP-RK2 averaging.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — ncells varies per mesh partition
      ! (see CLAUDE.md outer-shim + flat-impl pattern).
      real(wp), intent(in) :: hTr_layer(:, :)
      real(wp), intent(inout) :: hTr_layer0(:, :)  ! assumed-shape-ok: unstructured outer-shim; ncells varies per mesh partition

      integer :: cell, k, nz, n

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      !$omp parallel do collapse(2)
      do cell = 1, n
      do k = 1, nz
         hTr_layer0(k, cell) = hTr_layer(k, cell)
      end do
      end do
      !$omp end parallel do
   end subroutine ml_save_tracer_unstr

   subroutine ml_rk2_average_tracer_unstr(state, hTr_layer, hTr_layer0)
      !! RK2 blend for tracer fields: hTr := 0.5*(hTr + hTr^n).
      !!
      !! This element-wise blend is mass-conservative
      !! (sum_k 0.5*(hTr_layer0+hTr_layer) = 0.5*(M_0 + M_after_stage2),
      !! the SSP-RK2 of advected column mass) so the tracer total
      !! tracks the barotropic continuity exactly.  It is geometrically
      !! consistent when the layer layout doesn't change across the
      !! step, which is the steady-H regime.  Layout-mismatch issues
      !! at zsigma transition zones are addressed by initialising
      !! `h_layer` directly from vcoord (`ml_init_layers_unstr`) so
      !! step-to-step layouts already match.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — ncells varies per mesh partition
      ! (see CLAUDE.md outer-shim + flat-impl pattern).
      real(wp), intent(inout) :: hTr_layer(:, :)
      real(wp), intent(in) :: hTr_layer0(:, :)  ! assumed-shape-ok: unstructured outer-shim; ncells varies per mesh partition

      integer :: cell, k, nz, n

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      !$omp parallel do collapse(2)
      do cell = 1, n
      do k = 1, nz
         hTr_layer(k, cell) = 0.5_wp*(hTr_layer(k, cell) &
                                      + hTr_layer0(k, cell))
      end do
      end do
      !$omp end parallel do
   end subroutine ml_rk2_average_tracer_unstr

   subroutine ml_distribute_tracer_unstr(state, hTr_layer)
      !! Preserve tracer concentration across a barotropic depth change.
      !!
      !! Reads `h_layer` (pre-update, still carries `dsig_k * H_old`),
      !! computes `Tr_k = hTr_k / h_layer_k`, and writes
      !!   hTr_k := (dsig_k * H_new) * Tr_k
      !! so that when `ml_distribute_barotropic_unstr` overwrites
      !! `h_layer` with the new sigma thicknesses the concentration
      !! stays fixed.
      !!
      !! Must be called BEFORE `ml_distribute_barotropic_unstr` each
      !! RK stage, once per tracer.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — ncells varies per mesh partition
      ! (see CLAUDE.md outer-shim + flat-impl pattern).
      real(wp), intent(inout) :: hTr_layer(:, :)

      integer :: cell, k, nz, n
      real(wp) :: h_val, h_k_old, h_k_new, tr_k

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      !$omp parallel do private(k, h_val, h_k_old, h_k_new, tr_k)
      do cell = 1, n
         h_val = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
         if (h_val < 0.0_wp) h_val = 0.0_wp

         do k = 1, nz
            h_k_old = state%multilayer%h_layer(k, cell)
            h_k_new = state%multilayer%dsig_ml(k)*h_val
            if (h_k_old > DRY_TOLERANCE) then
               tr_k = hTr_layer(k, cell)/h_k_old
            else
               tr_k = 0.0_wp
            end if
            hTr_layer(k, cell) = h_k_new*tr_k
         end do
      end do
      !$omp end parallel do
   end subroutine ml_distribute_tracer_unstr

   subroutine ml_clamp_tracer_unstr(state, hTr_layer, tr_min, tr_max)
      !! Enforce physical tracer bounds [tr_min, tr_max] conservatively
      !! within each column.  Two-pass column-conservative limiter
      !! (mirror of the structured `ml_clamp_tracer`):
      !!
      !!   Pass 1 — clip violators, accumulate excess_above (above
      !!            tr_max, destroyed), deficit_below (below tr_min,
      !!            created), and per-cell spare capacity in non-
      !!            saturated layers.
      !!
      !!   Pass 2 — redistribute the net imbalance into available
      !!            capacity within the same column, proportional to
      !!            each non-saturated layer's room-to-fill (or room-
      !!            to-empty for net < 0).
      !!
      !! Per-cell tracer total preserved exactly when column-mean is in
      !! bounds.  Critical for ZSTAR_FULL deep-shallow scenarios where
      !! `ml_distribute_with_remap_full`'s `dz_old_col = remap_dz_old *
      !! scale_old` rescaling inflates per-layer S as a numerical
      !! artefact; the previous hard clip would destroy real tracer
      !! that the conservative remap had just preserved.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout), contiguous :: hTr_layer(:, :)
      real(wp), intent(in) :: tr_min, tr_max

      if (.not. state%use_multilayer) return

      call ml_clamp_tracer_unstr_impl(state%multilayer%nz_ml, state%ncells_total, &
                                      state%multilayer%h_layer, hTr_layer, tr_min, tr_max)
   end subroutine ml_clamp_tracer_unstr

   subroutine ml_clamp_tracer_unstr_impl(nz, n, h_layer, hTr_layer, tr_min, tr_max)
      !! Flat-impl of [[ml_clamp_tracer_unstr]]: explicit-shape dummies so the
      !! device kernel addresses raw pointers, not derived-type descriptors.
      integer, intent(in) :: nz, n
      real(wp), intent(in) :: h_layer(nz, n)
      real(wp), intent(inout) :: hTr_layer(nz, n)
      real(wp), intent(in) :: tr_min, tr_max

      integer :: cell, k
      real(wp) :: h_k, tr_k, room_k
      real(wp) :: excess_above, deficit_below, cap_to_add, cap_to_remove
      real(wp) :: net, frac
      real(wp) :: col_hTr, col_h
      real(wp), parameter :: TINY_CAP = 1.0e-30_wp

      !$omp parallel do private(k, h_k, tr_k, room_k, excess_above, deficit_below, cap_to_add, cap_to_remove, net, frac, col_hTr, col_h)
      do cell = 1, n

         ! Pass 0: column mean.  The two-pass clip+redistribute below is
         ! EXACTLY conservative if and only if the column-mean concentration
         ! lies within [tr_min, tr_max] (net <= cap_to_add follows from
         ! mean <= tr_max; see the structured twin `ml_clamp_tracer` for the
         ! algebra).  Outside that range no clamp can both hold the bound and
         ! conserve, and the old code silently destroyed the un-absorbable
         ! residual.  Conservation wins: skip the column.
         col_hTr = 0.0_wp
         col_h = 0.0_wp
         do k = 1, nz
            h_k = h_layer(k, cell)
            if (h_k > DRY_TOLERANCE) then
               col_hTr = col_hTr + hTr_layer(k, cell)
               col_h = col_h + h_k
            end if
         end do
         if (col_h > DRY_TOLERANCE) then
            tr_k = col_hTr/col_h
            if (tr_k > tr_max .or. tr_k < tr_min) cycle
         else
            cycle
         end if

         ! Pass 1: hard clip violators, accumulate excess/deficit and capacity.
         excess_above = 0.0_wp
         deficit_below = 0.0_wp
         cap_to_add = 0.0_wp
         cap_to_remove = 0.0_wp
         do k = 1, nz
            h_k = h_layer(k, cell)
            if (h_k > DRY_TOLERANCE) then
               tr_k = hTr_layer(k, cell)/h_k
               if (tr_k > tr_max) then
                  excess_above = excess_above + (tr_k - tr_max)*h_k
                  hTr_layer(k, cell) = tr_max*h_k
                  cap_to_remove = cap_to_remove + (tr_max - tr_min)*h_k
               else if (tr_k < tr_min) then
                  deficit_below = deficit_below + (tr_min - tr_k)*h_k
                  hTr_layer(k, cell) = tr_min*h_k
                  cap_to_add = cap_to_add + (tr_max - tr_min)*h_k
               else
                  cap_to_add = cap_to_add + (tr_max - tr_k)*h_k
                  cap_to_remove = cap_to_remove + (tr_k - tr_min)*h_k
               end if
            end if
         end do

         ! Pass 2: redistribute net imbalance into available capacity.
         net = excess_above - deficit_below
         if (net > 0.0_wp .and. cap_to_add > TINY_CAP) then
            frac = min(1.0_wp, net/cap_to_add)
            do k = 1, nz
               h_k = h_layer(k, cell)
               if (h_k > DRY_TOLERANCE) then
                  tr_k = hTr_layer(k, cell)/h_k
                  room_k = (tr_max - tr_k)*h_k
                  if (room_k > 0.0_wp) then
                     hTr_layer(k, cell) = hTr_layer(k, cell) + frac*room_k
                  end if
               end if
            end do
         else if (net < 0.0_wp .and. cap_to_remove > TINY_CAP) then
            frac = min(1.0_wp, (-net)/cap_to_remove)
            do k = 1, nz
               h_k = h_layer(k, cell)
               if (h_k > DRY_TOLERANCE) then
                  tr_k = hTr_layer(k, cell)/h_k
                  room_k = (tr_k - tr_min)*h_k
                  if (room_k > 0.0_wp) then
                     hTr_layer(k, cell) = hTr_layer(k, cell) - frac*room_k
                  end if
               end if
            end do
         end if
      end do
      !$omp end parallel do
   end subroutine ml_clamp_tracer_unstr_impl

   subroutine ml_advect_tracer_unstr(state, mesh, hTr_layer, dt, bc_tr_ext)
      !! First-order upwind horizontal tracer advection on triangles
      !! using a SHEAR-AWARE consistency-with-continuity (CWC) per-layer
      !! mass flux. Mirrors what the structured `ml_advect_tracer` does
      !! with its `0.5*(hu_layer left + right)` face flux, but reshaped
      !! for the cell-centric triangle pattern.
      !!
      !! Per-layer mass flux at an interior edge:
      !!   v_layer(k) = -(hu_edge_k * n1 + hv_edge_k * n2) * L_e
      !!     where hu_edge_k = 0.5 * (hu_layer(k, here) + hu_layer(k, nbr))
      !!   correction = mass_flux_e - sum_k v_layer(k)
      !!   f_layer(k) = v_layer(k) + dsig_k * correction
      !!
      !! The construction enforces `sum_k f_layer(k) = mass_flux_e`
      !! exactly, so per-cell tracer mass conservation matches the
      !! barotropic continuity update. The *distribution* between
      !! layers reflects the actual per-layer momentum, so a buoyant
      !! plume (surface seaward / bottom return flow) advects tracer
      !! in the right direction at every layer.
      !!
      !! For depth-uniform velocity (no baroclinic shear) `v_layer(k) =
      !! dsig_k * sum_k v_layer` and the correction term cancels the
      !! cell-avg-edge-vs-CU-flux mismatch exactly, so the scheme
      !! reduces to `f_layer(k) = dsig_k * mass_flux_e`. Uniform
      !! tracer is preserved bit-perfectly when there is no shear,
      !! which is the property `test_tracer_uniform_stays_uniform`
      !! checks (it gates on `beta_S = alpha_T = 0`).
      !!
      !! For boundary edges we fall back to `f_layer(k) = dsig_k *
      !! mass_flux_e` (no neighbour to average against; the BC kernel
      !! has already prescribed the barotropic mass flux there).
      !!
      !! Known limitation: per-cell mass conservation is exact, but
      !! per-layer hS evolution does not exactly track per-layer
      !! `dsig_k * H_new` evolution under shear. As a result, the
      !! scheme can produce small bounded over- or undershoots
      !! (typically << 1 PSU at sharp baroclinic fronts). The
      !! `S_min`/`S_max` clamp in `ml_clamp_tracer_unstr` is the
      !! safety net. The proper fix is an ALE-style remap of hS
      !! between the per-layer-advected mass and the sigma h after
      !! distribute -- not landing in this PR.
      !!
      !! Boundary handling:
      !!   - Wall edges: CU flux kernel writes mass_flux_e = 0, so no
      !!     tracer flux through walls automatically.
      !!   - Inflow edges (mass_flux_e > 0 at boundary): use bc_tr_ext
      !!     indexed by the edge's boundary group.
      !!   - Outflow edges: upwind is interior cell, handled normally.
      !!
      !! Cell-centric double-visit pattern preserved: each owned
      !! triangle writes only into its own `hTr_layer` column. The
      !! pair of visits on interior edges produces equal-and-opposite
      !! contributions because `mass_flux_e(e, K) = -mass_flux_e(e',
      !! N)` across a shared edge and both sides see the same upwind
      !! tracer value.
      !!
      !! Jacobi-vs-Gauss-Seidel note: we snapshot the input `hTr_layer`
      !! into `state%multilayer%tracer_advect_scratch` on entry and have the
      !! kernel body read from the scratch for BOTH the own-cell and
      !! neighbour-cell concentrations, writing only to `hTr_layer`.
      !! Without this, the enclosing `do concurrent` over owned cells
      !! produces correct results on GPU (launch-time Jacobi semantics)
      !! but biased results on CPU, where `do concurrent` collapses to
      !! a sequential / OMP loop and later cell iterations observe
      !! already-updated neighbour values. The result is a drift of
      !! ~3e-6 per step on a uniform tracer that compounds over time.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      ! assumed-shape-ok: unstructured outer-shim — ncells varies per mesh partition
      ! (see CLAUDE.md outer-shim + flat-impl pattern).
      real(wp), intent(inout) :: hTr_layer(:, :)
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: bc_tr_ext(MAX_BC_GROUPS)
         !! External tracer concentration for each boundary group.
         !! Used for inflow BCs. If not provided, inflow boundaries
         !! use zero tracer concentration.

      integer :: cell, k, nz, n
      real(wp) :: bc_ext_local(MAX_BC_GROUPS)

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      ! Copy bc_tr_ext to a local array so it can be passed to the
      ! pure per-cell kernel. Default to zero if not provided.
      if (present(bc_tr_ext)) then
         bc_ext_local = bc_tr_ext
      else
         bc_ext_local = 0.0_wp
      end if

      ! Snapshot the input tracer into the scratch buffer so the per-cell
      ! kernel can read neighbour concentrations from a frozen image and
      ! write the update into the live buffer without a race.
      !$omp parallel do collapse(2)
      do cell = 1, n
      do k = 1, nz
         state%multilayer%tracer_advect_scratch(k, cell) = hTr_layer(k, cell)
      end do
      end do
      !$omp end parallel do

      !$omp parallel do
      do cell = 1, mesh%ncells_owned
         call advect_tracer_one_cell(state, mesh, hTr_layer, nz, dt, cell, &
                                     bc_ext_local)
      end do
      !$omp end parallel do
   end subroutine ml_advect_tracer_unstr

   pure subroutine advect_tracer_one_cell(state, mesh, hTr_layer, nz, dt, cell, &
                                          bc_tr_ext)
      !! Per-cell body, marked pure to document the "this cell only"
      !! write pattern that makes the enclosing `do concurrent` safe.
      !!
      !! Per-layer mass flux uses the vanishing-layer-gated centred CWC
      !! split (mirrors the structured `ml_advect_tracer`):
      !!
      !!   ! For each edge (per-edge, both cells see the same value):
      !!   if h_here(k) > tol .and. h_nbr(k) > tol then
      !!      v_layer(k, e)   = -(hu_edge*n1 + hv_edge*n2) * L_e
      !!         where hu_edge = 0.5*(hu_layer(here) + hu_layer(nbr))
      !!      hf_layer(k, e)  = 0.5*(h_layer(here) + h_layer(nbr))
      !!   else
      !!      v_layer(k, e)   = 0     ! vanishing layer carries no flux
      !!      hf_layer(k, e)  = 0
      !!   end if
      !!   H_active(e) = sum_k hf_layer(k, e)
      !!   correction(e) = mass_flux_e - sum_k v_layer(k, e)
      !!   frac_e(k, e) = hf_layer(k, e) / H_active(e)   ! 0 if vanishing
      !!   f_layer(k, e) = v_layer(k, e) + frac_e(k, e) * correction(e)
      !!
      !! `v_layer`, `hf_layer`, and `frac_e` are all symmetric across
      !! the edge → per-layer face mass conservation holds (the
      !! contribution that this cell sees through edge e is equal-and-
      !! opposite to what its neighbour sees through the shared edge).
      !! `sum_k f_layer = mass_flux_e` (CWC).
      !!
      !! For sigma / zsigma / zstar-lite where every layer is non-
      !! vanishing on both sides, frac_e(k) reduces to dsig_k for
      !! uniform sigma, and the formula degenerates to the original
      !! centred CWC scheme — zero risk to existing tests.
      !!
      !! For ZSTAR_FULL deep-shallow edges, vanishing bed layers in
      !! the shallow cell are excluded from the per-layer split, so the
      !! column flux is redistributed across the (active) non-vanishing
      !! layers instead of being centred-averaged into the tiny
      !! vanishing receivers.  No bloat of S in vanishing layers, no
      !! clamp clipping, no salt-mass leak.
      !!
      !! Boundary edges have no neighbour to average against, so we
      !! fall back to `dsig_k * mass_flux_e` there (CWC preserved
      !! trivially, BC-prescribed mass flux honoured as-is).
      !!
      !! Neighbour and own-cell tracer reads go through
      !! `state%multilayer%tracer_advect_scratch` (the frozen copy of `hTr_layer`
      !! taken at the start of the advect step) so CPU serial loops
      !! observe the same Jacobi update GPU parallel launches do.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(inout) :: hTr_layer(:, :)
      integer, intent(in) :: nz, cell
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: bc_tr_ext(MAX_BC_GROUPS)
         !! External tracer for each boundary group

      integer :: k, e, nbr, bc_grp
      real(wp) :: inv_area
      real(wp) :: f_mass, f_layer
      real(wp) :: h_here, h_nbr, tr_up
      real(wp) :: accum
      real(wp) :: n1, n2, L_e
      real(wp) :: hu_edge, hv_edge
      real(wp) :: v_layer(BPG_NZ_MAX, 3)   !! Per-layer per-edge raw flux
      real(wp) :: hf_layer(BPG_NZ_MAX, 3)  !! Per-layer per-edge active thickness
      real(wp) :: H_active(3)               !! sum_k hf_layer(k, e)
      real(wp) :: inv_H_active(3)           !! 1 / H_active where positive
      real(wp) :: correction(3)             !! mass_flux_e - sum_k v_layer
      real(wp) :: frac(BPG_NZ_MAX)          !! Column-local layer fraction (boundary fallback)
      real(wp) :: h_total, inv_h_total

      inv_area = 1.0_wp/mesh%area(cell)

      ! Column-local layer fraction h_layer(k)/H_total — used for the
      ! BOUNDARY-edge fallback only.  For interior edges algorithm N
      ! uses per-edge `frac_e = hf_layer(k, e) / H_active(e)`.
      h_total = 0.0_wp
      do k = 1, nz
         h_total = h_total + state%multilayer%h_layer(k, cell)
      end do
      if (h_total > DRY_TOLERANCE) then
         inv_h_total = 1.0_wp/h_total
         do k = 1, nz
            frac(k) = state%multilayer%h_layer(k, cell)*inv_h_total
         end do
      else
         do k = 1, nz
            frac(k) = state%multilayer%dsig_ml(k)
         end do
      end if

      ! ---- Pass 1: build per-edge per-layer flux + active thickness +
      !              the CWC correction (vanishing-layer-gated centred). ----
      do e = 1, 3
         nbr = mesh%neighbours(e, cell)
         f_mass = state%barotropic%mass_flux_e(e, cell)

         if (nbr > 0) then
            n1 = mesh%edge_n1(e, cell)
            n2 = mesh%edge_n2(e, cell)
            L_e = mesh%edge_len(e, cell)
            correction(e) = f_mass
            H_active(e) = 0.0_wp
            do k = 1, nz
               h_here = state%multilayer%h_layer(k, cell)
               h_nbr = state%multilayer%h_layer(k, nbr)
               if (h_here > VANISHING_LAYER_TOL_UNSTR .and. h_nbr > VANISHING_LAYER_TOL_UNSTR) then
                  hu_edge = 0.5_wp*(state%multilayer%hu_layer(k, cell) + state%multilayer%hu_layer(k, nbr))
                  hv_edge = 0.5_wp*(state%multilayer%hv_layer(k, cell) + state%multilayer%hv_layer(k, nbr))
                  ! mass_flux_e convention: > 0 = inflow. (hu*n1 + hv*n2)*L
                  ! is OUTFLOW for outward normal n, so negate.
                  v_layer(k, e) = -(hu_edge*n1 + hv_edge*n2)*L_e
                  hf_layer(k, e) = 0.5_wp*(h_here + h_nbr)
               else
                  v_layer(k, e) = 0.0_wp
                  hf_layer(k, e) = 0.0_wp
               end if
               correction(e) = correction(e) - v_layer(k, e)
               H_active(e) = H_active(e) + hf_layer(k, e)
            end do
            if (H_active(e) > 0.0_wp) then
               inv_H_active(e) = 1.0_wp/H_active(e)
            else
               inv_H_active(e) = 0.0_wp
            end if
         else
            ! Boundary edge: no neighbour. Split mass flux by the local
            ! layer-thickness fraction so f_layer reflects the actual
            ! vertical layout (matches sigma when frac == dsig_ml).
            correction(e) = 0.0_wp
            H_active(e) = 0.0_wp
            inv_H_active(e) = 0.0_wp
            do k = 1, nz
               v_layer(k, e) = frac(k)*f_mass
               hf_layer(k, e) = 0.0_wp
            end do
         end if
      end do

      ! ---- Pass 2: per-layer upwind tracer flux accumulation ----
      do k = 1, nz
         accum = 0.0_wp
         h_here = state%multilayer%h_layer(k, cell)

         do e = 1, 3
            nbr = mesh%neighbours(e, cell)

            if (nbr > 0) then
               f_layer = v_layer(k, e) + hf_layer(k, e)*inv_H_active(e)*correction(e)
            else
               ! Boundary edge: v_layer already carries the full
               ! frac(k)*f_mass split (no correction term needed).
               f_layer = v_layer(k, e)
            end if

            ! Pick upwind concentration. If the layer is losing mass
            ! (`f_layer < 0` = outflow) the upwind is this cell;
            ! otherwise it's the neighbour.
            if (f_layer < 0.0_wp) then
               ! Outflow: upwind is this cell
               if (h_here > DRY_TOLERANCE) then
                  tr_up = state%multilayer%tracer_advect_scratch(k, cell)/h_here
               else
                  tr_up = 0.0_wp
               end if
            else if (f_layer > 0.0_wp) then
               ! Inflow: upwind is neighbour or boundary
               if (nbr > 0) then
                  ! Interior edge: use neighbour value
                  h_nbr = state%multilayer%h_layer(k, nbr)
                  if (h_nbr > DRY_TOLERANCE) then
                     tr_up = state%multilayer%tracer_advect_scratch(k, nbr)/h_nbr
                  else
                     tr_up = 0.0_wp
                  end if
               else
                  ! Boundary inflow: use external tracer from BC group
                  bc_grp = mesh%bc_group(e, cell)
                  if (bc_grp > 0 .and. bc_grp <= MAX_BC_GROUPS) then
                     tr_up = bc_tr_ext(bc_grp)
                  else
                     ! Fallback for undefined BC group
                     tr_up = 0.0_wp
                  end if
               end if
            else
               cycle
            end if

            ! hTr_K += dt/area * Tr_up * F_layer_inflow
            ! (positive f_layer = inflow, so sign of accum matches
            ! the mass change through this edge).
            accum = accum + tr_up*f_layer
         end do

         hTr_layer(k, cell) = hTr_layer(k, cell) + dt*inv_area*accum
      end do
   end subroutine advect_tracer_one_cell

   subroutine ml_advect_tracer_muscl_unstr(state, mesh, hTr_layer, dt, bc_tr_ext)
      !! Second-order MUSCL horizontal tracer advection on triangles.
      !!
      !! Uses the same CWC-consistent per-layer mass flux as the
      !! first-order scheme, but replaces the piecewise-constant upwind
      !! concentration with a linearly reconstructed value:
      !!
      !!   Tr_edge = Tr_cell + phi * grad(Tr) . (x_edge - x_centroid)
      !!
      !! where grad(Tr) is the Green-Gauss gradient on the triangle and
      !! phi is the Barth-Jespersen limiter that prevents new extrema.
      !!
      !! Two-phase approach:
      !!   Phase 1 — Gradient + limiter: compute limited gradient for
      !!             every cell (owned + halo) so that Phase 2 can read
      !!             neighbour gradients.
      !!   Phase 2 — MUSCL advection: same CWC flux construction as the
      !!             first-order kernel, but upwind concentration is the
      !!             limited edge-extrapolated value.
      !!
      !! Cell-centric double-visit pattern preserved: each owned cell
      !! writes only to its own hTr_layer column.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(inout), contiguous :: hTr_layer(:, :)
      real(wp), intent(in) :: dt
      real(wp), intent(in), optional :: bc_tr_ext(MAX_BC_GROUPS)

      real(wp) :: bc_ext_local(MAX_BC_GROUPS)

      if (.not. state%use_multilayer) return

      if (present(bc_tr_ext)) then
         bc_ext_local = bc_tr_ext
      else
         bc_ext_local = 0.0_wp
      end if

      ! Outer shim: forward the state/mesh array components to the flat _impl
      ! as explicit-shape dummies, so the device kernels address contiguous raw
      ! pointers instead of derived-type descriptors.  Descriptor access drags
      ! base+strides+bounds into registers on every reference (144 regs/thread
      ! -> 18.75% occupancy on this kernel); explicit-shape dummies cut that.
      call ml_advect_tracer_muscl_unstr_impl( &
         state%multilayer%nz_ml, state%ncells_total, mesh%ncells_owned, &
         state%multilayer%h_layer, state%multilayer%hu_layer, &
         state%multilayer%hv_layer, state%multilayer%dsig_ml, &
         state%multilayer%tracer_advect_scratch, &
         state%multilayer%ml_lgrad_x, state%multilayer%ml_lgrad_y, &
         state%barotropic%mass_flux_e, mesh%area, mesh%neighbours, &
         mesh%surrogate_neighbours, mesh%edge_n1, mesh%edge_n2, &
         mesh%edge_len, mesh%edge_mid_x, mesh%edge_mid_y, &
         mesh%centroid_x, mesh%centroid_y, mesh%bc_group, &
         hTr_layer, dt, bc_ext_local)
   end subroutine ml_advect_tracer_muscl_unstr

   subroutine ml_advect_tracer_muscl_unstr_impl(nz, n, ncells_owned, &
                                                     h_layer, hu_layer, hv_layer, &
                                                     dsig_ml, tracer_advect_scratch, &
                                                     lgrad_x, lgrad_y, mass_flux_e, &
                                                     area, neighbours, surrogate_neighbours, &
                                                     edge_n1, edge_n2, edge_len, &
                                                     edge_mid_x, edge_mid_y, &
                                                     centroid_x, centroid_y, bc_group, &
                                                     hTr_layer, dt, bc_ext_local)
      !! Flat-impl of [[ml_advect_tracer_muscl_unstr]]: explicit-shape dummies
      !! so the device kernels see contiguous raw pointers, not derived-type
      !! descriptors (descriptor access inflates per-thread registers and caps
      !! occupancy).  Body is identical to the historical kernel modulo the
      !! stripped `state%`/`mesh%` prefixes and the removed `associate`.
      integer, intent(in) :: nz, n, ncells_owned
      real(wp), intent(in) :: h_layer(nz, n), hu_layer(nz, n), hv_layer(nz, n)
      real(wp), intent(in) :: dsig_ml(nz)
      real(wp), intent(inout) :: tracer_advect_scratch(nz, n)
      real(wp), intent(inout) :: lgrad_x(nz, n), lgrad_y(nz, n)
      real(wp), intent(in) :: mass_flux_e(3, n)
      real(wp), intent(in) :: area(n)
      integer, intent(in) :: neighbours(3, n), surrogate_neighbours(3, n)
      real(wp), intent(in) :: edge_n1(3, n), edge_n2(3, n), edge_len(3, n)
      real(wp), intent(in) :: edge_mid_x(3, n), edge_mid_y(3, n)
      real(wp), intent(in) :: centroid_x(n), centroid_y(n)
      integer, intent(in) :: bc_group(3, n)
      real(wp), intent(inout) :: hTr_layer(nz, n)
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: bc_ext_local(MAX_BC_GROUPS)

      integer :: cell, k, e, nbr
      real(wp) :: tr_here, tr_nbr, tr_edge_val, tr_max, tr_min
      real(wp) :: gradx, grady, inv_area
      real(wp) :: dq_e, dx_e, dy_e, phi
      real(wp) :: n1, n2, L_e
      real(wp) :: f_mass, f_layer, h_here, h_nbr, tr_up, accum
      real(wp) :: hu_edge, hv_edge
      real(wp) :: v_layer(BPG_NZ_MAX, 3)
      real(wp) :: hf_layer(BPG_NZ_MAX, 3)   !! Per-layer per-edge active thickness
      real(wp) :: H_active(3)               !! sum_k hf_layer(k, e)
      real(wp) :: inv_H_active(3)           !! 1 / H_active where positive
      real(wp) :: correction(3)
      real(wp) :: frac(BPG_NZ_MAX)          !! Column-local layer fraction (boundary fallback)
      real(wp) :: h_total, inv_h_total
      integer :: bc_grp

      ! Snapshot the input tracer for Jacobi consistency.
      !$omp parallel do collapse(2)
      do cell = 1, n
      do k = 1, nz
         tracer_advect_scratch(k, cell) = hTr_layer(k, cell)
      end do
      end do
      !$omp end parallel do

      ! ---- Phase 1: Green-Gauss gradient + Barth-Jespersen limiter ----
      ! Computed for ALL cells (including halo) so Phase 2 can read
      ! the neighbour's gradient for inflow faces.
      !$omp parallel do private(k, e, nbr, tr_here, tr_nbr, tr_edge_val, tr_max, tr_min, gradx, grady, inv_area, n1, n2, L_e, dq_e, dx_e, dy_e, phi)
      do cell = 1, n
         inv_area = 1.0_wp/area(cell)

         do k = 1, nz
            h_here = h_layer(k, cell)
            if (h_here > DRY_TOLERANCE) then
               tr_here = tracer_advect_scratch(k, cell)/h_here
            else
               tr_here = 0.0_wp
            end if

            ! Green-Gauss gradient using surrogate neighbours
            gradx = 0.0_wp
            grady = 0.0_wp
            tr_max = tr_here
            tr_min = tr_here

            do e = 1, 3
               nbr = surrogate_neighbours(e, cell)
               h_nbr = h_layer(k, nbr)
               if (h_nbr > DRY_TOLERANCE) then
                  tr_nbr = tracer_advect_scratch(k, nbr)/h_nbr
               else
                  tr_nbr = 0.0_wp
               end if

               tr_edge_val = 0.5_wp*(tr_here + tr_nbr)
               n1 = edge_n1(e, cell)
               n2 = edge_n2(e, cell)
               L_e = edge_len(e, cell)
               gradx = gradx + tr_edge_val*n1*L_e
               grady = grady + tr_edge_val*n2*L_e

               if (tr_nbr > tr_max) tr_max = tr_nbr
               if (tr_nbr < tr_min) tr_min = tr_nbr
            end do

            gradx = gradx*inv_area
            grady = grady*inv_area

            ! Barth-Jespersen limiter: find the minimum phi over all
            ! edge midpoints such that the extrapolated value stays
            ! within [tr_min, tr_max].
            phi = 1.0_wp

            do e = 1, 3
               dx_e = edge_mid_x(e, cell) - centroid_x(cell)
               dy_e = edge_mid_y(e, cell) - centroid_y(cell)
               dq_e = gradx*dx_e + grady*dy_e

               if (dq_e > 1.0e-15_wp) then
                  phi = min(phi, (tr_max - tr_here)/dq_e)
               else if (dq_e < -1.0e-15_wp) then
                  phi = min(phi, (tr_min - tr_here)/dq_e)
               end if
            end do

            ! Clamp phi to [0, 1]
            if (phi < 0.0_wp) phi = 0.0_wp
            if (phi > 1.0_wp) phi = 1.0_wp

            ! Store limited gradient (pre-multiplied by phi)
            lgrad_x(k, cell) = phi*gradx
            lgrad_y(k, cell) = phi*grady
         end do
      end do
      !$omp end parallel do

      ! ---- Phase 2: MUSCL advection (same CWC flux + limited upwind) ----
      !$omp parallel do private(k, e, nbr, bc_grp, inv_area, f_mass, f_layer, h_here, h_nbr, tr_up, accum, n1, n2, L_e, hu_edge, hv_edge, dx_e, dy_e, tr_here, tr_nbr, v_layer, hf_layer, H_active, inv_H_active, correction, frac, h_total, inv_h_total)
      do cell = 1, ncells_owned
         inv_area = 1.0_wp/area(cell)

         ! Column-local layer fraction (used only for the boundary-edge
         ! fallback below).  Interior edges use per-edge active frac.
         h_total = 0.0_wp
         do k = 1, nz
            h_total = h_total + h_layer(k, cell)
         end do
         if (h_total > DRY_TOLERANCE) then
            inv_h_total = 1.0_wp/h_total
            do k = 1, nz
               frac(k) = h_layer(k, cell)*inv_h_total
            end do
         else
            do k = 1, nz
               frac(k) = dsig_ml(k)
            end do
         end if

         ! ---- Pass 1: vanishing-layer-gated centred CWC flux per edge ----
         do e = 1, 3
            nbr = neighbours(e, cell)
            f_mass = mass_flux_e(e, cell)

            if (nbr > 0) then
               n1 = edge_n1(e, cell)
               n2 = edge_n2(e, cell)
               L_e = edge_len(e, cell)
               correction(e) = f_mass
               H_active(e) = 0.0_wp
               do k = 1, nz
                  h_here = h_layer(k, cell)
                  h_nbr = h_layer(k, nbr)
                  if (h_here > VANISHING_LAYER_TOL_UNSTR .and. h_nbr > VANISHING_LAYER_TOL_UNSTR) then
                     hu_edge = 0.5_wp*(hu_layer(k, cell) + hu_layer(k, nbr))
                     hv_edge = 0.5_wp*(hv_layer(k, cell) + hv_layer(k, nbr))
                     v_layer(k, e) = -(hu_edge*n1 + hv_edge*n2)*L_e
                     hf_layer(k, e) = 0.5_wp*(h_here + h_nbr)
                  else
                     v_layer(k, e) = 0.0_wp
                     hf_layer(k, e) = 0.0_wp
                  end if
                  correction(e) = correction(e) - v_layer(k, e)
                  H_active(e) = H_active(e) + hf_layer(k, e)
               end do
               if (H_active(e) > 0.0_wp) then
                  inv_H_active(e) = 1.0_wp/H_active(e)
               else
                  inv_H_active(e) = 0.0_wp
               end if
            else
               correction(e) = 0.0_wp
               H_active(e) = 0.0_wp
               inv_H_active(e) = 0.0_wp
               do k = 1, nz
                  v_layer(k, e) = frac(k)*f_mass
                  hf_layer(k, e) = 0.0_wp
               end do
            end if
         end do

         ! ---- Pass 2: per-layer MUSCL upwind flux accumulation ----
         do k = 1, nz
            accum = 0.0_wp
            h_here = h_layer(k, cell)
            if (h_here > DRY_TOLERANCE) then
               tr_here = tracer_advect_scratch(k, cell)/h_here
            else
               tr_here = 0.0_wp
            end if

            do e = 1, 3
               nbr = neighbours(e, cell)
               if (nbr > 0) then
                  f_layer = v_layer(k, e) + hf_layer(k, e)*inv_H_active(e)*correction(e)
               else
                  f_layer = v_layer(k, e)
               end if

               if (f_layer < 0.0_wp) then
                  ! Outflow: upwind is this cell, use limited extrapolation
                  if (h_here > DRY_TOLERANCE) then
                     dx_e = edge_mid_x(e, cell) - centroid_x(cell)
                     dy_e = edge_mid_y(e, cell) - centroid_y(cell)
                     tr_up = tr_here + lgrad_x(k, cell)*dx_e + lgrad_y(k, cell)*dy_e
                  else
                     tr_up = 0.0_wp
                  end if
               else if (f_layer > 0.0_wp) then
                  ! Inflow: upwind is neighbour or boundary
                  if (nbr > 0) then
                     h_nbr = h_layer(k, nbr)
                     if (h_nbr > DRY_TOLERANCE) then
                        tr_nbr = tracer_advect_scratch(k, nbr)/h_nbr
                        ! Extrapolate from neighbour's centroid to this edge midpoint
                        dx_e = edge_mid_x(e, cell) - centroid_x(nbr)
                        dy_e = edge_mid_y(e, cell) - centroid_y(nbr)
                        tr_up = tr_nbr + lgrad_x(k, nbr)*dx_e + lgrad_y(k, nbr)*dy_e
                     else
                        tr_up = 0.0_wp
                     end if
                  else
                     bc_grp = bc_group(e, cell)
                     if (bc_grp > 0 .and. bc_grp <= MAX_BC_GROUPS) then
                        tr_up = bc_ext_local(bc_grp)
                     else
                        tr_up = 0.0_wp
                     end if
                  end if
               else
                  cycle
               end if

               accum = accum + tr_up*f_layer
            end do

            hTr_layer(k, cell) = hTr_layer(k, cell) + dt*inv_area*accum
         end do
      end do
      !$omp end parallel do
   end subroutine ml_advect_tracer_muscl_unstr_impl

   subroutine ml_apply_tracer_bcs_unstr(state, mesh, hTr_layer, bc_tr_ext)
      !! Apply tracer boundary conditions for inflow edges.
      !!
      !! For boundary edges where mass is flowing into the domain (inflow),
      !! this routine sets the tracer concentration to bc_tr_ext. The value
      !! is written into `tracer_advect_scratch` so that `ml_advect_tracer_unstr`
      !! sees the correct upwind value at inflow boundaries.
      !!
      !! Must be called AFTER `ml_advect_tracer_unstr` has snapshotted
      !! `hTr_layer` into `tracer_advect_scratch`, but BEFORE the per-cell
      !! advection loop runs. Since `ml_advect_tracer_unstr` does both in
      !! sequence, this routine should be called between the snapshot and
      !! the advection, or the advection kernel needs to be modified to
      !! call this internally.
      !!
      !! Alternative approach (implemented here): For each boundary edge
      !! with inflow (mass_flux_e > 0), we don't need to modify the scratch
      !! — instead we can post-process the advection to add the missing
      !! boundary contribution. But this is tricky with the current kernel.
      !!
      !! Current implementation: mark boundary edges in state and let the
      !! advection kernel check `bc_tr_ext` directly. This requires modifying
      !! `advect_tracer_one_cell` to accept bc_tr_ext and check for boundary
      !! inflow.
      !!
      !! For now, this routine is a placeholder that documents the interface.
      !! The actual BC application is integrated into the advection kernel.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(inout) :: hTr_layer(:, :)
      real(wp), intent(in) :: bc_tr_ext(MAX_BC_GROUPS)

      ! Placeholder — BCs are handled inline in the advection kernel
      ! by checking mesh%neighbours(e, cell) < 0 and state%barotropic%mass_flux_e > 0
      ! for inflow boundaries.
      !
      ! This routine exists so the solver can call a uniform interface.
      ! Future: implement open-boundary handling here if needed.

   end subroutine ml_apply_tracer_bcs_unstr

   subroutine ml_vertical_exchange_tracer_unstr(state, hTr_layer, dt)
      !! Vertical tracer exchange using upwind advection at interfaces.
      !!
      !! Same pattern as momentum exchange: upwind selection based on w,
      !! 50% limiter on donor layer.
      !!
      !! Uses two-pass approach: first compute all interface fluxes from
      !! the current state, then apply them. This avoids the race condition
      !! where modifying layer k affects the flux computation at layer k+1.
      !!
      !! Requires w_interface to be populated by ml_compute_vertical_velocity_unstr.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout), contiguous :: hTr_layer(:, :)
      real(wp), intent(in) :: dt

      if (.not. state%use_multilayer) return
      ! Single-layer case: no vertical exchange
      if (state%multilayer%nz_ml < 2) return

      call ml_vertical_exchange_tracer_unstr_impl( &
         state%multilayer%nz_ml, state%ncells_total, state%ncells_owned, &
         state%barotropic%height_cv, state%multilayer%h_layer, &
         state%multilayer%w_interface, hTr_layer, dt)
   end subroutine ml_vertical_exchange_tracer_unstr

   subroutine ml_vertical_exchange_tracer_unstr_impl(nz, n, ncells_owned, height_cv, &
                                                          h_layer, w_interface, hTr_layer, dt)
      !! Flat-impl of [[ml_vertical_exchange_tracer_unstr]]: explicit-shape dummies
      !! so the device kernel addresses raw pointers, not derived-type descriptors.
      integer, intent(in) :: nz, n, ncells_owned
      real(wp), intent(in) :: height_cv(n)
      real(wp), intent(in) :: h_layer(nz, n)
      real(wp), intent(in) :: w_interface(nz + 1, n)
      real(wp), intent(inout) :: hTr_layer(nz, n)
      real(wp), intent(in) :: dt

      integer :: cell, k
      real(wp) :: w_face, tr_up, flux_tr, h_val
      real(wp) :: flux_iface(BPG_NZ_MAX)

      !$omp parallel do private(k, w_face, tr_up, flux_tr, h_val, flux_iface)
      do cell = 1, ncells_owned
         ! Skip thin columns
         if (height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         ! Pass 1: compute all interface fluxes from current state
         do k = 1, nz - 1
            w_face = w_interface(k + 1, cell)

            ! Upwind tracer selection
            if (w_face > 0.0_wp) then
               ! Upward flow: donor is layer k
               h_val = h_layer(k, cell)
               if (h_val > DRY_TOLERANCE) then
                  tr_up = hTr_layer(k, cell)/h_val
               else
                  tr_up = 0.0_wp
               end if
            else
               ! Downward flow: donor is layer k+1
               h_val = h_layer(k + 1, cell)
               if (h_val > DRY_TOLERANCE) then
                  tr_up = hTr_layer(k + 1, cell)/h_val
               else
                  tr_up = 0.0_wp
               end if
            end if

            ! Flux = w * tracer_concentration * dt
            flux_tr = w_face*tr_up*dt

            ! 50% limiter
            if (w_face > 0.0_wp) then
               if (abs(flux_tr) > 0.5_wp*abs(hTr_layer(k, cell))) then
                  flux_tr = sign(0.5_wp*abs(hTr_layer(k, cell)), flux_tr)
               end if
            else
               if (abs(flux_tr) > 0.5_wp*abs(hTr_layer(k + 1, cell))) then
                  flux_tr = sign(0.5_wp*abs(hTr_layer(k + 1, cell)), flux_tr)
               end if
            end if

            flux_iface(k) = flux_tr
         end do

         ! Pass 2: apply fluxes to layers
         ! Bottom layer (k=1): only has the interface ABOVE it
         hTr_layer(1, cell) = hTr_layer(1, cell) - flux_iface(1)
         ! Interior layers: interface below and above
         do k = 2, nz - 1
            hTr_layer(k, cell) = hTr_layer(k, cell) + flux_iface(k - 1) - flux_iface(k)
         end do
         ! Top / surface layer (k=nz): only has the interface BELOW it
         hTr_layer(nz, cell) = hTr_layer(nz, cell) + flux_iface(nz - 1)
      end do
      !$omp end parallel do
   end subroutine ml_vertical_exchange_tracer_unstr_impl

   subroutine ml_horizontal_diffusion_tracer_unstr(state, mesh, hTr_layer, kappa, dt)
      !! Horizontal Laplacian diffusion of a scalar tracer on triangles
      !! in flux form, conservative per-edge under non-uniform layer
      !! thicknesses.
      !!
      !! Per edge the diffusive flux is
      !!   flux_edge = h_face * kappa_eff * (Tr_nbr - Tr_self) * L_e / d_cn
      !! where `h_face = 0.5 * (h_self + h_nbr)` is symmetric across the
      !! edge.  Accumulated as `d(hTr)/dt = (1/A) * sum_e flux_edge`, so
      !! what cell A gains from edge A-B equals what cell B loses
      !! through the shared edge — mass conservation holds even when
      !! the per-cell layer thicknesses differ across the edge.
      !!
      !! History: the previous formulation computed `lap_tr = (1/A) *
      !! sum_e (Tr_nbr - Tr_self) * L_e / d_cn` (per-edge concentration
      !! gradient with no `h_face` weighting) and updated `hTr_new =
      !! h_self * (Tr_self + kappa·dt·lap_tr)`.  Each cell used its
      !! OWN `h` on the RHS, so adjacent cells with different
      !! `h_layer(k)` (zsigma / zstar at bathymetric transitions) saw
      !! different mass changes through the same shared edge, leaking
      !! total mass.  Currently dormant in production runs
      !! (`hdiff_kappa = 0` in the LE / Sydney Harbour configs) but
      !! would bite the moment any user enables it.
      !!
      !! Boundary edges (surrogate_neighbours == self) and dry-neighbour
      !! edges contribute zero (closed face), preserving the symmetric
      !! per-face conservation property.
      !!
      !! Jacobi snapshot: `hTr_pre` freezes the input so CPU-sequential
      !! and GPU-parallel builds give bit-identical results.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(inout), contiguous :: hTr_layer(:, :)
      real(wp), intent(in) :: kappa
      real(wp), intent(in) :: dt

      if (kappa <= 0.0_wp) return
      if (.not. state%use_multilayer) return

      call ml_horizontal_diffusion_tracer_unstr_impl( &
         state%multilayer%nz_ml, state%ncells_total, state%ncells_owned, &
         state%barotropic%height_cv, state%multilayer%h_layer, &
         state%multilayer%ml_hdiff_hTr_pre, mesh%area, mesh%radius, &
         mesh%centroid_x, mesh%centroid_y, mesh%edge_len, &
         mesh%surrogate_neighbours, hTr_layer, kappa, dt)
   end subroutine ml_horizontal_diffusion_tracer_unstr

   subroutine ml_horizontal_diffusion_tracer_unstr_impl(nz, n, ncells_owned, height_cv, &
                                                             h_layer, hTr_pre, area, radius, &
                                                             centroid_x, centroid_y, edge_len, &
                                                             surrogate_neighbours, hTr_layer, &
                                                             kappa, dt)
      !! Flat-impl of [[ml_horizontal_diffusion_tracer_unstr]]: explicit-shape
      !! dummies so the device kernels address raw pointers, not derived-type
      !! descriptors.  `hTr_pre` is the persistent Jacobi-snapshot workspace.
      integer, intent(in) :: nz, n, ncells_owned
      real(wp), intent(in) :: height_cv(n)
      real(wp), intent(in) :: h_layer(nz, n)
      real(wp), intent(inout) :: hTr_pre(nz, n)
      real(wp), intent(in) :: area(n), radius(n)
      real(wp), intent(in) :: centroid_x(n), centroid_y(n)
      real(wp), intent(in) :: edge_len(3, n)
      integer, intent(in) :: surrogate_neighbours(3, n)
      real(wp), intent(inout) :: hTr_layer(nz, n)
      real(wp), intent(in) :: kappa, dt

      integer :: cell, k, e, nbr
      real(wp) :: kappa_eff, tr_c, tr_nbr, h_k, h_nbr, inv_area
      real(wp) :: dx_cn, dy_cn, d_cn, flux_accum

      !$omp parallel do collapse(2)
      do k = 1, nz
      do cell = 1, n
         hTr_pre(k, cell) = hTr_layer(k, cell)
      end do
      end do
      !$omp end parallel do

      !$omp parallel do private(k, e, nbr, kappa_eff, tr_c, tr_nbr, h_k, h_nbr, inv_area, dx_cn, dy_cn, d_cn, flux_accum)
      do cell = 1, ncells_owned

         if (height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         inv_area = 1.0_wp/area(cell)
         ! Per-cell stability limit: kappa * dt / radius^2 <= 0.125
         kappa_eff = min(kappa, 0.125_wp*radius(cell)**2/(dt + 1.0e-30_wp))

         do k = 1, nz
            h_k = h_layer(k, cell)
            if (h_k < THIN_LAYER_THRESHOLD) cycle

            tr_c = hTr_pre(k, cell)/h_k

            flux_accum = 0.0_wp
            do e = 1, 3
               nbr = surrogate_neighbours(e, cell)
               if (nbr == cell) cycle  ! boundary edge — zero-flux wall

               h_nbr = h_layer(k, nbr)
               if (h_nbr < THIN_LAYER_THRESHOLD) cycle  ! dry-nbr — zero-flux wall

               tr_nbr = hTr_pre(k, nbr)/h_nbr

               dx_cn = centroid_x(nbr) - centroid_x(cell)
               dy_cn = centroid_y(nbr) - centroid_y(cell)
               d_cn = sqrt(dx_cn*dx_cn + dy_cn*dy_cn)

               flux_accum = flux_accum + 0.5_wp*(h_k + h_nbr) &
                            *(tr_nbr - tr_c)*edge_len(e, cell)/d_cn
            end do

            hTr_layer(k, cell) = hTr_layer(k, cell) &
                                 + kappa_eff*dt*flux_accum*inv_area
         end do
      end do
      !$omp end parallel do
   end subroutine ml_horizontal_diffusion_tracer_unstr_impl

   subroutine ml_vertical_diffusion_tracer_unstr(state, hTr_layer, kappa_bg, dt)
      !! Ri-dependent vertical diffusion of a scalar tracer.
      !! Uses stronger Ri suppression for tracers:
      !!   kappa_t = nu_max / (1 + alpha*Ri)^3 + kappa_bg.
      !!
      !! Column-local — no neighbour access, safe for `do concurrent`.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout), contiguous :: hTr_layer(:, :)
      real(wp), intent(in) :: kappa_bg
      real(wp), intent(in) :: dt

      if (.not. state%use_multilayer) return
      if (state%multilayer%nz_ml < 2) return

      call ml_vertical_diffusion_tracer_unstr_impl( &
         state%multilayer%nz_ml, state%ncells_total, state%ncells_owned, &
         state%barotropic%height_cv, state%multilayer%h_layer, &
         state%multilayer%hu_layer, state%multilayer%hv_layer, &
         state%multilayer%rho_layer, state%multilayer%rho_0, &
         state%multilayer%vmix_ri_alpha, state%multilayer%vmix_nu_max, &
         hTr_layer, kappa_bg, dt)
   end subroutine ml_vertical_diffusion_tracer_unstr

   subroutine ml_vertical_diffusion_tracer_unstr_impl(nz, n, ncells_owned, height_cv, &
                                                           h_layer, hu_layer, hv_layer, &
                                                           rho_layer, rho_0, vmix_ri_alpha, &
                                                           vmix_nu_max, hTr_layer, kappa_bg, dt)
      !! Flat-impl of [[ml_vertical_diffusion_tracer_unstr]]: explicit-shape
      !! dummies so the device kernel addresses raw pointers, not derived-type
      !! descriptors.  Column-local; no neighbour access.
      integer, intent(in) :: nz, n, ncells_owned
      real(wp), intent(in) :: height_cv(n)
      real(wp), intent(in) :: h_layer(nz, n), hu_layer(nz, n), hv_layer(nz, n)
      real(wp), intent(in) :: rho_layer(nz, n)
      real(wp), intent(in) :: rho_0, vmix_ri_alpha, vmix_nu_max
      real(wp), intent(inout) :: hTr_layer(nz, n)
      real(wp), intent(in) :: kappa_bg, dt

      integer :: cell, k
      real(wp) :: h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1
      real(wp) :: du, dv, drho, Ri, denom, kappa_t
      real(wp) :: tr_k, tr_k1, dtr, flux_tr
      real(wp) :: flux_iface(BPG_NZ_MAX)

      real(wp), parameter :: EPS = 1.0e-20_wp

      !$omp parallel do private(k, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, du, dv, drho, Ri, denom, kappa_t, tr_k, tr_k1, dtr, flux_tr, flux_iface)
      do cell = 1, ncells_owned

         if (height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         ! Pass 1: compute interface fluxes
         do k = 1, nz - 1
            h_k = h_layer(k, cell)
            h_k1 = h_layer(k + 1, cell)
            h_avg = 0.5_wp*(h_k + h_k1)

            flux_iface(k) = 0.0_wp

            ! VCOORD_ZSTAR_FULL: no-flux through vanishing-layer interfaces.
            if (h_avg > VANISHING_LAYER_TOL_UNSTR .and. &
                h_k > VANISHING_LAYER_TOL_UNSTR .and. &
                h_k1 > VANISHING_LAYER_TOL_UNSTR) then
               if (h_k > DRY_TOLERANCE) then
                  u_k = hu_layer(k, cell)/h_k
                  v_k = hv_layer(k, cell)/h_k
                  tr_k = hTr_layer(k, cell)/h_k
               else
                  u_k = 0.0_wp
                  v_k = 0.0_wp
                  tr_k = 0.0_wp
               end if
               if (h_k1 > DRY_TOLERANCE) then
                  u_k1 = hu_layer(k + 1, cell)/h_k1
                  v_k1 = hv_layer(k + 1, cell)/h_k1
                  tr_k1 = hTr_layer(k + 1, cell)/h_k1
               else
                  u_k1 = 0.0_wp
                  v_k1 = 0.0_wp
                  tr_k1 = 0.0_wp
               end if

               du = u_k1 - u_k
               dv = v_k1 - v_k
               drho = rho_layer(k, cell) - rho_layer(k + 1, cell)

               Ri = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + EPS)
               if (Ri < 0.0_wp) Ri = 0.0_wp

               denom = 1.0_wp + vmix_ri_alpha*Ri
               kappa_t = vmix_nu_max/(denom*denom*denom) + kappa_bg

               ! Stability: kappa_t * dt / h_avg^2 <= 0.5
               if (kappa_t*dt > 0.5_wp*h_avg*h_avg) then
                  kappa_t = 0.5_wp*h_avg*h_avg/dt
               end if

               dtr = tr_k1 - tr_k
               flux_iface(k) = kappa_t*dtr/h_avg*dt
            end if
         end do

         ! Pass 2: apply fluxes.  flux_iface(k) is the diffusive transfer
         ! across interface k+1/2 (units [tr]*m); the same quantity adds to
         ! layer k and subtracts from layer k+1.  Multiplying each side by
         ! h_layer (previous form) was dimensionally wrong and non-conservative
         ! under non-uniform dz (e.g. VCOORD_ZSTAR_FULL).
         hTr_layer(1, cell) = hTr_layer(1, cell) + flux_iface(1)
         do k = 2, nz - 1
            hTr_layer(k, cell) = hTr_layer(k, cell) &
                                 + (flux_iface(k) - flux_iface(k - 1))
         end do
         hTr_layer(nz, cell) = hTr_layer(nz, cell) - flux_iface(nz - 1)
      end do
      !$omp end parallel do
   end subroutine ml_vertical_diffusion_tracer_unstr_impl

   subroutine ml_apply_surface_tracer_fluxes_unstr(state, hS_layer, hT_layer)
      !! Apply prescribed surface fluxes to the top layer (k=nz) of every
      !! cell.  Two evaporation modes, mirroring the structured
      !! `ml_apply_surface_tracer_fluxes`:
      !!
      !! **Default (evap_mass_sink = .false.) — virtual salt flux:**
      !!   hT(nz, c) += dt * heat_flux(c) / (rho_0 * cp)    ! thermal
      !!   hS(nz, c) += dt * (E - P)(c) * S_top             ! virtual salt
      !!
      !! **Mass-sink mode (evap_mass_sink = .true.) — real water removal:**
      !! net evaporation dh = E*dt is removed COLUMN-PROPORTIONALLY — every
      !! layer's h, hu, hv, hT is scaled by fac = 1 - dh/h_col while hS is
      !! held (salt concentrates, T/u/v invariant).  The removal is clamped
      !! to the water above ML_MIN_DEPTH_UNSTR; the unsatisfied fraction is
      !! dropped honestly.  Precipitation (E - P < 0) is NOT added here.
      !! The barotropic state is refreshed from the scaled layers by the
      !! `ml_compute_barotropic_momentum_unstr` call that follows every
      !! surface-flux call in the unstructured RK2 pipeline, so the depth
      !! drawdown propagates to `height_cv`/`stage_cv`.  This is a
      !! per-cell kernel writing only its OWN column (no neighbour scatter).
      !!
      !! **evap_latent_cooling** (needs evap_mass_sink): the evaporated mass
      !! additionally withdraws its latent heat of vaporization from the top
      !! layer.  hT(nz) is fac-scaled by the sensible sink FIRST; the latent
      !! term is the additive increment applied AFTER that scaling (NOT the
      !! full hT(nz) update): hT(nz) += -(rho_water*L_v/(rho_0*cp))*dh_eff —
      !! an ADDITIVE surface cooling, distinct from the sensible heat that
      !! leaves with the removed mass (no double-count).
      !!
      !! Layer-index convention: the unstructured solver indexes layers
      !! bottom-up (k=1 is the bottom sigma layer, k=nz is the surface
      !! layer just below the free surface) — see the header of
      !! `rki_ml_dynamics_unstr`.  Surface forcings therefore land at
      !! k=nz, mirroring the structured side after both adopted the
      !! same ROMS convention.
      !!
      !! Heat flux is in W/m^2, positive into the ocean.  E - P is in
      !! m/s, positive for net evaporation (concentrates salinity in the
      !! top layer).  Neither term changes h_layer -- this is pure
      !! tracer forcing in the ROMS / MOM convention.
      !!
      !! Mirrors `ml_apply_surface_tracer_fluxes` on the structured side
      !! (`rki_ml_tracers.F90`) modulo the axis-order swap (unstructured
      !! is `(nz, ncells_total)`).  Must run on every RK2 stage to
      !! preserve second-order accuracy in time.
      !!
      !! Reads the forcing fields directly off `state%forcing`; both
      !! arrays are zero-filled at init so the kernel is a bit-exact
      !! no-op until a setter pushes non-zero values.
      !!
      !! Only owned-cell tracer columns are written.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — ncells varies per mesh partition;
      ! caller host-dereferences S/T tracer slots before passing
      ! (see CLAUDE.md outer-shim + flat-impl pattern).
      real(wp), intent(inout) :: hS_layer(:, :)
      real(wp), intent(inout) :: hT_layer(:, :)  ! assumed-shape-ok: unstructured outer-shim; ncells varies per mesh partition

      integer  :: cell, k, nz, n
      real(wp) :: inv_rho0cp, dt, h_top, h_col, S_top, kin_q
      real(wp) :: dh_req, dh_eff, fac, latent_gain
      logical  :: evap_mass_sink

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total
      dt = state%dt
      inv_rho0cp = 1.0_wp/(state%multilayer%rho_0*CP_WATER)
      ! Hoist the scalar flags to host locals (avoid per-iteration struct
      ! deref on device; latent_gain == 0 when off -> bit-identical default).
      evap_mass_sink = state%multilayer%evap_mass_sink
      latent_gain = 0.0_wp
      if (state%multilayer%evap_latent_cooling) then
         latent_gain = RHO_WATER*LATENT_HEAT_VAPORIZATION*inv_rho0cp
      end if

      !$omp parallel do private(k, h_top, h_col, S_top, kin_q, dh_req, dh_eff, fac)
      do cell = 1, n
         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + state%multilayer%h_layer(k, cell)
         end do
         if (h_col > ML_MIN_DEPTH_UNSTR) then
            h_top = state%multilayer%h_layer(nz, cell)
            if (h_top > THIN_LAYER_THRESHOLD) then
               ! Thermal (both modes): dt * Q_net / (rho_0 * cp)  [m * K]
               kin_q = state%forcing%heat_flux(cell)*inv_rho0cp
               hT_layer(nz, cell) = hT_layer(nz, cell) + dt*kin_q

               if (evap_mass_sink) then
                  ! Column-proportional water-mass sink: hS held; h, hu, hv,
                  ! hT scale by fac so S concentrates while T, u, v stay
                  ! invariant.  Clamp to the water above the dry floor.
                  dh_req = dt*max(state%forcing%evap_minus_precip(cell), 0.0_wp)
                  dh_eff = min(dh_req, max(h_col - ML_MIN_DEPTH_UNSTR, 0.0_wp))
                  if (dh_eff > 0.0_wp) then
                     fac = 1.0_wp - dh_eff/h_col
                     do k = 1, nz
                        state%multilayer%h_layer(k, cell) = state%multilayer%h_layer(k, cell)*fac
                        state%multilayer%hu_layer(k, cell) = state%multilayer%hu_layer(k, cell)*fac
                        state%multilayer%hv_layer(k, cell) = state%multilayer%hv_layer(k, cell)*fac
                        hT_layer(k, cell) = hT_layer(k, cell)*fac
                        ! hS_layer held -> concentration rises naturally.
                     end do
                     ! Latent-heat cooling (additive, top layer; not scaled
                     ! by fac).  latent_gain == 0 when off -> bit-identical.
                     hT_layer(nz, cell) = hT_layer(nz, cell) - latent_gain*dh_eff
                  end if
               else
                  ! Legacy virtual salt flux (bit-identical default).
                  S_top = hS_layer(nz, cell)/h_top
                  hS_layer(nz, cell) = hS_layer(nz, cell) &
                                       + dt*state%forcing%evap_minus_precip(cell)*S_top
               end if
            end if
         end if
      end do
      !$omp end parallel do
   end subroutine ml_apply_surface_tracer_fluxes_unstr

end module rki_ml_tracers_unstr
