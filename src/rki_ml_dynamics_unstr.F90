!! Barotropic-baroclinic coupling and layer lifecycle for the unstructured
!! multilayer solver.
module rki_ml_dynamics_unstr
   !! Continuity-side machinery: layer initialisation from the barotropic
   !! state, distribute / distribute-with-remap (sigma + zsigma ALE),
   !! barotropic recompute from layer momenta, save / RK2-average for
   !! layer fields, diagnostic vertical velocity, vertical exchange of
   !! momentum.
   !!
   !! **Layer index convention**: `k = 1` is the BOTTOM sigma layer,
   !! `k = nz_ml` is the TOP / surface layer. `w_interface` is indexed
   !! `1..nz+1` with `w_interface(1) = 0` at the bed and `w_interface(nz+1)`
   !! at the free surface.
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            VCOORD_ZSTAR_FULL
#else
   use rki_constants, only: NZ_STACK_MAX, wp, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            VCOORD_ZSTAR_FULL
#endif
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_vcoord, only: vcoord_target_dz_column, &
                         vcoord_target_dz_column_zstar_full, &
                         zstar_full_build_column
   use rki_remap_column, only: remap_column
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR, BPG_NZ_MAX, &
                                  VANISHING_LAYER_TOL_UNSTR
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ml_init_layers_unstr
   public :: ml_save_barotropic_pre_unstr
   public :: ml_distribute_barotropic_unstr
   public :: ml_distribute_with_remap_unstr
   public :: ml_distribute_with_remap_full_unstr
   public :: ml_remap_tracer_unstr
   public :: ml_compute_barotropic_momentum_unstr
   public :: ml_compute_total_depth_unstr
   public :: ml_save_layer_state_unstr
   public :: ml_rk2_average_layers_unstr
   public :: ml_compute_vertical_velocity_unstr
   public :: ml_vertical_exchange_unstr
   public :: ml_build_vcoord_zstar_full_unstr

contains

   subroutine ml_init_layers_unstr(state)
      !! Fill the per-layer arrays from the barotropic state. Pure sigma:
      !!   h_k = dsig_k*H, hu_k = dsig_k*HU, hv_k = dsig_k*HV.
      !! For zsigma / hybrid coords `vcoord_target_dz_column` lays out the
      !! per-column thicknesses (matching distribute_with_remap), so the
      !! first step's tracer RK2 average stays geometrically coherent.
      !!
      !! Called once after bathymetry + ICs are set. Runs on the device
      !! (assumes the state went through `state_unstr_enter_data`).
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n
      integer :: coord_type
      real(wp) :: depth_transition, blend_width
      real(wp) :: H_col
      real(wp) :: dz_col(NZ_STACK_MAX)

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      coord_type = state%multilayer%vcoord%coord_type
      depth_transition = state%multilayer%vcoord%depth_transition
      blend_width = state%multilayer%vcoord%blend_width

      ! `height_cv` is populated lazily by the extrapolation kernel (zero
      ! at IC time); compute it from `stage_cv - bed_cv` for a consistent H.
      do concurrent(cell=1:n)
         state%barotropic%height_cv(cell) = max(state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell), 0.0_wp)
      end do

      ! Per-column thicknesses from vcoord, momentum split uniformly.
      ! `vcoord_target_dz_column` outputs dz_col in ROMS order
      ! (dz_col(1)=bottom, dz_col(nz)=surface), matching h_layer's
      ! bottom-up indexing — no reversal required.
      if (coord_type == VCOORD_ZSTAR_FULL) then
         block
            real(wp) :: h_min
            real(wp) :: z_ref_local(0:NZ_STACK_MAX)
            integer :: kk
            h_min = state%multilayer%vcoord%zstar_h_min
            !$acc parallel loop present(state) private(k, kk, H_col, dz_col, z_ref_local)
            do cell = 1, n
               H_col = state%barotropic%height_cv(cell)
               do kk = 0, nz
                  z_ref_local(kk) = state%multilayer%vcoord_z_ref_col(kk, cell)
               end do
               call vcoord_target_dz_column_zstar_full(nz, H_col, &
                                                       z_ref_local(0:nz), &
                                                       h_min, dz_col)
               do k = 1, nz
                  state%multilayer%h_layer(k, cell) = dz_col(k)
                  state%multilayer%hu_layer(k, cell) = state%multilayer%dsig_ml(k)*state%barotropic%xmom_cv(cell)
                  state%multilayer%hv_layer(k, cell) = state%multilayer%dsig_ml(k)*state%barotropic%ymom_cv(cell)
               end do
            end do
         end block
      else
         !$acc parallel loop present(state) private(k, H_col, dz_col)
         do cell = 1, n
            H_col = state%barotropic%height_cv(cell)
            call vcoord_target_dz_column(coord_type, nz, H_col, &
                                         state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref, &
                                         depth_transition, blend_width, dz_col)
            do k = 1, nz
               state%multilayer%h_layer(k, cell) = dz_col(k)
               state%multilayer%hu_layer(k, cell) = state%multilayer%dsig_ml(k)*state%barotropic%xmom_cv(cell)
               state%multilayer%hv_layer(k, cell) = state%multilayer%dsig_ml(k)*state%barotropic%ymom_cv(cell)
            end do
         end do
      end if
   end subroutine ml_init_layers_unstr

   subroutine ml_build_vcoord_zstar_full_unstr(state)
      !! Fill the per-cell reference z-level array from the bathymetry.
      !! Runs once at solver init (host-side loop); pushes
      !! `vcoord_z_ref_col` to the device when done. `bed_cv` is bed
      !! elevation (negative when submerged); its magnitude is `h_bed`.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, n, nz, n_surf, stretching
      real(wp) :: h_bed, h_surf_target

      if (.not. state%use_multilayer) return
      if (state%multilayer%vcoord%coord_type /= VCOORD_ZSTAR_FULL) return
      if (.not. allocated(state%multilayer%vcoord_z_ref_col)) return

      n = state%ncells_total
      nz = state%multilayer%nz_ml
      h_surf_target = state%multilayer%vcoord%zstar_h_surf_target
      n_surf = state%multilayer%vcoord%zstar_n_surf
      stretching = state%multilayer%vcoord%zstar_stretching

      !$acc update self(state%barotropic%bed_cv) if_present

      do cell = 1, n
         ! Clamp to non-negative bed depth: bed_cv > 0 (dry land) gives a
         ! zero-depth column rather than negative `z_ref` interfaces.
         h_bed = max(-state%barotropic%bed_cv(cell), 0.0_wp)
         call zstar_full_build_column(nz, h_bed, h_surf_target, &
                                      n_surf, stretching, &
                                      state%multilayer%vcoord_z_ref_col(:, cell))
      end do

      !$acc update device(state%multilayer%vcoord_z_ref_col) if_present
   end subroutine ml_build_vcoord_zstar_full_unstr

   pure subroutine ml_save_barotropic_pre_unstr(state)
      !! Snapshot the barotropic momentum before the HLL update so the
      !! distribute kernel can compute the acceleration increment
      !! `du = (hu_new - hu_pre) / H_new` at cell level.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, n

      if (.not. state%use_multilayer) return

      n = state%ncells_total
      do concurrent(cell=1:n)
         state%multilayer%xmom_pre(cell) = state%barotropic%xmom_cv(cell)
         state%multilayer%ymom_pre(cell) = state%barotropic%ymom_cv(cell)
      end do
   end subroutine ml_save_barotropic_pre_unstr

   pure subroutine ml_distribute_barotropic_unstr(state)
      !! Distribute the barotropic HLL acceleration across layers (sigma):
      !!   du = (hu_new - hu_pre)/H_new; u_k_old = hu_k_old/h_k_old;
      !!   h_k = dsig_k*H_new; hu_k = h_k*(u_k_old + du)
      !! so the acceleration lands uniformly while preserving shear.
      !! Requires `xmom_pre`/`ymom_pre` from `ml_save_barotropic_pre_unstr`.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n
      real(wp) :: h_val, du, dv, h_k_old, u_k, v_k, h_k_new

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      do concurrent(cell=1:n) local(k, h_val, du, dv, h_k_old, u_k, v_k, h_k_new)
         ! `height_cv` is stale right after the barotropic update;
         ! recompute H_new from `stage_cv - bed_cv`, clamped at zero so
         ! emptied cells get a clean dry layer profile.
         h_val = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
         if (h_val < 0.0_wp) h_val = 0.0_wp
         if (h_val > DRY_TOLERANCE) then
            du = (state%barotropic%xmom_cv(cell) - state%multilayer%xmom_pre(cell))/h_val
            dv = (state%barotropic%ymom_cv(cell) - state%multilayer%ymom_pre(cell))/h_val
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         do k = 1, nz
            h_k_old = state%multilayer%h_layer(k, cell)
            if (h_k_old > DRY_TOLERANCE) then
               u_k = state%multilayer%hu_layer(k, cell)/h_k_old
               v_k = state%multilayer%hv_layer(k, cell)/h_k_old
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_k_new = state%multilayer%dsig_ml(k)*h_val
            state%multilayer%h_layer(k, cell) = h_k_new

            if (h_k_new < THIN_LAYER_THRESHOLD) then
               state%multilayer%hu_layer(k, cell) = 0.0_wp
               state%multilayer%hv_layer(k, cell) = 0.0_wp
            else
               state%multilayer%hu_layer(k, cell) = h_k_new*(u_k + du)
               state%multilayer%hv_layer(k, cell) = h_k_new*(v_k + dv)
            end if
         end do
      end do
   end subroutine ml_distribute_barotropic_unstr

   subroutine ml_distribute_with_remap_unstr(state)
      !! Barotropic-to-layer distribution with vcoord reshaping for
      !! non-sigma coords (z-sigma hybrid, zstar-lite) — momentum +
      !! thickness only:
      !! 1. Save old thicknesses into `remap_dz_old` (read by
      !!    ml_remap_tracer_unstr for the per-tracer remap).
      !! 2. Compute new target thicknesses via vcoord, write h_layer.
      !! 3. Distribute the barotropic acceleration to the new layers.
      !! Per-tracer remap is done by ml_remap_tracer_unstr (looped over
      !! the registry by the solver); touches only top-level allocatables
      !! so no outer shim needed.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n
      real(wp) :: h_val, du, dv, h_k_old, u_k, v_k, h_k_new
      real(wp) :: dz_new(NZ_STACK_MAX)
      integer :: coord_type
      real(wp) :: depth_transition, blend_width

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total
      ! Extract vcoord params to local scalars so the GPU loop avoids
      ! nested derived-type access on device.
      coord_type = state%multilayer%vcoord%coord_type
      depth_transition = state%multilayer%vcoord%depth_transition
      blend_width = state%multilayer%vcoord%blend_width

      ! Step 1: save old layer thicknesses for the per-tracer remap.
      do concurrent(cell=1:n, k=1:nz)
         state%multilayer%remap_dz_old(k, cell) = state%multilayer%h_layer(k, cell)
      end do

      ! Steps 2-3: target thicknesses + momentum distribution.
      ! `!$acc parallel loop` (not do concurrent): the column body calls
      ! vcoord_target_dz_column cross-module, which NVHPC stdpar won't
      ! lower from do concurrent.
      !$acc parallel loop present(state) &
      !$acc&   private(k, h_val, du, dv, h_k_old, u_k, v_k, h_k_new, dz_new)
      do cell = 1, n

         ! `height_cv` is stale after the barotropic update; recompute
         ! H_new from `stage_cv - bed_cv`.
         h_val = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
         if (h_val < 0.0_wp) h_val = 0.0_wp

         if (h_val > DRY_TOLERANCE) then
            du = (state%barotropic%xmom_cv(cell) - state%multilayer%xmom_pre(cell))/h_val
            dv = (state%barotropic%ymom_cv(cell) - state%multilayer%ymom_pre(cell))/h_val
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         ! New target layer thicknesses.  Top-level dsig_ml / vcoord_z_ref
         ! (device-mapped) rather than the nested vcoord pointers.
         ! Output is ROMS order (dz_new(1)=bottom, dz_new(nz)=surface),
         ! matching h_layer's bottom-up indexing — no reversal required.
         call vcoord_target_dz_column(coord_type, nz, h_val, &
                                      state%multilayer%dsig_ml, state%multilayer%vcoord_z_ref, &
                                      depth_transition, blend_width, dz_new)

         do k = 1, nz
            h_k_old = state%multilayer%remap_dz_old(k, cell)
            if (h_k_old > DRY_TOLERANCE) then
               u_k = state%multilayer%hu_layer(k, cell)/h_k_old
               v_k = state%multilayer%hv_layer(k, cell)/h_k_old
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_k_new = dz_new(k)
            state%multilayer%h_layer(k, cell) = h_k_new

            if (h_k_new < THIN_LAYER_THRESHOLD) then
               state%multilayer%hu_layer(k, cell) = 0.0_wp
               state%multilayer%hv_layer(k, cell) = 0.0_wp
            else
               state%multilayer%hu_layer(k, cell) = h_k_new*(u_k + du)
               state%multilayer%hv_layer(k, cell) = h_k_new*(v_k + dv)
            end if
         end do
      end do
   end subroutine ml_distribute_with_remap_unstr

   subroutine ml_remap_tracer_unstr(state, hTr_layer)
      !! Conservatively remap ONE tracer from the saved pre-update layer
      !! thicknesses (state%multilayer%remap_dz_old, written by the
      !! distribute-with-remap routines) onto the current h_layer (the new
      !! vcoord target).  Vcoord-agnostic: shared by the sigma / zsigma /
      !! zstar-lite and the ZSTAR_FULL paths, since both write the new
      !! target into h_layer before this runs.  Split out of the fused
      !! remap so the solver can loop it over the whole tracer registry;
      !! the per-tracer arithmetic is exactly what the fused impls used
      !! for S/T (scale dz_old to sum to H_new, depth-weighted column-mean
      !! fallback at vanishing layers, conservative remap_column), so
      !! salinity / temperature results are unchanged.
      !!
      !! Outer shim: the solver passes state%multilayer%tracers(it)%hTr as
      !! the flat 2D actual argument; the array-of-derived-types
      !! indirection is resolved host-side.  Must run AFTER the distribute
      !! (which saves remap_dz_old and writes the new h_layer).
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(inout) :: hTr_layer(state%multilayer%nz_ml, state%ncells_total)

      integer :: cell, k, nz, n, remap_method
      real(wp) :: h_val, h_old_total, scale_old, tr_col_sum, tr_col_mean
      real(wp) :: dz_old_col(NZ_STACK_MAX), dz_new_col(NZ_STACK_MAX)
      real(wp) :: tr_old(NZ_STACK_MAX), tr_new(NZ_STACK_MAX)

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total
      remap_method = state%multilayer%vcoord%remap_method

      !$acc parallel loop present(state, hTr_layer) &
      !$acc&   private(k, h_val, h_old_total, scale_old, tr_col_sum, &
      !$acc&           tr_col_mean, dz_old_col, dz_new_col, tr_old, tr_new)
      do cell = 1, n

         ! Recompute H_new from stage_cv - bed_cv (height_cv is stale right
         ! after the barotropic update), matching the distribute routines.
         h_val = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
         if (h_val < 0.0_wp) h_val = 0.0_wp

         ! Advection added flux mass using H_new while remap_dz_old still
         ! sums to H_old; scale dz_old uniformly to H_new so remap_column's
         ! `sum dz_old == sum dz_new` conservation guarantee holds.
         h_old_total = 0.0_wp
         do k = 1, nz
            h_old_total = h_old_total + state%multilayer%remap_dz_old(k, cell)
         end do
         if (h_old_total > DRY_TOLERANCE) then
            scale_old = h_val/h_old_total
         else
            scale_old = 1.0_wp
         end if

         ! Depth-weighted column-mean concentration as the vanishing-layer
         ! fallback (falling back to 0 leaks tracer through PPM
         ! redistribution over pinch-recreate cycles; column-mean keeps the
         ! remap conservative).
         tr_col_sum = 0.0_wp
         do k = 1, nz
            tr_col_sum = tr_col_sum + hTr_layer(k, cell)
         end do
         if (h_val > DRY_TOLERANCE) then
            tr_col_mean = tr_col_sum/h_val
         else
            tr_col_mean = 0.0_wp
         end if

         ! dz_new_col is the new vcoord target the distribute already wrote
         ! into h_layer, so this routine is vcoord-agnostic.
         do k = 1, nz
            dz_old_col(k) = state%multilayer%remap_dz_old(k, cell)*scale_old
            dz_new_col(k) = state%multilayer%h_layer(k, cell)
            if (dz_old_col(k) > DRY_TOLERANCE) then
               tr_old(k) = hTr_layer(k, cell)/dz_old_col(k)
            else
               tr_old(k) = tr_col_mean
            end if
         end do

         call remap_column(remap_method, nz, dz_old_col, dz_new_col, tr_old, tr_new)

         do k = 1, nz
            hTr_layer(k, cell) = dz_new_col(k)*tr_new(k)
         end do
      end do
   end subroutine ml_remap_tracer_unstr

   subroutine ml_distribute_with_remap_full_unstr(state)
      !! VCOORD_ZSTAR_FULL variant of `ml_distribute_with_remap_unstr`
      !! (momentum + thickness only).  Builds the per-cell target dz_new
      !! via `vcoord_target_dz_column_zstar_full` from the cell's local
      !! `vcoord_z_ref_col(:, cell)` slice.  Per-tracer remap is done by
      !! ml_remap_tracer_unstr (looped over the registry by the solver),
      !! which is vcoord-agnostic — it reads the new h_layer this routine
      !! writes.  Shares the remap_dz_old workspace with the sigma path.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n, kk
      real(wp) :: h_val, du, dv, h_k_old, u_k, v_k, h_k_new, h_min
      real(wp) :: dz_new(NZ_STACK_MAX)
      real(wp) :: z_ref_local(0:NZ_STACK_MAX)

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total
      h_min = state%multilayer%vcoord%zstar_h_min

      do concurrent(cell=1:n, k=1:nz)
         state%multilayer%remap_dz_old(k, cell) = state%multilayer%h_layer(k, cell)
      end do

      !$acc parallel loop present(state) &
      !$acc&   private(k, kk, h_val, du, dv, h_k_old, u_k, v_k, h_k_new, &
      !$acc&           dz_new, z_ref_local)
      do cell = 1, n

         h_val = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
         if (h_val < 0.0_wp) h_val = 0.0_wp

         if (h_val > DRY_TOLERANCE) then
            du = (state%barotropic%xmom_cv(cell) - state%multilayer%xmom_pre(cell))/h_val
            dv = (state%barotropic%ymom_cv(cell) - state%multilayer%ymom_pre(cell))/h_val
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         ! Per-cell z_ref slice into a private stack array (explicit-shape
         ! slice for the pure callee).
         do kk = 0, nz
            z_ref_local(kk) = state%multilayer%vcoord_z_ref_col(kk, cell)
         end do

         call vcoord_target_dz_column_zstar_full(nz, h_val, &
                                                 z_ref_local(0:nz), &
                                                 h_min, dz_new)

         do k = 1, nz
            h_k_old = state%multilayer%remap_dz_old(k, cell)
            if (h_k_old > DRY_TOLERANCE) then
               u_k = state%multilayer%hu_layer(k, cell)/h_k_old
               v_k = state%multilayer%hv_layer(k, cell)/h_k_old
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_k_new = dz_new(k)
            state%multilayer%h_layer(k, cell) = h_k_new

            if (h_k_new < THIN_LAYER_THRESHOLD) then
               state%multilayer%hu_layer(k, cell) = 0.0_wp
               state%multilayer%hv_layer(k, cell) = 0.0_wp
            else
               state%multilayer%hu_layer(k, cell) = h_k_new*(u_k + du)
               state%multilayer%hv_layer(k, cell) = h_k_new*(v_k + dv)
            end if
         end do
      end do
   end subroutine ml_distribute_with_remap_full_unstr

   pure subroutine ml_compute_barotropic_momentum_unstr(state)
      !! Depth-integrate layer momentum back onto the barotropic fields:
      !!   HU = sum_k hu_k,  HV = sum_k hv_k.
      !! Also refreshes `height_cv = sum_k h_k` and `stage_cv` from the
      !! updated height so the 2D barotropic state stays in sync with
      !! the layers. Called between RK stages and at the end of the
      !! step after the layer RK2 average.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n
      real(wp) :: su, sv, sh

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      do concurrent(cell=1:n) local(k, su, sv, sh)
         su = 0.0_wp
         sv = 0.0_wp
         sh = 0.0_wp
         do k = 1, nz
            su = su + state%multilayer%hu_layer(k, cell)
            sv = sv + state%multilayer%hv_layer(k, cell)
            sh = sh + state%multilayer%h_layer(k, cell)
         end do
         state%barotropic%xmom_cv(cell) = su
         state%barotropic%ymom_cv(cell) = sv
         state%barotropic%height_cv(cell) = sh
         state%barotropic%stage_cv(cell) = sh + state%barotropic%bed_cv(cell)
      end do
   end subroutine ml_compute_barotropic_momentum_unstr

   pure subroutine ml_compute_total_depth_unstr(state)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Refresh `height_cv = sum_k h_k` without touching momentum.
      !! Useful when only the layer thicknesses have changed.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n
      real(wp) :: sh

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      do concurrent(cell=1:n) local(k, sh)
         sh = 0.0_wp
         do k = 1, nz
            sh = sh + state%multilayer%h_layer(k, cell)
         end do
         state%barotropic%height_cv(cell) = sh
         state%barotropic%stage_cv(cell) = sh + state%barotropic%bed_cv(cell)
      end do
   end subroutine ml_compute_total_depth_unstr

   pure subroutine ml_save_layer_state_unstr(state)
      !! Snapshot the layer state at the start of the step for the
      !! SSP-RK2 average.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      do concurrent(cell=1:n, k=1:nz)
         state%multilayer%h_layer0(k, cell) = state%multilayer%h_layer(k, cell)
         state%multilayer%hu_layer0(k, cell) = state%multilayer%hu_layer(k, cell)
         state%multilayer%hv_layer0(k, cell) = state%multilayer%hv_layer(k, cell)
      end do
   end subroutine ml_save_layer_state_unstr

   pure subroutine ml_rk2_average_layers_unstr(state)
      !! Public only for the unit-test suite (no production module imports it);
      !! ignore when developing production code in other modules.
      !! Final RK2 blend:  Q := 0.5 * Q + 0.5 * Q^n
      !! The barotropic blend is done in `fused_update_stage2_unstr`;
      !! this routine handles the layer fields.
      type(unstr_state_t), intent(inout) :: state

      integer :: cell, k, nz, n

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_total

      do concurrent(cell=1:n, k=1:nz)
         state%multilayer%h_layer(k, cell) = 0.5_wp*(state%multilayer%h_layer(k, cell) &
                                                     + state%multilayer%h_layer0(k, cell))
         state%multilayer%hu_layer(k, cell) = 0.5_wp*(state%multilayer%hu_layer(k, cell) &
                                                      + state%multilayer%hu_layer0(k, cell))
         state%multilayer%hv_layer(k, cell) = 0.5_wp*(state%multilayer%hv_layer(k, cell) &
                                                      + state%multilayer%hv_layer0(k, cell))
      end do
   end subroutine ml_rk2_average_layers_unstr

   ! ======================================================================
   ! Vertical velocity and momentum exchange
   !
   ! The diagnostic vertical velocity w is computed from the horizontal
   ! divergence of layer momentum using Green-Gauss on triangles:
   !   div(hu_k, hv_k) = (1/A) * sum_e (hu_edge * n1 + hv_edge * n2) * L
   ! with edge values taken as the average of the cell and its neighbour.
   ! Then w is integrated bottom-up:
   !   w(k+1) = w(k) - div(hu_k, hv_k)
   !
   ! Vertical momentum exchange uses upwind advection at layer interfaces
   ! based on the sign of w. A 50% limiter prevents draining more than
   ! half the donor layer in a single step.
   ! ======================================================================

   pure subroutine ml_compute_vertical_velocity_unstr(state, mesh)
      !! Compute diagnostic vertical velocity from horizontal divergence.
      !!
      !! Uses the barotropic mass flux stashed by the flux kernel to ensure
      !! consistency with the barotropic continuity equation:
      !!   div(HU, HV) = -(1/A) * sum_e mass_flux_e(e, cell)
      !!
      !! For sigma coordinates with fixed fractions:
      !!   div(hu_k, hv_k) = dsig_k * div(HU, HV)
      !!
      !! Then integrates bottom-up:
      !!   w(1)   = 0     (no penetration at bottom)
      !!   w(k+1) = w(k) - div(hu_k, hv_k)
      !!
      !! This ensures the computed w is exactly consistent with the
      !! barotropic mass fluxes, so uniform ICs stay uniform.
      !!
      !! Writes to state%multilayer%w_interface(nz+1, ncells_total). Run after the
      !! barotropic distribute step has updated layer velocities.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: cell, k, e, nz, n
      real(wp) :: inv_area, div_barotropic, div_layer

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = mesh%ncells_owned

      do concurrent(cell=1:n) local(k, e, inv_area, div_barotropic, div_layer)
         inv_area = 1.0_wp/mesh%area(cell)

         ! No penetration at bottom
         state%multilayer%w_interface(1, cell) = 0.0_wp

         ! Skip thin columns — vertical exchange is meaningless
         if (state%barotropic%height_cv(cell) < ML_MIN_DEPTH_UNSTR) then
            do k = 2, nz + 1
               state%multilayer%w_interface(k, cell) = 0.0_wp
            end do
         else
            ! Compute barotropic divergence from the mass flux:
            ! mass_flux_e > 0 means inflow (negative divergence)
            div_barotropic = 0.0_wp
            do e = 1, 3
               div_barotropic = div_barotropic - state%barotropic%mass_flux_e(e, cell)
            end do
            div_barotropic = div_barotropic*inv_area

            ! Integrate bottom-up: w(k+1) = w(k) - div(hu_k, hv_k)
            ! For sigma: div(hu_k) = dsig_k * div(HU)
            do k = 1, nz
               div_layer = state%multilayer%dsig_ml(k)*div_barotropic
               state%multilayer%w_interface(k + 1, cell) = state%multilayer%w_interface(k, cell) - div_layer
            end do
         end if
      end do
   end subroutine ml_compute_vertical_velocity_unstr

   pure subroutine ml_vertical_exchange_unstr(state, dt)
      !! Vertical momentum exchange using upwind advection at interfaces.
      !!
      !! For each interior interface k+1/2 (between layers k and k+1):
      !!   - If w > 0 (upward): flux = w * u_k (donor is layer k)
      !!   - If w < 0 (downward): flux = w * u_{k+1} (donor is layer k+1)
      !!
      !! A 50% limiter prevents draining more than half the donor layer
      !! momentum in a single step.
      !!
      !! Uses two-pass approach: first compute all interface fluxes from
      !! the current state, then apply them. This avoids the race condition
      !! where modifying layer k affects the flux computation at layer k+1.
      !!
      !! Requires w_interface to be populated by ml_compute_vertical_velocity_unstr.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt

      integer :: cell, k, nz, n
      real(wp) :: w_face, u_up, v_up, flux_u, flux_v, h_val
      real(wp) :: flux_u_iface(BPG_NZ_MAX), flux_v_iface(BPG_NZ_MAX)

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_owned

      ! Single-layer case: no vertical exchange
      if (nz < 2) return

      do concurrent(cell=1:n) local(k, w_face, u_up, v_up, flux_u, flux_v, h_val, &
                                    flux_u_iface, flux_v_iface)
         ! Skip thin columns
         if (state%barotropic%height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         ! Pass 1: compute all interface fluxes from current state
         do k = 1, nz - 1
            w_face = state%multilayer%w_interface(k + 1, cell)

            ! Upwind velocity selection
            if (w_face > 0.0_wp) then
               ! Upward flow: donor is layer k (below)
               h_val = state%multilayer%h_layer(k, cell)
               if (h_val > DRY_TOLERANCE) then
                  u_up = state%multilayer%hu_layer(k, cell)/h_val
                  v_up = state%multilayer%hv_layer(k, cell)/h_val
               else
                  u_up = 0.0_wp
                  v_up = 0.0_wp
               end if
            else
               ! Downward flow: donor is layer k+1 (above)
               h_val = state%multilayer%h_layer(k + 1, cell)
               if (h_val > DRY_TOLERANCE) then
                  u_up = state%multilayer%hu_layer(k + 1, cell)/h_val
                  v_up = state%multilayer%hv_layer(k + 1, cell)/h_val
               else
                  u_up = 0.0_wp
                  v_up = 0.0_wp
               end if
            end if

            ! Flux = w * velocity * dt
            flux_u = w_face*u_up*dt
            flux_v = w_face*v_up*dt

            ! 50% limiter: don't extract more than half donor layer momentum
            if (w_face > 0.0_wp) then
               ! Donor is layer k
               if (abs(flux_u) > 0.5_wp*abs(state%multilayer%hu_layer(k, cell))) then
                  flux_u = sign(0.5_wp*abs(state%multilayer%hu_layer(k, cell)), flux_u)
               end if
               if (abs(flux_v) > 0.5_wp*abs(state%multilayer%hv_layer(k, cell))) then
                  flux_v = sign(0.5_wp*abs(state%multilayer%hv_layer(k, cell)), flux_v)
               end if
            else
               ! Donor is layer k+1
               if (abs(flux_u) > 0.5_wp*abs(state%multilayer%hu_layer(k + 1, cell))) then
                  flux_u = sign(0.5_wp*abs(state%multilayer%hu_layer(k + 1, cell)), flux_u)
               end if
               if (abs(flux_v) > 0.5_wp*abs(state%multilayer%hv_layer(k + 1, cell))) then
                  flux_v = sign(0.5_wp*abs(state%multilayer%hv_layer(k + 1, cell)), flux_v)
               end if
            end if

            flux_u_iface(k) = flux_u
            flux_v_iface(k) = flux_v
         end do

         ! Pass 2: apply fluxes to layers
         ! Bottom layer (k=1): only has the interface ABOVE it
         state%multilayer%hu_layer(1, cell) = state%multilayer%hu_layer(1, cell) - flux_u_iface(1)
         state%multilayer%hv_layer(1, cell) = state%multilayer%hv_layer(1, cell) - flux_v_iface(1)
         ! Interior layers: interface below and above
         do k = 2, nz - 1
            state%multilayer%hu_layer(k, cell) = state%multilayer%hu_layer(k, cell) + flux_u_iface(k - 1) - flux_u_iface(k)
            state%multilayer%hv_layer(k, cell) = state%multilayer%hv_layer(k, cell) + flux_v_iface(k - 1) - flux_v_iface(k)
         end do
         ! Top / surface layer (k=nz): only has the interface BELOW it
         state%multilayer%hu_layer(nz, cell) = state%multilayer%hu_layer(nz, cell) + flux_u_iface(nz - 1)
         state%multilayer%hv_layer(nz, cell) = state%multilayer%hv_layer(nz, cell) + flux_v_iface(nz - 1)
      end do
   end subroutine ml_vertical_exchange_unstr

end module rki_ml_dynamics_unstr
