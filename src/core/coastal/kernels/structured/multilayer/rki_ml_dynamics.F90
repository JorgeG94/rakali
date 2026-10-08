!! Layer-integrated dynamics for the structured multilayer solver.
module rki_ml_dynamics
   !! Barotropic-baroclinic coupling and vertical structure: distribute
   !! barotropic acceleration to layers (preserving shear), maintain layer
   !! thickness (sigma constraint or vcoord remap), diagnostic w + interface
   !! momentum exchange, layer-state RK2 plumbing and BC distribution.
   !! Owns the persistent vcoord-remap workspace; teardown via
   !! `ml_dynamics_workspaces_cleanup`.
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            RHO_WATER, VCOORD_ZSTAR_FULL
#else
   use rki_constants, only: NZ_STACK_MAX, wp, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            RHO_WATER, VCOORD_ZSTAR_FULL
#endif
   use rki_grid, only: hgrid_t
   use rki_boundary_types, only: BC_WALL
   use rki_vcoord, only: vcoord_target_dz_column, &
                         vcoord_target_dz_column_zstar_full, &
                         zstar_full_build_column
   use rki_remap_column, only: remap_column
   use rki_ml_common, only: ML_MIN_DEPTH
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ml_distribute_barotropic
   public :: ml_distribute_with_remap
   public :: ml_remap_tracer
   public :: ml_compute_total_depth
   public :: ml_compute_barotropic_momentum
   public :: ml_compute_vertical_velocity
   public :: ml_compute_vertical_exchange
   public :: ml_rk2_average_layers
   public :: ml_init_layers_from_barotropic
   public :: ml_save_layer_state
   public :: ml_apply_layer_bcs
   public :: ml_apply_layer_bcs_full
   public :: ml_build_vcoord_zstar_full
   public :: ml_distribute_with_remap_full
   public :: ml_apply_wind_to_surface_layer

   public :: ml_dynamics_workspaces_cleanup

   ! ---- Persistent per-kernel workspaces ----

   ! ml_distribute_with_remap
   real(wp), allocatable :: remap_dz_old(:, :, :)

contains

   pure subroutine ml_distribute_barotropic(h, hu, hv, hu_pre, hv_pre, &
                                            dsig, h_layer, hu_layer, hv_layer, &
                                            grid, nz)
      !! Distribute barotropic HLL acceleration to layers, preserving
      !! per-layer velocity anomalies (shear). h_k = dsig_k*H (sigma
      !! constraint); hu_k = h_k*(u_k_old + du), du = (HU_new-HU_pre)/H.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
         !! Total depth AFTER barotropic update
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
         !! Barotropic x-momentum AFTER update
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
         !! Barotropic y-momentum AFTER update
      real(wp), intent(in) :: hu_pre(grid%nx_total, grid%ny_total)
         !! Barotropic x-momentum BEFORE update
      real(wp), intent(in) :: hv_pre(grid%nx_total, grid%ny_total)
         !! Barotropic y-momentum BEFORE update
      real(wp), intent(in) :: dsig(nz)
         !! Sigma layer fractions
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
         !! Layer thickness (old on entry, dsig_k*H on exit)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
         !! Layer x-momentum (updated)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
         !! Layer y-momentum (updated)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_val, du, dv, h_k_old, u_k, v_k

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, h_val, du, dv, h_k_old, u_k, v_k)
         h_val = h(i, j)
         if (h_val > DRY_TOLERANCE) then
            du = (hu(i, j) - hu_pre(i, j))/h_val
            dv = (hv(i, j) - hv_pre(i, j))/h_val
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         do k = 1, nz
            h_k_old = h_layer(i, j, k)
            if (h_k_old > DRY_TOLERANCE) then
               u_k = hu_layer(i, j, k)/h_k_old
               v_k = hv_layer(i, j, k)/h_k_old
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_layer(i, j, k) = dsig(k)*h_val

            ! Thin-layer zeroing consistent with the barotropic update
            if (h_layer(i, j, k) < THIN_LAYER_THRESHOLD) then
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            else
               hu_layer(i, j, k) = h_layer(i, j, k)*(u_k + du)
               hv_layer(i, j, k) = h_layer(i, j, k)*(v_k + dv)
            end if
         end do
      end do

   end subroutine ml_distribute_barotropic

   subroutine ml_distribute_with_remap(h, hu, hv, hu_pre, hv_pre, &
                                       coord_type, &
                                       dsig, z_ref, &
                                       depth_transition, blend_width, &
                                       h_layer, hu_layer, hv_layer, &
                                       grid, nz)
      !! Barotropic-to-layer distribution with vertical-coordinate
      !! reshaping — momentum + thickness only. Saves old thicknesses into
      !! `remap_dz_old`, sets new target `h_layer` via vcoord, distributes
      !! the barotropic acceleration. Per-tracer remap is done separately by
      !! `ml_remap_tracer` (looped over the registry by the solver).
      !! All array args must be device-resident under GPU offload; vcoord_t
      !! is NOT passed (raw arrays/scalars instead) to avoid
      !! allocatable-component resolution issues on the device.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu_pre(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv_pre(grid%nx_total, grid%ny_total)
      integer, intent(in) :: coord_type
         !! VCOORD_SIGMA, VCOORD_ZSIGMA, etc.
      real(wp), intent(in) :: dsig(nz)
         !! Target sigma fractions (device-resident)
      real(wp), intent(in) :: z_ref(0:nz)
         !! Reference z-level interfaces (device-resident)
      real(wp), intent(in) :: depth_transition
      real(wp), intent(in) :: blend_width
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_val, du, dv, h_k_old, u_k, v_k
      real(wp) :: dz_new(NZ_STACK_MAX)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      call remap_workspace_ensure(nx, ny, nz)

      ! Save old layer thicknesses for the per-tracer remap.
      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         remap_dz_old(i, j, k) = h_layer(i, j, k)
      end do

      ! !$acc parallel loop (not do concurrent): the body calls
      ! vcoord_target_dz_column from another module (rki_vcoord), which
      ! NVHPC stdpar GPU offload does not lower from inside do concurrent.
      !$acc parallel loop collapse(2) &
      !$acc&   present(h, hu, hv, hu_pre, hv_pre, dsig, z_ref, &
      !$acc&           h_layer, hu_layer, hv_layer, remap_dz_old) &
      !$acc&   private(k, h_val, du, dv, h_k_old, u_k, v_k, dz_new)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost

         h_val = h(i, j)

         if (h_val > DRY_TOLERANCE) then
            du = (hu(i, j) - hu_pre(i, j))/h_val
            dv = (hv(i, j) - hv_pre(i, j))/h_val
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         call vcoord_target_dz_column(coord_type, nz, h_val, &
                                      dsig, z_ref, &
                                      depth_transition, blend_width, &
                                      dz_new)

         do k = 1, nz
            h_k_old = remap_dz_old(i, j, k)
            if (h_k_old > DRY_TOLERANCE) then
               u_k = hu_layer(i, j, k)/h_k_old
               v_k = hv_layer(i, j, k)/h_k_old
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_layer(i, j, k) = dz_new(k)

            if (h_layer(i, j, k) < THIN_LAYER_THRESHOLD) then
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            else
               hu_layer(i, j, k) = h_layer(i, j, k)*(u_k + du)
               hv_layer(i, j, k) = h_layer(i, j, k)*(v_k + dv)
            end if
         end do
      end do
      end do

   end subroutine ml_distribute_with_remap

   pure subroutine ml_apply_wind_to_surface_layer(hu_layer, hv_layer, h_layer, &
                                                  wind_x, wind_y, grid, nz, dt)
      !! Apply the surface wind stress to the TOP layer (k=nz) momentum only,
      !! so the wind is surface-concentrated and generates vertical shear that
      !! the vertical mixing then spreads downward — instead of the
      !! depth-uniform barotropic forcing (`apply_friction`).  tau =
      !! (wind_x, wind_y) in Pa; adds dt*tau/RHO_WATER to hu/hv_layer(:,:,nz).
      !! Gated by `wind_surface_layer`; the caller suppresses the barotropic
      !! wind so the stress is not applied twice.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt

      integer :: i, j, nx, ny, nghost

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost)
         if (h_layer(i, j, nz) > THIN_LAYER_THRESHOLD) then
            hu_layer(i, j, nz) = hu_layer(i, j, nz) + dt*wind_x(i, j)/RHO_WATER
            hv_layer(i, j, nz) = hv_layer(i, j, nz) + dt*wind_y(i, j)/RHO_WATER
         end if
      end do

   end subroutine ml_apply_wind_to_surface_layer

   subroutine ml_remap_tracer(hTr, h, remap_method, h_layer, grid, nz)
      !! Conservatively remap ONE tracer from the saved pre-update layer
      !! thicknesses (`remap_dz_old`) onto the current `h_layer` (new vcoord
      !! target): scale dz_old to sum to H_new, column-mean fallback at
      !! vanishing layers, conservative `remap_column`. Looped over the
      !! tracer registry by the solver. Must run AFTER the dynamics remap
      !! (which writes `remap_dz_old` and the new `h_layer`).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hTr(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      integer, intent(in) :: remap_method
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_val, h_old_total, scale_old, tr_col_sum, tr_col_mean
      real(wp) :: dz_old_col(NZ_STACK_MAX), dz_new_col(NZ_STACK_MAX)
      real(wp) :: tr_old(NZ_STACK_MAX), tr_new(NZ_STACK_MAX)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$acc parallel loop collapse(2) &
      !$acc&   present(hTr, h, h_layer, remap_dz_old) &
      !$acc&   private(k, h_val, h_old_total, scale_old, tr_col_sum, &
      !$acc&           tr_col_mean, dz_old_col, dz_new_col, tr_old, tr_new)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost

         h_val = h(i, j)

         ! Scale dz_old uniformly to H_new so remap_column's
         ! `sum dz_old == sum dz_new` conservation guarantee holds.
         h_old_total = 0.0_wp
         do k = 1, nz
            h_old_total = h_old_total + remap_dz_old(i, j, k)
         end do
         if (h_old_total > DRY_TOLERANCE) then
            scale_old = h_val/h_old_total
         else
            scale_old = 1.0_wp
         end if

         ! Column-mean concentration as the vanishing-layer fallback
         ! (falling back to 0 leaks tracer through PPM over pinch cycles).
         tr_col_sum = 0.0_wp
         do k = 1, nz
            tr_col_sum = tr_col_sum + hTr(i, j, k)
         end do
         if (h_val > DRY_TOLERANCE) then
            tr_col_mean = tr_col_sum/h_val
         else
            tr_col_mean = 0.0_wp
         end if

         do k = 1, nz
            dz_old_col(k) = remap_dz_old(i, j, k)*scale_old
            dz_new_col(k) = h_layer(i, j, k)
            if (dz_old_col(k) > DRY_TOLERANCE) then
               tr_old(k) = hTr(i, j, k)/dz_old_col(k)
            else
               tr_old(k) = tr_col_mean
            end if
         end do

         call remap_column(remap_method, nz, dz_old_col, dz_new_col, tr_old, tr_new)

         do k = 1, nz
            hTr(i, j, k) = dz_new_col(k)*tr_new(k)
         end do
      end do
      end do

   end subroutine ml_remap_tracer

   subroutine ml_distribute_with_remap_full(h, hu, hv, hu_pre, hv_pre, &
                                            z_ref_col, h_min, &
                                            h_layer, hu_layer, hv_layer, &
                                            grid, nz)
      !! VCOORD_ZSTAR_FULL variant of `ml_distribute_with_remap` (momentum +
      !! thickness only): builds per-cell dz_new via
      !! `vcoord_target_dz_column_zstar_full` from the per-column reference
      !! profile `z_ref_col(0:nz,i,j)`. Per-tracer remap by `ml_remap_tracer`.
      !! Shares the `remap_dz_old` persistent workspace.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu_pre(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv_pre(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: z_ref_col(0:nz, grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: h_min
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_val, du, dv, h_k_old, u_k, v_k
      real(wp) :: dz_new(NZ_STACK_MAX)
      real(wp) :: z_ref_local(0:NZ_STACK_MAX)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      call remap_workspace_ensure(nx, ny, nz)

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         remap_dz_old(i, j, k) = h_layer(i, j, k)
      end do

      !$acc parallel loop collapse(2) &
      !$acc&   present(h, hu, hv, hu_pre, hv_pre, z_ref_col, &
      !$acc&           h_layer, hu_layer, hv_layer, remap_dz_old) &
      !$acc&   private(k, h_val, du, dv, h_k_old, u_k, v_k, &
      !$acc&           dz_new, z_ref_local)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost

         h_val = h(i, j)

         if (h_val > DRY_TOLERANCE) then
            du = (hu(i, j) - hu_pre(i, j))/h_val
            dv = (hv(i, j) - hv_pre(i, j))/h_val
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         ! Copy this column's reference profile into a private stack
         ! slice so the pure callee sees an explicit-shape array slice.
         do k = 0, nz
            z_ref_local(k) = z_ref_col(k, i, j)
         end do

         call vcoord_target_dz_column_zstar_full(nz, h_val, &
                                                 z_ref_local(0:nz), &
                                                 h_min, dz_new)

         do k = 1, nz
            h_k_old = remap_dz_old(i, j, k)
            if (h_k_old > DRY_TOLERANCE) then
               u_k = hu_layer(i, j, k)/h_k_old
               v_k = hv_layer(i, j, k)/h_k_old
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_layer(i, j, k) = dz_new(k)

            if (h_layer(i, j, k) < THIN_LAYER_THRESHOLD) then
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            else
               hu_layer(i, j, k) = h_layer(i, j, k)*(u_k + du)
               hv_layer(i, j, k) = h_layer(i, j, k)*(v_k + dv)
            end if
         end do
      end do
      end do

   end subroutine ml_distribute_with_remap_full

   pure subroutine ml_compute_total_depth(h_layer, h, grid, nz)
      !! Compute total water depth: H = sum_k h_k
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: h(grid%nx_total, grid%ny_total)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: total

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, total)
         total = 0.0_wp
         do k = 1, nz
            total = total + h_layer(i, j, k)
         end do
         h(i, j) = total
      end do

   end subroutine ml_compute_total_depth

   pure subroutine ml_compute_barotropic_momentum(h_layer, hu_layer, hv_layer, &
                                                  hu, hv, grid, nz)
      !! Compute depth-integrated momentum: hu = sum_k hu_k, hv = sum_k hv_k
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hv(grid%nx_total, grid%ny_total)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: su, sv

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, su, sv)
         su = 0.0_wp
         sv = 0.0_wp
         do k = 1, nz
            su = su + hu_layer(i, j, k)
            sv = sv + hv_layer(i, j, k)
         end do
         hu(i, j) = su
         hv(i, j) = sv
      end do

   end subroutine ml_compute_barotropic_momentum

   pure subroutine ml_compute_vertical_velocity(h_layer, hu_layer, hv_layer, &
                                                dsig, w, grid, nz)
      !! DIASURFACE (sigma-relative) vertical mass-transport velocity from
      !! per-layer continuity, integrated bottom-up.
      !!
      !! `w_ml` is the volume flux PER UNIT AREA across each moving sigma
      !! interface — i.e. the sigma-relative (diasurface) velocity that the
      !! vertical exchange kernels advect momentum/tracers with, NOT the
      !! physical Cartesian vertical velocity `w_phys`.
      !!
      !! In sigma coordinates the interfaces move with the free surface,
      !! `∂h_k/∂t = dsig_k·∂H/∂t`. Layer continuity is
      !!   ∂h_k/∂t + ∂(hu_k)/∂x + ∂(hv_k)/∂y = w(k+1) − w(k),
      !! so the flux across the MOVING interface is
      !!   w(k+1) = w(k) − [div(hu_k) + dsig_k·∂H/∂t],
      !! with the column mass balance ∂H/∂t = −Σ_m div(hu_m) = −div_col.
      !! Substituting gives the implemented increment
      !!   w(k+1) = w(k) − div(hu_k) + dsig_k·div_col.
      !!
      !! Two consequences vs. the former `w(k+1)=w(k)−div(hu_k)` form (which
      !! omitted the moving-coordinate term and so reported `w_phys`):
      !!   * Under pure barotropic sigma flow (hu_k = dsig_k·hu_bt) every
      !!     interface flux is 0 to round-off — no spurious diapycnal
      !!     transport (was ~1e-3, drove a salt leak in the exchange).
      !!   * The free surface w(:,:,nz+1) ≡ 0 by construction
      !!     (Σ_k dsig_k·div_col = div_col cancels div_col) — no diasurface
      !!     flux through the moving top, the correct kinematic BC.
      !! Baroclinic (sheared) flow still yields sensible nonzero interior w.
      !!
      !! Boundary condition: w(:,:,1) = 0 (no-penetration at the bed).
      !!
      !! NOTE: the `dsig_k·∂H/∂t` interface-motion model is exact for pure
      !! sigma layering (the coastal default `VCOORD_SIGMA`); for hybrid
      !! z/z* coordinates the interface velocity is only approximated here.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(out) :: w(grid%nx_total, grid%ny_total, nz + 1)

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: div_hu, div_hv, div_col, inv_2dx, inv_2dy, dx, dy
      real(wp) :: h_col, h_col_im1, h_col_ip1, h_col_jm1, h_col_jp1

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      inv_2dx = 0.5_wp/dx
      inv_2dy = 0.5_wp/dy

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, kk, div_hu, div_hv, div_col, &
               h_col, h_col_im1, h_col_ip1, h_col_jm1, h_col_jp1)
         w(i, j, 1) = 0.0_wp
         ! Compute total column depth for this cell and its neighbours
         h_col = 0.0_wp
         h_col_im1 = 0.0_wp
         h_col_ip1 = 0.0_wp
         h_col_jm1 = 0.0_wp
         h_col_jp1 = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
            h_col_im1 = h_col_im1 + h_layer(i - 1, j, kk)
            h_col_ip1 = h_col_ip1 + h_layer(i + 1, j, kk)
            h_col_jm1 = h_col_jm1 + h_layer(i, j - 1, kk)
            h_col_jp1 = h_col_jp1 + h_layer(i, j + 1, kk)
         end do
         ! Skip dry columns and columns adjacent to dry cells where
         ! central differences would produce spurious divergence
         if (all_neighbours_wet(h_col, h_col_ip1, h_col_im1, &
                                h_col_jp1, h_col_jm1)) then
            ! Column-integrated horizontal transport divergence
            ! div_col = Σ_m [div(hu_m) + div(hv_m)] = −∂H/∂t.
            div_col = 0.0_wp
            do kk = 1, nz
               div_col = div_col &
                         + (hu_layer(i + 1, j, kk) - hu_layer(i - 1, j, kk))*inv_2dx &
                         + (hv_layer(i, j + 1, kk) - hv_layer(i, j - 1, kk))*inv_2dy
            end do
            do k = 1, nz
               div_hu = (hu_layer(i + 1, j, k) - hu_layer(i - 1, j, k))*inv_2dx
               div_hv = (hv_layer(i, j + 1, k) - hv_layer(i, j - 1, k))*inv_2dy
               ! Diasurface increment: subtract the moving-interface term
               ! dsig_k·∂H/∂t = −dsig_k·div_col.
               w(i, j, k + 1) = w(i, j, k) - (div_hu + div_hv) + dsig(k)*div_col
            end do
         else
            do k = 1, nz
               w(i, j, k + 1) = 0.0_wp
            end do
         end if
      end do

   end subroutine ml_compute_vertical_velocity

   pure subroutine ml_compute_vertical_exchange(h_layer, hu_layer, hv_layer, &
                                                w, grid, nz, dt)
      !! Vertical momentum exchange using upwind advection at layer interfaces.
      !! Transfers momentum between adjacent layers based on diagnostic w.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: w_face, u_up, v_up, flux_u, flux_v, h_val, h_col

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      ! Single-layer case: no vertical exchange needed
      if (nz < 2) return

      ! Each thread processes its full column — no cross-column dependency
      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, kk, w_face, u_up, v_up, flux_u, flux_v, h_val, h_col)
         ! Compute total column depth
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         ! Skip thin columns entirely — vertical exchange is meaningless
         ! when the water column is near-dry
         if (h_col > ML_MIN_DEPTH) then
            do k = 1, nz - 1
               w_face = w(i, j, k + 1)

               ! Upwind velocity selection
               if (w_face > 0.0_wp) then
                  h_val = h_layer(i, j, k)
                  if (h_val > DRY_TOLERANCE) then
                     u_up = hu_layer(i, j, k)/h_val
                     v_up = hv_layer(i, j, k)/h_val
                  else
                     u_up = 0.0_wp
                     v_up = 0.0_wp
                  end if
               else
                  h_val = h_layer(i, j, k + 1)
                  if (h_val > DRY_TOLERANCE) then
                     u_up = hu_layer(i, j, k + 1)/h_val
                     v_up = hv_layer(i, j, k + 1)/h_val
                  else
                     u_up = 0.0_wp
                     v_up = 0.0_wp
                  end if
               end if

               ! Flux = w * u_up * dt (mass flux * velocity)
               ! Limit: don't extract more than 50% of donor layer momentum
               flux_u = w_face*u_up*dt
               flux_v = w_face*v_up*dt
               if (w_face > 0.0_wp) then
                  ! Donor is layer k
                  if (abs(flux_u) > 0.5_wp*abs(hu_layer(i, j, k))) then
                     flux_u = sign(0.5_wp*abs(hu_layer(i, j, k)), flux_u)
                  end if
                  if (abs(flux_v) > 0.5_wp*abs(hv_layer(i, j, k))) then
                     flux_v = sign(0.5_wp*abs(hv_layer(i, j, k)), flux_v)
                  end if
               else
                  ! Donor is layer k+1
                  if (abs(flux_u) > 0.5_wp*abs(hu_layer(i, j, k + 1))) then
                     flux_u = sign(0.5_wp*abs(hu_layer(i, j, k + 1)), flux_u)
                  end if
                  if (abs(flux_v) > 0.5_wp*abs(hv_layer(i, j, k + 1))) then
                     flux_v = sign(0.5_wp*abs(hv_layer(i, j, k + 1)), flux_v)
                  end if
               end if

               hu_layer(i, j, k) = hu_layer(i, j, k) - flux_u
               hu_layer(i, j, k + 1) = hu_layer(i, j, k + 1) + flux_u
               hv_layer(i, j, k) = hv_layer(i, j, k) - flux_v
               hv_layer(i, j, k + 1) = hv_layer(i, j, k + 1) + flux_v
            end do
         end if
      end do

   end subroutine ml_compute_vertical_exchange

   pure subroutine ml_rk2_average_layers(h_layer, hu_layer, hv_layer, &
                                         h_layer0, hu_layer0, hv_layer0, &
                                         grid, nz)
      !! SSP-RK2 averaging: U^{n+1} = 0.5*(U^n + U**)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer0(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny, nghost

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(k=1:nz, j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost)
         h_layer(i, j, k) = 0.5_wp*(h_layer0(i, j, k) + h_layer(i, j, k))
         hu_layer(i, j, k) = 0.5_wp*(hu_layer0(i, j, k) + hu_layer(i, j, k))
         hv_layer(i, j, k) = 0.5_wp*(hv_layer0(i, j, k) + hv_layer(i, j, k))

         if (h_layer(i, j, k) < 0.0_wp) then
            h_layer(i, j, k) = 0.0_wp
            hu_layer(i, j, k) = 0.0_wp
            hv_layer(i, j, k) = 0.0_wp
         else if (h_layer(i, j, k) < THIN_LAYER_THRESHOLD) then
            hu_layer(i, j, k) = 0.0_wp
            hv_layer(i, j, k) = 0.0_wp
         end if
      end do

   end subroutine ml_rk2_average_layers

   pure subroutine ml_init_layers_from_barotropic(h, hu, hv, dsig, &
                                                  h_layer, hu_layer, hv_layer, &
                                                  grid, nz)
      !! Initialise layer variables from barotropic state.
      !! h_k = dsig_k * H, hu_k = dsig_k * hu, hv_k = dsig_k * hv
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(out) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hv_layer(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         h_layer(i, j, k) = dsig(k)*h(i, j)
         hu_layer(i, j, k) = dsig(k)*hu(i, j)
         hv_layer(i, j, k) = dsig(k)*hv(i, j)
      end do

   end subroutine ml_init_layers_from_barotropic

   pure subroutine ml_save_layer_state(h_layer, hu_layer, hv_layer, &
                                       h_layer0, hu_layer0, hv_layer0, &
                                       grid, nz)
      !! Save layer state for RK2 averaging
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: h_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hu_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hv_layer0(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      do concurrent(k=1:nz, j=1:ny, i=1:nx)
         h_layer0(i, j, k) = h_layer(i, j, k)
         hu_layer0(i, j, k) = hu_layer(i, j, k)
         hv_layer0(i, j, k) = hv_layer(i, j, k)
      end do

   end subroutine ml_save_layer_state

   pure subroutine ml_apply_layer_bcs(h_layer, hu_layer, hv_layer, &
                                      h, hu, hv, dsig, &
                                      grid, nz, &
                                      bc_west, bc_east, bc_south, bc_north)
      !! Enforce layer consistency in ghost cells.
      !! Call after applying barotropic BCs and before flux computation.
      !!
      !! Two regimes:
      !!
      !! * **Wall edges** (BC type supplied and equal to `BC_WALL`) mirror the
      !!   interior PER-LAYER momentum profile, normal component negated and
      !!   tangential mirrored — matching the barotropic wall convention in
      !!   `rki_boundary`.  This makes a closed wall tracer-TIGHT.
      !! * **Every other edge** distributes the barotropic ghost state
      !!   proportionally to the sigma fractions, as before.  There is no
      !!   interior profile to mirror meaningfully for inflow/tidal/open.
      !!
      !! ## Why the wall case needs the per-layer mirror
      !!
      !! The old code always did the `dsig` distribution, i.e.
      !! `hu_layer(g,j,k) = dsig(k)*hu(g,j)` — a SHEAR-FREE ghost.  At a wall
      !! the barotropic ghost carries `hu(g,j) = -hu(mirror,j)`, so the layer
      !! momenta sum correctly to zero net barotropic flux and WATER mass is
      !! conserved to machine precision.  But the per-layer face velocity
      !! `v_face(k) = 0.5*(hu_ghost(k) + hu_int(k))` is then NON-ZERO in every
      !! layer whenever the interior column is sheared: the ghost's uniform
      !! profile cannot cancel the interior's sheared one layer by layer.
      !! Each layer therefore trades tracer with the ghost at its own
      !! concentration, and the net tracer flux through a CLOSED wall is
      !! `sum_k v(k)*Tr(k) /= 0` whenever Tr varies with depth.
      !!
      !! Measured before the fix, in a single `ml_advect_tracer` call with a
      !! sheared interior column (+q/-q) and Tr = 40/35 by layer beside a
      !! wall: water-mass change exactly 0.0, tracer change 2.000000000044e-2
      !! against the analytic prediction `ny*dt/dx*(q/2)*(40-35)` = 2.0e-2 —
      !! an exact match, confirming the mechanism.  With the per-layer mirror
      !! `v_face(k)` vanishes identically and the wall is tight layer by layer.
      !!
      !! The BC arguments are OPTIONAL so existing call sites that do not (or
      !! cannot) supply them keep the previous behaviour; callers that know
      !! their edge types should pass them.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dsig(nz)
      integer, intent(in), optional :: bc_west, bc_east, bc_south, bc_north

      integer :: i, j, k, g, m, nx, ny, nghost
      logical :: wall_w, wall_e, wall_s, wall_n

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      wall_w = .false.
      wall_e = .false.
      wall_s = .false.
      wall_n = .false.
      if (present(bc_west)) wall_w = (bc_west == BC_WALL)
      if (present(bc_east)) wall_e = (bc_east == BC_WALL)
      if (present(bc_south)) wall_s = (bc_south == BC_WALL)
      if (present(bc_north)) wall_n = (bc_north == BC_WALL)

      ! West ghost cells.  Mirror pairing matches `rki_boundary`:
      ! ghost (nghost + 1 - g)  <->  interior (nghost + g),  g = 1..nghost.
      do g = 1, nghost
         m = 2*nghost + 1 - g
         if (wall_w) then
            do concurrent(j=1:ny, k=1:nz)
               h_layer(g, j, k) = h_layer(m, j, k)
               hu_layer(g, j, k) = -hu_layer(m, j, k)
               hv_layer(g, j, k) = hv_layer(m, j, k)
            end do
         else
            do concurrent(j=1:ny, k=1:nz)
               h_layer(g, j, k) = dsig(k)*h(g, j)
               hu_layer(g, j, k) = dsig(k)*hu(g, j)
               hv_layer(g, j, k) = dsig(k)*hv(g, j)
            end do
         end if
      end do

      ! East ghost cells
      do g = nx - nghost + 1, nx
         m = 2*(nx - nghost) + 1 - g
         if (wall_e) then
            do concurrent(j=1:ny, k=1:nz)
               h_layer(g, j, k) = h_layer(m, j, k)
               hu_layer(g, j, k) = -hu_layer(m, j, k)
               hv_layer(g, j, k) = hv_layer(m, j, k)
            end do
         else
            do concurrent(j=1:ny, k=1:nz)
               h_layer(g, j, k) = dsig(k)*h(g, j)
               hu_layer(g, j, k) = dsig(k)*hu(g, j)
               hv_layer(g, j, k) = dsig(k)*hv(g, j)
            end do
         end if
      end do

      ! South ghost cells
      do g = 1, nghost
         m = 2*nghost + 1 - g
         if (wall_s) then
            do concurrent(i=1:nx, k=1:nz)
               h_layer(i, g, k) = h_layer(i, m, k)
               hu_layer(i, g, k) = hu_layer(i, m, k)
               hv_layer(i, g, k) = -hv_layer(i, m, k)
            end do
         else
            do concurrent(i=1:nx, k=1:nz)
               h_layer(i, g, k) = dsig(k)*h(i, g)
               hu_layer(i, g, k) = dsig(k)*hu(i, g)
               hv_layer(i, g, k) = dsig(k)*hv(i, g)
            end do
         end if
      end do

      ! North ghost cells
      do g = ny - nghost + 1, ny
         m = 2*(ny - nghost) + 1 - g
         if (wall_n) then
            do concurrent(i=1:nx, k=1:nz)
               h_layer(i, g, k) = h_layer(i, m, k)
               hu_layer(i, g, k) = hu_layer(i, m, k)
               hv_layer(i, g, k) = -hv_layer(i, m, k)
            end do
         else
            do concurrent(i=1:nx, k=1:nz)
               h_layer(i, g, k) = dsig(k)*h(i, g)
               hu_layer(i, g, k) = dsig(k)*hu(i, g)
               hv_layer(i, g, k) = dsig(k)*hv(i, g)
            end do
         end if
      end do

   end subroutine ml_apply_layer_bcs

   subroutine ml_apply_layer_bcs_full(h_layer, hu_layer, hv_layer, &
                                      h, hu, hv, dsig, &
                                      vcoord_z_ref_col, h_min, &
                                      grid, nz)
      !! VCOORD_ZSTAR_FULL variant of `ml_apply_layer_bcs`. Rebuilds
      !! ghost-cell `h_layer` from the per-column `vcoord_z_ref_col` profile
      !! (as the interior `ml_distribute_with_remap_full`) so the layout is
      !! consistent across the boundary — a sigma fallback would mismatch and
      !! make CWC advection inject tracer at the wrong per-layer ratio.
      !! Momentum is still split by `dsig(k)` (ghost barotropic momentum
      !! carries no internal shear).
      use rki_vcoord, only: vcoord_target_dz_column_zstar_full
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: vcoord_z_ref_col(0:nz, grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: h_min

      integer :: i, j, k, g, nx, ny, nghost
      real(wp) :: dz_col(NZ_STACK_MAX)
      real(wp) :: z_ref_local(0:NZ_STACK_MAX)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      ! West ghost cells
      !$acc parallel loop collapse(2) &
      !$acc&   present(h_layer, hu_layer, hv_layer, h, hu, hv, dsig, &
      !$acc&           vcoord_z_ref_col) &
      !$acc&   private(k, dz_col, z_ref_local)
      do g = 1, nghost
         do j = 1, ny
            do k = 0, nz
               z_ref_local(k) = vcoord_z_ref_col(k, g, j)
            end do
            call vcoord_target_dz_column_zstar_full(nz, h(g, j), &
                                                    z_ref_local(0:nz), &
                                                    h_min, dz_col)
            do k = 1, nz
               h_layer(g, j, k) = dz_col(k)
               hu_layer(g, j, k) = dsig(k)*hu(g, j)
               hv_layer(g, j, k) = dsig(k)*hv(g, j)
            end do
         end do
      end do

      ! East ghost cells
      !$acc parallel loop collapse(2) &
      !$acc&   present(h_layer, hu_layer, hv_layer, h, hu, hv, dsig, &
      !$acc&           vcoord_z_ref_col) &
      !$acc&   private(k, dz_col, z_ref_local)
      do g = nx - nghost + 1, nx
         do j = 1, ny
            do k = 0, nz
               z_ref_local(k) = vcoord_z_ref_col(k, g, j)
            end do
            call vcoord_target_dz_column_zstar_full(nz, h(g, j), &
                                                    z_ref_local(0:nz), &
                                                    h_min, dz_col)
            do k = 1, nz
               h_layer(g, j, k) = dz_col(k)
               hu_layer(g, j, k) = dsig(k)*hu(g, j)
               hv_layer(g, j, k) = dsig(k)*hv(g, j)
            end do
         end do
      end do

      ! South ghost cells
      !$acc parallel loop collapse(2) &
      !$acc&   present(h_layer, hu_layer, hv_layer, h, hu, hv, dsig, &
      !$acc&           vcoord_z_ref_col) &
      !$acc&   private(k, dz_col, z_ref_local)
      do g = 1, nghost
         do i = 1, nx
            do k = 0, nz
               z_ref_local(k) = vcoord_z_ref_col(k, i, g)
            end do
            call vcoord_target_dz_column_zstar_full(nz, h(i, g), &
                                                    z_ref_local(0:nz), &
                                                    h_min, dz_col)
            do k = 1, nz
               h_layer(i, g, k) = dz_col(k)
               hu_layer(i, g, k) = dsig(k)*hu(i, g)
               hv_layer(i, g, k) = dsig(k)*hv(i, g)
            end do
         end do
      end do

      ! North ghost cells
      !$acc parallel loop collapse(2) &
      !$acc&   present(h_layer, hu_layer, hv_layer, h, hu, hv, dsig, &
      !$acc&           vcoord_z_ref_col) &
      !$acc&   private(k, dz_col, z_ref_local)
      do g = ny - nghost + 1, ny
         do i = 1, nx
            do k = 0, nz
               z_ref_local(k) = vcoord_z_ref_col(k, i, g)
            end do
            call vcoord_target_dz_column_zstar_full(nz, h(i, g), &
                                                    z_ref_local(0:nz), &
                                                    h_min, dz_col)
            do k = 1, nz
               h_layer(i, g, k) = dz_col(k)
               hu_layer(i, g, k) = dsig(k)*hu(i, g)
               hv_layer(i, g, k) = dsig(k)*hv(i, g)
            end do
         end do
      end do

   end subroutine ml_apply_layer_bcs_full

   subroutine remap_workspace_ensure(nx, ny, nz)
      integer, intent(in) :: nx, ny, nz
      logical :: need_alloc

      need_alloc = .false.
      if (.not. allocated(remap_dz_old)) then
         need_alloc = .true.
      else if (size(remap_dz_old, 1) /= nx .or. &
               size(remap_dz_old, 2) /= ny .or. &
               size(remap_dz_old, 3) /= nz) then
         call remap_workspace_cleanup()
         need_alloc = .true.
      end if

      if (need_alloc) then
         allocate (remap_dz_old(nx, ny, nz))
         !$acc enter data create(remap_dz_old)
      end if
   end subroutine remap_workspace_ensure

   subroutine remap_workspace_cleanup()
      if (.not. allocated(remap_dz_old)) return
      !$acc exit data delete(remap_dz_old)
      deallocate (remap_dz_old)
   end subroutine remap_workspace_cleanup

   pure logical function all_neighbours_wet(h_col, h_col_ip1, h_col_im1, &
                                            h_col_jp1, h_col_jm1) result(all_wet)
      !$acc routine seq
      !! Returns `.true.` if the centre cell and its four face-neighbours all
      !! have total column depth >= `ML_MIN_DEPTH`. Guards central-difference
      !! operators (BPG, divergence) against spurious gradients at dry cells.
      real(wp), intent(in) :: h_col     !! column depth at (i,j)
      real(wp), intent(in) :: h_col_ip1  !! column depth at (i+1,j)
      real(wp), intent(in) :: h_col_im1  !! column depth at (i-1,j)
      real(wp), intent(in) :: h_col_jp1  !! column depth at (i,j+1)
      real(wp), intent(in) :: h_col_jm1  !! column depth at (i,j-1)
      all_wet = h_col >= ML_MIN_DEPTH .and. &
                h_col_ip1 >= ML_MIN_DEPTH .and. &
                h_col_im1 >= ML_MIN_DEPTH .and. &
                h_col_jp1 >= ML_MIN_DEPTH .and. &
                h_col_jm1 >= ML_MIN_DEPTH
   end function all_neighbours_wet

   subroutine ml_build_vcoord_zstar_full(b, vcoord_z_ref_col, grid, nz, &
                                         h_surf_target, n_surf, stretching)
      !! Fill the per-column reference z-level array from the bathymetry.
      !! Host-side loop (runs at init / bathymetry change); updates the
      !! device copy after writing.
      type(hgrid_t), intent(in)  :: grid
      integer, intent(in)        :: nz
      real(wp), intent(in)       :: b(grid%nx_total, grid%ny_total)
         !! Bed ELEVATION (`state%barotropic%b`): negative for submerged
         !! cells; negated here to the POSITIVE bed depth that
         !! `zstar_full_build_column` expects.
      real(wp), intent(inout)    :: vcoord_z_ref_col(0:nz, grid%nx_total, grid%ny_total)
      real(wp), intent(in)       :: h_surf_target
         !! Target surface-layer thickness (m), see zstar_full_build_column
      integer, intent(in)        :: n_surf
      integer, intent(in)        :: stretching

      integer :: i, j

      !$acc update self(b) if_present

      do j = 1, grid%ny_total
         do i = 1, grid%nx_total
            ! Negate to bed depth and clamp non-negative: LAND cells (b > 0)
            ! get a zero-depth degenerate column instead of negative z_ref
            ! interfaces, which would propagate to negative dz_new and NaN CFL.
            call zstar_full_build_column(nz, max(-b(i, j), 0.0_wp), &
                                         h_surf_target, &
                                         n_surf, stretching, &
                                         vcoord_z_ref_col(:, i, j))
         end do
      end do

      !$acc update device(vcoord_z_ref_col) if_present
   end subroutine ml_build_vcoord_zstar_full

   subroutine ml_dynamics_workspaces_cleanup()
      !! Release the vcoord-remap persistent workspace.  Called by the
      !! solver at simulation end.
      call remap_workspace_cleanup()
   end subroutine ml_dynamics_workspaces_cleanup

end module rki_ml_dynamics
