!! Fused per-cell post-update kernels for the multilayer structured solver.
!! Collapse the per-stage post-flux kernel chain (update → distribute →
!! bottom friction → truncate → barotropic momentum → wind/pressure/precip
!! → coriolis) into one launch per stage; no horizontal neighbour reads
!! except the atmospheric pressure gradient.
!! Caller must run tracer advection BEFORE these kernels (CWC needs OLD h_layer).
module rki_kernel_fused_ml
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            DESING_EPS, FROUDE_BACKSTOP, RHO_WATER, &
                            LATENT_HEAT_VAPORIZATION
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: ml_fused_post_update
   public :: ml_fused_rk2_finalise
   public :: ml_fused_tracer_post

contains

   pure subroutine ml_fused_tracer_post(hS_layer, hT_layer, h_layer, &
                                        hu_layer, hv_layer, &
                                        h_bt, hu_bt, hv_bt, &
                                        heat_flux, evap_minus_precip, &
                                        evap_mass_sink, evap_latent_cooling, &
                                        S_min, S_max, T_min, T_max, &
                                        rho_0, cp, min_depth, &
                                        grid, nz, dt)
      !! Fuses surface tracer fluxes + per-layer clamp (S, T) into one
      !! per-(i,j)-column kernel (no horizontal neighbour reads).
      !!
      !! Surface flux (top layer k=nz, wet column only):
      !! - evap_mass_sink=.false. (default): virtual salt flux
      !!     hT(top) += dt * heat_flux / (rho_0 * cp);  hS(top) += dt * (E-P) * S_top
      !! - evap_mass_sink=.true.: real water-mass sink — dh = E*dt removed
      !!   COLUMN-PROPORTIONALLY (h, hu, hv, hT scale by 1 - dh/h_col; hS
      !!   held), so S concentrates while T, u, v are exactly invariant.
      !!   With evap_latent_cooling, the evaporated mass ALSO withdraws its
      !!   latent heat from the top layer (additive surface cooling; distinct
      !!   from the sensible heat carried away by the fac-scaled hT — no
      !!   double-count).  Removal clamps to the water available above
      !!   min_depth; the unsatisfied fraction is dropped honestly.  Precipitation
      !!   (evap_minus_precip < 0) is NOT added — that is the barotropic
      !!   precipitation field's job.  See ml_apply_surface_tracer_fluxes
      !!   for the full derivation.
      !!
      !! Clamp: bounds [tr_min, tr_max] per layer where h_layer > DRY_TOLERANCE,
      !! else zeroes the column tracer mass.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: h_bt(grid%nx_total, grid%ny_total)
         !! Barotropic depth — the nz=1 fast path maintains h_bt == h_layer
         !! BY CONSTRUCTION (no recompute-from-layers), so the mass sink
         !! must update both or the next stage rebuilds the pre-sink depth.
      real(wp), intent(inout) :: hu_bt(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv_bt(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: heat_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: evap_minus_precip(grid%nx_total, grid%ny_total)
      logical, intent(in) :: evap_mass_sink
      logical, intent(in) :: evap_latent_cooling
      real(wp), intent(in) :: S_min, S_max, T_min, T_max
      real(wp), intent(in) :: rho_0, cp, dt, min_depth

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_col, h_top, S_top, kin_heat_flux, inv_rho0cp
      real(wp) :: h_k, hs_k, ht_k, s_val, t_val
      real(wp) :: dh_req, dh_eff, fac, latent_gain
      logical :: column_wet

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      inv_rho0cp = 1.0_wp/(rho_0*cp)
      ! Kinematic latent-heat withdrawal per metre of evaporated water
      ! (0 unless the knob is on -> default bit-identical).
      latent_gain = 0.0_wp
      if (evap_latent_cooling) latent_gain = RHO_WATER*LATENT_HEAT_VAPORIZATION*inv_rho0cp

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, h_col, h_top, S_top, kin_heat_flux, h_k, hs_k, ht_k, &
               s_val, t_val, column_wet, dh_req, dh_eff, fac)

         ! ---- Total column depth (surface-flux wet gate) ----
         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + h_layer(i, j, k)
         end do
         column_wet = (h_col > min_depth)

         ! ---- Surface flux (top layer, k=nz under bottom-up ROMS indexing) ----
         if (column_wet) then
            h_top = h_layer(i, j, nz)
            if (h_top > THIN_LAYER_THRESHOLD) then
               ! Thermal (unchanged in both modes): dt * Q_net / (rho_0 * cp)
               kin_heat_flux = heat_flux(i, j)*inv_rho0cp
               hT_layer(i, j, nz) = hT_layer(i, j, nz) + dt*kin_heat_flux

               if (evap_mass_sink) then
                  ! Column-proportional water-mass sink (see
                  ! ml_apply_surface_tracer_fluxes).  hS held; h, hu, hv,
                  ! hT scale so S concentrates while T, u, v stay invariant.
                  dh_req = dt*max(evap_minus_precip(i, j), 0.0_wp)
                  dh_eff = min(dh_req, max(h_col - min_depth, 0.0_wp))
                  if (dh_eff > 0.0_wp) then
                     fac = 1.0_wp - dh_eff/h_col
                     do k = 1, nz
                        h_layer(i, j, k) = h_layer(i, j, k)*fac
                        hu_layer(i, j, k) = hu_layer(i, j, k)*fac
                        hv_layer(i, j, k) = hv_layer(i, j, k)*fac
                        hT_layer(i, j, k) = hT_layer(i, j, k)*fac
                        ! hS_layer held -> concentration rises naturally.
                     end do
                     ! Keep the barotropic state consistent with the layers
                     ! (the nz=1 path never recomputes it from the layers).
                     h_bt(i, j) = h_bt(i, j)*fac
                     hu_bt(i, j) = hu_bt(i, j)*fac
                     hv_bt(i, j) = hv_bt(i, j)*fac
                     ! Latent-heat cooling (additive, top layer; not scaled
                     ! by fac).  latent_gain == 0 when off -> bit-identical.
                     hT_layer(i, j, nz) = hT_layer(i, j, nz) - latent_gain*dh_eff
                  end if
               else
                  ! Legacy virtual salt flux (bit-identical default).
                  S_top = hS_layer(i, j, nz)/h_top
                  hS_layer(i, j, nz) = hS_layer(i, j, nz) &
                                       + dt*evap_minus_precip(i, j)*S_top
               end if
            end if
         end if

         ! ---- Clamp every layer to physical bounds ----
         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k > DRY_TOLERANCE) then
               hs_k = hS_layer(i, j, k)
               ht_k = hT_layer(i, j, k)

               s_val = hs_k/h_k
               s_val = max(S_min, min(S_max, s_val))
               hS_layer(i, j, k) = h_k*s_val

               t_val = ht_k/h_k
               t_val = max(T_min, min(T_max, t_val))
               hT_layer(i, j, k) = h_k*t_val
            else
               hS_layer(i, j, k) = 0.0_wp
               hT_layer(i, j, k) = 0.0_wp
            end if
         end do
      end do

   end subroutine ml_fused_tracer_post

   pure subroutine ml_fused_rk2_finalise(h, hu, hv, &
                                         h_layer, hu_layer, hv_layer, &
                                         h_layer0, hu_layer0, hv_layer0, &
                                         hS_layer, hT_layer, &
                                         hS_layer0, hT_layer0, &
                                         grid, nz)
      !! End-of-step fused kernel: RK2-averages all layer arrays with their
      !! `*_0` counterparts and accumulates the depth-integrated `h, hu, hv`
      !! in one inner k-loop. No neighbour reads.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(out) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hS_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hT_layer0(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_k, hu_k, hv_k
      real(wp) :: hs_k, ht_k
      real(wp) :: h_sum, hu_sum, hv_sum

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, h_k, hu_k, hv_k, hs_k, ht_k, h_sum, hu_sum, hv_sum)

         h_sum = 0.0_wp
         hu_sum = 0.0_wp
         hv_sum = 0.0_wp

         do k = 1, nz
            ! ---- RK2 average of layer state ----
            h_k = 0.5_wp*(h_layer0(i, j, k) + h_layer(i, j, k))
            hu_k = 0.5_wp*(hu_layer0(i, j, k) + hu_layer(i, j, k))
            hv_k = 0.5_wp*(hv_layer0(i, j, k) + hv_layer(i, j, k))

            ! Wet-cell guard, mirrors ml_rk2_average_layers
            if (h_k < 0.0_wp) then
               h_k = 0.0_wp
               hu_k = 0.0_wp
               hv_k = 0.0_wp
            else if (h_k < THIN_LAYER_THRESHOLD) then
               hu_k = 0.0_wp
               hv_k = 0.0_wp
            end if

            h_layer(i, j, k) = h_k
            hu_layer(i, j, k) = hu_k
            hv_layer(i, j, k) = hv_k

            ! Accumulate barotropic quantities
            h_sum = h_sum + h_k
            hu_sum = hu_sum + hu_k
            hv_sum = hv_sum + hv_k

            ! ---- Tracer RK2 average ----
            hs_k = 0.5_wp*(hS_layer0(i, j, k) + hS_layer(i, j, k))
            ht_k = 0.5_wp*(hT_layer0(i, j, k) + hT_layer(i, j, k))
            if (hs_k < 0.0_wp) hs_k = 0.0_wp
            if (ht_k < 0.0_wp) ht_k = 0.0_wp
            hS_layer(i, j, k) = hs_k
            hT_layer(i, j, k) = ht_k
         end do

         h(i, j) = h_sum
         hu(i, j) = hu_sum
         hv(i, j) = hv_sum
      end do

   end subroutine ml_fused_rk2_finalise

   pure subroutine ml_fused_post_update(h, hu, hv, &
                                        h_layer, hu_layer, hv_layer, &
                                        h0, hu0, hv0, &
                                        h_layer0, hu_layer0, hv_layer0, &
                                        flux_h, flux_hu, flux_hv, &
                                        dsig, &
                                        wind_x, wind_y, &
                                        atm_pressure, precipitation, &
                                        grid, nz, dt, &
                                        manning_n, use_linear_drag, drag_r, coriolis_f, &
                                        forcing_active, do_save)
      !! Single fused kernel covering compute_update through coriolis.
      !!
      !! Manning friction is applied *only* to the bottom layer (k=nz),
      !! mirroring `ml_apply_bottom_friction`. The depth-averaged Manning
      !! pass that the un-fused path does inside `apply_friction` is
      !! omitted intentionally: per-layer friction already accounts for
      !! it. Wind, pressure-gradient, and precipitation forcing are still
      !! applied to the depth-averaged state.
      !!
      !! Pre-update barotropic momentum is captured in registers at the
      !! top of each (i,j) iteration, so the un-fused path's separate
      !! `hu_ml_pre` / `hv_ml_pre` copies are not needed.
      !!
      !! When `do_save` is true (RK stage 1), absorbs the work that
      !! `solver_save_state` and `ml_save_layer_state` would have done:
      !! writes pre-update `h, hu, hv` to `h0, hu0, hv0` and pre-update
      !! `h_layer, hu_layer, hv_layer` to the corresponding `*0` arrays.
      !! Tracer save is still a separate kernel (it has to run before
      !! the tracer pipeline, which precedes this kernel).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      ! intent(inout) so that when do_save=.false. we can leave these
      ! arrays untouched (intent(out) would require full assignment).
      real(wp), intent(inout) :: h0(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu0(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv0(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: h_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: flux_h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: flux_hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: atm_pressure(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: precipitation(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt, manning_n
      logical, intent(in) :: use_linear_drag
      real(wp), intent(in) :: drag_r
      real(wp), intent(in) :: coriolis_f
      logical, intent(in) :: forcing_active, do_save

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: h_new, hu_new, hv_new, hu_pre, hv_pre
      real(wp) :: eps4, h4, u_ds, v_ds, speed_sq, fr_max_sq, fr_scale
      real(wp) :: du, dv, h_k, h_k_new, u_k, v_k
      real(wp) :: u_b, v_b, speed_b, cf, denom
      real(wp) :: u_abs, v_abs, vel_max_x, vel_max_y
      real(wp) :: hu_sum, hv_sum
      real(wp) :: dpdx, dpdy, hu_old, hv_old
      real(wp) :: n_sq, cosf, sinf
      logical :: has_friction, has_coriolis

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      n_sq = manning_n*manning_n
      has_friction = (manning_n > 0.0_wp) .or. (use_linear_drag .and. drag_r > 0.0_wp)
      has_coriolis = (abs(coriolis_f) > 1.0e-20_wp)
      if (has_coriolis) then
         cosf = cos(coriolis_f*dt)
         sinf = sin(coriolis_f*dt)
      end if

      vel_max_x = 0.9_wp*0.8_wp*dx/max(dt, 1.0e-10_wp)
      vel_max_y = 0.9_wp*0.8_wp*dy/max(dt, 1.0e-10_wp)

      ! KP07 desingularization length for the barotropic clamp, shared
      ! convention with compute_update / the CFL timestep estimator.
      eps4 = DESING_EPS**4

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, h_new, hu_new, hv_new, hu_pre, hv_pre, &
               h4, u_ds, v_ds, speed_sq, fr_max_sq, fr_scale, &
               du, dv, h_k, h_k_new, u_k, v_k, &
               u_b, v_b, speed_b, cf, denom, &
               u_abs, v_abs, hu_sum, hv_sum, &
               dpdx, dpdy, hu_old, hv_old)

         ! ---- Capture pre-update barotropic momentum in registers ----
         hu_pre = hu(i, j)
         hv_pre = hv(i, j)

         ! ---- Save pre-update barotropic state (RK stage 1 only) ----
         if (do_save) then
            h0(i, j) = h(i, j)
            hu0(i, j) = hu_pre
            hv0(i, j) = hv_pre
         end if

         ! ---- Euler update for the barotropic state ----
         h_new = h(i, j) + dt*flux_h(i, j)
         hu_new = hu_pre + dt*flux_hu(i, j)
         hv_new = hv_pre + dt*flux_hv(i, j)

         ! ---- Negative-depth clamp + KP07 desingularization ----
         ! (mirrors compute_update)
         if (h_new < 0.0_wp) then
            h_new = 0.0_wp
            hu_new = 0.0_wp
            hv_new = 0.0_wp
         else if (h_new < DRY_TOLERANCE) then
            hu_new = 0.0_wp
            hv_new = 0.0_wp
         else
            h4 = h_new**4
            if (h4 < eps4) then
               ! u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4)); hu := h*u.
               ! Identity for h >= eps, so applied only below eps.
               u_ds = sqrt(2.0_wp)*h_new*hu_new/sqrt(h4 + eps4)
               v_ds = sqrt(2.0_wp)*h_new*hv_new/sqrt(h4 + eps4)
               hu_new = h_new*u_ds
               hv_new = h_new*v_ds
            end if
            ! Froude backstop: blowup arrester only (triggers on a TRUE
            ! comparison, so NaN passes through and fails loudly).
            speed_sq = (hu_new**2 + hv_new**2)/h_new**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h_new
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu_new = hu_new*fr_scale
               hv_new = hv_new*fr_scale
            end if
         end if

         ! ---- Distribute barotropic acceleration to layers ----
         ! du = (hu_new - hu_pre) / h_new is the depth-averaged acceleration
         ! Per-layer u_k retains the existing baroclinic anomaly.
         if (h_new > DRY_TOLERANCE) then
            du = (hu_new - hu_pre)/h_new
            dv = (hv_new - hv_pre)/h_new
         else
            du = 0.0_wp
            dv = 0.0_wp
         end if

         do k = 1, nz
            h_k = h_layer(i, j, k)
            ! Save OLD layer state into *_0 arrays before we overwrite
            ! (RK stage 1 only — absorbs `ml_save_layer_state`).
            if (do_save) then
               h_layer0(i, j, k) = h_k
               hu_layer0(i, j, k) = hu_layer(i, j, k)
               hv_layer0(i, j, k) = hv_layer(i, j, k)
            end if
            if (h_k > DRY_TOLERANCE) then
               u_k = hu_layer(i, j, k)/h_k
               v_k = hv_layer(i, j, k)/h_k
            else
               u_k = 0.0_wp
               v_k = 0.0_wp
            end if

            h_k_new = dsig(k)*h_new
            h_layer(i, j, k) = h_k_new

            if (h_k_new < THIN_LAYER_THRESHOLD) then
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            else
               hu_layer(i, j, k) = h_k_new*(u_k + du)
               hv_layer(i, j, k) = h_k_new*(v_k + dv)
            end if
         end do

         ! ---- Bottom friction on layer 1 (semi-implicit) ----
         ! k=1 is the bottom layer under bottom-up ROMS indexing (matches
         ! BPG / rho_ref / NH init / standalone ml_apply_bottom_friction).
         ! Linear (Rayleigh) or quadratic Manning, mirroring the un-fused path.
         if (has_friction) then
            h_k = h_layer(i, j, 1)
            if (h_k >= THIN_LAYER_THRESHOLD .and. h_new >= THIN_LAYER_THRESHOLD) then
               if (use_linear_drag) then
                  denom = 1.0_wp + dt*drag_r
               else
                  u_b = hu_layer(i, j, 1)/h_k
                  v_b = hv_layer(i, j, 1)/h_k
                  speed_b = sqrt(u_b*u_b + v_b*v_b)
                  cf = GRAVITY*n_sq/(h_new**(1.0_wp/3.0_wp))
                  denom = 1.0_wp + dt*cf*speed_b/h_k
               end if
               hu_layer(i, j, 1) = hu_layer(i, j, 1)/denom
               hv_layer(i, j, 1) = hv_layer(i, j, 1)/denom
            end if
         end if

         ! ---- CFL-based velocity truncation (MOM6-style, cfl_trunc=0.8) ----
         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k >= THIN_LAYER_THRESHOLD) then
               u_abs = abs(hu_layer(i, j, k))/h_k
               v_abs = abs(hv_layer(i, j, k))/h_k
               if (u_abs > vel_max_x) then
                  hu_layer(i, j, k) = sign(vel_max_x*h_k, hu_layer(i, j, k))
               end if
               if (v_abs > vel_max_y) then
                  hv_layer(i, j, k) = sign(vel_max_y*h_k, hv_layer(i, j, k))
               end if
            end if
         end do

         ! ---- Recompute depth-integrated momentum from layers ----
         hu_sum = 0.0_wp
         hv_sum = 0.0_wp
         do k = 1, nz
            hu_sum = hu_sum + hu_layer(i, j, k)
            hv_sum = hv_sum + hv_layer(i, j, k)
         end do
         hu_new = hu_sum
         hv_new = hv_sum

         ! ---- Wind / atmospheric pressure / precipitation on barotropic ----
         ! (Manning was applied per-layer above so we skip it here.)
         if (forcing_active .and. h_new > DRY_TOLERANCE) then
            hu_new = hu_new + dt*wind_x(i, j)/RHO_WATER
            hv_new = hv_new + dt*wind_y(i, j)/RHO_WATER
            dpdx = (atm_pressure(i + 1, j) - atm_pressure(i - 1, j))/(2.0_wp*dx)
            dpdy = (atm_pressure(i, j + 1) - atm_pressure(i, j - 1))/(2.0_wp*dy)
            hu_new = hu_new - dt*h_new*dpdx/RHO_WATER
            hv_new = hv_new - dt*h_new*dpdy/RHO_WATER
         end if

         ! ---- Coriolis rotation ----
         if (has_coriolis .and. h_new > DRY_TOLERANCE) then
            hu_old = hu_new
            hv_old = hv_new
            hu_new = hu_old*cosf + hv_old*sinf
            hv_new = -hu_old*sinf + hv_old*cosf
         end if

         ! ---- Precipitation on barotropic h (all cells) ----
         if (forcing_active) then
            h_new = h_new + dt*precipitation(i, j)
         end if

         ! ---- Write back ----
         h(i, j) = h_new
         hu(i, j) = hu_new
         hv(i, j) = hv_new
      end do

   end subroutine ml_fused_post_update

end module rki_kernel_fused_ml
