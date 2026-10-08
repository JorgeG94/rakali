!! KPP (Large, McWilliams & Doney 1994) surface-boundary-layer mixing
!! for the unstructured (triangular) coastal multilayer path.
module rki_ml_kpp_unstr
   !! Unstructured port of the structured A-grid KPP overlay
   !! (`rki_ml_kpp`).  Same physics, adapted to the cell-centric
   !! `(nz, ncells)` layout and the two-pass (compute-then-apply)
   !! flux style the unstructured PP81 kernels use.  Like the structured
   !! path, the closure is *fused* with the explicit-flux apply because
   !! the unstructured solver has no implicit vertical-diffusion solve:
   !!
   !!   * `ml_kpp_boundary_layer_unstr` — per cell, computes u_* (from
   !!     the barotropic wind stress), the surface buoyancy flux B_0,
   !!     and the bulk-Richardson boundary-layer depth h_b (with the
   !!     LMD94 eq-23 V_t^2 term).  Stores all three on the multilayer
   !!     state; `kpp_bl_depth` persists across steps to seed V_t^2.
   !!
   !!   * `ml_kpp_viscosity_unstr` — momentum: nu_t = max(PP81 interior,
   !!     KPP shape-function overlay), CFL-clamped explicit flux to
   !!     hu/hv.  Drop-in for `ml_vertical_viscosity_unstr`.
   !!
   !!   * `ml_kpp_diffusion_tracers_unstr` — S + T: same overlay on the
   !!     PP81 cubic diffusivity, plus the non-local counter-gradient
   !!     transport folded into the interface flux.  Drop-in for the two
   !!     `ml_vertical_diffusion_tracer_unstr` calls.
   !!
   !! When the overlay is inactive (no wind + no destabilizing flux ⇒
   !! w_s = 0), every interface falls through to the bit-identical PP81
   !! arithmetic, so `vmix_use_kpp` off leaves existing results
   !! unchanged.  u_* comes from the spatially-uniform
   !! `state%barotropic%wind_tau_x/y` (the unstructured wind-stress
   !! coupling); a per-cell wind field is a future enhancement.
   !!
   !! Layer convention is bottom-up (k=1 bed, k=nz surface), matching
   !! the rest of the unstructured multilayer path.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, CP_WATER
   use rki_state_unstr, only: unstr_state_t
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR, BPG_NZ_MAX, &
                                  VANISHING_LAYER_TOL_UNSTR
   implicit none
   private

   public :: ml_kpp_boundary_layer_unstr
   public :: ml_kpp_viscosity_unstr
   public :: ml_kpp_diffusion_tracers_unstr

   real(wp), parameter :: EPS = 1.0e-20_wp
      !! Shear² floor in the interior gradient-Ri (matches PP81 kernels).
   real(wp), parameter :: SHEAR2_FLOOR = 1.0e-10_wp
      !! Velocity-difference floor in the bulk-Ri denominator (matches
      !! the structured KPP); stops Ri_bulk blowing up in quiescent water.

contains

   pure subroutine ml_kpp_boundary_layer_unstr(state, hS_layer)
      !! Pass 1: per-cell u_*, B_0, and the bulk-Ri boundary-layer depth.
      !! Mirrors `ml_kpp_boundary_layer` (structured) on the cell-centric
      !! layout.  u_* uses the uniform barotropic wind stress; B_0 uses
      !! the per-cell surface heat / virtual-salt fluxes.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — caller host-dereferences
      ! state%multilayer%tracers(idx)%hTr_layer; ncells varies per mesh partition.
      real(wp), intent(in) :: hS_layer(:, :)

      integer :: cell, k, nz, n
      real(wp) :: h_col, h_top, s_top, f_t_kin, f_s_kin, b_0, us
      real(wp) :: tau_x, tau_y, rho_0, alpha_T, beta_S, ri_crit, c_vt2
      real(wp) :: u_ref, v_ref, b_ref, d_centre_ref
      real(wp) :: h_b_lagged, wstar3_lagged, w_s_col
      real(wp) :: u_k, v_k, b_k, d_centre_k, delta_d, delta_b
      real(wp) :: n_brunt, vt2, shear2, ri_bulk, ri_prev, d_prev
      real(wp) :: denom, frac, h_b, h_k, h_k1
      logical :: crossing_found

      if (.not. state%use_multilayer) return
      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      if (nz < 2) return

      rho_0 = state%multilayer%rho_0
      alpha_T = state%multilayer%alpha_T
      beta_S = state%multilayer%beta_S
      ri_crit = state%multilayer%kpp_ri_crit
      c_vt2 = state%multilayer%kpp_c_vt2
      tau_x = state%barotropic%wind_tau_x
      tau_y = state%barotropic%wind_tau_y

      do concurrent(cell=1:n) &
         local(k, h_col, h_top, s_top, f_t_kin, f_s_kin, b_0, us, &
               u_ref, v_ref, b_ref, d_centre_ref, h_b_lagged, &
               wstar3_lagged, w_s_col, u_k, v_k, b_k, d_centre_k, &
               delta_d, delta_b, n_brunt, vt2, shear2, ri_bulk, ri_prev, &
               d_prev, denom, frac, h_b, h_k, h_k1, crossing_found)
         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + state%multilayer%h_layer(k, cell)
         end do
         if (h_col <= ML_MIN_DEPTH_UNSTR) then
            state%multilayer%kpp_u_star(cell) = 0.0_wp
            state%multilayer%kpp_buoy_flux(cell) = 0.0_wp
            state%multilayer%kpp_bl_depth(cell) = 0.0_wp
            cycle
         end if

         ! Surface friction velocity from the uniform wind stress (Pa).
         us = sqrt(sqrt(tau_x*tau_x + tau_y*tau_y)/rho_0)

         ! Surface buoyancy flux (dimensional EOS coefficients ⇒ the
         ! 1/rho_0 factor is explicit; see the structured kernel).
         h_top = state%multilayer%h_layer(nz, cell)
         if (h_top > DRY_TOLERANCE) then
            s_top = hS_layer(nz, cell)/h_top
            u_ref = state%multilayer%hu_layer(nz, cell)/h_top
            v_ref = state%multilayer%hv_layer(nz, cell)/h_top
         else
            s_top = 0.0_wp
            u_ref = 0.0_wp
            v_ref = 0.0_wp
         end if
         f_t_kin = state%forcing%heat_flux(cell)/(rho_0*CP_WATER)
         f_s_kin = state%forcing%evap_minus_precip(cell)*s_top
         b_0 = (GRAVITY/rho_0)*(alpha_T*f_t_kin - beta_S*f_s_kin)

         state%multilayer%kpp_u_star(cell) = us
         state%multilayer%kpp_buoy_flux(cell) = b_0

         b_ref = -GRAVITY*state%multilayer%rho_layer(nz, cell)/rho_0
         d_centre_ref = 0.5_wp*h_top

         h_b_lagged = state%multilayer%kpp_bl_depth(cell)
         wstar3_lagged = max(0.0_wp, -b_0)*h_b_lagged
         w_s_col = sqrt(us*us + wstar3_lagged**(2.0_wp/3.0_wp))

         ri_prev = 0.0_wp
         d_prev = d_centre_ref
         h_b = 0.0_wp
         crossing_found = .false.
         d_centre_k = d_centre_ref
         do k = nz - 1, 1, -1
            if (.not. crossing_found) then
               h_k = state%multilayer%h_layer(k, cell)
               h_k1 = state%multilayer%h_layer(k + 1, cell)
               d_centre_k = d_centre_k + 0.5_wp*(h_k1 + h_k)
               if (h_k > DRY_TOLERANCE) then
                  u_k = state%multilayer%hu_layer(k, cell)/h_k
                  v_k = state%multilayer%hv_layer(k, cell)/h_k
               else
                  u_k = 0.0_wp
                  v_k = 0.0_wp
               end if
               b_k = -GRAVITY*state%multilayer%rho_layer(k, cell)/rho_0

               delta_d = d_centre_k - d_centre_ref
               delta_b = b_ref - b_k
               n_brunt = sqrt(max(0.0_wp, delta_b/delta_d))
               vt2 = (c_vt2/ri_crit)*d_centre_k*n_brunt*w_s_col
               shear2 = max((u_k - u_ref)*(u_k - u_ref) + &
                            (v_k - v_ref)*(v_k - v_ref) + vt2, SHEAR2_FLOOR)
               ri_bulk = delta_b*delta_d/shear2

               if (ri_bulk >= ri_crit) then
                  denom = ri_bulk - ri_prev
                  if (abs(denom) > 1.0e-12_wp) then
                     frac = (ri_crit - ri_prev)/denom
                  else
                     frac = 0.5_wp
                  end if
                  h_b = d_prev + frac*(d_centre_k - d_prev)
                  crossing_found = .true.
               end if

               ri_prev = ri_bulk
               d_prev = d_centre_k
            end if
         end do

         if (.not. crossing_found) h_b = h_col
         state%multilayer%kpp_bl_depth(cell) = h_b
      end do
   end subroutine ml_kpp_boundary_layer_unstr

   pure subroutine ml_kpp_viscosity_unstr(state, dt)
      !! Pass 2 (momentum): PP81 interior viscosity with the KPP shape-
      !! function overlay max'd in, then the explicit CFL-clamped
      !! diffusive flux applied to hu/hv (two-pass, mirroring
      !! `ml_vertical_viscosity_unstr`).
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt

      integer :: cell, k, nz, n
      real(wp) :: h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1
      real(wp) :: du, dv, drho, Ri, denom, nu_t
      real(wp) :: h_col, h_b, b_0, us, wstar3, w_star, w_s
      real(wp) :: d_below, d_iface, sigma, g_shape, nu_kpp
      real(wp) :: flux_u_iface(BPG_NZ_MAX), flux_v_iface(BPG_NZ_MAX)

      if (.not. state%use_multilayer) return
      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      if (nz < 2) return

      do concurrent(cell=1:n) &
         local(k, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, du, dv, drho, &
               Ri, denom, nu_t, h_col, h_b, b_0, us, wstar3, w_star, w_s, &
               d_below, d_iface, sigma, g_shape, nu_kpp, &
               flux_u_iface, flux_v_iface)
         if (state%barotropic%height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + state%multilayer%h_layer(k, cell)
         end do
         h_b = state%multilayer%kpp_bl_depth(cell)
         b_0 = state%multilayer%kpp_buoy_flux(cell)
         us = state%multilayer%kpp_u_star(cell)
         wstar3 = max(0.0_wp, -b_0)*h_b
         w_star = wstar3**(1.0_wp/3.0_wp)
         w_s = sqrt(us*us + w_star*w_star)

         ! Pass 1: interface fluxes (PP81 + KPP overlay)
         d_below = 0.0_wp
         do k = 1, nz - 1
            h_k = state%multilayer%h_layer(k, cell)
            h_k1 = state%multilayer%h_layer(k + 1, cell)
            h_avg = 0.5_wp*(h_k + h_k1)
            d_below = d_below + h_k
            d_iface = h_col - d_below

            flux_u_iface(k) = 0.0_wp
            flux_v_iface(k) = 0.0_wp

            if (h_avg > VANISHING_LAYER_TOL_UNSTR .and. &
                h_k > VANISHING_LAYER_TOL_UNSTR .and. &
                h_k1 > VANISHING_LAYER_TOL_UNSTR) then
               if (h_k > DRY_TOLERANCE) then
                  u_k = state%multilayer%hu_layer(k, cell)/h_k
                  v_k = state%multilayer%hv_layer(k, cell)/h_k
               else
                  u_k = 0.0_wp
                  v_k = 0.0_wp
               end if
               if (h_k1 > DRY_TOLERANCE) then
                  u_k1 = state%multilayer%hu_layer(k + 1, cell)/h_k1
                  v_k1 = state%multilayer%hv_layer(k + 1, cell)/h_k1
               else
                  u_k1 = 0.0_wp
                  v_k1 = 0.0_wp
               end if

               du = u_k1 - u_k
               dv = v_k1 - v_k
               drho = state%multilayer%rho_layer(k, cell) - state%multilayer%rho_layer(k + 1, cell)
               Ri = GRAVITY*(drho/state%multilayer%rho_0)*h_avg/(du*du + dv*dv + EPS)
               if (Ri < 0.0_wp) Ri = 0.0_wp

               denom = 1.0_wp + state%multilayer%vmix_ri_alpha*Ri
               nu_t = state%multilayer%vmix_nu_max/(denom*denom) + state%multilayer%vmix_nu_bg

               ! KPP overlay (g_shape = 0 outside the BL ⇒ no-op).
               g_shape = 0.0_wp
               if (w_s > 0.0_wp .and. h_b > 0.0_wp .and. d_iface < h_b) then
                  sigma = d_iface/h_b
                  g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
               end if
               nu_kpp = h_b*w_s*g_shape
               nu_t = max(nu_t, nu_kpp)

               if (nu_t*dt > 0.5_wp*h_avg*h_avg) nu_t = 0.5_wp*h_avg*h_avg/dt

               flux_u_iface(k) = nu_t*du/h_avg*dt
               flux_v_iface(k) = nu_t*dv/h_avg*dt
            end if
         end do

         ! Pass 2: apply interface fluxes to layers.
         state%multilayer%hu_layer(1, cell) = state%multilayer%hu_layer(1, cell) + flux_u_iface(1)
         state%multilayer%hv_layer(1, cell) = state%multilayer%hv_layer(1, cell) + flux_v_iface(1)
         do k = 2, nz - 1
            state%multilayer%hu_layer(k, cell) = state%multilayer%hu_layer(k, cell) &
                                                 + (flux_u_iface(k) - flux_u_iface(k - 1))
            state%multilayer%hv_layer(k, cell) = state%multilayer%hv_layer(k, cell) &
                                                 + (flux_v_iface(k) - flux_v_iface(k - 1))
         end do
         state%multilayer%hu_layer(nz, cell) = state%multilayer%hu_layer(nz, cell) - flux_u_iface(nz - 1)
         state%multilayer%hv_layer(nz, cell) = state%multilayer%hv_layer(nz, cell) - flux_v_iface(nz - 1)
      end do
   end subroutine ml_kpp_viscosity_unstr

   pure subroutine ml_kpp_diffusion_tracers_unstr(state, hS_layer, hT_layer, &
                                                  kappa_S_bg, kappa_T_bg, dt)
      !! Pass 2 (tracers): PP81 cubic interior diffusivity with the KPP
      !! overlay max'd in, plus the non-local counter-gradient transport
      !! folded into the interface flux, applied two-pass on salinity +
      !! temperature.  Fuses the two `ml_vertical_diffusion_tracer_unstr`
      !! calls; bit-identical to them where the overlay + non-local
      !! terms vanish.
      !!
      !! The non-local term gamma(iface) = cs·G(sigma)·F_kin is a
      !! downward flux; adding dt·gamma to the diffusive interface
      !! transfer makes the shared apply pass produce its divergence.
      !! gamma = 0 at the bed and free surface (no iface there), so the
      !! column integral is conserved.
      type(unstr_state_t), intent(inout) :: state
      ! assumed-shape-ok: unstructured outer-shim — caller host-dereferences
      ! state%multilayer%tracers(idx)%hTr_layer; ncells varies per mesh partition.
      real(wp), intent(inout) :: hS_layer(:, :)
      real(wp), intent(inout) :: hT_layer(:, :)  ! assumed-shape-ok: unstructured outer-shim; ncells varies per mesh partition
      real(wp), intent(in) :: kappa_S_bg, kappa_T_bg, dt

      integer :: cell, k, nz, n
      real(wp) :: h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1
      real(wp) :: s_k, s_k1, t_k, t_k1, du, dv, drho, Ri, denom
      real(wp) :: kappa_S, kappa_T, h_col, h_b, b_0, us, wstar3, w_star, w_s
      real(wp) :: h_top, s_top, f_t_kin, f_s_kin, cs_nl
      real(wp) :: d_below, d_iface, sigma, g_shape, nu_kpp
      logical :: destabilizing
      real(wp) :: flux_s_iface(BPG_NZ_MAX), flux_t_iface(BPG_NZ_MAX)

      if (.not. state%use_multilayer) return
      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      if (nz < 2) return

      cs_nl = state%multilayer%kpp_cs_nonlocal

      do concurrent(cell=1:n) &
         local(k, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, &
               s_k, s_k1, t_k, t_k1, du, dv, drho, Ri, denom, &
               kappa_S, kappa_T, h_col, h_b, b_0, us, wstar3, w_star, w_s, &
               h_top, s_top, f_t_kin, f_s_kin, d_below, d_iface, sigma, &
               g_shape, nu_kpp, destabilizing, flux_s_iface, flux_t_iface)
         if (state%barotropic%height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + state%multilayer%h_layer(k, cell)
         end do
         h_b = state%multilayer%kpp_bl_depth(cell)
         b_0 = state%multilayer%kpp_buoy_flux(cell)
         us = state%multilayer%kpp_u_star(cell)
         wstar3 = max(0.0_wp, -b_0)*h_b
         w_star = wstar3**(1.0_wp/3.0_wp)
         w_s = sqrt(us*us + w_star*w_star)
         destabilizing = (b_0 < 0.0_wp)

         h_top = state%multilayer%h_layer(nz, cell)
         if (h_top > DRY_TOLERANCE) then
            s_top = hS_layer(nz, cell)/h_top
         else
            s_top = 0.0_wp
         end if
         f_t_kin = state%forcing%heat_flux(cell)/(state%multilayer%rho_0*CP_WATER)
         f_s_kin = state%forcing%evap_minus_precip(cell)*s_top

         ! Pass 1: interface fluxes (PP81 cubic + overlay + non-local)
         d_below = 0.0_wp
         do k = 1, nz - 1
            h_k = state%multilayer%h_layer(k, cell)
            h_k1 = state%multilayer%h_layer(k + 1, cell)
            h_avg = 0.5_wp*(h_k + h_k1)
            d_below = d_below + h_k
            d_iface = h_col - d_below

            flux_s_iface(k) = 0.0_wp
            flux_t_iface(k) = 0.0_wp

            if (h_avg > VANISHING_LAYER_TOL_UNSTR .and. &
                h_k > VANISHING_LAYER_TOL_UNSTR .and. &
                h_k1 > VANISHING_LAYER_TOL_UNSTR) then
               if (h_k > DRY_TOLERANCE) then
                  u_k = state%multilayer%hu_layer(k, cell)/h_k
                  v_k = state%multilayer%hv_layer(k, cell)/h_k
                  s_k = hS_layer(k, cell)/h_k
                  t_k = hT_layer(k, cell)/h_k
               else
                  u_k = 0.0_wp
                  v_k = 0.0_wp
                  s_k = 0.0_wp
                  t_k = 0.0_wp
               end if
               if (h_k1 > DRY_TOLERANCE) then
                  u_k1 = state%multilayer%hu_layer(k + 1, cell)/h_k1
                  v_k1 = state%multilayer%hv_layer(k + 1, cell)/h_k1
                  s_k1 = hS_layer(k + 1, cell)/h_k1
                  t_k1 = hT_layer(k + 1, cell)/h_k1
               else
                  u_k1 = 0.0_wp
                  v_k1 = 0.0_wp
                  s_k1 = 0.0_wp
                  t_k1 = 0.0_wp
               end if

               du = u_k1 - u_k
               dv = v_k1 - v_k
               drho = state%multilayer%rho_layer(k, cell) - state%multilayer%rho_layer(k + 1, cell)
               Ri = GRAVITY*(drho/state%multilayer%rho_0)*h_avg/(du*du + dv*dv + EPS)
               if (Ri < 0.0_wp) Ri = 0.0_wp

               denom = 1.0_wp + state%multilayer%vmix_ri_alpha*Ri
               kappa_S = state%multilayer%vmix_nu_max/(denom*denom*denom) + kappa_S_bg
               kappa_T = state%multilayer%vmix_nu_max/(denom*denom*denom) + kappa_T_bg

               g_shape = 0.0_wp
               if (w_s > 0.0_wp .and. h_b > 0.0_wp .and. d_iface < h_b) then
                  sigma = d_iface/h_b
                  g_shape = sigma*(1.0_wp - sigma)*(1.0_wp - sigma)
               end if
               nu_kpp = h_b*w_s*g_shape
               kappa_S = max(kappa_S, nu_kpp)
               kappa_T = max(kappa_T, nu_kpp)

               if (kappa_S*dt > 0.5_wp*h_avg*h_avg) kappa_S = 0.5_wp*h_avg*h_avg/dt
               if (kappa_T*dt > 0.5_wp*h_avg*h_avg) kappa_T = 0.5_wp*h_avg*h_avg/dt

               flux_s_iface(k) = kappa_S*(s_k1 - s_k)/h_avg*dt
               flux_t_iface(k) = kappa_T*(t_k1 - t_k)/h_avg*dt

               ! Non-local counter-gradient transport (g_shape = 0 outside
               ! the BL ⇒ no-op; skipped entirely under stabilizing flux).
               if (destabilizing) then
                  flux_s_iface(k) = flux_s_iface(k) + dt*cs_nl*g_shape*f_s_kin
                  flux_t_iface(k) = flux_t_iface(k) + dt*cs_nl*g_shape*f_t_kin
               end if
            end if
         end do

         ! Pass 2: apply interface fluxes to layers.
         hS_layer(1, cell) = hS_layer(1, cell) + flux_s_iface(1)
         hT_layer(1, cell) = hT_layer(1, cell) + flux_t_iface(1)
         do k = 2, nz - 1
            hS_layer(k, cell) = hS_layer(k, cell) + (flux_s_iface(k) - flux_s_iface(k - 1))
            hT_layer(k, cell) = hT_layer(k, cell) + (flux_t_iface(k) - flux_t_iface(k - 1))
         end do
         hS_layer(nz, cell) = hS_layer(nz, cell) - flux_s_iface(nz - 1)
         hT_layer(nz, cell) = hT_layer(nz, cell) - flux_t_iface(nz - 1)
      end do
   end subroutine ml_kpp_diffusion_tracers_unstr

end module rki_ml_kpp_unstr
