!! Vertical viscosity for layer momentum on the unstructured grid.
module rki_ml_vertical_viscosity_unstr
   !! Currently implements the gradient-Richardson-number closure of
   !! Pacanowski & Philander (1981).  Computes a turbulent viscosity at
   !! each interior interface and applies a diffusive flux to layer
   !! momentum:
   !!
   !!   nu_t = nu_max / (1 + alpha*Ri)^2 + nu_bg
   !!
   !! Column-local, no neighbour reads — safe for `do concurrent`.
   !! When an ocean-grade KPP closure lands, it will live alongside this
   !! one and the solver will pick between them via a config flag.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_state_unstr, only: unstr_state_t
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR, BPG_NZ_MAX, &
                                  VANISHING_LAYER_TOL_UNSTR
   implicit none
   private

   public :: ml_vertical_viscosity_unstr

contains

   subroutine ml_vertical_viscosity_unstr(state, dt)
      !! Ri-dependent vertical viscosity for momentum.
      !! Computes gradient Richardson number at each interface and applies
      !! turbulent diffusive flux: nu_t = nu_max / (1 + alpha*Ri)^2 + nu_bg.
      !!
      !! Column-local — no neighbour access, safe for `do concurrent`.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt

      integer :: cell, k, nz, n
      real(wp) :: h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1
      real(wp) :: du, dv, drho, Ri, denom, nu_t, flux_u, flux_v
      real(wp) :: flux_u_iface(BPG_NZ_MAX), flux_v_iface(BPG_NZ_MAX)

      real(wp), parameter :: EPS = 1.0e-20_wp

      if (.not. state%use_multilayer) return
      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      if (nz < 2) return

      !$omp parallel do private(k, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, du, dv, drho, Ri, denom, nu_t, flux_u, flux_v, flux_u_iface, flux_v_iface)
      do cell = 1, n

         if (state%barotropic%height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

         ! Pass 1: compute interface fluxes
         do k = 1, nz - 1
            h_k = state%multilayer%h_layer(k, cell)
            h_k1 = state%multilayer%h_layer(k + 1, cell)
            h_avg = 0.5_wp*(h_k + h_k1)

            flux_u_iface(k) = 0.0_wp
            flux_v_iface(k) = 0.0_wp

            ! VCOORD_ZSTAR_FULL: no-flux through interfaces adjacent to
            ! vanishing layers (otherwise nu_t * du / h_avg blows up
            ! when h_avg = h_min / 2).
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

               ! Stability: nu_t * dt / h_avg^2 <= 0.5
               if (nu_t*dt > 0.5_wp*h_avg*h_avg) then
                  nu_t = 0.5_wp*h_avg*h_avg/dt
               end if

               flux_u_iface(k) = nu_t*du/h_avg*dt
               flux_v_iface(k) = nu_t*dv/h_avg*dt
            end if
         end do

         ! Pass 2: apply fluxes to layers.  flux_*_iface(k) is the momentum
         ! transfer across interface k+1/2 (units m^2/s); same quantity adds
         ! to layer k and subtracts from layer k+1.  Previous h_layer
         ! multipliers were dimensionally wrong and non-conservative under
         ! non-uniform dz (e.g. VCOORD_ZSTAR_FULL).
         state%multilayer%hu_layer(1, cell) = state%multilayer%hu_layer(1, cell) &
                                              + flux_u_iface(1)
         state%multilayer%hv_layer(1, cell) = state%multilayer%hv_layer(1, cell) &
                                              + flux_v_iface(1)
         do k = 2, nz - 1
            state%multilayer%hu_layer(k, cell) = state%multilayer%hu_layer(k, cell) &
                                                 + (flux_u_iface(k) - flux_u_iface(k - 1))
            state%multilayer%hv_layer(k, cell) = state%multilayer%hv_layer(k, cell) &
                                                 + (flux_v_iface(k) - flux_v_iface(k - 1))
         end do
         state%multilayer%hu_layer(nz, cell) = state%multilayer%hu_layer(nz, cell) &
                                               - flux_u_iface(nz - 1)
         state%multilayer%hv_layer(nz, cell) = state%multilayer%hv_layer(nz, cell) &
                                               - flux_v_iface(nz - 1)
      end do
      !$omp end parallel do

   end subroutine ml_vertical_viscosity_unstr

end module rki_ml_vertical_viscosity_unstr
