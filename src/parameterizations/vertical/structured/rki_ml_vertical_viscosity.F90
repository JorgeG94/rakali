!! Vertical viscosity for layer momentum.
module rki_ml_vertical_viscosity
   !! Currently implements the gradient-Richardson-number closure of
   !! Pacanowski & Philander (1981).  Computes a turbulent viscosity at
   !! each interior interface and applies a diffusive flux to layer
   !! momentum:
   !!
   !!   nu_t = nu_max / (1 + alpha*Ri)^2 + nu_bg
   !!
   !! When an ocean-grade KPP closure lands, it will live alongside
   !! this one (e.g. `vertical/rki_ml_kpp.F90`) and the solver will
   !! pick between them via a config flag.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH, VANISHING_LAYER_TOL
   implicit none
   private

   public :: ml_vertical_viscosity

contains

   subroutine ml_vertical_viscosity(h_layer, hu_layer, hv_layer, rho_layer, &
                                         rho_0, nu_max, nu_bg, ri_alpha, &
                                         grid, nz, dt)
      !! Ri-dependent vertical viscosity for momentum (Pacanowski-Philander 1981).
      !! Computes gradient Richardson number at each interface and applies
      !! turbulent diffusive flux: nu_t = nu_max / (1 + alpha*Ri)^2 + nu_bg.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0
         !! Reference density (kg/m^3)
      real(wp), intent(in) :: nu_max
         !! Maximum turbulent viscosity (m^2/s)
      real(wp), intent(in) :: nu_bg
         !! Background viscosity (m^2/s)
      real(wp), intent(in) :: ri_alpha
         !! Richardson number scaling coefficient
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: h_avg, h_k, h_k1, u_k, v_k, u_k1, v_k1, h_col
      real(wp) :: du, dv, drho, Ri, denom, nu_t, flux_u, flux_v

      real(wp), parameter :: EPS = 1.0e-20_wp

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (nz < 2) return

      !$omp target teams distribute parallel do collapse(2) private(k, kk, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, du, dv, drho, Ri, denom, nu_t, flux_u, flux_v, h_col)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then
            do k = 1, nz - 1
               h_k = h_layer(i, j, k)
               h_k1 = h_layer(i, j, k + 1)
               h_avg = 0.5_wp*(h_k + h_k1)

               ! VCOORD_ZSTAR_FULL: vanishing layers (h ≈ h_min) on either
               ! side of the interface produce a tiny h_avg, which makes
               ! the diffusive flux nu_t*du/h_avg blow up.  Treat as
               ! no-flux instead.  VANISHING_LAYER_TOL ≥ DRY_TOLERANCE,
               ! so this is strictly more conservative than the legacy guard.
               if (h_avg > VANISHING_LAYER_TOL .and. &
                   h_k > VANISHING_LAYER_TOL .and. &
                   h_k1 > VANISHING_LAYER_TOL) then
                  if (h_k > DRY_TOLERANCE) then
                     u_k = hu_layer(i, j, k)/h_k
                     v_k = hv_layer(i, j, k)/h_k
                  else
                     u_k = 0.0_wp
                     v_k = 0.0_wp
                  end if
                  if (h_k1 > DRY_TOLERANCE) then
                     u_k1 = hu_layer(i, j, k + 1)/h_k1
                     v_k1 = hv_layer(i, j, k + 1)/h_k1
                  else
                     u_k1 = 0.0_wp
                     v_k1 = 0.0_wp
                  end if

                  du = u_k1 - u_k
                  dv = v_k1 - v_k
                  drho = rho_layer(i, j, k) - rho_layer(i, j, k + 1)

                  Ri = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + EPS)
                  if (Ri < 0.0_wp) Ri = 0.0_wp

                  denom = 1.0_wp + ri_alpha*Ri
                  nu_t = nu_max/(denom*denom) + nu_bg

                  ! Limit diffusive flux: stability requires nu_t*dt/h_avg^2 < 0.5
                  if (nu_t*dt > 0.5_wp*h_avg*h_avg) then
                     nu_t = 0.5_wp*h_avg*h_avg/dt
                  end if

                  flux_u = nu_t*du/h_avg*dt
                  flux_v = nu_t*dv/h_avg*dt

                  ! Conservative + dimensionally correct: same momentum flux
                  ! both sides of the interface (previous h_k / h_k1 multipliers
                  ! were dimensionally wrong and non-conservative under
                  ! non-uniform layer thicknesses, e.g. VCOORD_ZSTAR_FULL).
                  hu_layer(i, j, k) = hu_layer(i, j, k) + flux_u
                  hu_layer(i, j, k + 1) = hu_layer(i, j, k + 1) - flux_u
                  hv_layer(i, j, k) = hv_layer(i, j, k) + flux_v
                  hv_layer(i, j, k + 1) = hv_layer(i, j, k + 1) - flux_v
               end if
            end do
         end if
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine ml_vertical_viscosity

end module rki_ml_vertical_viscosity
