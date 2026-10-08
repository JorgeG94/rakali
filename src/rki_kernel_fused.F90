!! Fused update kernels for the 2D barotropic structured solver.
!!
!! Merge save + update + clamp/desingularization + friction + wind + pressure +
!! precipitation + coriolis into one kernel launch per RK stage, keeping
!! state in registers across the sequence. 2D barotropic path ONLY (the
!! NH/multilayer paths insert layer ops between update and friction).
!!   fused_update_stage1 — save U^n + Euler update + source terms
!!   fused_update_stage2 — Euler update + source terms + RK2 average
module rki_kernel_fused
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, DESING_EPS, &
                            FROUDE_BACKSTOP, RHO_WATER
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: fused_update_stage1
   public :: fused_update_stage2
   public :: fused_update_stage_generic

contains

   pure subroutine fused_update_stage_generic(h, hu, hv, &
                                              flux_h, flux_hu, flux_hv, &
                                              grid, dt, &
                                              manning_n, wind_x, wind_y, &
                                              atm_pressure, precipitation, &
                                              coriolis_f, forcing_active)
      !! One forward-Euler stage for the selectable-integrator path:
      !! Euler update + clamp/desingularization + friction + wind + pressure +
      !! precipitation + coriolis, all scaled by the per-stage `dt`
      !! (`beta(i)*dt_step`).  No state save and no RK2 average — the
      !! start-of-step backup is taken once by solver_save_state and the
      !! convex blend is a separate compute_ssp_combine call.  The body is
      !! identical to fused_update_stage1 minus the h0/hu0/hv0 writes so
      !! every scheme shares one source-term implementation.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_h(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt, manning_n
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: atm_pressure(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: precipitation(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: coriolis_f
      logical, intent(in) :: forcing_active

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: h_val, hu_val, hv_val
      real(wp) :: eps4, h4, speed_sq, fr_max_sq, fr_scale
      real(wp) :: u_val, v_val, speed, cf, denom
      real(wp) :: dpdx, dpdy, hu_old, hv_old
      real(wp) :: n_sq, cosf, sinf
      logical :: has_friction, has_coriolis

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      ! KP07 desingularization length, shared with the CFL timestep
      ! estimator (see DESING_EPS).
      eps4 = DESING_EPS**4

      n_sq = manning_n*manning_n
      has_friction = (manning_n > 0.0_wp)
      has_coriolis = (abs(coriolis_f) > 1.0e-20_wp)
      if (has_coriolis) then
         cosf = cos(coriolis_f*dt)
         sinf = sin(coriolis_f*dt)
      end if

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(h_val, hu_val, hv_val, h4, speed_sq, fr_max_sq, fr_scale, &
               u_val, v_val, speed, cf, denom, dpdx, dpdy, hu_old, hv_old)

         ! ---- Read state ----
         h_val = h(i, j)
         hu_val = hu(i, j)
         hv_val = hv(i, j)

         ! ---- Euler update ----
         h_val = h_val + dt*flux_h(i, j)
         hu_val = hu_val + dt*flux_hu(i, j)
         hv_val = hv_val + dt*flux_hv(i, j)

         ! ---- Clamp + thin-cell desingularization (KP07) ----
         if (h_val < 0.0_wp) then
            h_val = 0.0_wp
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else if (h_val < DRY_TOLERANCE) then
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else
            h4 = h_val**4
            if (h4 < eps4) then
               ! u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4)); hu := h*u.
               ! Identity for h >= eps, so applied only below eps.
               u_val = sqrt(2.0_wp)*h_val*hu_val/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h_val*hv_val/sqrt(h4 + eps4)
               hu_val = h_val*u_val
               hv_val = h_val*v_val
            end if
            ! Froude backstop: blowup arrester only (triggers on a TRUE
            ! comparison, so NaN passes through and fails loudly).
            speed_sq = (hu_val**2 + hv_val**2)/h_val**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h_val
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu_val = hu_val*fr_scale
               hv_val = hv_val*fr_scale
            end if
         end if

         ! ---- Friction + forcing ----
         if (h_val > DRY_TOLERANCE) then
            if (has_friction) then
               u_val = hu_val/h_val
               v_val = hv_val/h_val
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h_val**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h_val
               hu_val = hu_val/denom
               hv_val = hv_val/denom
            end if

            if (forcing_active) then
               hu_val = hu_val + dt*wind_x(i, j)/RHO_WATER
               hv_val = hv_val + dt*wind_y(i, j)/RHO_WATER
               dpdx = (atm_pressure(i + 1, j) - atm_pressure(i - 1, j))/(2.0_wp*dx)
               dpdy = (atm_pressure(i, j + 1) - atm_pressure(i, j - 1))/(2.0_wp*dy)
               hu_val = hu_val - dt*h_val*dpdx/RHO_WATER
               hv_val = hv_val - dt*h_val*dpdy/RHO_WATER
            end if

            ! ---- Coriolis ----
            if (has_coriolis) then
               hu_old = hu_val
               hv_old = hv_val
               hu_val = hu_old*cosf + hv_old*sinf
               hv_val = -hu_old*sinf + hv_old*cosf
            end if
         end if

         ! ---- Precipitation (all cells) ----
         if (forcing_active) then
            h_val = h_val + dt*precipitation(i, j)
         end if

         ! ---- Write back ----
         h(i, j) = h_val
         hu(i, j) = hu_val
         hv(i, j) = hv_val
      end do

   end subroutine fused_update_stage_generic

   pure subroutine fused_update_stage1(h, hu, hv, h0, hu0, hv0, &
                                       flux_h, flux_hu, flux_hv, &
                                       grid, dt, &
                                       manning_n, wind_x, wind_y, &
                                       atm_pressure, precipitation, &
                                       coriolis_f, forcing_active)
      !! Stage 1: save + Euler update + clamp + friction + coriolis
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(out)   :: h0(grid%nx_total, grid%ny_total)
      real(wp), intent(out)   :: hu0(grid%nx_total, grid%ny_total)
      real(wp), intent(out)   :: hv0(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_h(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt, manning_n
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: atm_pressure(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: precipitation(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: coriolis_f
      logical, intent(in) :: forcing_active

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: h_val, hu_val, hv_val
      real(wp) :: eps4, h4, speed_sq, fr_max_sq, fr_scale
      real(wp) :: u_val, v_val, speed, cf, denom
      real(wp) :: dpdx, dpdy, hu_old, hv_old
      real(wp) :: n_sq, cosf, sinf
      logical :: has_friction, has_coriolis

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      ! KP07 desingularization length, shared with the CFL timestep
      ! estimator (see DESING_EPS).
      eps4 = DESING_EPS**4

      n_sq = manning_n*manning_n
      has_friction = (manning_n > 0.0_wp)
      has_coriolis = (abs(coriolis_f) > 1.0e-20_wp)
      if (has_coriolis) then
         cosf = cos(coriolis_f*dt)
         sinf = sin(coriolis_f*dt)
      end if

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(h_val, hu_val, hv_val, h4, speed_sq, fr_max_sq, fr_scale, &
               u_val, v_val, speed, cf, denom, dpdx, dpdy, hu_old, hv_old)

         ! ---- Read + save backup ----
         h_val = h(i, j)
         hu_val = hu(i, j)
         hv_val = hv(i, j)
         h0(i, j) = h_val
         hu0(i, j) = hu_val
         hv0(i, j) = hv_val

         ! ---- Euler update ----
         h_val = h_val + dt*flux_h(i, j)
         hu_val = hu_val + dt*flux_hu(i, j)
         hv_val = hv_val + dt*flux_hv(i, j)

         ! ---- Clamp + thin-cell desingularization (KP07) ----
         if (h_val < 0.0_wp) then
            h_val = 0.0_wp
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else if (h_val < DRY_TOLERANCE) then
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else
            h4 = h_val**4
            if (h4 < eps4) then
               ! u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4)); hu := h*u.
               ! Identity for h >= eps, so applied only below eps.
               u_val = sqrt(2.0_wp)*h_val*hu_val/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h_val*hv_val/sqrt(h4 + eps4)
               hu_val = h_val*u_val
               hv_val = h_val*v_val
            end if
            ! Froude backstop: blowup arrester only (triggers on a TRUE
            ! comparison, so NaN passes through and fails loudly).
            speed_sq = (hu_val**2 + hv_val**2)/h_val**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h_val
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu_val = hu_val*fr_scale
               hv_val = hv_val*fr_scale
            end if
         end if

         ! ---- Friction + forcing ----
         if (h_val > DRY_TOLERANCE) then
            if (has_friction) then
               u_val = hu_val/h_val
               v_val = hv_val/h_val
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h_val**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h_val
               hu_val = hu_val/denom
               hv_val = hv_val/denom
            end if

            if (forcing_active) then
               hu_val = hu_val + dt*wind_x(i, j)/RHO_WATER
               hv_val = hv_val + dt*wind_y(i, j)/RHO_WATER
               dpdx = (atm_pressure(i + 1, j) - atm_pressure(i - 1, j))/(2.0_wp*dx)
               dpdy = (atm_pressure(i, j + 1) - atm_pressure(i, j - 1))/(2.0_wp*dy)
               hu_val = hu_val - dt*h_val*dpdx/RHO_WATER
               hv_val = hv_val - dt*h_val*dpdy/RHO_WATER
            end if

            ! ---- Coriolis ----
            if (has_coriolis) then
               hu_old = hu_val
               hv_old = hv_val
               hu_val = hu_old*cosf + hv_old*sinf
               hv_val = -hu_old*sinf + hv_old*cosf
            end if
         end if

         ! ---- Precipitation (all cells) ----
         if (forcing_active) then
            h_val = h_val + dt*precipitation(i, j)
         end if

         ! ---- Write back ----
         h(i, j) = h_val
         hu(i, j) = hu_val
         hv(i, j) = hv_val
      end do

   end subroutine fused_update_stage1

   pure subroutine fused_update_stage2(h, hu, hv, h0, hu0, hv0, &
                                       flux_h, flux_hu, flux_hv, &
                                       grid, dt, &
                                       manning_n, wind_x, wind_y, &
                                       atm_pressure, precipitation, &
                                       coriolis_f, forcing_active)
      !! Stage 2: Euler update + clamp + friction + coriolis + RK2 average + clamp
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: h0(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: hu0(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: hv0(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_h(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in)    :: flux_hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt, manning_n
      real(wp), intent(in) :: wind_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: wind_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: atm_pressure(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: precipitation(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: coriolis_f
      logical, intent(in) :: forcing_active

      integer :: i, j, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: h_val, hu_val, hv_val
      real(wp) :: eps4, h4, speed_sq, fr_max_sq, fr_scale
      real(wp) :: u_val, v_val, speed, cf, denom
      real(wp) :: dpdx, dpdy, hu_old, hv_old
      real(wp) :: n_sq, cosf, sinf
      logical :: has_friction, has_coriolis

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      ! KP07 desingularization length, shared with the CFL timestep
      ! estimator (see DESING_EPS).
      eps4 = DESING_EPS**4

      n_sq = manning_n*manning_n
      has_friction = (manning_n > 0.0_wp)
      has_coriolis = (abs(coriolis_f) > 1.0e-20_wp)
      if (has_coriolis) then
         cosf = cos(coriolis_f*dt)
         sinf = sin(coriolis_f*dt)
      end if

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(h_val, hu_val, hv_val, h4, speed_sq, fr_max_sq, fr_scale, &
               u_val, v_val, speed, cf, denom, dpdx, dpdy, hu_old, hv_old)

         ! ---- Read state ----
         h_val = h(i, j)
         hu_val = hu(i, j)
         hv_val = hv(i, j)

         ! ---- Euler update ----
         h_val = h_val + dt*flux_h(i, j)
         hu_val = hu_val + dt*flux_hu(i, j)
         hv_val = hv_val + dt*flux_hv(i, j)

         ! ---- Clamp + thin-cell desingularization (KP07) ----
         if (h_val < 0.0_wp) then
            h_val = 0.0_wp
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else if (h_val < DRY_TOLERANCE) then
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else
            h4 = h_val**4
            if (h4 < eps4) then
               ! u = sqrt(2) h hu / sqrt(h^4 + max(h^4, eps^4)); hu := h*u.
               ! Identity for h >= eps, so applied only below eps.
               u_val = sqrt(2.0_wp)*h_val*hu_val/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h_val*hv_val/sqrt(h4 + eps4)
               hu_val = h_val*u_val
               hv_val = h_val*v_val
            end if
            ! Froude backstop: blowup arrester only (triggers on a TRUE
            ! comparison, so NaN passes through and fails loudly).
            speed_sq = (hu_val**2 + hv_val**2)/h_val**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h_val
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu_val = hu_val*fr_scale
               hv_val = hv_val*fr_scale
            end if
         end if

         ! ---- Friction + forcing ----
         if (h_val > DRY_TOLERANCE) then
            if (has_friction) then
               u_val = hu_val/h_val
               v_val = hv_val/h_val
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h_val**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h_val
               hu_val = hu_val/denom
               hv_val = hv_val/denom
            end if

            if (forcing_active) then
               hu_val = hu_val + dt*wind_x(i, j)/RHO_WATER
               hv_val = hv_val + dt*wind_y(i, j)/RHO_WATER
               dpdx = (atm_pressure(i + 1, j) - atm_pressure(i - 1, j))/(2.0_wp*dx)
               dpdy = (atm_pressure(i, j + 1) - atm_pressure(i, j - 1))/(2.0_wp*dy)
               hu_val = hu_val - dt*h_val*dpdx/RHO_WATER
               hv_val = hv_val - dt*h_val*dpdy/RHO_WATER
            end if

            ! ---- Coriolis ----
            if (has_coriolis) then
               hu_old = hu_val
               hv_old = hv_val
               hu_val = hu_old*cosf + hv_old*sinf
               hv_val = -hu_old*sinf + hv_old*cosf
            end if
         end if

         ! ---- Precipitation (all cells) ----
         if (forcing_active) then
            h_val = h_val + dt*precipitation(i, j)
         end if

         ! ---- RK2 average: 0.5*(U^n + U**) ----
         h_val = 0.5_wp*(h0(i, j) + h_val)
         hu_val = 0.5_wp*(hu0(i, j) + hu_val)
         hv_val = 0.5_wp*(hv0(i, j) + hv_val)

         ! ---- Final clamp + thin-cell desingularization (post-average) ----
         if (h_val < 0.0_wp) then
            h_val = 0.0_wp
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else if (h_val < DRY_TOLERANCE) then
            hu_val = 0.0_wp
            hv_val = 0.0_wp
         else
            h4 = h_val**4
            if (h4 < eps4) then
               u_val = sqrt(2.0_wp)*h_val*hu_val/sqrt(h4 + eps4)
               v_val = sqrt(2.0_wp)*h_val*hv_val/sqrt(h4 + eps4)
               hu_val = h_val*u_val
               hv_val = h_val*v_val
            end if
            speed_sq = (hu_val**2 + hv_val**2)/h_val**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h_val
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               hu_val = hu_val*fr_scale
               hv_val = hv_val*fr_scale
            end if
         end if

         ! ---- Write back ----
         h(i, j) = h_val
         hu(i, j) = hu_val
         hv(i, j) = hv_val
      end do

   end subroutine fused_update_stage2

end module rki_kernel_fused
