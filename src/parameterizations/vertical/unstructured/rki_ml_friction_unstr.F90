!! Per-layer bottom friction and Coriolis rotation on the unstructured
!! multilayer column.
module rki_ml_friction_unstr
   !! Provides `ml_apply_coriolis_unstr` (closed-form per-layer Coriolis
   !! rotation, mirrors structured `ml_apply_coriolis`) and
   !! `ml_apply_bottom_friction_unstr` (semi-implicit Manning quadratic
   !! drag on the lowest active layer, mirrors structured
   !! `ml_apply_bottom_friction`).
   !!
   !! Layer-index convention: bottom-up (`k=1` bed, `k=nz` surface).
   !! Coriolis is uniform in the vertical (scalar `f`), applied to every
   !! layer's `(hu_k, hv_k)`. Bottom friction acts on the lowest
   !! non-vanishing layer (usually `k=1`; the deepest active layer for
   !! ZSTAR_FULL where the bed layer can be parked at h_min). Both are
   !! column-local — safe for `do concurrent`.
   use rki_constants, only: wp, GRAVITY, THIN_LAYER_THRESHOLD
   use rki_state_unstr, only: unstr_state_t
   use rki_ml_common_unstr, only: VANISHING_LAYER_TOL_UNSTR
   implicit none
   private

   public :: ml_apply_bottom_friction_unstr
   public :: ml_apply_coriolis_unstr

contains

   subroutine ml_apply_bottom_friction_unstr(state, dt)
      !! Apply semi-implicit Manning bottom friction to the lowest active
      !! layer (bottom-up `k=1` bed) using that layer's velocity:
      !!   hu_new = hu / (1 + dt*C_d*|u_b|/h_k),  C_d = g*n^2/H^{1/3},
      !! where H is the total column depth. No-op when manning_n <= 0.
      !! Mirrors structured `ml_apply_bottom_friction` (Manning path); the
      !! unstructured barotropic state carries no linear-drag rate.
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt

      integer :: cell, k, kk, k_bot, nz, n
      real(wp) :: manning_n, n_sq, h_k, h_total, u_b, v_b, speed_b, cf, denom

      if (.not. state%use_multilayer) return
      manning_n = state%barotropic%manning_n
      if (manning_n <= 0.0_wp) return
      if (dt <= 0.0_wp) return

      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      n_sq = manning_n*manning_n

      !$omp target teams distribute parallel do private(k, kk, k_bot, h_k, h_total, u_b, v_b, speed_b, cf, denom)
      do cell = 1, n

         ! Lowest active layer: k=1 for SIGMA/ZSIGMA/ZSTAR-lite; the first
         ! non-vanishing layer for ZSTAR_FULL (bed layer may be parked at
         ! h_min in shallow columns).
         k_bot = 1
         do kk = 1, nz
            if (state%multilayer%h_layer(kk, cell) > VANISHING_LAYER_TOL_UNSTR) then
               k_bot = kk
               exit
            end if
         end do

         h_k = state%multilayer%h_layer(k_bot, cell)
         if (h_k >= THIN_LAYER_THRESHOLD) then
            h_total = 0.0_wp
            do kk = 1, nz
               h_total = h_total + state%multilayer%h_layer(kk, cell)
            end do
            if (h_total >= THIN_LAYER_THRESHOLD) then
               u_b = state%multilayer%hu_layer(k_bot, cell)/h_k
               v_b = state%multilayer%hv_layer(k_bot, cell)/h_k
               speed_b = sqrt(u_b*u_b + v_b*v_b)
               cf = GRAVITY*n_sq/(h_total**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed_b/h_k
               state%multilayer%hu_layer(k_bot, cell) = &
                  state%multilayer%hu_layer(k_bot, cell)/denom
               state%multilayer%hv_layer(k_bot, cell) = &
                  state%multilayer%hv_layer(k_bot, cell)/denom
            end if
         end if
      end do
      !$omp end target teams distribute parallel do

   end subroutine ml_apply_bottom_friction_unstr

   subroutine ml_apply_coriolis_unstr(state, dt)
      !! Rotate every wet layer's `(hu_k, hv_k)` through angle `f*dt` with
      !! the closed-form rotation matrix (energy-conserving, det = 1).
      !! Skips layers thinner than THIN_LAYER_THRESHOLD (already zeroed by
      !! the distribute / friction kernels).
      type(unstr_state_t), intent(inout) :: state
      real(wp), intent(in) :: dt

      integer :: cell, k, nz, n
      real(wp) :: f, cosf, sinf, hu_old, hv_old

      if (.not. state%use_multilayer) return
      f = state%barotropic%coriolis_f
      if (f == 0.0_wp) return

      nz = state%multilayer%nz_ml
      n = state%ncells_owned

      cosf = cos(f*dt)
      sinf = sin(f*dt)

      !$omp target teams distribute parallel do private(k, hu_old, hv_old)
      do cell = 1, n
         do k = 1, nz
            if (state%multilayer%h_layer(k, cell) > THIN_LAYER_THRESHOLD) then
               hu_old = state%multilayer%hu_layer(k, cell)
               hv_old = state%multilayer%hv_layer(k, cell)
               state%multilayer%hu_layer(k, cell) = hu_old*cosf + hv_old*sinf
               state%multilayer%hv_layer(k, cell) = -hu_old*sinf + hv_old*cosf
            end if
         end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine ml_apply_coriolis_unstr

end module rki_ml_friction_unstr
