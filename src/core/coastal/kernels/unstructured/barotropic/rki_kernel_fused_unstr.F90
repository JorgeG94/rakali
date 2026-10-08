!! Fused update kernels for the unstructured solver: merge save +
!! update + friction + wind + protect into one kernel launch per RK
!! stage, keeping state in registers across the sequence. Entry points:
!!   fused_update_stage1_unstr — save U^n + Euler + friction + protect
!!   fused_update_stage2_unstr — Euler + friction + protect + RK2 avg + protect
module rki_kernel_fused_unstr
   use rki_constants, only: wp, GRAVITY, RHO_WATER, DRY_TOLERANCE, &
                            DESING_EPS, FROUDE_BACKSTOP
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   implicit none
   private

   public :: fused_update_stage1_unstr
   public :: fused_update_stage2_unstr
   public :: fused_save_state_unstr
   public :: fused_update_stage_generic_unstr
   public :: ssp_combine_unstr

contains

   subroutine fused_save_state_unstr(state, mesh)
      !! Snapshot U^n into the `_backup` slots once at step start, for the
      !! selectable-integrator generic path (the convex blend reads them).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k, n

      n = mesh%ncells_owned

      !$omp target teams distribute parallel do
      do k = 1, n
         state%barotropic%stage_backup(k) = state%barotropic%stage_cv(k)
         state%barotropic%xmom_backup(k) = state%barotropic%xmom_cv(k)
         state%barotropic%ymom_backup(k) = state%barotropic%ymom_cv(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine fused_save_state_unstr

   subroutine fused_update_stage_generic_unstr(state, mesh, dt)
      !! One forward-Euler stage for the selectable-integrator path:
      !! Euler update + friction/wind + protect, scaled by the per-stage
      !! `dt` (`beta(i)*dt_step`).  Identical to fused_update_stage1_unstr
      !! minus the start-of-step backup writes (taken once by
      !! fused_save_state_unstr) so every scheme shares one source-term body.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: k, n
      real(wp) :: stage, xmom, ymom, h, eps4, h4
      real(wp) :: speed_sq, fr_max_sq, fr_scale
      real(wp) :: u_val, v_val, speed, cf, denom
      real(wp) :: n_sq, min_h, tau_x, tau_y
      real(wp) :: mass_error_local
      logical :: has_friction, has_wind

      n = mesh%ncells_owned
      has_friction = (state%barotropic%manning_n > 0.0_wp)
      has_wind = (state%barotropic%wind_tau_x /= 0.0_wp) .or. &
                 (state%barotropic%wind_tau_y /= 0.0_wp) .or. &
                 state%forcing%wind_active
      n_sq = state%barotropic%manning_n*state%barotropic%manning_n
      min_h = state%minimum_allowed_height
      eps4 = DESING_EPS**4
      tau_x = state%barotropic%wind_tau_x/RHO_WATER
      tau_y = state%barotropic%wind_tau_y/RHO_WATER
      mass_error_local = 0.0_wp

      !$omp target teams distribute parallel do private(stage, xmom, ymom, h, h4, speed_sq, fr_max_sq, fr_scale, u_val, v_val, speed, cf, denom)
      do k = 1, n
         ! ---- Read state ----
         stage = state%barotropic%stage_cv(k)
         xmom = state%barotropic%xmom_cv(k)
         ymom = state%barotropic%ymom_cv(k)

         ! ---- Euler update ----
         stage = stage + dt*state%barotropic%stage_eu(k)
         xmom = xmom + dt*state%barotropic%xmom_eu(k)
         ymom = ymom + dt*state%barotropic%ymom_eu(k)

         ! ---- Friction + wind ----
         h = stage - state%barotropic%bed_cv(k)
         if (h > min_h) then
            if (has_friction) then
               u_val = xmom/h
               v_val = ymom/h
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h
               xmom = xmom/denom
               ymom = ymom/denom
            end if
            if (has_wind) then
               xmom = xmom + dt*(tau_x + state%forcing%wind_x(k)/RHO_WATER)
               ymom = ymom + dt*(tau_y + state%forcing%wind_y(k)/RHO_WATER)
            end if
         end if

         ! ---- Protect ----
         if (h < DRY_TOLERANCE) then
            ! Dry (or negative) -- zero momentum.
            xmom = 0.0_wp
            ymom = 0.0_wp
         else
            ! Thin-but-wet: KP07 desingularized velocity
            ! u = sqrt(2) h uh / sqrt(h^4 + max(h^4, eps^4)), uh := h*u
            ! (identity for h >= eps, so applied only below eps), then
            ! the Froude backstop — blowup arrester only (triggers on a
            ! TRUE comparison, so NaN passes through and fails loudly).
            h4 = h**4
            if (h4 < eps4) then
               xmom = h*(sqrt(2.0_wp)*h*xmom/sqrt(h4 + eps4))
               ymom = h*(sqrt(2.0_wp)*h*ymom/sqrt(h4 + eps4))
            end if
            speed_sq = (xmom**2 + ymom**2)/(h*h)
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               xmom = xmom*fr_scale
               ymom = ymom*fr_scale
            end if
         end if
         if (h < 0.0_wp) then
            mass_error_local = mass_error_local + (-h)*mesh%area(k)
            stage = state%barotropic%bed_cv(k)
         end if

         ! ---- Write back ----
         state%barotropic%stage_cv(k) = stage
         state%barotropic%xmom_cv(k) = xmom
         state%barotropic%ymom_cv(k) = ymom
      end do
      !$omp end target teams distribute parallel do

      state%mass_error = state%mass_error + mass_error_local

   end subroutine fused_update_stage_generic_unstr

   subroutine ssp_combine_unstr(state, mesh, a_n)
      !! Generalized SSP Shu-Osher combine: U = a_n*U^n + (1 - a_n)*U,
      !! where U is the post-Euler-stage iterate and `_backup` = U^n.
      !! Reduces exactly to the 0.5 average in fused_update_stage2_unstr for
      !! a_n = 0.5.  Applies the same post-blend protect (mass-error bucket
      !! included) so each stage stays wet/dry-safe.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: a_n

      integer :: k, n
      real(wp) :: stage, xmom, ymom, h, eps4, h4, a_prev
      real(wp) :: speed_sq, fr_max_sq, fr_scale
      real(wp) :: mass_error_local

      n = mesh%ncells_owned
      a_prev = 1.0_wp - a_n
      mass_error_local = 0.0_wp
      eps4 = DESING_EPS**4

      !$omp target teams distribute parallel do private(stage, xmom, ymom, h, h4, speed_sq, fr_max_sq, fr_scale)
      do k = 1, n
         stage = a_n*state%barotropic%stage_backup(k) + a_prev*state%barotropic%stage_cv(k)
         xmom = a_n*state%barotropic%xmom_backup(k) + a_prev*state%barotropic%xmom_cv(k)
         ymom = a_n*state%barotropic%ymom_backup(k) + a_prev*state%barotropic%ymom_cv(k)

         h = stage - state%barotropic%bed_cv(k)
         if (h < DRY_TOLERANCE) then
            ! Dry (or negative) -- zero momentum.
            xmom = 0.0_wp
            ymom = 0.0_wp
         else
            ! Thin-but-wet: KP07 desingularized velocity
            ! u = sqrt(2) h uh / sqrt(h^4 + max(h^4, eps^4)), uh := h*u
            ! (identity for h >= eps, so applied only below eps), then
            ! the Froude backstop — blowup arrester only (triggers on a
            ! TRUE comparison, so NaN passes through and fails loudly).
            h4 = h**4
            if (h4 < eps4) then
               xmom = h*(sqrt(2.0_wp)*h*xmom/sqrt(h4 + eps4))
               ymom = h*(sqrt(2.0_wp)*h*ymom/sqrt(h4 + eps4))
            end if
            speed_sq = (xmom**2 + ymom**2)/(h*h)
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               xmom = xmom*fr_scale
               ymom = ymom*fr_scale
            end if
         end if
         if (h < 0.0_wp) then
            mass_error_local = mass_error_local + (-h)*mesh%area(k)
            stage = state%barotropic%bed_cv(k)
         end if

         state%barotropic%stage_cv(k) = stage
         state%barotropic%xmom_cv(k) = xmom
         state%barotropic%ymom_cv(k) = ymom
      end do
      !$omp end target teams distribute parallel do

      state%mass_error = state%mass_error + mass_error_local

   end subroutine ssp_combine_unstr

   subroutine fused_update_stage1_unstr(state, mesh, dt)
      !! Stage 1: save + Euler update + friction/wind + protect
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: k, n
      real(wp) :: stage, xmom, ymom, h, eps4, h4
      real(wp) :: speed_sq, fr_max_sq, fr_scale
      real(wp) :: u_val, v_val, speed, cf, denom
      real(wp) :: n_sq, min_h, tau_x, tau_y
      real(wp) :: mass_error_local
      logical :: has_friction, has_wind

      n = mesh%ncells_owned
      has_friction = (state%barotropic%manning_n > 0.0_wp)
      has_wind = (state%barotropic%wind_tau_x /= 0.0_wp) .or. &
                 (state%barotropic%wind_tau_y /= 0.0_wp) .or. &
                 state%forcing%wind_active
      n_sq = state%barotropic%manning_n*state%barotropic%manning_n
      min_h = state%minimum_allowed_height
      eps4 = DESING_EPS**4
      tau_x = state%barotropic%wind_tau_x/RHO_WATER
      tau_y = state%barotropic%wind_tau_y/RHO_WATER
      mass_error_local = 0.0_wp

      !$omp target teams distribute parallel do private(stage, xmom, ymom, h, h4, speed_sq, fr_max_sq, fr_scale, u_val, v_val, speed, cf, denom)
      do k = 1, n
         ! ---- Read + save backup ----
         stage = state%barotropic%stage_cv(k)
         xmom = state%barotropic%xmom_cv(k)
         ymom = state%barotropic%ymom_cv(k)
         state%barotropic%stage_backup(k) = stage
         state%barotropic%xmom_backup(k) = xmom
         state%barotropic%ymom_backup(k) = ymom

         ! ---- Euler update ----
         stage = stage + dt*state%barotropic%stage_eu(k)
         xmom = xmom + dt*state%barotropic%xmom_eu(k)
         ymom = ymom + dt*state%barotropic%ymom_eu(k)

         ! ---- Friction + wind ----
         h = stage - state%barotropic%bed_cv(k)
         if (h > min_h) then
            if (has_friction) then
               u_val = xmom/h
               v_val = ymom/h
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h
               xmom = xmom/denom
               ymom = ymom/denom
            end if
            if (has_wind) then
               xmom = xmom + dt*(tau_x + state%forcing%wind_x(k)/RHO_WATER)
               ymom = ymom + dt*(tau_y + state%forcing%wind_y(k)/RHO_WATER)
            end if
         end if

         ! ---- Protect ----
         if (h < DRY_TOLERANCE) then
            ! Dry (or negative) -- zero momentum.
            xmom = 0.0_wp
            ymom = 0.0_wp
         else
            ! Thin-but-wet: KP07 desingularized velocity
            ! u = sqrt(2) h uh / sqrt(h^4 + max(h^4, eps^4)), uh := h*u
            ! (identity for h >= eps, so applied only below eps), then
            ! the Froude backstop — blowup arrester only (triggers on a
            ! TRUE comparison, so NaN passes through and fails loudly).
            h4 = h**4
            if (h4 < eps4) then
               xmom = h*(sqrt(2.0_wp)*h*xmom/sqrt(h4 + eps4))
               ymom = h*(sqrt(2.0_wp)*h*ymom/sqrt(h4 + eps4))
            end if
            speed_sq = (xmom**2 + ymom**2)/(h*h)
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               xmom = xmom*fr_scale
               ymom = ymom*fr_scale
            end if
         end if
         if (h < 0.0_wp) then
            mass_error_local = mass_error_local + (-h)*mesh%area(k)
            stage = state%barotropic%bed_cv(k)
         end if

         ! ---- Write back ----
         state%barotropic%stage_cv(k) = stage
         state%barotropic%xmom_cv(k) = xmom
         state%barotropic%ymom_cv(k) = ymom
      end do
      !$omp end target teams distribute parallel do

      state%mass_error = state%mass_error + mass_error_local

   end subroutine fused_update_stage1_unstr

   subroutine fused_update_stage2_unstr(state, mesh, dt)
      !! Stage 2: Euler update + friction/wind + protect + RK2 average + protect
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: k, n
      real(wp) :: stage, xmom, ymom, h, eps4, h4
      real(wp) :: speed_sq, fr_max_sq, fr_scale
      real(wp) :: u_val, v_val, speed, cf, denom
      real(wp) :: n_sq, min_h, tau_x, tau_y
      real(wp) :: mass_error_local
      logical :: has_friction, has_wind

      n = mesh%ncells_owned
      has_friction = (state%barotropic%manning_n > 0.0_wp)
      has_wind = (state%barotropic%wind_tau_x /= 0.0_wp) .or. &
                 (state%barotropic%wind_tau_y /= 0.0_wp) .or. &
                 state%forcing%wind_active
      n_sq = state%barotropic%manning_n*state%barotropic%manning_n
      min_h = state%minimum_allowed_height
      eps4 = DESING_EPS**4
      tau_x = state%barotropic%wind_tau_x/RHO_WATER
      tau_y = state%barotropic%wind_tau_y/RHO_WATER
      mass_error_local = 0.0_wp

      !$omp target teams distribute parallel do private(stage, xmom, ymom, h, h4, speed_sq, fr_max_sq, fr_scale, u_val, v_val, speed, cf, denom)
      do k = 1, n
         ! ---- Read state ----
         stage = state%barotropic%stage_cv(k)
         xmom = state%barotropic%xmom_cv(k)
         ymom = state%barotropic%ymom_cv(k)

         ! ---- Euler update ----
         stage = stage + dt*state%barotropic%stage_eu(k)
         xmom = xmom + dt*state%barotropic%xmom_eu(k)
         ymom = ymom + dt*state%barotropic%ymom_eu(k)

         ! ---- Friction + wind ----
         h = stage - state%barotropic%bed_cv(k)
         if (h > min_h) then
            if (has_friction) then
               u_val = xmom/h
               v_val = ymom/h
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h
               xmom = xmom/denom
               ymom = ymom/denom
            end if
            if (has_wind) then
               xmom = xmom + dt*(tau_x + state%forcing%wind_x(k)/RHO_WATER)
               ymom = ymom + dt*(tau_y + state%forcing%wind_y(k)/RHO_WATER)
            end if
         end if

         ! ---- Protect (substage 2) ----
         if (h < DRY_TOLERANCE) then
            ! Dry (or negative) -- zero momentum.
            xmom = 0.0_wp
            ymom = 0.0_wp
         else
            ! Thin-but-wet: KP07 desingularized velocity
            ! u = sqrt(2) h uh / sqrt(h^4 + max(h^4, eps^4)), uh := h*u
            ! (identity for h >= eps, so applied only below eps), then
            ! the Froude backstop — blowup arrester only (triggers on a
            ! TRUE comparison, so NaN passes through and fails loudly).
            h4 = h**4
            if (h4 < eps4) then
               xmom = h*(sqrt(2.0_wp)*h*xmom/sqrt(h4 + eps4))
               ymom = h*(sqrt(2.0_wp)*h*ymom/sqrt(h4 + eps4))
            end if
            speed_sq = (xmom**2 + ymom**2)/(h*h)
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               xmom = xmom*fr_scale
               ymom = ymom*fr_scale
            end if
         end if
         if (h < 0.0_wp) then
            stage = state%barotropic%bed_cv(k)
         end if

         ! ---- RK2 average: 0.5*(U^n + U**) ----
         stage = 0.5_wp*stage + 0.5_wp*state%barotropic%stage_backup(k)
         xmom = 0.5_wp*xmom + 0.5_wp*state%barotropic%xmom_backup(k)
         ymom = 0.5_wp*ymom + 0.5_wp*state%barotropic%ymom_backup(k)

         ! ---- Final protect (post-average) ----
         h = stage - state%barotropic%bed_cv(k)
         if (h < DRY_TOLERANCE) then
            ! Dry (or negative) -- zero momentum.
            xmom = 0.0_wp
            ymom = 0.0_wp
         else
            ! Thin-but-wet: KP07 desingularized velocity
            ! u = sqrt(2) h uh / sqrt(h^4 + max(h^4, eps^4)), uh := h*u
            ! (identity for h >= eps, so applied only below eps), then
            ! the Froude backstop — blowup arrester only (triggers on a
            ! TRUE comparison, so NaN passes through and fails loudly).
            h4 = h**4
            if (h4 < eps4) then
               xmom = h*(sqrt(2.0_wp)*h*xmom/sqrt(h4 + eps4))
               ymom = h*(sqrt(2.0_wp)*h*ymom/sqrt(h4 + eps4))
            end if
            speed_sq = (xmom**2 + ymom**2)/(h*h)
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               xmom = xmom*fr_scale
               ymom = ymom*fr_scale
            end if
         end if
         if (h < 0.0_wp) then
            mass_error_local = mass_error_local + (-h)*mesh%area(k)
            stage = state%barotropic%bed_cv(k)
         end if

         ! ---- Write back ----
         state%barotropic%stage_cv(k) = stage
         state%barotropic%xmom_cv(k) = xmom
         state%barotropic%ymom_cv(k) = ymom
      end do
      !$omp end target teams distribute parallel do

      state%mass_error = state%mass_error + mass_error_local

   end subroutine fused_update_stage2_unstr

end module rki_kernel_fused_unstr
