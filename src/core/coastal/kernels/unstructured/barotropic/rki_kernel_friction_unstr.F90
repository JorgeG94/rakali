!! Bottom friction + wind stress kernel for the unstructured solver
module rki_kernel_friction_unstr
   !! Apply Manning bottom friction and surface wind stress to the
   !! centroid momenta (bundled to mirror the structured solver).
   !!
   !! Friction (semi-implicit, unconditionally stable):
   !!   c_f = g*n^2/h^{1/3}; denom = 1 + dt*c_f*|u|/h; M /= denom.
   !! Wind stress (explicit): M += dt*tau/RHO_WATER, with the uniform
   !! `state%wind_tau_{x,y}` in Pa.
   !!
   !! No-ops by default: friction skips when `manning_n == 0`, wind when
   !! both components are zero; both inactive returns immediately. Dry
   !! cells (depth below `state%minimum_allowed_height`) are skipped to
   !! avoid the `1/h^{4/3}` singularity.
   use rki_constants, only: wp, GRAVITY, RHO_WATER
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   implicit none
   private

   public :: apply_manning_friction_unstr

contains

   subroutine apply_manning_friction_unstr(state, mesh, dt)
      !! Semi-implicit Manning friction + explicit wind stress over
      !! every owned cell.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: k
      real(wp) :: h_val, u_val, v_val, speed, cf, denom
      real(wp) :: n_sq, min_h, tau_x_over_rho, tau_y_over_rho
      logical :: has_friction, has_wind

      if (dt <= 0.0_wp) return

      has_friction = state%barotropic%manning_n > 0.0_wp
      has_wind = (state%barotropic%wind_tau_x /= 0.0_wp) .or. (state%barotropic%wind_tau_y /= 0.0_wp)
      if (.not. has_friction .and. .not. has_wind) return

      n_sq = state%barotropic%manning_n*state%barotropic%manning_n
      min_h = state%minimum_allowed_height
      tau_x_over_rho = state%barotropic%wind_tau_x/RHO_WATER
      tau_y_over_rho = state%barotropic%wind_tau_y/RHO_WATER

      !$omp target teams distribute parallel do private(h_val, u_val, v_val, speed, cf, denom)
      do k = 1, mesh%ncells_owned
         h_val = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)
         if (h_val > min_h) then
            if (has_friction) then
               u_val = state%barotropic%xmom_cv(k)/h_val
               v_val = state%barotropic%ymom_cv(k)/h_val
               speed = sqrt(u_val*u_val + v_val*v_val)
               cf = GRAVITY*n_sq/(h_val**(1.0_wp/3.0_wp))
               denom = 1.0_wp + dt*cf*speed/h_val
               state%barotropic%xmom_cv(k) = state%barotropic%xmom_cv(k)/denom
               state%barotropic%ymom_cv(k) = state%barotropic%ymom_cv(k)/denom
            end if
            if (has_wind) then
               state%barotropic%xmom_cv(k) = state%barotropic%xmom_cv(k) + dt*tau_x_over_rho
               state%barotropic%ymom_cv(k) = state%barotropic%ymom_cv(k) + dt*tau_y_over_rho
            end if
         end if
      end do
      !$omp end target teams distribute parallel do
   end subroutine apply_manning_friction_unstr

end module rki_kernel_friction_unstr
