!! CFL-based velocity truncation for layer momentum (unstructured twin).
module rki_ml_truncate_unstr
   !! Per-cell-per-layer limiter (not a closure) that clips runaway
   !! velocities from any upstream operator (BPG, vertical exchange,
   !! viscosity) to `0.9*cfl_trunc*L/dt` so a few bad cells don't tank the
   !! global timestep. Unstructured twin of structured `ml_truncate_velocity`.
   !!
   !! The isotropic length scale on a triangle is the inscribed-circle CFL
   !! radius `mesh%radius(cell)` (the same measure the barotropic CFL uses),
   !! replacing the structured `dx`/`dy` pair. Column-local — safe for
   !! `do concurrent`. Layer index is bottom-up (`k=1` bed, `k=nz` surface).
   use rki_constants, only: wp, THIN_LAYER_THRESHOLD
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   implicit none
   private

   public :: ml_truncate_velocity_unstr

contains

   subroutine ml_truncate_velocity_unstr(state, mesh, dt, cfl_trunc)
      !! Clip any layer velocity with CFL > cfl_trunc to
      !! `0.9*cfl_trunc*radius/dt`, per owned cell and layer. Skips layers
      !! thinner than THIN_LAYER_THRESHOLD (already zeroed upstream).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt
      real(wp), intent(in) :: cfl_trunc

      integer :: cell, k, nz, n
      real(wp) :: vel_max, h_k, u_abs, v_abs, inv_dt, coef

      if (.not. state%use_multilayer) return
      if (cfl_trunc <= 0.0_wp) return

      nz = state%multilayer%nz_ml
      n = state%ncells_owned
      inv_dt = 1.0_wp/max(dt, 1.0e-10_wp)
      coef = 0.9_wp*cfl_trunc*inv_dt

      !$omp parallel do private(k, vel_max, h_k, u_abs, v_abs)
      do cell = 1, n
         vel_max = coef*mesh%radius(cell)
         do k = 1, nz
            h_k = state%multilayer%h_layer(k, cell)
            if (h_k < THIN_LAYER_THRESHOLD) cycle

            u_abs = abs(state%multilayer%hu_layer(k, cell))/h_k
            v_abs = abs(state%multilayer%hv_layer(k, cell))/h_k

            if (u_abs > vel_max) then
               state%multilayer%hu_layer(k, cell) = &
                  sign(vel_max*h_k, state%multilayer%hu_layer(k, cell))
            end if
            if (v_abs > vel_max) then
               state%multilayer%hv_layer(k, cell) = &
                  sign(vel_max*h_k, state%multilayer%hv_layer(k, cell))
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine ml_truncate_velocity_unstr

end module rki_ml_truncate_unstr
