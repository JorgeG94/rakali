!! Explicit update and SSP-RK2 helpers for the unstructured solver
module rki_kernel_update_unstr
   !! Cell-centric kernels that combine the flux-pass divergence arrays
   !! `(stage_eu, xmom_eu, ymom_eu)` into the next time level:
   !!   - backup/restore conserved centroid arrays into `*_backup`,
   !!   - forward Euler `Q := Q + dt*Q_eu`,
   !!   - SSP-RK2 average `Q := a*Q + b*Q_backup`.
   !! All loops are `do concurrent` over owned cells.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   implicit none
   private

   public :: backup_conserved_unstr
   public :: restore_conserved_unstr
   public :: apply_explicit_update_unstr
   public :: saxpy_conserved_unstr

contains

   subroutine backup_conserved_unstr(state, mesh)
      !! Save (stage_cv, xmom_cv, ymom_cv) into the *_backup arrays.
      !! Called once at the start of each SSP-RK2 step.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k

      !$omp target teams distribute parallel do
      do k = 1, mesh%ncells_owned
         state%barotropic%stage_backup(k) = state%barotropic%stage_cv(k)
         state%barotropic%xmom_backup(k) = state%barotropic%xmom_cv(k)
         state%barotropic%ymom_backup(k) = state%barotropic%ymom_cv(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine backup_conserved_unstr

   subroutine restore_conserved_unstr(state, mesh)
      !! Inverse of backup_conserved_unstr (restores from `*_backup`).
      !! Public only for the unit-test suite; no production module uses it.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k

      !$omp target teams distribute parallel do
      do k = 1, mesh%ncells_owned
         state%barotropic%stage_cv(k) = state%barotropic%stage_backup(k)
         state%barotropic%xmom_cv(k) = state%barotropic%xmom_backup(k)
         state%barotropic%ymom_cv(k) = state%barotropic%ymom_backup(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine restore_conserved_unstr

   subroutine apply_explicit_update_unstr(state, mesh, dt)
      !! Forward Euler update `Q := Q + dt * Q_eu` over owned cells.
      !! Requires the flux kernel to have written the explicit-update
      !! arrays first. `dt` is THIS RK substage's step, from the driver.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: k

      !$omp target teams distribute parallel do
      do k = 1, mesh%ncells_owned
         state%barotropic%stage_cv(k) = state%barotropic%stage_cv(k) + dt*state%barotropic%stage_eu(k)
         state%barotropic%xmom_cv(k) = state%barotropic%xmom_cv(k) + dt*state%barotropic%xmom_eu(k)
         state%barotropic%ymom_cv(k) = state%barotropic%ymom_cv(k) + dt*state%barotropic%ymom_eu(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine apply_explicit_update_unstr

   subroutine saxpy_conserved_unstr(state, mesh, a, b)
      !! Linear combination `Q := a*Q + b*Q_backup` over the conserved
      !! centroid arrays. SSP-RK2 average uses `(a, b) = (0.5, 0.5)`.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: a
      real(wp), intent(in) :: b

      integer :: k

      !$omp target teams distribute parallel do
      do k = 1, mesh%ncells_owned
         state%barotropic%stage_cv(k) = a*state%barotropic%stage_cv(k) + b*state%barotropic%stage_backup(k)
         state%barotropic%xmom_cv(k) = a*state%barotropic%xmom_cv(k) + b*state%barotropic%xmom_backup(k)
         state%barotropic%ymom_cv(k) = a*state%barotropic%ymom_cv(k) + b*state%barotropic%ymom_backup(k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine saxpy_conserved_unstr

end module rki_kernel_update_unstr
