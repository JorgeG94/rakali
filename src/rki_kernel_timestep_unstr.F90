!! CFL timestep kernel for the unstructured solver
module rki_kernel_timestep_unstr
   !! Global CFL-limited timestep from the per-cell wave speeds the flux
   !! kernel writes into `state%barotropic%max_speed(k)`. Per cell
   !!   dt_k = mesh%radius(k) / max_speed(k)
   !! where `radius(k)` is the min centroid-to-edge-midpoint distance
   !! (ANUGA `radii`); global `dt = cfl * min_k dt_k`. Cells with
   !! `max_speed` below `epsilon` are skipped. Reduction via
   !! `!$acc parallel loop reduction(min:..)`.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   implicit none
   private

   public :: compute_timestep_unstr
   public :: TIMESTEP_LARGE

   real(wp), parameter :: TIMESTEP_EPSILON = 1.0e-15_wp
      !! Cells with `max_speed` below this contribute nothing to the
      !! timestep reduction (no characteristic to constrain).
   real(wp), parameter :: TIMESTEP_LARGE = huge(1.0_wp)
      !! Initial running minimum (effectively infinity); `huge` fits
      !! either precision (`1.0e+100_wp` overflows `real32`).

contains

   subroutine compute_timestep_unstr(state, mesh)
      !! Compute the CFL-limited timestep and store it in `state%dt`.
      !! Reads `max_speed`, `state%cfl`, `mesh%radius` only. Must run
      !! after the flux kernel; before any flux pass it returns
      !! `cfl * TIMESTEP_LARGE` (all max_speed still 0).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k
      real(wp) :: dt_local, cell_dt

      dt_local = TIMESTEP_LARGE

      do concurrent(k=1:mesh%ncells_owned) reduce(min:dt_local)
         if (state%barotropic%max_speed(k) > TIMESTEP_EPSILON) then
            cell_dt = mesh%radius(k)/state%barotropic%max_speed(k)
            if (cell_dt < dt_local) dt_local = cell_dt
         end if
      end do

      ! Clamp +inf back to TIMESTEP_LARGE: NVHPC's reduction(min:..) can
      ! return +inf when the loop never updates the accumulator (all
      ! max_speed = 0), so the "no characteristics => dt = cfl *
      ! TIMESTEP_LARGE" contract holds bit-exactly.
      if (.not. (dt_local <= TIMESTEP_LARGE)) dt_local = TIMESTEP_LARGE

      state%dt = state%cfl*dt_local
   end subroutine compute_timestep_unstr

end module rki_kernel_timestep_unstr
