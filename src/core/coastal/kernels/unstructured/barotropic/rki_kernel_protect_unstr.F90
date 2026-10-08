!! Protect against negative depths in the unstructured solver
module rki_kernel_protect_unstr
   !! Negative-depth guard after each flux + update step. Cells with
   !! depth `h = stage - bed` below `DRY_TOLERANCE` have their centroid
   !! momenta zeroed; thin-but-wet cells (`h < eps = DESING_EPS`,
   !! the Kurganov & Petrova 2007 desingularization length shared with the
   !! structured kernels) get the desingularized velocity
   !! `u = sqrt(2) h uh / sqrt(h^4 + eps^4)` and `uh := h*u`, which
   !! bounds thin-cell velocity smoothly instead of freezing thin water
   !! (the old hard zero below `minimum_allowed_height` pinned wet/dry
   !! fronts); a `FROUDE_BACKSTOP` rescale arrests thin-film velocity
   !! runaway (see the constant's docstring for the measurements).
   !! Strictly-negative cells clip stage back to bed (h=0) and
   !! accumulate the removed mass `(-h)*area` into `state%mass_error`.
   !! Uses `!$acc parallel loop reduction(+:..)` (array write + scalar
   !! reduction don't fit `do concurrent`).
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, DESING_EPS, &
                            FROUDE_BACKSTOP
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   implicit none
   private

   public :: protect_unstr

contains

   subroutine protect_unstr(state, mesh)
      !! Apply the negative-depth and shallow-cell guards. Adds the
      !! clipped negative mass to `state%mass_error` (does not reset).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k
      real(wp) :: h
      real(wp) :: eps4, h4, u_ds, v_ds, speed_sq, fr_max_sq, fr_scale
      real(wp) :: mass_error_local

      mass_error_local = 0.0_wp
      eps4 = DESING_EPS**4

      !$omp target teams distribute parallel do private(h, h4, u_ds, v_ds, speed_sq, fr_max_sq, fr_scale)
      do k = 1, mesh%ncells_owned
         h = state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k)

         if (h < DRY_TOLERANCE) then
            ! Dry (or negative) -- no meaningful velocity; zero momentum.
            state%barotropic%xmom_cv(k) = 0.0_wp
            state%barotropic%ymom_cv(k) = 0.0_wp
         else
            h4 = h**4
            if (h4 < eps4) then
               ! u = sqrt(2) h uh / sqrt(h^4 + max(h^4, eps^4)); uh := h*u.
               ! Identity for h >= eps, so applied only below eps.
               u_ds = sqrt(2.0_wp)*h*state%barotropic%xmom_cv(k)/sqrt(h4 + eps4)
               v_ds = sqrt(2.0_wp)*h*state%barotropic%ymom_cv(k)/sqrt(h4 + eps4)
               state%barotropic%xmom_cv(k) = h*u_ds
               state%barotropic%ymom_cv(k) = h*v_ds
            end if
            ! Froude backstop: blowup arrester only (triggers on a TRUE
            ! comparison, so NaN passes through and fails loudly).
            speed_sq = (state%barotropic%xmom_cv(k)**2 &
                        + state%barotropic%ymom_cv(k)**2)/h**2
            fr_max_sq = FROUDE_BACKSTOP**2*GRAVITY*h
            if (speed_sq > fr_max_sq) then
               fr_scale = sqrt(fr_max_sq/speed_sq)
               state%barotropic%xmom_cv(k) = state%barotropic%xmom_cv(k)*fr_scale
               state%barotropic%ymom_cv(k) = state%barotropic%ymom_cv(k)*fr_scale
            end if
         end if

         if (h < 0.0_wp) then
            ! Strictly negative -- clip stage to bed and remember the
            ! unphysical mass we just removed.
            mass_error_local = mass_error_local + (-h)*mesh%area(k)
            state%barotropic%stage_cv(k) = state%barotropic%bed_cv(k)
         end if
      end do
      !$omp end target teams distribute parallel do

      ! Accumulate into the running mass-error tally on the state.
      state%mass_error = state%mass_error + mass_error_local
   end subroutine protect_unstr

end module rki_kernel_protect_unstr
