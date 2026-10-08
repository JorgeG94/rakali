!! Sponge layer relaxation kernel for the unstructured solver
module rki_kernel_sponge_unstr
   !! Explicit volumetric damping of h and momentum toward a reference
   !! state at cells adjacent to boundary edges tagged `BC_SPONGE`.
   !!
   !! Algorithm: for each owned cell, scan the three edges.  If any edge
   !! is a boundary edge (`mesh%neighbours(e, cell) < 0`) belonging to a
   !! BC_SPONGE group, collect the group's `strength` [1/s] and `h_ref`
   !! [m]. When multiple sponge edges exist (a corner cell), use the
   !! maximum strength and its paired `h_ref`.  Then apply one forward-
   !! Euler relaxation step:
   !!
   !!   h <- h - dt*alpha*(h - h_ref)
   !!   xmom <- xmom*(1 - dt*alpha)
   !!   ymom <- ymom*(1 - dt*alpha)
   !!
   !! where `h = stage_cv - bed_cv` and `alpha = bc_sponge_strength(g)`.
   !!
   !! The kernel self-gates: cells not adjacent to any active BC_SPONGE
   !! edge are skipped implicitly (alpha stays at 0).  No-op when
   !! `strength <= 0` for all sponge groups.
   !!
   !! v1 limitation: the sponge band is one cell deep — only cells with
   !! a direct boundary edge are affected.  Graded-width profiles require
   !! cell-distance information that is not available on triangular meshes
   !! without additional mesh topology data.  On the structured side the
   !! distance is implicit in the (i,j) indices; no such shortcut exists
   !! here.  Multi-cell sponge widths are deferred.
   use rki_constants, only: wp
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_boundary_unstr, only: BC_SPONGE
   implicit none
   private

   public :: apply_sponge_unstr

contains

   pure subroutine apply_sponge_unstr(state, mesh, dt)
      !! Apply sponge relaxation to all owned cells adjacent to a BC_SPONGE edge.
      !!
      !! Relaxes `h = stage_cv - bed_cv` toward `bc_sponge_h_ref(g)` and
      !! `xmom_cv`/`ymom_cv` toward 0.  No-op when no sponge edges exist or
      !! all sponge strengths are 0.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in)    :: mesh
      real(wp), intent(in)    :: dt

      integer :: cell, e, g
      real(wp) :: alpha, h_ref_best, h_water, h_new

      do concurrent(cell=1:state%ncells_owned) &
         local(e, g, alpha, h_ref_best, h_water, h_new)

         ! Scan all three edges for the highest-strength sponge group.
         alpha = 0.0_wp
         h_ref_best = 0.0_wp
         do e = 1, 3
            if (mesh%neighbours(e, cell) < 0) then
               g = mesh%bc_group(e, cell)
               if (g > 0) then
                  if (state%bc_type(g) == BC_SPONGE .and. &
                      state%bc_sponge_strength(g) > alpha) then
                     alpha = state%bc_sponge_strength(g)
                     h_ref_best = state%bc_sponge_h_ref(g)
                  end if
               end if
            end if
         end do

         if (alpha > 0.0_wp) then
            h_water = state%barotropic%stage_cv(cell) - state%barotropic%bed_cv(cell)
            h_new = h_water - dt*alpha*(h_water - h_ref_best)
            state%barotropic%stage_cv(cell) = h_new + state%barotropic%bed_cv(cell)
            state%barotropic%xmom_cv(cell) = state%barotropic%xmom_cv(cell)*(1.0_wp - dt*alpha)
            state%barotropic%ymom_cv(cell) = state%barotropic%ymom_cv(cell)*(1.0_wp - dt*alpha)
         end if

      end do

   end subroutine apply_sponge_unstr

end module rki_kernel_sponge_unstr
