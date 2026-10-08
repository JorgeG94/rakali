!! Horizontal viscosity for layer momentum on the unstructured grid.
module rki_ml_horizontal_viscosity_unstr
   !! Smagorinsky (1963) eddy-viscosity closure on triangles.
   !! Jacobi snapshots (`ml_horvisc_hu_pre/hv_pre`) live on `unstr_state_t`
   !! (not module-level), per the "every allocatable component must be
   !! attached" rule — see `state_unstr_init` / `state_unstr_enter_data`.
   use rki_constants, only: wp, DRY_TOLERANCE, THIN_LAYER_THRESHOLD
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR
   implicit none
   private

   public :: ml_horizontal_viscosity_smag_unstr

contains

   subroutine ml_horizontal_viscosity_smag_unstr(state, mesh, dt)
      !! Smagorinsky (1963) horizontal viscosity for layer momentum on triangles.
      !! Velocity gradients via Green-Gauss; nu_h = (Cs^2*A)*|S|, stability-limited;
      !! FVM Laplacian. Jacobi snapshot on hu/hv for CPU/GPU consistency.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: cell, k, e, nbr, nz, n
      real(wp) :: cs, cs2_area, inv_area, nu_max_cell
      real(wp) :: h_k, u_c, v_c, u_nbr, v_nbr, u_edge, v_edge
      real(wp) :: dudx, dvdy, dudy, dvdx, strain, nu_h
      real(wp) :: dx_cn, dy_cn, d_cn, lap_u, lap_v

      cs = state%multilayer%smagorinsky_cs
      if (cs <= 0.0_wp) return
      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = state%ncells_owned

      ! Jacobi snapshots live on the state struct (persistent, GPU-resident).
      associate (hu_pre => state%multilayer%ml_horvisc_hu_pre, &
                 hv_pre => state%multilayer%ml_horvisc_hv_pre)

         !$omp target teams distribute parallel do collapse(2)
         do k = 1, nz
         do cell = 1, mesh%ncells_total
            hu_pre(k, cell) = state%multilayer%hu_layer(k, cell)
            hv_pre(k, cell) = state%multilayer%hv_layer(k, cell)
         end do
         end do
         !$omp end target teams distribute parallel do

         !$omp target teams distribute parallel do private(k, e, nbr, cs2_area, inv_area, nu_max_cell, h_k, u_c, v_c, u_nbr, v_nbr, u_edge, v_edge, dudx, dvdy, dudy, dvdx, strain, nu_h, dx_cn, dy_cn, d_cn, lap_u, lap_v)
         do cell = 1, n

            if (state%barotropic%height_cv(cell) < ML_MIN_DEPTH_UNSTR) cycle

            inv_area = 1.0_wp/mesh%area(cell)
            cs2_area = cs*cs*mesh%area(cell)
            nu_max_cell = 0.125_wp*mesh%radius(cell)**2/(dt + 1.0e-30_wp)

            do k = 1, nz
               h_k = state%multilayer%h_layer(k, cell)
               if (h_k < THIN_LAYER_THRESHOLD) cycle

               u_c = hu_pre(k, cell)/h_k
               v_c = hv_pre(k, cell)/h_k

               ! ----- Velocity gradient via Green-Gauss -----
               dudx = 0.0_wp
               dvdy = 0.0_wp
               dudy = 0.0_wp
               dvdx = 0.0_wp

               do e = 1, 3
                  nbr = mesh%surrogate_neighbours(e, cell)
                  if (nbr == cell) then
                     u_edge = u_c
                     v_edge = v_c
                  else if (state%multilayer%h_layer(k, nbr) < THIN_LAYER_THRESHOLD) then
                     u_edge = u_c
                     v_edge = v_c
                  else
                     u_nbr = hu_pre(k, nbr)/state%multilayer%h_layer(k, nbr)
                     v_nbr = hv_pre(k, nbr)/state%multilayer%h_layer(k, nbr)
                     u_edge = 0.5_wp*(u_c + u_nbr)
                     v_edge = 0.5_wp*(v_c + v_nbr)
                  end if

                  dudx = dudx + u_edge*mesh%edge_n1(e, cell)*mesh%edge_len(e, cell)
                  dudy = dudy + u_edge*mesh%edge_n2(e, cell)*mesh%edge_len(e, cell)
                  dvdx = dvdx + v_edge*mesh%edge_n1(e, cell)*mesh%edge_len(e, cell)
                  dvdy = dvdy + v_edge*mesh%edge_n2(e, cell)*mesh%edge_len(e, cell)
               end do
               dudx = dudx*inv_area
               dudy = dudy*inv_area
               dvdx = dvdx*inv_area
               dvdy = dvdy*inv_area

               ! Strain rate magnitude
               strain = sqrt((dudx - dvdy)**2 + (dudy + dvdx)**2)

               ! Smagorinsky eddy viscosity with stability limiter
               nu_h = min(cs2_area*strain, nu_max_cell)

               ! ----- FVM Laplacian of velocity -----
               lap_u = 0.0_wp
               lap_v = 0.0_wp

               do e = 1, 3
                  nbr = mesh%surrogate_neighbours(e, cell)
                  if (nbr == cell) cycle  ! boundary — zero gradient

                  if (state%multilayer%h_layer(k, nbr) < THIN_LAYER_THRESHOLD) cycle

                  u_nbr = hu_pre(k, nbr)/state%multilayer%h_layer(k, nbr)
                  v_nbr = hv_pre(k, nbr)/state%multilayer%h_layer(k, nbr)

                  dx_cn = mesh%centroid_x(nbr) - mesh%centroid_x(cell)
                  dy_cn = mesh%centroid_y(nbr) - mesh%centroid_y(cell)
                  d_cn = sqrt(dx_cn*dx_cn + dy_cn*dy_cn)

                  lap_u = lap_u + (u_nbr - u_c)*mesh%edge_len(e, cell)/d_cn
                  lap_v = lap_v + (v_nbr - v_c)*mesh%edge_len(e, cell)/d_cn
               end do
               lap_u = lap_u*inv_area
               lap_v = lap_v*inv_area

               ! Update momentum
               state%multilayer%hu_layer(k, cell) = hu_pre(k, cell) + dt*h_k*nu_h*lap_u
               state%multilayer%hv_layer(k, cell) = hv_pre(k, cell) + dt*h_k*nu_h*lap_v
            end do
         end do
         !$omp end target teams distribute parallel do

      end associate

   end subroutine ml_horizontal_viscosity_smag_unstr

end module rki_ml_horizontal_viscosity_unstr
