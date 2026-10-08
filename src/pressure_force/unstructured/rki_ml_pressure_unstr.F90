!! Baroclinic pressure gradient for the unstructured multilayer solver.
module rki_ml_pressure_unstr
   !! Two BPG variants:
   !! * `ml_baroclinic_pressure_gradient_unstr` (sigma method): integrates
   !!   the density anomaly bottom-up at the centre column, Green-Gauss
   !!   horizontal gradient at constant sigma. Carries the sigma PGE.
   !! * `ml_baroclinic_pressure_gradient_jacobian_unstr` (Shchepetkin &
   !!   McWilliams 2003): interpolates each neighbour column's pressure to
   !!   the physical depth of the centre column's layer centre, then
   !!   differences at constant z. Uniform density ⇒ zero gradient,
   !!   eliminating the sigma PGE. Default for the regional path.
   !! Both subtract `rho_ref(z)` (from `ml_compute_reference_density_unstr`)
   !! before differencing to cut the residual sigma PGE in stratified water.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR, BPG_NZ_MAX, &
                                  VANISHING_LAYER_TOL_UNSTR
   use rki_ml_eos_unstr, only: rho_ref_at
   implicit none
   private

   public :: ml_baroclinic_pressure_gradient_unstr
   public :: ml_baroclinic_pressure_gradient_jacobian_unstr

contains

   subroutine ml_baroclinic_pressure_gradient_unstr(state, mesh, dt)
      !! Apply baroclinic pressure gradient as explicit source term.
      !!
      !! The hydrostatic internal pressure anomaly at layer k center is:
      !!   p_int(k) = g * [sum_{j=k+1}^{nz} (rho_j - rho_0) * h_j
      !!                   + 0.5 * (rho_k - rho_0) * h_k]
      !!
      !! The baroclinic PGF for layer k is:
      !!   BPG_k = -(1/rho_0) * grad(p_int(k))
      !!
      !! Applied as:
      !!   hu_k += dt * h_k * BPG_k^x
      !!   hv_k += dt * h_k * BPG_k^y
      !!
      !! Uses Green-Gauss for the gradient on triangles:
      !!   grad(p) = (1/A) * sum_e p_edge * n_e * L_e
      !!
      !! Includes a reduced-gravity acceleration cap: sqrt(g'*H)/dt
      !! to prevent sigma PGE from creating supersonic velocities.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: cell, k, kn, e, nbr, nz, n, nz_ref
      real(wp) :: inv_area, inv_rho0, rho_ref_dz, rho_0
      real(wp) :: acc, p_here, p_nbr, p_edge
      real(wp) :: dp_dx, dp_dy, n1, n2, L_e
      real(wp) :: g_prime, H_col, bpg_max, bpg_mag, bpg_scale
      real(wp) :: z_top, z_bot, z_mid_lay
      ! Stack workspace for internal pressure (nz values per cell)
      real(wp) :: p_int(BPG_NZ_MAX)
      real(wp) :: p_int_nbr(BPG_NZ_MAX)

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = mesh%ncells_owned

      if (nz > BPG_NZ_MAX) then
         ! Defence-in-depth only: validate_config refuses unstructured
         ! multilayer configs with nz + 1 > BPG_NZ_MAX up front
         ! (bpg_nz_is_sufficient in rki_ml_common_unstr), so this cannot
         ! fire on a validated run.  Kept so a caller that bypasses
         ! config validation overruns nothing — but note the cost of
         ! reaching it is silently skipping the whole baroclinic
         ! pressure gradient.
         return
      end if

      inv_rho0 = 1.0_wp/state%multilayer%rho_0
      rho_0 = state%multilayer%rho_0
      rho_ref_dz = state%multilayer%rho_ref_dz
      nz_ref = state%multilayer%nz_ref

      !$omp parallel do private(k, kn, e, nbr, inv_area, acc, p_here, p_nbr, p_edge, dp_dx, dp_dy, n1, n2, L_e, g_prime, H_col, bpg_max, bpg_mag, bpg_scale, z_top, z_bot, z_mid_lay, p_int, p_int_nbr)
      do cell = 1, n
         ! Skip thin columns
         H_col = state%barotropic%height_cv(cell)
         if (H_col < ML_MIN_DEPTH_UNSTR) cycle

         inv_area = 1.0_wp/mesh%area(cell)

         ! ---- Compute internal pressure at this cell's layer centers ----
         ! Integrate from surface (k = nz) downward to bed (k = 1) using
         ! the density anomaly `(rho - rho_ref(z))` (Phase B2). When the
         ! reference profile is unset (rho_ref_dz==0) the helper returns
         ! `rho_0`, recovering the pre-B2 anomaly.
         acc = 0.0_wp
         z_top = 0.0_wp
         do k = nz, 1, -1
            z_bot = z_top + state%multilayer%h_layer(k, cell)
            z_mid_lay = 0.5_wp*(z_top + z_bot)
            acc = acc + (state%multilayer%rho_layer(k, cell) &
                         - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, rho_ref_dz, &
                                      nz_ref, rho_0)) &
                  *state%multilayer%h_layer(k, cell)
            p_int(k) = acc - 0.5_wp*(state%multilayer%rho_layer(k, cell) &
                                     - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, &
                                                  rho_ref_dz, nz_ref, rho_0)) &
                       *state%multilayer%h_layer(k, cell)
            z_top = z_bot
         end do

         ! ---- Compute gradient using Green-Gauss ----
         ! dp/dx = (1/A) * sum_e p_edge * n1 * L
         ! dp/dy = (1/A) * sum_e p_edge * n2 * L

         do k = 1, nz
            dp_dx = 0.0_wp
            dp_dy = 0.0_wp
            p_here = p_int(k)

            do e = 1, 3
               nbr = mesh%neighbours(e, cell)
               n1 = mesh%edge_n1(e, cell)
               n2 = mesh%edge_n2(e, cell)
               L_e = mesh%edge_len(e, cell)

               if (nbr > 0) then
                  ! Interior edge: compute neighbour's internal pressure
                  ! We need p_int at layer k for the neighbour
                  if (state%barotropic%height_cv(nbr) >= ML_MIN_DEPTH_UNSTR) then
                     acc = 0.0_wp
                     z_top = 0.0_wp
                     do kn = nz, 1, -1
                        z_bot = z_top + state%multilayer%h_layer(kn, nbr)
                        z_mid_lay = 0.5_wp*(z_top + z_bot)
                        acc = acc + (state%multilayer%rho_layer(kn, nbr) &
                                     - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, &
                                                  rho_ref_dz, nz_ref, rho_0)) &
                              *state%multilayer%h_layer(kn, nbr)
                        if (kn == k) then
                           p_nbr = acc - 0.5_wp*(state%multilayer%rho_layer(k, nbr) &
                                                 - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, &
                                                              rho_ref_dz, nz_ref, rho_0)) &
                                   *state%multilayer%h_layer(k, nbr)
                        end if
                        z_top = z_bot
                     end do
                  else
                     p_nbr = p_here  ! Dry neighbour: no gradient
                  end if
                  p_edge = 0.5_wp*(p_here + p_nbr)
               else
                  ! Boundary edge: extrapolate (use cell value)
                  p_edge = p_here
               end if

               ! Green-Gauss contribution
               dp_dx = dp_dx + p_edge*n1*L_e
               dp_dy = dp_dy + p_edge*n2*L_e
            end do

            dp_dx = dp_dx*inv_area
            dp_dy = dp_dy*inv_area

            ! Apply acceleration: a = -(g/rho_0) * dp/dx
            ! Momentum: hu += dt * h * a = -dt * h * (g/rho_0) * dp/dx
            ! But p_int already has g factored in? No - p_int = (rho-rho_0)*h
            ! So we need to multiply by g: a = -(g/rho_0) * g * dp/dx = -g²/rho_0 * dp/dx
            ! Wait, let me check the structured code...
            ! p_int(k) = acc - 0.5*(rho-rho_0)*h  (units: kg/m² = pressure/g)
            ! dp_dx has units kg/m³
            ! Then: GRAVITY * inv_rho0 * dp_dx gives acceleration (m/s²)
            ! So we multiply dp by GRAVITY to get pressure gradient, then by 1/rho_0
            dp_dx = GRAVITY*inv_rho0*dp_dx
            dp_dy = GRAVITY*inv_rho0*dp_dy

            ! Limiter: cap BPG acceleration to sqrt(g'*H)/dt
            ! g' = g * (rho_bot - rho_top) / rho_0
            g_prime = GRAVITY*abs(state%multilayer%rho_layer(nz, cell) - state%multilayer%rho_layer(1, cell)) &
                      *inv_rho0
            g_prime = max(g_prime, 1.0e-6_wp)  ! Floor to avoid sqrt(0)
            bpg_max = sqrt(g_prime*H_col)/max(dt, 1.0e-10_wp)

            bpg_mag = sqrt(dp_dx*dp_dx + dp_dy*dp_dy)
            if (bpg_mag > bpg_max .and. bpg_mag > 0.0_wp) then
               bpg_scale = bpg_max/bpg_mag
               dp_dx = dp_dx*bpg_scale
               dp_dy = dp_dy*bpg_scale
            end if

            ! Apply to layer momentum.
            ! VCOORD_ZSTAR_FULL: vanishing layers carry no momentum.
            if (state%multilayer%h_layer(k, cell) > VANISHING_LAYER_TOL_UNSTR) then
               state%multilayer%hu_layer(k, cell) = state%multilayer%hu_layer(k, cell) &
                                                    - dt*state%multilayer%h_layer(k, cell)*dp_dx
               state%multilayer%hv_layer(k, cell) = state%multilayer%hv_layer(k, cell) &
                                                    - dt*state%multilayer%h_layer(k, cell)*dp_dy
            else
               state%multilayer%hu_layer(k, cell) = 0.0_wp
               state%multilayer%hv_layer(k, cell) = 0.0_wp
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine ml_baroclinic_pressure_gradient_unstr

   subroutine ml_baroclinic_pressure_gradient_jacobian_unstr(state, mesh, dt)
      !! Apply the Jacobian (Shchepetkin-McWilliams 2003) baroclinic
      !! pressure gradient as an explicit momentum source.
      !!
      !! Per cell, per layer k:
      !!   1. z_mid_c = depth of layer k centre in the centre column
      !!   2. p_here  = integrated density anomaly from surface to z_mid_c
      !!   3. For each of 3 edges:
      !!        - interior with wet neighbour whose column reaches z_mid_c:
      !!          build the neighbour's (zn_if, pn_if) column and
      !!          linear-interpolate p to z_mid_c.
      !!        - else (dry/missing neighbour, OR z_mid_c below the
      !!          neighbour's bed at a shelf break): p_nbr = p_here
      !!          (reflective wall) -- avoids the spurious horizontal
      !!          gradient from clamping to neighbour-bed pressure.
      !!   4. p_edge = 0.5 * (p_here + p_nbr); accumulate
      !!        dp_dx += p_edge * n1 * L_e
      !!        dp_dy += p_edge * n2 * L_e
      !!   5. (dp_dx, dp_dy) *= 1/area; multiply by g/rho_0; cap by
      !!      sqrt(g'*H)/dt; subtract dt*h*acc from layer momentum.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: dt

      integer :: cell, k, kn, e, nbr, nz, n, nz_ref
      real(wp) :: inv_area, inv_rho0, H_col, rho_ref_dz, rho_0
      real(wp) :: z_if(BPG_NZ_MAX + 1), p_if(BPG_NZ_MAX + 1)
      ! Per-edge neighbour columns. Built once per cell in the prepass
      ! and then reused for every layer's search/interpolate. Hoisting
      ! these out of the inner `k` loop is what cuts BPG runtime in
      ! half: the column data is independent of `z_mid`, only the
      ! query into it depends on `k`.
      real(wp) :: zn_if_e(BPG_NZ_MAX + 1, 3), pn_if_e(BPG_NZ_MAX + 1, 3)
      integer  :: nbr_idx(3)
      logical  :: nbr_valid(3)
      real(wp) :: z_mid, z_mid_lay, p_here, p_nbr, p_edge, frac
      real(wp) :: dp_dx, dp_dy, n1, n2, L_e
      real(wp) :: g_prime, bpg_max, bpg_mag, bpg_scale

      if (.not. state%use_multilayer) return

      nz = state%multilayer%nz_ml
      n = mesh%ncells_owned

      ! Defence-in-depth only: validate_config (bpg_nz_is_sufficient in
      ! rki_ml_common_unstr) refuses nz + 1 > BPG_NZ_MAX up front, so
      ! this cannot fire on a validated run; reaching it would silently
      ! skip the whole baroclinic pressure gradient.
      if (nz > BPG_NZ_MAX) return

      inv_rho0 = 1.0_wp/state%multilayer%rho_0
      rho_0 = state%multilayer%rho_0
      rho_ref_dz = state%multilayer%rho_ref_dz
      nz_ref = state%multilayer%nz_ref

      !$omp parallel do private(k, kn, e, nbr, inv_area, H_col, z_if, p_if, zn_if_e, pn_if_e, nbr_idx, nbr_valid, z_mid, z_mid_lay, p_here, p_nbr, p_edge, frac, dp_dx, dp_dy, n1, n2, L_e, g_prime, bpg_max, bpg_mag, bpg_scale)
      do cell = 1, n
         H_col = state%barotropic%height_cv(cell)
         if (H_col < ML_MIN_DEPTH_UNSTR) cycle

         inv_area = 1.0_wp/mesh%area(cell)

         ! --- Centre column: interface depths and integrated pressure ---
         ! z_if(k) = depth from surface to the LOWER interface of layer k
         ! (z_if(nz+1) = 0 at surface, z_if(1) = bed depth).
         ! p_if(k) = integral of (rho - rho_ref(z)) * dz from surface to
         ! z_if(k). rho_ref(z) is the horizontally-uniform reference
         ! profile (Phase B2); when unset (rho_ref_dz==0) the helper
         ! returns rho_0, recovering the pre-B2 anomaly.
         z_if(nz + 1) = 0.0_wp
         p_if(nz + 1) = 0.0_wp
         do k = nz, 1, -1
            z_if(k) = z_if(k + 1) + state%multilayer%h_layer(k, cell)
            z_mid_lay = 0.5_wp*(z_if(k) + z_if(k + 1))
            p_if(k) = p_if(k + 1) &
                      + (state%multilayer%rho_layer(k, cell) &
                         - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, rho_ref_dz, &
                                      nz_ref, rho_0)) &
                      *state%multilayer%h_layer(k, cell)
         end do

         ! --- Build the 3 neighbour columns ONCE per cell (was: nz times) ---
         ! Each column is independent of `k`; the per-layer work in the
         ! inner loop is just a search/interpolate against these columns.
         ! Edges that are boundaries or face a thin neighbour are
         ! flagged invalid so the inner loop falls through to the
         ! reflective-wall fallback (`p_nbr = p_here`) without searching.
         do e = 1, 3
            nbr = mesh%neighbours(e, cell)
            nbr_idx(e) = nbr
            if (nbr > 0) then
               if (state%barotropic%height_cv(nbr) >= ML_MIN_DEPTH_UNSTR) then
                  nbr_valid(e) = .true.
                  zn_if_e(nz + 1, e) = 0.0_wp
                  pn_if_e(nz + 1, e) = 0.0_wp
                  do kn = nz, 1, -1
                     zn_if_e(kn, e) = zn_if_e(kn + 1, e) + state%multilayer%h_layer(kn, nbr)
                     z_mid_lay = 0.5_wp*(zn_if_e(kn, e) + zn_if_e(kn + 1, e))
                     pn_if_e(kn, e) = pn_if_e(kn + 1, e) &
                                      + (state%multilayer%rho_layer(kn, nbr) &
                                         - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, &
                                                      rho_ref_dz, nz_ref, rho_0)) &
                                      *state%multilayer%h_layer(kn, nbr)
                  end do
               else
                  nbr_valid(e) = .false.
               end if
            else
               nbr_valid(e) = .false.
            end if
         end do

         do k = 1, nz
            ! Centre column's layer k centre depth and p there.
            z_mid = 0.5_wp*(z_if(k) + z_if(k + 1))
            p_here = p_if(k + 1) &
                     + 0.5_wp*(state%multilayer%rho_layer(k, cell) &
                               - rho_ref_at(z_mid, state%multilayer%rho_ref, rho_ref_dz, &
                                            nz_ref, rho_0)) &
                     *state%multilayer%h_layer(k, cell)

            dp_dx = 0.0_wp
            dp_dy = 0.0_wp

            do e = 1, 3
               n1 = mesh%edge_n1(e, cell)
               n2 = mesh%edge_n2(e, cell)
               L_e = mesh%edge_len(e, cell)

               if (nbr_valid(e)) then
                  nbr = nbr_idx(e)
                  ! Default: z_mid below neighbour's bed -> shelf-break wall
                  ! (treat the edge as reflective so no spurious gradient is
                  ! generated by clamping to neighbour-bed pressure).
                  p_nbr = p_here
                  do kn = 1, nz
                     if (z_mid >= zn_if_e(kn + 1, e) .and. z_mid <= zn_if_e(kn, e)) then
                        if (state%multilayer%h_layer(kn, nbr) > DRY_TOLERANCE) then
                           frac = (z_mid - zn_if_e(kn + 1, e))/state%multilayer%h_layer(kn, nbr)
                        else
                           frac = 0.0_wp
                        end if
                        z_mid_lay = 0.5_wp*(zn_if_e(kn, e) + zn_if_e(kn + 1, e))
                        p_nbr = pn_if_e(kn + 1, e) + frac &
                                *(state%multilayer%rho_layer(kn, nbr) &
                                  - rho_ref_at(z_mid_lay, state%multilayer%rho_ref, &
                                               rho_ref_dz, nz_ref, rho_0)) &
                                *state%multilayer%h_layer(kn, nbr)
                        exit
                     end if
                  end do
               else
                  ! Boundary edge or thin neighbour column: reflective wall.
                  p_nbr = p_here
               end if

               p_edge = 0.5_wp*(p_here + p_nbr)
               dp_dx = dp_dx + p_edge*n1*L_e
               dp_dy = dp_dy + p_edge*n2*L_e
            end do

            dp_dx = dp_dx*inv_area
            dp_dy = dp_dy*inv_area

            ! Convert from integrated-density-anomaly to acceleration:
            !   a = -(g/rho_0) * dp/dx  (and likewise y)
            dp_dx = GRAVITY*inv_rho0*dp_dx
            dp_dy = GRAVITY*inv_rho0*dp_dy

            ! Reduced-gravity CFL cap (same as the sigma kernel).
            g_prime = GRAVITY*abs(state%multilayer%rho_layer(nz, cell) - state%multilayer%rho_layer(1, cell)) &
                      *inv_rho0
            g_prime = max(g_prime, 1.0e-6_wp)
            bpg_max = sqrt(g_prime*H_col)/max(dt, 1.0e-10_wp)

            bpg_mag = sqrt(dp_dx*dp_dx + dp_dy*dp_dy)
            if (bpg_mag > bpg_max .and. bpg_mag > 0.0_wp) then
               bpg_scale = bpg_max/bpg_mag
               dp_dx = dp_dx*bpg_scale
               dp_dy = dp_dy*bpg_scale
            end if

            ! VCOORD_ZSTAR_FULL: vanishing layers carry no momentum.
            if (state%multilayer%h_layer(k, cell) > VANISHING_LAYER_TOL_UNSTR) then
               state%multilayer%hu_layer(k, cell) = state%multilayer%hu_layer(k, cell) &
                                                    - dt*state%multilayer%h_layer(k, cell)*dp_dx
               state%multilayer%hv_layer(k, cell) = state%multilayer%hv_layer(k, cell) &
                                                    - dt*state%multilayer%h_layer(k, cell)*dp_dy
            else
               state%multilayer%hu_layer(k, cell) = 0.0_wp
               state%multilayer%hv_layer(k, cell) = 0.0_wp
            end if
         end do
      end do
      !$omp end parallel do
   end subroutine ml_baroclinic_pressure_gradient_jacobian_unstr

end module rki_ml_pressure_unstr
