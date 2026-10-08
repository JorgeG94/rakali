!! Baroclinic pressure gradient (BPG) for the structured multilayer solver.
module rki_ml_pressure
   !! Two BPG implementations:
   !!   - `ml_baroclinic_pressure_gradient` (sigma method): differences
   !!     internal pressure at the same sigma level. Cheap; PGE ~ dH/dx
   !!     on steep bathymetry.
   !!   - `ml_baroclinic_pressure_gradient_jacobian` (Shchepetkin &
   !!     McWilliams 2003): interpolates neighbour-column pressure to the
   !!     centre column's z-mid level. Eliminates the leading sigma PGE;
   !!     ~3x more expensive per cell.
   !! Both use the density anomaly `rho - rho_ref(z)` to keep the
   !! integrand small. `rho_ref_at` and `all_neighbours_wet` are private
   !! same-module helpers because the BPG kernels call them inside
   !! `do concurrent`, where NVHPC stdpar can't lower cross-module calls.
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
#else
   use rki_constants, only: BPG_NZ_MAX, wp, GRAVITY, DRY_TOLERANCE
#endif
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH, VANISHING_LAYER_TOL
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: BPG_NZ_MAX = 128
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=128).
#endif

   public :: ml_baroclinic_pressure_gradient
   public :: ml_baroclinic_pressure_gradient_jacobian

contains

   pure function rho_ref_at(z, rho_ref, dz, nz_ref, rho_0) result(rho)
      !$omp declare target
      !! Linear interpolation of the reference-density profile at depth z.
      !! Out-of-range falls back to the nearest endpoint. When the profile
      !! is unset (`dz == 0`) returns rho_0 so the kernel reduces to the
      !! pre-B2 anomaly.
      real(wp), intent(in) :: z, dz, rho_0
      integer, intent(in) :: nz_ref
      real(wp), intent(in) :: rho_ref(0:nz_ref)
      real(wp) :: rho, zr, frac
      integer :: i

      if (dz <= 0.0_wp .or. nz_ref < 1) then
         rho = rho_0
         return
      end if

      zr = max(0.0_wp, min(z, dz*real(nz_ref, wp)))
      i = int(zr/dz)
      if (i >= nz_ref) then
         rho = rho_ref(nz_ref)
      else
         frac = zr/dz - real(i, wp)
         rho = (1.0_wp - frac)*rho_ref(i) + frac*rho_ref(i + 1)
      end if
   end function rho_ref_at
   pure subroutine ml_baroclinic_pressure_gradient(h_layer, rho_layer, &
                                                   hu_layer, hv_layer, &
                                                   rho_0, rho_ref, rho_ref_dz, nz_ref, &
                                                   grid, nz, dt)
      !! Apply baroclinic pressure gradient as explicit source term —
      !! only the density-driven part; the barotropic PGF (g·deta/dx) is
      !! handled by the HLL step.
      !!
      !! Hydrostatic internal pressure at the centre of layer k:
      !!   p_k = g * [ 0.5*rho_k*h_k + sum_{j=k+1}^{nz} rho_j*h_j ]
      !! Baroclinic PGF for layer k:
      !!   BPG_k^x = -(g/rho_0) * h_k * dp_internal/dx
      !! Accumulated top-down per column, then central-differenced.
      !!
      !! OBC: the stencil reads ghost columns. OBC_WALL ghosts keep their
      !! IC state ⇒ zero-gradient PGF at the wall. For OPEN/TIDAL/CLAMPED,
      !! callers wanting PGF to see the OBC value must populate the ghost
      !! `h_layer` first; per-edge dispatch is a deferred v2.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz, nz_ref
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0
      real(wp), intent(in) :: rho_ref(0:nz_ref)
      real(wp), intent(in) :: rho_ref_dz
      real(wp), intent(in) :: dt

      ! Internal pressure at layer centres.
      ! For small nz (2-5) these fit comfortably in registers/stack per thread.
      ! Workspace capped by BPG_NZ_MAX (see module-level comment);
      ! validate_config rejects namelists that exceed this at startup.
      real(wp) :: p_int(BPG_NZ_MAX)

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: inv_2dx, inv_2dy, inv_rho0
      real(wp) :: acc, h_col, h_col_im1, h_col_ip1, h_col_jm1, h_col_jp1
      real(wp) :: dp_dx, dp_dy
      real(wp) :: z_top, z_bot, z_mid_lay
      real(wp) :: p_int_ip1(BPG_NZ_MAX), p_int_im1(BPG_NZ_MAX)
      real(wp) :: p_int_jp1(BPG_NZ_MAX), p_int_jm1(BPG_NZ_MAX)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      inv_2dx = 0.5_wp/dx
      inv_2dy = 0.5_wp/dy
      inv_rho0 = 1.0_wp/rho_0

      ! Process each interior cell
      ! We need the internal pressure at (i,j), (i±1,j), (i,j±1) for
      ! central differences. Compute all five columns per cell.
      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, kk, acc, p_int, p_int_ip1, p_int_im1, p_int_jp1, p_int_jm1, &
               dp_dx, dp_dy, h_col, h_col_im1, h_col_ip1, h_col_jm1, h_col_jp1, &
               z_top, z_bot, z_mid_lay)

         ! Compute total column depth for this cell and its neighbours
         h_col = 0.0_wp
         h_col_im1 = 0.0_wp
         h_col_ip1 = 0.0_wp
         h_col_jm1 = 0.0_wp
         h_col_jp1 = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
            h_col_im1 = h_col_im1 + h_layer(i - 1, j, kk)
            h_col_ip1 = h_col_ip1 + h_layer(i + 1, j, kk)
            h_col_jm1 = h_col_jm1 + h_layer(i, j - 1, kk)
            h_col_jp1 = h_col_jp1 + h_layer(i, j + 1, kk)
         end do

         ! Skip cells near dry columns — BPG central differences across
         ! wet/dry fronts produce spurious pressure gradients
         if (.not. all_neighbours_wet(h_col, h_col_ip1, h_col_im1, &
                                      h_col_jp1, h_col_jm1)) cycle

         ! Compute baroclinic pressure anomaly using
         !   rho' = rho - rho_ref(z)
         ! where rho_ref(z) is the horizontally-uniform reference profile
         ! (Phase B2). When unset (rho_ref_dz==0) the helper returns
         ! rho_0, recovering the pre-B2 anomaly.

         ! Column (i,j)
         acc = 0.0_wp
         z_top = 0.0_wp
         do k = nz, 1, -1
            z_bot = z_top + h_layer(i, j, k)
            z_mid_lay = 0.5_wp*(z_top + z_bot)
            acc = acc + (rho_layer(i, j, k) &
                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                  *h_layer(i, j, k)
            p_int(k) = acc - 0.5_wp*(rho_layer(i, j, k) &
                                     - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                       *h_layer(i, j, k)
            z_top = z_bot
         end do

         ! Column (i+1,j)
         acc = 0.0_wp
         z_top = 0.0_wp
         do k = nz, 1, -1
            z_bot = z_top + h_layer(i + 1, j, k)
            z_mid_lay = 0.5_wp*(z_top + z_bot)
            acc = acc + (rho_layer(i + 1, j, k) &
                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                  *h_layer(i + 1, j, k)
            p_int_ip1(k) = acc - 0.5_wp*(rho_layer(i + 1, j, k) &
                                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                           *h_layer(i + 1, j, k)
            z_top = z_bot
         end do

         ! Column (i-1,j)
         acc = 0.0_wp
         z_top = 0.0_wp
         do k = nz, 1, -1
            z_bot = z_top + h_layer(i - 1, j, k)
            z_mid_lay = 0.5_wp*(z_top + z_bot)
            acc = acc + (rho_layer(i - 1, j, k) &
                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                  *h_layer(i - 1, j, k)
            p_int_im1(k) = acc - 0.5_wp*(rho_layer(i - 1, j, k) &
                                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                           *h_layer(i - 1, j, k)
            z_top = z_bot
         end do

         ! Column (i,j+1)
         acc = 0.0_wp
         z_top = 0.0_wp
         do k = nz, 1, -1
            z_bot = z_top + h_layer(i, j + 1, k)
            z_mid_lay = 0.5_wp*(z_top + z_bot)
            acc = acc + (rho_layer(i, j + 1, k) &
                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                  *h_layer(i, j + 1, k)
            p_int_jp1(k) = acc - 0.5_wp*(rho_layer(i, j + 1, k) &
                                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                           *h_layer(i, j + 1, k)
            z_top = z_bot
         end do

         ! Column (i,j-1)
         acc = 0.0_wp
         z_top = 0.0_wp
         do k = nz, 1, -1
            z_bot = z_top + h_layer(i, j - 1, k)
            z_mid_lay = 0.5_wp*(z_top + z_bot)
            acc = acc + (rho_layer(i, j - 1, k) &
                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                  *h_layer(i, j - 1, k)
            p_int_jm1(k) = acc - 0.5_wp*(rho_layer(i, j - 1, k) &
                                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                           *h_layer(i, j - 1, k)
            z_top = z_bot
         end do

         ! Apply pressure gradient to each layer.
         ! Cap the BPG acceleration to prevent the sigma-coordinate PGE
         ! from creating velocities exceeding the local wave speed.
         do k = 1, nz
            dp_dx = (p_int_ip1(k) - p_int_im1(k))*inv_2dx
            dp_dy = (p_int_jp1(k) - p_int_jm1(k))*inv_2dy

            dp_dx = GRAVITY*inv_rho0*dp_dx
            dp_dy = GRAVITY*inv_rho0*dp_dy

            ! Limit: BPG acceleration <= sqrt(g'*H)/dt per step
            acc = abs(rho_layer(i, j, 1) - rho_layer(i, j, nz))*inv_rho0
            acc = sqrt(max(acc, 1.0e-6_wp)*GRAVITY &
                       *max(h_col, ML_MIN_DEPTH)) &
                  /max(dt, 1.0e-10_wp)
            if (abs(dp_dx) > acc) dp_dx = sign(acc, dp_dx)
            if (abs(dp_dy) > acc) dp_dy = sign(acc, dp_dy)

            if (h_layer(i, j, k) > VANISHING_LAYER_TOL) then
               hu_layer(i, j, k) = hu_layer(i, j, k) - dt*h_layer(i, j, k)*dp_dx
               hv_layer(i, j, k) = hv_layer(i, j, k) - dt*h_layer(i, j, k)*dp_dy
            else
               ! VCOORD_ZSTAR_FULL invariant: vanishing layers carry no
               ! momentum.  Keeps hu/h velocity bounded in vanishing rows
               ! where h_layer ≈ h_min and any residual momentum would
               ! show up as an unphysical velocity.
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            end if
         end do
      end do

   end subroutine ml_baroclinic_pressure_gradient
   pure subroutine ml_baroclinic_pressure_gradient_jacobian(h_layer, rho_layer, &
                                                            hu_layer, hv_layer, &
                                                            rho_0, rho_ref, rho_ref_dz, nz_ref, &
                                                            grid, nz, dt)
      !! Baroclinic pressure gradient via the pressure Jacobian method
      !! (Shchepetkin & McWilliams 2003): the pressure gradient along
      !! constant-z surfaces, interpolating the neighbour column's pressure
      !! to the centre column's physical depth. Uniform density ⇒ zero BPG.
      !!
      !! Per cell (i,j): z_mid(k) = layer-k centre depth from surface;
      !! p(z) = integral of (rho-rho_ref)*g from surface; interpolate
      !! neighbour p(z) to z_mid(k); dp/dx = [p_nbr - p_centre]/dx.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz, nz_ref
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0
      real(wp), intent(in) :: rho_ref(0:nz_ref)
      real(wp), intent(in) :: rho_ref_dz
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, kn, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: inv_2dx, inv_2dy, inv_rho0
      real(wp) :: h_col, h_col_im1, h_col_ip1, h_col_jm1, h_col_jp1
      ! Centre column: depth of layer interfaces and pressure at interfaces.
      ! Workspace capped by BPG_NZ_MAX (see module-level comment);
      ! validate_config rejects namelists that exceed this at startup.
      real(wp) :: z_if(BPG_NZ_MAX + 1), p_if(BPG_NZ_MAX + 1)
      ! One neighbour column per direction. Built once per cell and
      ! reused for every layer's search/interpolate. Hoisting these out
      ! of the inner `k` loop is what cuts BPG runtime in half -- the
      ! columns are independent of `z_mid`, only the queries into them
      ! depend on `k`.
      real(wp) :: zn_if_ip1(BPG_NZ_MAX + 1), pn_if_ip1(BPG_NZ_MAX + 1)
      real(wp) :: zn_if_im1(BPG_NZ_MAX + 1), pn_if_im1(BPG_NZ_MAX + 1)
      real(wp) :: zn_if_jp1(BPG_NZ_MAX + 1), pn_if_jp1(BPG_NZ_MAX + 1)
      real(wp) :: zn_if_jm1(BPG_NZ_MAX + 1), pn_if_jm1(BPG_NZ_MAX + 1)
      real(wp) :: z_mid, z_mid_lay, p_here, p_nbr, dp_dx, dp_dy, frac, acc

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      inv_2dx = 0.5_wp/dx
      inv_2dy = 0.5_wp/dy
      inv_rho0 = 1.0_wp/rho_0

      do concurrent(j=nghost + 1:ny - nghost, i=nghost + 1:nx - nghost) &
         local(k, kk, kn, z_if, p_if, &
               zn_if_ip1, pn_if_ip1, zn_if_im1, pn_if_im1, &
               zn_if_jp1, pn_if_jp1, zn_if_jm1, pn_if_jm1, &
               z_mid, z_mid_lay, p_here, p_nbr, &
               dp_dx, dp_dy, frac, acc, &
               h_col, h_col_im1, h_col_ip1, h_col_jm1, h_col_jp1)

         ! Compute total column depth for this cell and its neighbours
         h_col = 0.0_wp
         h_col_im1 = 0.0_wp
         h_col_ip1 = 0.0_wp
         h_col_jm1 = 0.0_wp
         h_col_jp1 = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
            h_col_im1 = h_col_im1 + h_layer(i - 1, j, kk)
            h_col_ip1 = h_col_ip1 + h_layer(i + 1, j, kk)
            h_col_jm1 = h_col_jm1 + h_layer(i, j - 1, kk)
            h_col_jp1 = h_col_jp1 + h_layer(i, j + 1, kk)
         end do

         ! Skip cells near dry columns
         if (.not. all_neighbours_wet(h_col, h_col_ip1, h_col_im1, &
                                      h_col_jp1, h_col_jm1)) cycle

         ! --- Centre column: interface depths (from surface) and pressure ---
         ! z_if(nz+1) = 0 (surface), z_if(k) = depth of interface below layer k.
         ! Pressure uses the density anomaly `(rho - rho_ref(z))` (Phase B2);
         ! when rho_ref is unset (rho_ref_dz==0) the helper returns rho_0.
         z_if(nz + 1) = 0.0_wp
         p_if(nz + 1) = 0.0_wp
         do k = nz, 1, -1
            z_if(k) = z_if(k + 1) + h_layer(i, j, k)
            z_mid_lay = 0.5_wp*(z_if(k) + z_if(k + 1))
            p_if(k) = p_if(k + 1) &
                      + (rho_layer(i, j, k) &
                         - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                      *h_layer(i, j, k)
         end do

         ! --- Build neighbour columns ONCE per cell (was: nz times) ---
         ! Each (zn_if_*, pn_if_*) is independent of `k`; the per-layer
         ! work in the inner loop is just a search/interpolate against
         ! these prebuilt columns.

         ! Column (i+1, j)
         zn_if_ip1(nz + 1) = 0.0_wp
         pn_if_ip1(nz + 1) = 0.0_wp
         do kn = nz, 1, -1
            zn_if_ip1(kn) = zn_if_ip1(kn + 1) + h_layer(i + 1, j, kn)
            z_mid_lay = 0.5_wp*(zn_if_ip1(kn) + zn_if_ip1(kn + 1))
            pn_if_ip1(kn) = pn_if_ip1(kn + 1) &
                            + (rho_layer(i + 1, j, kn) &
                               - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                            *h_layer(i + 1, j, kn)
         end do

         ! Column (i-1, j)
         zn_if_im1(nz + 1) = 0.0_wp
         pn_if_im1(nz + 1) = 0.0_wp
         do kn = nz, 1, -1
            zn_if_im1(kn) = zn_if_im1(kn + 1) + h_layer(i - 1, j, kn)
            z_mid_lay = 0.5_wp*(zn_if_im1(kn) + zn_if_im1(kn + 1))
            pn_if_im1(kn) = pn_if_im1(kn + 1) &
                            + (rho_layer(i - 1, j, kn) &
                               - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                            *h_layer(i - 1, j, kn)
         end do

         ! Column (i, j+1)
         zn_if_jp1(nz + 1) = 0.0_wp
         pn_if_jp1(nz + 1) = 0.0_wp
         do kn = nz, 1, -1
            zn_if_jp1(kn) = zn_if_jp1(kn + 1) + h_layer(i, j + 1, kn)
            z_mid_lay = 0.5_wp*(zn_if_jp1(kn) + zn_if_jp1(kn + 1))
            pn_if_jp1(kn) = pn_if_jp1(kn + 1) &
                            + (rho_layer(i, j + 1, kn) &
                               - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                            *h_layer(i, j + 1, kn)
         end do

         ! Column (i, j-1)
         zn_if_jm1(nz + 1) = 0.0_wp
         pn_if_jm1(nz + 1) = 0.0_wp
         do kn = nz, 1, -1
            zn_if_jm1(kn) = zn_if_jm1(kn + 1) + h_layer(i, j - 1, kn)
            z_mid_lay = 0.5_wp*(zn_if_jm1(kn) + zn_if_jm1(kn + 1))
            pn_if_jm1(kn) = pn_if_jm1(kn + 1) &
                            + (rho_layer(i, j - 1, kn) &
                               - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                            *h_layer(i, j - 1, kn)
         end do

         do k = 1, nz
            z_mid = 0.5_wp*(z_if(k) + z_if(k + 1))
            ! Centre-column pressure at z_mid (used as the shelf-break wall
            ! fallback whenever a neighbour column is shallower than z_mid).
            p_here = p_if(k + 1) &
                     + 0.5_wp*(rho_layer(i, j, k) &
                               - rho_ref_at(z_mid, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                     *h_layer(i, j, k)

            ! --- x-direction: search & interpolate p at z_mid ---

            ! Column (i+1, j)
            p_nbr = p_here  ! default: z_mid below neighbour bed -> wall
            do kn = 1, nz
               if (z_mid >= zn_if_ip1(kn + 1) .and. z_mid <= zn_if_ip1(kn)) then
                  if (h_layer(i + 1, j, kn) > DRY_TOLERANCE) then
                     frac = (z_mid - zn_if_ip1(kn + 1))/h_layer(i + 1, j, kn)
                  else
                     frac = 0.0_wp
                  end if
                  z_mid_lay = 0.5_wp*(zn_if_ip1(kn) + zn_if_ip1(kn + 1))
                  p_nbr = pn_if_ip1(kn + 1) + frac &
                          *(rho_layer(i + 1, j, kn) &
                            - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                          *h_layer(i + 1, j, kn)
                  exit
               end if
            end do
            dp_dx = p_nbr  ! save p at z_mid in i+1 column

            ! Column (i-1, j)
            p_nbr = p_here
            do kn = 1, nz
               if (z_mid >= zn_if_im1(kn + 1) .and. z_mid <= zn_if_im1(kn)) then
                  if (h_layer(i - 1, j, kn) > DRY_TOLERANCE) then
                     frac = (z_mid - zn_if_im1(kn + 1))/h_layer(i - 1, j, kn)
                  else
                     frac = 0.0_wp
                  end if
                  z_mid_lay = 0.5_wp*(zn_if_im1(kn) + zn_if_im1(kn + 1))
                  p_nbr = pn_if_im1(kn + 1) + frac &
                          *(rho_layer(i - 1, j, kn) &
                            - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                          *h_layer(i - 1, j, kn)
                  exit
               end if
            end do
            dp_dx = (dp_dx - p_nbr)*inv_2dx  ! central: (p_ip1 - p_im1)/(2*dx)

            ! --- y-direction ---

            ! Column (i, j+1)
            p_nbr = p_here
            do kn = 1, nz
               if (z_mid >= zn_if_jp1(kn + 1) .and. z_mid <= zn_if_jp1(kn)) then
                  if (h_layer(i, j + 1, kn) > DRY_TOLERANCE) then
                     frac = (z_mid - zn_if_jp1(kn + 1))/h_layer(i, j + 1, kn)
                  else
                     frac = 0.0_wp
                  end if
                  z_mid_lay = 0.5_wp*(zn_if_jp1(kn) + zn_if_jp1(kn + 1))
                  p_nbr = pn_if_jp1(kn + 1) + frac &
                          *(rho_layer(i, j + 1, kn) &
                            - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                          *h_layer(i, j + 1, kn)
                  exit
               end if
            end do
            dp_dy = p_nbr  ! save p at z_mid in j+1 column

            ! Column (i, j-1)
            p_nbr = p_here
            do kn = 1, nz
               if (z_mid >= zn_if_jm1(kn + 1) .and. z_mid <= zn_if_jm1(kn)) then
                  if (h_layer(i, j - 1, kn) > DRY_TOLERANCE) then
                     frac = (z_mid - zn_if_jm1(kn + 1))/h_layer(i, j - 1, kn)
                  else
                     frac = 0.0_wp
                  end if
                  z_mid_lay = 0.5_wp*(zn_if_jm1(kn) + zn_if_jm1(kn + 1))
                  p_nbr = pn_if_jm1(kn + 1) + frac &
                          *(rho_layer(i, j - 1, kn) &
                            - rho_ref_at(z_mid_lay, rho_ref, rho_ref_dz, nz_ref, rho_0)) &
                          *h_layer(i, j - 1, kn)
                  exit
               end if
            end do
            dp_dy = (dp_dy - p_nbr)*inv_2dy  ! central: (p_jp1 - p_jm1)/(2*dy)

            ! Apply
            dp_dx = GRAVITY*inv_rho0*dp_dx
            dp_dy = GRAVITY*inv_rho0*dp_dy

            ! Limit: BPG acceleration <= sqrt(g'*H)/dt per step
            ! Use reduced gravity g' = g*|drho_max|/rho_0 from this column
            ! so the cap scales with the actual baroclinic forcing, not the
            ! barotropic wave speed (which is orders of magnitude too loose).
            acc = abs(rho_layer(i, j, 1) - rho_layer(i, j, nz))*inv_rho0
            acc = sqrt(max(acc, 1.0e-6_wp)*GRAVITY &
                       *max(h_col, ML_MIN_DEPTH)) &
                  /max(dt, 1.0e-10_wp)
            if (abs(dp_dx) > acc) dp_dx = sign(acc, dp_dx)
            if (abs(dp_dy) > acc) dp_dy = sign(acc, dp_dy)

            if (h_layer(i, j, k) > VANISHING_LAYER_TOL) then
               hu_layer(i, j, k) = hu_layer(i, j, k) - dt*h_layer(i, j, k)*dp_dx
               hv_layer(i, j, k) = hv_layer(i, j, k) - dt*h_layer(i, j, k)*dp_dy
            else
               ! VCOORD_ZSTAR_FULL invariant: vanishing layers carry no momentum
               hu_layer(i, j, k) = 0.0_wp
               hv_layer(i, j, k) = 0.0_wp
            end if
         end do
      end do

   end subroutine ml_baroclinic_pressure_gradient_jacobian
   pure logical function all_neighbours_wet(h_col, h_col_ip1, h_col_im1, &
                                            h_col_jp1, h_col_jm1) result(all_wet)
      !$omp declare target
      !! Returns `.true.` if the centre cell and its four face-neighbours all
      !! have a total column depth at or above `ML_MIN_DEPTH`.
      !!
      !! Used to guard central-difference operators (BPG, divergence, etc.)
      !! that produce spurious gradients when any of the five stencil
      !! columns is dry.
      real(wp), intent(in) :: h_col     !! column depth at (i,j)
      real(wp), intent(in) :: h_col_ip1  !! column depth at (i+1,j)
      real(wp), intent(in) :: h_col_im1  !! column depth at (i-1,j)
      real(wp), intent(in) :: h_col_jp1  !! column depth at (i,j+1)
      real(wp), intent(in) :: h_col_jm1  !! column depth at (i,j-1)
      all_wet = h_col >= ML_MIN_DEPTH .and. &
                h_col_ip1 >= ML_MIN_DEPTH .and. &
                h_col_im1 >= ML_MIN_DEPTH .and. &
                h_col_jp1 >= ML_MIN_DEPTH .and. &
                h_col_jm1 >= ML_MIN_DEPTH
   end function all_neighbours_wet

end module rki_ml_pressure
