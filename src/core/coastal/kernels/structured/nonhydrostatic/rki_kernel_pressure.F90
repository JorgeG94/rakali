!! Non-hydrostatic pressure correction kernels
module rki_kernel_pressure
   !! Fractional-step NH pressure correction (Stelling & Zijlema, 2003).
   !! Multi-layer sigma coordinates with a CG Poisson solver.
   use rki_constants, only: wp, GRAVITY, NH_MIN_DEPTH
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: compute_depth_averaged_momentum
   public :: save_barotropic_state
   public :: compute_layer_from_barotropic
   public :: compute_vertical_velocity
   public :: compute_nh_rhs
   public :: apply_poisson_operator
   public :: compute_vert_coefficients
   public :: apply_poisson_operator_fast
   public :: apply_nh_correction
   public :: compute_rk2_average_layers
   public :: solve_nh_pressure
   public :: apply_vertical_preconditioner
   public :: compute_nh_baro_tendency
   public :: apply_baro_nh_tendency
   public :: reconcile_barotropic

contains

   subroutine compute_depth_averaged_momentum(h, u_layer, v_layer, dsig, &
                                                   hu, hv, grid, nz)
      !! Compute depth-integrated momentum from per-layer velocities
      !! hu = h * sum(dsig_k * u_k), hv = h * sum(dsig_k * v_k)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(out) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hv(grid%nx_total, grid%ny_total)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: u_avg, v_avg

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(2) private(u_avg, v_avg, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         u_avg = 0.0_wp
         v_avg = 0.0_wp
         do k = 1, nz
            u_avg = u_avg + dsig(k)*u_layer(i, j, k)
            v_avg = v_avg + dsig(k)*v_layer(i, j, k)
         end do
         hu(i, j) = h(i, j)*u_avg
         hv(i, j) = h(i, j)*v_avg
      end do
      end do
      !$omp end parallel do

   end subroutine compute_depth_averaged_momentum

   subroutine save_barotropic_state(hu, hv, hu_pre, hv_pre, grid)
      !! Save barotropic momentum before HLL update
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hu_pre(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hv_pre(grid%nx_total, grid%ny_total)

      integer :: i, j, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp parallel do collapse(2)
      do j = 1, ny
      do i = 1, nx
         hu_pre(i, j) = hu(i, j)
         hv_pre(i, j) = hv(i, j)
      end do
      end do
      !$omp end parallel do

   end subroutine save_barotropic_state

   subroutine compute_layer_from_barotropic(h, hu, hv, hu_pre, hv_pre, &
                                                 u_layer, v_layer, dsig, &
                                                 grid, nz, dt)
      !! Distribute barotropic acceleration (hu_new - hu_pre)/h uniformly
      !! to all layer velocities.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total), hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu_pre(grid%nx_total, grid%ny_total), hv_pre(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: h_val, du, dv

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(2) private(h_val, du, dv, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_val = h(i, j)
         if (h_val > NH_MIN_DEPTH) then
            du = (hu(i, j) - hu_pre(i, j))/h_val
            dv = (hv(i, j) - hv_pre(i, j))/h_val
            do k = 1, nz
               u_layer(i, j, k) = u_layer(i, j, k) + du
               v_layer(i, j, k) = v_layer(i, j, k) + dv
            end do
         else
            do k = 1, nz
               u_layer(i, j, k) = 0.0_wp
               v_layer(i, j, k) = 0.0_wp
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine compute_layer_from_barotropic

   subroutine compute_vertical_velocity(h, u_layer, v_layer, b, dsig, w, &
                                             grid, nz)
      !! Integrate continuity equation bottom-up to get vertical velocity
      !!
      !! w(:,:,1) = 0 at bottom
      !! w(:,:,k+1) = w(:,:,k) - dsig_k * [d(h*u_k)/dx + d(h*v_k)/dy]
      !! using central differences for horizontal divergence
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(out) :: w(grid%nx_total, grid%ny_total, nz + 1)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: div_hu, div_hv

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      !$omp parallel do collapse(2) private(div_hu, div_hv, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         w(i, j, 1) = 0.0_wp
         if (h(i, j) > NH_MIN_DEPTH) then
            do k = 1, nz
               div_hu = (h(i + 1, j)*u_layer(i + 1, j, k) - &
                         h(i - 1, j)*u_layer(i - 1, j, k))/(2.0_wp*dx)
               div_hv = (h(i, j + 1)*v_layer(i, j + 1, k) - &
                         h(i, j - 1)*v_layer(i, j - 1, k))/(2.0_wp*dy)
               w(i, j, k + 1) = w(i, j, k) - dsig(k)*(div_hu + div_hv)
            end do
         else
            do k = 1, nz
               w(i, j, k + 1) = 0.0_wp
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine compute_vertical_velocity

   subroutine compute_nh_rhs(h, w, dsig, rhs, nx, ny, nz, nghost, dx, dy, dt)
      !! Public only for the unit-test suite (no production module imports it).
      !! Compute RHS of NH Poisson equation from velocity divergence
      !!
      !! rhs_k = -(w_{k+1} - w_k) / (dsig_k * h * dt)
      integer, intent(in) :: nx, ny, nz, nghost
      real(wp), intent(in) :: h(nx, ny)
      real(wp), intent(in) :: w(nx, ny, nz + 1)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(out) :: rhs(nx, ny, nz)
      real(wp), intent(in) :: dx, dy, dt

      integer :: i, j, k
      real(wp) :: h_val

      !$omp parallel do collapse(2) private(h_val, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_val = h(i, j)
         if (h_val > NH_MIN_DEPTH) then
            do k = 1, nz
               rhs(i, j, k) = -(w(i, j, k + 1) - w(i, j, k)) &
                              /(dsig(k)*h_val*dt)
            end do
         else
            do k = 1, nz
               rhs(i, j, k) = 0.0_wp
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine compute_nh_rhs

   subroutine apply_poisson_operator(h, dsig, p, Ap, nx, ny, nz, nghost, dx, dy, dt)
      !! Public only for the unit-test suite (no production module imports it).
      !! Matrix-free Poisson operator: Ap = A * p
      !!
      !! Horizontal: 5-point Laplacian per layer
      !! Vertical: tridiagonal coupling between layers
      integer, intent(in) :: nx, ny, nz, nghost
      real(wp), intent(in) :: h(nx, ny)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: p(nx, ny, nz)
      real(wp), intent(out) :: Ap(nx, ny, nz)
      real(wp), intent(in) :: dx, dy, dt

      integer :: i, j, k
      real(wp) :: h_val, horiz, vert_lo, vert_hi
      real(wp) :: inv_dx2, inv_dy2, inv_dt2

      inv_dx2 = 1.0_wp/(dx*dx)
      inv_dy2 = 1.0_wp/(dy*dy)
      inv_dt2 = 1.0_wp/(dt*dt)

      !$omp parallel do collapse(2) private(h_val, horiz, vert_lo, vert_hi, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_val = h(i, j)
         if (h_val > NH_MIN_DEPTH) then
            do k = 1, nz
               ! Horizontal Laplacian
               horiz = h_val*inv_dx2*(p(i + 1, j, k) - 2.0_wp*p(i, j, k) + p(i - 1, j, k)) &
                       + h_val*inv_dy2*(p(i, j + 1, k) - 2.0_wp*p(i, j, k) + p(i, j - 1, k))

               ! Vertical coupling: second derivative in sigma
               vert_lo = 0.0_wp
               vert_hi = 0.0_wp
               if (k > 1) then
                  vert_lo = 2.0_wp*inv_dt2/(h_val*(dsig(k) + dsig(k - 1))*dsig(k))* &
                            (p(i, j, k - 1) - p(i, j, k))
               end if
               if (k < nz) then
                  vert_hi = 2.0_wp*inv_dt2/(h_val*(dsig(k) + dsig(k + 1))*dsig(k))* &
                            (p(i, j, k + 1) - p(i, j, k))
               end if

               Ap(i, j, k) = horiz + vert_lo + vert_hi
            end do
         else
            do k = 1, nz
               Ap(i, j, k) = 0.0_wp
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine apply_poisson_operator

   pure subroutine compute_vert_coefficients(dsig, nz, dt, vert_lo, vert_hi)
      !! Pre-compute vertical coupling coefficients (dsig-dependent part only)
      !!
      !! vert_lo(k) = 2 / (dt^2 * (dsig(k) + dsig(k-1)) * dsig(k))  for k=2..nz
      !! vert_hi(k) = 2 / (dt^2 * (dsig(k) + dsig(k+1)) * dsig(k))  for k=1..nz-1
      !! In the operator, multiply by 1/h(i,j) to get the full coefficient.
      integer, intent(in) :: nz
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: dt
      real(wp), intent(out) :: vert_lo(nz)
      real(wp), intent(out) :: vert_hi(nz)

      integer :: k
      real(wp) :: inv_dt2

      inv_dt2 = 1.0_wp/(dt*dt)

      vert_lo(1) = 0.0_wp
      do k = 2, nz
         vert_lo(k) = 2.0_wp*inv_dt2/((dsig(k) + dsig(k - 1))*dsig(k))
      end do

      do k = 1, nz - 1
         vert_hi(k) = 2.0_wp*inv_dt2/((dsig(k) + dsig(k + 1))*dsig(k))
      end do
      vert_hi(nz) = 0.0_wp

   end subroutine compute_vert_coefficients

   subroutine apply_poisson_operator_fast(h, p, Ap, nx, ny, nz, nghost, &
                                               inv_dx2, inv_dy2, vert_lo, vert_hi)
      !! Poisson operator (Ap = A*p) using pre-computed dsig-dependent
      !! vertical coefficients to avoid divisions; horizontal 5-point
      !! Laplacian + vertical tridiagonal coupling (coefficients / h).
      integer, intent(in) :: nx, ny, nz, nghost
      real(wp), intent(in) :: h(nx, ny)
      real(wp), intent(in) :: p(nx, ny, nz)
      real(wp), intent(out) :: Ap(nx, ny, nz)
      real(wp), intent(in) :: inv_dx2, inv_dy2
      real(wp), intent(in) :: vert_lo(nz)
      real(wp), intent(in) :: vert_hi(nz)

      integer :: i, j, k
      real(wp) :: h_val, horiz, vl, vh

      !$omp parallel do collapse(2) private(h_val, horiz, vl, vh, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_val = h(i, j)
         if (h_val > NH_MIN_DEPTH) then
            do k = 1, nz
               ! Horizontal Laplacian
               horiz = h_val*inv_dx2*(p(i + 1, j, k) - 2.0_wp*p(i, j, k) + p(i - 1, j, k)) &
                       + h_val*inv_dy2*(p(i, j + 1, k) - 2.0_wp*p(i, j, k) + p(i, j - 1, k))

               ! Vertical coupling: pre-computed coefficients / h
               vl = 0.0_wp
               vh = 0.0_wp
               if (k > 1) then
                  vl = vert_lo(k)/h_val*(p(i, j, k - 1) - p(i, j, k))
               end if
               if (k < nz) then
                  vh = vert_hi(k)/h_val*(p(i, j, k + 1) - p(i, j, k))
               end if

               Ap(i, j, k) = horiz + vl + vh
            end do
         else
            do k = 1, nz
               Ap(i, j, k) = 0.0_wp
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine apply_poisson_operator_fast

   subroutine apply_nh_correction(h, u_layer, v_layer, w, q_nh, dsig, &
                                       grid, nz, dt)
      !! Apply NH pressure gradient to correct layer velocities
      !!
      !! u_k -= dt * dq/dx
      !! v_k -= dt * dq/dy
      !! w corrected vertically
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: q_nh(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: dqdx, dqdy, h_val

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      !$omp parallel do collapse(2) private(dqdx, dqdy, h_val, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_val = h(i, j)
         if (h_val > NH_MIN_DEPTH) then
            do k = 1, nz
               ! Horizontal NH pressure gradient (central difference)
               dqdx = (q_nh(i + 1, j, k) - q_nh(i - 1, j, k))/(2.0_wp*dx)
               dqdy = (q_nh(i, j + 1, k) - q_nh(i, j - 1, k))/(2.0_wp*dy)
               u_layer(i, j, k) = u_layer(i, j, k) - dt*dqdx
               v_layer(i, j, k) = v_layer(i, j, k) - dt*dqdy
            end do

            ! Vertical velocity correction
            ! w_1 = 0 (bottom), already set
            do k = 1, nz
               if (k < nz) then
                  w(i, j, k + 1) = w(i, j, k + 1) + &
                                   dt*(q_nh(i, j, k + 1) - q_nh(i, j, k)) &
                                   *2.0_wp/(h_val*(dsig(k) + dsig(k + 1)))
               else
                  ! Top interface: q above surface = 0
                  w(i, j, k + 1) = w(i, j, k + 1) - &
                                   dt*q_nh(i, j, k)*2.0_wp/(h_val*dsig(k))
               end if
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine apply_nh_correction

   subroutine compute_rk2_average_layers(u_layer, v_layer, w, &
                                              u_layer0, v_layer0, w0, &
                                              grid, nz)
      !! RK2 averaging for layer fields: 0.5 * (old + new)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: u_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer0(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: w0(grid%nx_total, grid%ny_total, nz + 1)

      integer :: i, j, k, nx, ny, nghost

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         u_layer(i, j, k) = 0.5_wp*(u_layer0(i, j, k) + u_layer(i, j, k))
         v_layer(i, j, k) = 0.5_wp*(v_layer0(i, j, k) + v_layer(i, j, k))
      end do
      end do
      end do
      !$omp end parallel do

      !$omp parallel do collapse(3)
      do k = 1, nz + 1
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         w(i, j, k) = 0.5_wp*(w0(i, j, k) + w(i, j, k))
      end do
      end do
      end do
      !$omp end parallel do

   end subroutine compute_rk2_average_layers

   subroutine apply_vertical_preconditioner(h, dsig, r, z, nx, ny, nz, nghost, dx, dy, dt)
      !! Public only for the unit-test suite (no production module imports it).
      !! Vertical-line preconditioner: solve M*z = r per (i,j) column
      !!
      !! M is the vertical tridiagonal part of the Poisson operator plus
      !! the horizontal diagonal contribution. Thomas algorithm per column
      !! (serial in k, parallel over i,j via do concurrent).
      integer, intent(in) :: nx, ny, nz, nghost
      real(wp), intent(in) :: h(nx, ny)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: r(nx, ny, nz)
      real(wp), intent(out) :: z(nx, ny, nz)
      real(wp), intent(in) :: dx, dy, dt

      integer, parameter :: NZ_MAX = 64
      integer :: i, j, k
      real(wp) :: inv_dx2, inv_dy2, inv_dt2
      real(wp) :: h_val, h_diag, a_k, c_k, d_k
      real(wp) :: cp(NZ_MAX), dp(NZ_MAX)

      inv_dx2 = 1.0_wp/(dx*dx)
      inv_dy2 = 1.0_wp/(dy*dy)
      inv_dt2 = 1.0_wp/(dt*dt)

      ! Horizontal diagonal contribution (same for all layers)
      h_diag = 2.0_wp*inv_dx2 + 2.0_wp*inv_dy2

      !$omp parallel do collapse(2) private(h_val, a_k, c_k, d_k, cp, dp, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_val = h(i, j)
         if (h_val > NH_MIN_DEPTH) then
            ! Forward sweep (Thomas algorithm)
            ! k = 1: a_k = 0 (no layer below)
            if (nz > 1) then
               c_k = 2.0_wp*inv_dt2/(h_val*(dsig(1) + dsig(2))*dsig(1))
            else
               c_k = 0.0_wp
            end if
            d_k = -(c_k + h_val*h_diag)
            cp(1) = c_k/d_k
            dp(1) = r(i, j, 1)/d_k

            do k = 2, nz
               ! Sub-diagonal: coupling to k-1
               a_k = 2.0_wp*inv_dt2/(h_val*(dsig(k) + dsig(k - 1))*dsig(k))
               ! Super-diagonal: coupling to k+1
               if (k < nz) then
                  c_k = 2.0_wp*inv_dt2/(h_val*(dsig(k) + dsig(k + 1))*dsig(k))
               else
                  c_k = 0.0_wp
               end if
               ! Modified diagonal and RHS
               d_k = -(a_k + c_k + h_val*h_diag) - a_k*cp(k - 1)
               cp(k) = c_k/d_k
               dp(k) = (r(i, j, k) - a_k*dp(k - 1))/d_k
            end do

            ! Back substitution
            z(i, j, nz) = dp(nz)
            do k = nz - 1, 1, -1
               z(i, j, k) = dp(k) - cp(k)*z(i, j, k + 1)
            end do
         else
            do k = 1, nz
               z(i, j, k) = 0.0_wp
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine apply_vertical_preconditioner

   subroutine solve_nh_pressure(h, u_layer, v_layer, w, q_nh, dsig, &
                                cg_r, cg_p, cg_Ap, cg_z, &
                                grid, nz, dt, &
                                max_iter, tol, iter_out, &
                                use_preconditioner)
      !! Conjugate Gradient solver for the NH Poisson equation
      !!
      !! Warm-started from previous q_nh.
      !! Optionally uses vertical-line preconditioner (Thomas algorithm per column).
      !! Solves A * q = rhs where A is the Poisson operator and
      !! rhs comes from the velocity divergence.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(inout) :: q_nh(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(inout) :: cg_r(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: cg_p(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: cg_Ap(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: cg_z(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dt
      integer, intent(in) :: max_iter
      real(wp), intent(in) :: tol
      integer, intent(out) :: iter_out
      logical, intent(in), optional :: use_preconditioner

      logical :: use_pc
      integer :: i, j, k, iter, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: rz, rz_new, rr, pAp, alpha, beta
      real(wp) :: tol2, inv_dx2, inv_dy2
      real(wp) :: vlo(nz), vhi(nz)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      use_pc = .false.
      if (present(use_preconditioner)) use_pc = use_preconditioner

      tol2 = tol*tol
      iter_out = 0

      ! Pre-compute vertical coefficients and horizontal inverses once
      call compute_vert_coefficients(dsig, nz, dt, vlo, vhi)
      inv_dx2 = 1.0_wp/(dx*dx)
      inv_dy2 = 1.0_wp/(dy*dy)

      ! Compute RHS into cg_r
      call compute_nh_rhs(h, w, dsig, cg_r, nx, ny, nz, nghost, dx, dy, dt)

      ! r = rhs - A*q_nh (warm start: q_nh retains previous solution)
      call apply_poisson_operator_fast(h, q_nh, cg_Ap, nx, ny, nz, nghost, &
                                       inv_dx2, inv_dy2, vlo, vhi)

      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         cg_r(i, j, k) = cg_r(i, j, k) - cg_Ap(i, j, k)
      end do
      end do
      end do
      !$omp end parallel do

      ! Check if already converged (unpreconditioned residual)
      rr = 0.0_wp
      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         rr = rr + cg_r(i, j, k)*cg_r(i, j, k)
      end do
      end do
      end do
      !$omp end parallel do

      if (rr < tol2) return

      if (use_pc) then
         ! PCG: z = M^{-1} r, p = z
         call apply_vertical_preconditioner(h, dsig, cg_r, cg_z, nx, ny, nz, nghost, dx, dy, dt)
         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            cg_p(i, j, k) = cg_z(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do
         ! rz = dot(r, z)
         rz = 0.0_wp
         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            rz = rz + cg_r(i, j, k)*cg_z(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do
      else
         ! Plain CG: p = r, rz = rr
         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            cg_p(i, j, k) = cg_r(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do
         rz = rr
      end if

      do iter = 1, max_iter
         ! Ap = A * p
         call apply_poisson_operator_fast(h, cg_p, cg_Ap, nx, ny, nz, nghost, &
                                          inv_dx2, inv_dy2, vlo, vhi)

         ! pAp = dot(p, Ap)
         pAp = 0.0_wp
         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            pAp = pAp + cg_p(i, j, k)*cg_Ap(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do

         if (abs(pAp) < 1.0e-30_wp) exit

         alpha = rz/pAp

         ! q_nh += alpha * p
         ! r -= alpha * Ap
         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            q_nh(i, j, k) = q_nh(i, j, k) + alpha*cg_p(i, j, k)
            cg_r(i, j, k) = cg_r(i, j, k) - alpha*cg_Ap(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do

         ! Convergence check on unpreconditioned residual
         rr = 0.0_wp
         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            rr = rr + cg_r(i, j, k)*cg_r(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do

         iter_out = iter

         if (rr < tol2) exit

         if (use_pc) then
            ! PCG: z = M^{-1} r, rz_new = dot(r, z)
            call apply_vertical_preconditioner(h, dsig, cg_r, cg_z, nx, ny, nz, nghost, dx, dy, dt)
            rz_new = 0.0_wp
            !$omp parallel do collapse(3)
            do k = 1, nz
            do j = nghost + 1, ny - nghost
            do i = nghost + 1, nx - nghost
               rz_new = rz_new + cg_r(i, j, k)*cg_z(i, j, k)
            end do
            end do
            end do
            !$omp end parallel do
            beta = rz_new/rz
            rz = rz_new
            ! p = z + beta * p
            !$omp parallel do collapse(3)
            do k = 1, nz
            do j = nghost + 1, ny - nghost
            do i = nghost + 1, nx - nghost
               cg_p(i, j, k) = cg_z(i, j, k) + beta*cg_p(i, j, k)
            end do
            end do
            end do
            !$omp end parallel do
         else
            ! Plain CG: rz_new = rr
            beta = rr/rz
            rz = rr
            ! p = r + beta * p
            !$omp parallel do collapse(3)
            do k = 1, nz
            do j = nghost + 1, ny - nghost
            do i = nghost + 1, nx - nghost
               cg_p(i, j, k) = cg_r(i, j, k) + beta*cg_p(i, j, k)
            end do
            end do
            end do
            !$omp end parallel do
         end if
      end do

   end subroutine solve_nh_pressure

   subroutine compute_nh_baro_tendency(h, q_nh, dsig, tend_x, tend_y, &
                                            grid, nz)
      !! Compute depth-averaged NH pressure gradient for barotropic forcing
      !!
      !! tend_x(i,j) = -sum_k(dsig_k * dq_nh_k/dx)
      !! tend_y(i,j) = -sum_k(dsig_k * dq_nh_k/dy)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: q_nh(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(out) :: tend_x(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: tend_y(grid%nx_total, grid%ny_total)

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: tx, ty, dqdx, dqdy, inv_2dx, inv_2dy

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      inv_2dx = 1.0_wp/(2.0_wp*dx)
      inv_2dy = 1.0_wp/(2.0_wp*dy)

      !$omp parallel do collapse(2) private(tx, ty, dqdx, dqdy, k)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         tx = 0.0_wp
         ty = 0.0_wp
         if (h(i, j) > NH_MIN_DEPTH) then
            do k = 1, nz
               dqdx = (q_nh(i + 1, j, k) - q_nh(i - 1, j, k))*inv_2dx
               dqdy = (q_nh(i, j + 1, k) - q_nh(i, j - 1, k))*inv_2dy
               tx = tx - dsig(k)*dqdx
               ty = ty - dsig(k)*dqdy
            end do
         end if
         tend_x(i, j) = tx
         tend_y(i, j) = ty
      end do
      end do
      !$omp end parallel do

   end subroutine compute_nh_baro_tendency

   subroutine apply_baro_nh_tendency(hu, hv, h, tend_x, tend_y, &
                                          grid, dt)
      !! Add NH barotropic tendency as forcing during barotropic subcycling
      !!
      !! hu += dt * h * tend_x
      !! hv += dt * h * tend_y
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: tend_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: tend_y(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: dt

      integer :: i, j, nx, ny, nghost

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(2)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         hu(i, j) = hu(i, j) + dt*h(i, j)*tend_x(i, j)
         hv(i, j) = hv(i, j) + dt*h(i, j)*tend_y(i, j)
      end do
      end do
      !$omp end parallel do

   end subroutine apply_baro_nh_tendency

   subroutine reconcile_barotropic(h, hu, hv, u_layer, v_layer, dsig, &
                                   grid, nz)
      !! Recompute barotropic momentum from corrected layer velocities
      !!
      !! hu = h * sum(dsig_k * u_layer_k)
      !! hv = h * sum(dsig_k * v_layer_k)
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(out) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dsig(nz)

      ! This is identical to compute_depth_averaged_momentum
      call compute_depth_averaged_momentum(h, u_layer, v_layer, dsig, &
                                           hu, hv, grid, nz)

   end subroutine reconcile_barotropic

end module rki_kernel_pressure
