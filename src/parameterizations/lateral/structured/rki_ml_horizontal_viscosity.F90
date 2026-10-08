!! Horizontal viscosity for layer momentum.
module rki_ml_horizontal_viscosity
   !! Smagorinsky (1963) eddy-viscosity + biharmonic closures on the
   !! structured multilayer grid.  Owns persistent per-step
   !! momentum-snapshot workspaces, torn down via
   !! `ml_horizontal_viscosity_cleanup` at simulation end.
   use rki_constants, only: wp, THIN_LAYER_THRESHOLD
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH
   implicit none
   private

   public :: ml_horizontal_viscosity_smag
   public :: ml_horizontal_viscosity_biharmonic
   public :: ml_horizontal_viscosity_cleanup

   ! Persistent per-kernel workspaces: lazy-allocated, sized (nx,ny,nz),
   ! reused every call to avoid per-call allocate/deallocate GPU stalls.
   real(wp), allocatable :: horvisc_hu_pre(:, :, :)
   real(wp), allocatable :: horvisc_hv_pre(:, :, :)
   ! Intermediate velocity Laplacian for the biharmonic first del^2 pass.
   real(wp), allocatable :: bih_lap_u(:, :, :)
   real(wp), allocatable :: bih_lap_v(:, :, :)

contains

   subroutine horvisc_workspace_ensure(nx, ny, nz)
      !! Lazy-allocate the Smagorinsky momentum-snapshot scratch.
      integer, intent(in) :: nx, ny, nz
      logical :: need_alloc

      need_alloc = .false.
      if (.not. allocated(horvisc_hu_pre)) then
         need_alloc = .true.
      else if (size(horvisc_hu_pre, 1) /= nx .or. &
               size(horvisc_hu_pre, 2) /= ny .or. &
               size(horvisc_hu_pre, 3) /= nz) then
         call ml_horizontal_viscosity_cleanup()
         need_alloc = .true.
      end if

      if (need_alloc) then
         allocate (horvisc_hu_pre(nx, ny, nz), horvisc_hv_pre(nx, ny, nz))
         !$omp target enter data map(alloc: horvisc_hu_pre, horvisc_hv_pre)
      end if
   end subroutine horvisc_workspace_ensure

   subroutine bih_workspace_ensure(nx, ny, nz)
      !! Lazy-allocate the biharmonic intermediate-Laplacian scratch.
      integer, intent(in) :: nx, ny, nz
      logical :: need_alloc

      need_alloc = .false.
      if (.not. allocated(bih_lap_u)) then
         need_alloc = .true.
      else if (size(bih_lap_u, 1) /= nx .or. &
               size(bih_lap_u, 2) /= ny .or. &
               size(bih_lap_u, 3) /= nz) then
         !$omp target exit data map(delete: bih_lap_u, bih_lap_v)
         deallocate (bih_lap_u, bih_lap_v)
         need_alloc = .true.
      end if

      if (need_alloc) then
         allocate (bih_lap_u(nx, ny, nz), bih_lap_v(nx, ny, nz))
         !$omp target enter data map(alloc: bih_lap_u, bih_lap_v)
      end if
   end subroutine bih_workspace_ensure

   subroutine ml_horizontal_viscosity_cleanup()
      !! Release the lateral-viscosity persistent workspaces.  Called by the
      !! solver at simulation end.
      if (allocated(horvisc_hu_pre)) then
         !$omp target exit data map(delete: horvisc_hu_pre, horvisc_hv_pre)
         deallocate (horvisc_hu_pre, horvisc_hv_pre)
      end if
      if (allocated(bih_lap_u)) then
         !$omp target exit data map(delete: bih_lap_u, bih_lap_v)
         deallocate (bih_lap_u, bih_lap_v)
      end if
   end subroutine ml_horizontal_viscosity_cleanup

   subroutine ml_horizontal_viscosity_smag(h_layer, hu_layer, hv_layer, &
                                           cs, grid, nz, dt)
      !! Smagorinsky (1963) horizontal viscosity for layer momentum:
      !! nu_h = (Cs*delta)^2*|S| Laplacian diffusion, delta=sqrt(dx*dy).
      !!
      !! Read-write hazard: the update reads neighbour velocities while
      !! writing the own cell.  `do concurrent` is Jacobi on GPU but
      !! sequential (Gauss-Seidel) on CPU, so we snapshot hu/hv_layer into
      !! frozen scratch and read neighbours from it for bit-identical
      !! results across backends.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: cs, dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: delta2, inv_dx, inv_dy, inv_2dx, inv_2dy, nu_max, h_col
      real(wp) :: h_k, u_c, v_c
      real(wp) :: u_ip, u_im, u_jp, u_jm, v_ip, v_im, v_jp, v_jm
      real(wp) :: dudx, dvdy, dudy, dvdx
      real(wp) :: strain, nu_h, lap_u, lap_v

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      if (cs <= 0.0_wp) return

      delta2 = cs*cs*dx*dy
      inv_dx = 1.0_wp/dx
      inv_dy = 1.0_wp/dy
      inv_2dx = 0.5_wp*inv_dx
      inv_2dy = 0.5_wp*inv_dy
      ! Stability: nu_h * dt * (1/dx^2 + 1/dy^2) <= 0.5
      nu_max = 0.5_wp/(dt*(inv_dx*inv_dx + inv_dy*inv_dy) + 1.0e-30_wp)

      call horvisc_workspace_ensure(nx, ny, nz)
      associate (hu_pre => horvisc_hu_pre, hv_pre => horvisc_hv_pre)

         ! Snapshot the input momentum so neighbour reads come from the
         ! pre-update image on every backend.
         !$omp target teams distribute parallel do collapse(3)
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
            hu_pre(i, j, k) = hu_layer(i, j, k)
            hv_pre(i, j, k) = hv_layer(i, j, k)
         end do
         end do
         end do
         !$omp end target teams distribute parallel do

         !$omp target teams distribute parallel do collapse(2) private(k, kk, h_k, u_c, v_c, u_ip, u_im, u_jp, u_jm, v_ip, v_im, v_jp, v_jm, dudx, dvdy, dudy, dvdx, strain, nu_h, lap_u, lap_v, h_col)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost

            h_col = 0.0_wp
            do kk = 1, nz
               h_col = h_col + h_layer(i, j, kk)
            end do
            if (h_col < ML_MIN_DEPTH) cycle

            do k = 1, nz
               h_k = h_layer(i, j, k)
               if (h_k < THIN_LAYER_THRESHOLD) cycle

               ! Centre velocity
               u_c = hu_pre(i, j, k)/h_k
               v_c = hv_pre(i, j, k)/h_k

               ! Neighbour velocities (zero-gradient at dry cells)
               if (h_layer(i + 1, j, k) > THIN_LAYER_THRESHOLD) then
                  u_ip = hu_pre(i + 1, j, k)/h_layer(i + 1, j, k)
                  v_ip = hv_pre(i + 1, j, k)/h_layer(i + 1, j, k)
               else
                  u_ip = u_c
                  v_ip = v_c
               end if
               if (h_layer(i - 1, j, k) > THIN_LAYER_THRESHOLD) then
                  u_im = hu_pre(i - 1, j, k)/h_layer(i - 1, j, k)
                  v_im = hv_pre(i - 1, j, k)/h_layer(i - 1, j, k)
               else
                  u_im = u_c
                  v_im = v_c
               end if
               if (h_layer(i, j + 1, k) > THIN_LAYER_THRESHOLD) then
                  u_jp = hu_pre(i, j + 1, k)/h_layer(i, j + 1, k)
                  v_jp = hv_pre(i, j + 1, k)/h_layer(i, j + 1, k)
               else
                  u_jp = u_c
                  v_jp = v_c
               end if
               if (h_layer(i, j - 1, k) > THIN_LAYER_THRESHOLD) then
                  u_jm = hu_pre(i, j - 1, k)/h_layer(i, j - 1, k)
                  v_jm = hv_pre(i, j - 1, k)/h_layer(i, j - 1, k)
               else
                  u_jm = u_c
                  v_jm = v_c
               end if

               ! Strain rate components (central differences)
               dudx = (u_ip - u_im)*inv_2dx
               dvdy = (v_jp - v_jm)*inv_2dy
               dudy = (u_jp - u_jm)*inv_2dy
               dvdx = (v_ip - v_im)*inv_2dx

               ! |S| = sqrt( (du/dx - dv/dy)^2 + (du/dy + dv/dx)^2 )
               strain = sqrt((dudx - dvdy)**2 + (dudy + dvdx)**2)

               ! Smagorinsky viscosity with stability limiter
               nu_h = min(delta2*strain, nu_max)

               ! Laplacian of velocity
               lap_u = (u_ip - 2.0_wp*u_c + u_im)*inv_dx*inv_dx &
                       + (u_jp - 2.0_wp*u_c + u_jm)*inv_dy*inv_dy
               lap_v = (v_ip - 2.0_wp*v_c + v_im)*inv_dx*inv_dx &
                       + (v_jp - 2.0_wp*v_c + v_jm)*inv_dy*inv_dy

               ! Update momentum -- read from the frozen pre-image, write
               ! to the live array so later iterations see the old state.
               hu_layer(i, j, k) = hu_pre(i, j, k) + dt*h_k*nu_h*lap_u
               hv_layer(i, j, k) = hv_pre(i, j, k) + dt*h_k*nu_h*lap_v
            end do
         end do
         end do
         !$omp end target teams distribute parallel do

      end associate

   end subroutine ml_horizontal_viscosity_smag

   subroutine ml_horizontal_viscosity_biharmonic(h_layer, hu_layer, hv_layer, &
                                                 nu_4, grid, nz, dt)
      !! Biharmonic (del^4) lateral viscosity on layer momentum:
      !! d(hu)/dt = -nu_4*h*del^2(del^2 u) in two Laplacian passes (pass 1
      !! over interior + one ghost ring; pass 2 over the interior).
      !! Scale-selective vs the harmonic Smagorinsky operator.  `nu_4`
      !! capped at the explicit stability limit
      !! nu_4*dt*((pi/dx)^2+(pi/dy)^2)^2 <= 2.  Assumes nghost >= 2.
      real(wp), parameter :: PI = 3.14159265358979323846_wp
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: nu_4, dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: inv_dx2, inv_dy2, nu4_max, nu4_eff
      real(wp) :: h_k, u_c, v_c, u_ip, u_im, u_jp, u_jm, v_ip, v_im, v_jp, v_jm
      real(wp) :: lu_c, lu_ip, lu_im, lu_jp, lu_jm
      real(wp) :: lv_c, lv_ip, lv_im, lv_jp, lv_jm
      real(wp) :: lap2_u, lap2_v

      if (nu_4 <= 0.0_wp) return

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      inv_dx2 = 1.0_wp/(grid%dx*grid%dx)
      inv_dy2 = 1.0_wp/(grid%dy*grid%dy)
      ! Explicit forward-step stability cap on nu_4.
      nu4_max = 2.0_wp/(dt*(PI*PI*inv_dx2 + PI*PI*inv_dy2)**2 + 1.0e-30_wp)
      nu4_eff = min(nu_4, nu4_max)

      call bih_workspace_ensure(nx, ny, nz)

      ! ---- Pass 1: L = del^2(u), over interior + one ghost ring ----
      !$omp target teams distribute parallel do collapse(2) private(k, h_k, u_c, v_c, u_ip, u_im, u_jp, u_jm, v_ip, v_im, v_jp, v_jm)
      do j = nghost, ny - nghost + 1
      do i = nghost, nx - nghost + 1
         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k < THIN_LAYER_THRESHOLD) then
               bih_lap_u(i, j, k) = 0.0_wp
               bih_lap_v(i, j, k) = 0.0_wp
               cycle
            end if
            u_c = hu_layer(i, j, k)/h_k
            v_c = hv_layer(i, j, k)/h_k
            if (h_layer(i + 1, j, k) > THIN_LAYER_THRESHOLD) then
               u_ip = hu_layer(i + 1, j, k)/h_layer(i + 1, j, k)
               v_ip = hv_layer(i + 1, j, k)/h_layer(i + 1, j, k)
            else
               u_ip = u_c
               v_ip = v_c
            end if
            if (h_layer(i - 1, j, k) > THIN_LAYER_THRESHOLD) then
               u_im = hu_layer(i - 1, j, k)/h_layer(i - 1, j, k)
               v_im = hv_layer(i - 1, j, k)/h_layer(i - 1, j, k)
            else
               u_im = u_c
               v_im = v_c
            end if
            if (h_layer(i, j + 1, k) > THIN_LAYER_THRESHOLD) then
               u_jp = hu_layer(i, j + 1, k)/h_layer(i, j + 1, k)
               v_jp = hv_layer(i, j + 1, k)/h_layer(i, j + 1, k)
            else
               u_jp = u_c
               v_jp = v_c
            end if
            if (h_layer(i, j - 1, k) > THIN_LAYER_THRESHOLD) then
               u_jm = hu_layer(i, j - 1, k)/h_layer(i, j - 1, k)
               v_jm = hv_layer(i, j - 1, k)/h_layer(i, j - 1, k)
            else
               u_jm = u_c
               v_jm = v_c
            end if
            bih_lap_u(i, j, k) = (u_ip - 2.0_wp*u_c + u_im)*inv_dx2 &
                                 + (u_jp - 2.0_wp*u_c + u_jm)*inv_dy2
            bih_lap_v(i, j, k) = (v_ip - 2.0_wp*v_c + v_im)*inv_dx2 &
                                 + (v_jp - 2.0_wp*v_c + v_jm)*inv_dy2
         end do
      end do
      end do
      !$omp end target teams distribute parallel do

      ! ---- Pass 2: hu -= dt*h*nu_4*del^2(L), over the interior ----
      !$omp target teams distribute parallel do collapse(2) private(k, h_k, lu_c, lu_ip, lu_im, lu_jp, lu_jm, lv_c, lv_ip, lv_im, lv_jp, lv_jm, lap2_u, lap2_v)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k < THIN_LAYER_THRESHOLD) cycle
            lu_c = bih_lap_u(i, j, k)
            lv_c = bih_lap_v(i, j, k)
            if (h_layer(i + 1, j, k) > THIN_LAYER_THRESHOLD) then
               lu_ip = bih_lap_u(i + 1, j, k)
               lv_ip = bih_lap_v(i + 1, j, k)
            else
               lu_ip = lu_c
               lv_ip = lv_c
            end if
            if (h_layer(i - 1, j, k) > THIN_LAYER_THRESHOLD) then
               lu_im = bih_lap_u(i - 1, j, k)
               lv_im = bih_lap_v(i - 1, j, k)
            else
               lu_im = lu_c
               lv_im = lv_c
            end if
            if (h_layer(i, j + 1, k) > THIN_LAYER_THRESHOLD) then
               lu_jp = bih_lap_u(i, j + 1, k)
               lv_jp = bih_lap_v(i, j + 1, k)
            else
               lu_jp = lu_c
               lv_jp = lv_c
            end if
            if (h_layer(i, j - 1, k) > THIN_LAYER_THRESHOLD) then
               lu_jm = bih_lap_u(i, j - 1, k)
               lv_jm = bih_lap_v(i, j - 1, k)
            else
               lu_jm = lu_c
               lv_jm = lv_c
            end if
            lap2_u = (lu_ip - 2.0_wp*lu_c + lu_im)*inv_dx2 &
                     + (lu_jp - 2.0_wp*lu_c + lu_jm)*inv_dy2
            lap2_v = (lv_ip - 2.0_wp*lv_c + lv_im)*inv_dx2 &
                     + (lv_jp - 2.0_wp*lv_c + lv_jm)*inv_dy2
            hu_layer(i, j, k) = hu_layer(i, j, k) - dt*h_k*nu4_eff*lap2_u
            hv_layer(i, j, k) = hv_layer(i, j, k) - dt*h_k*nu4_eff*lap2_v
         end do
      end do
      end do
      !$omp end target teams distribute parallel do

   end subroutine ml_horizontal_viscosity_biharmonic

end module rki_ml_horizontal_viscosity
