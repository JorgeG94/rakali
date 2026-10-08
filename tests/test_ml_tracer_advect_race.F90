!! Unit test for the do-concurrent read-write race in `ml_advect_tracer`.
!!
!! The race: inside the original `do concurrent(j,i)` the kernel read
!! neighbour cells `hTr_layer(i±1,j,k)` / `hTr_layer(i,j±1,k)` while
!! writing its own cell `hTr_layer(i,j,k)`.  On sequential backends
!! (gfortran CPU, nvfortran serial) later iterations see already-updated
!! neighbour values; the two visits of a shared face disagree and
!! ~1e-7 tracer mass leaks per step (UB).  The fix freezes a
!! `adv_hTr_snap` copy before the kernel and redirects all neighbour
!! reads to the snapshot (Jacobi semantics).
!!
!! Test strategy: periodic closed-basin (ghost = wrap-around).  Under
!! periodic BCs total interior tracer mass is exactly conserved by the
!! first-order upwind scheme; the racy kernel violates this on a
!! sequential backend.
module test_ml_tracer_advect_race
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_ml_tracers, only: ml_advect_tracer, tracer_scratch_t, &
                             tracer_scratch_ensure, tracer_scratch_destroy
   implicit none
   private

   public :: collect_ml_tracer_advect_race_tests

   integer, parameter :: NX_PHYS = 8
   integer, parameter :: NY_PHYS = 8
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 3
   real(wp), parameter :: DX_VAL = 100.0_wp
   real(wp), parameter :: DY_VAL = 100.0_wp
   real(wp), parameter :: H_VAL = 10.0_wp   !! uniform layer thickness
   real(wp), parameter :: U_FACE = 5.0_wp   !! uniform face velocity (m/s)
   real(wp), parameter :: DT_VAL = 5.0_wp   !! time step (s)

contains

   subroutine collect_ml_tracer_advect_race_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ml_advect_tracer_periodic_conservation", &
                               test_advect_tracer_periodic_conservation) &
                  ]
   end subroutine collect_ml_tracer_advect_race_tests

   subroutine test_advect_tracer_periodic_conservation(error)
      !! Periodic closed-basin conservation test for `ml_advect_tracer`.
      !!
      !! Set up an 8x8x3 grid with uniform layer thickness H_VAL = 10 m,
      !! a strong checkerboard tracer gradient (conc = 10 or 30 per cell),
      !! and a uniform non-zero advective velocity so every face sees
      !! real flux.  Fill ghost cells as periodic wrap-arounds so the
      !! interior is a closed periodic domain and total tracer mass must
      !! be conserved to round-off after one kernel call.
      !!
      !! On the racy (pre-fix) kernel sequential backends see ~1e-7 relative
      !! mass change; on the fixed (Jacobi-snapshot) kernel the change is
      !! below 1e-12 relative.
      type(error_type), allocatable, intent(out) :: error

      type(hgrid_t) :: grid
      integer :: i, j, k, nx, ny, ng

      real(wp), allocatable :: h_layer(:, :, :)
      real(wp), allocatable :: hu_layer(:, :, :)
      real(wp), allocatable :: hv_layer(:, :, :)
      real(wp), allocatable :: hTr_layer(:, :, :)
      real(wp), allocatable :: mass_flux_x(:, :)
      real(wp), allocatable :: mass_flux_y(:, :)
      type(tracer_scratch_t) :: scratch

      real(wp) :: m0, m1, conc_ij
      real(wp), parameter :: TOL_CONS = 1.0e-9_wp

      ! ---- build a minimal hgrid_t directly ----
      grid%nx_total = NX_PHYS + 2*NGHOST
      grid%ny_total = NY_PHYS + 2*NGHOST
      grid%nghost = NGHOST
      grid%dx = DX_VAL
      grid%dy = DY_VAL

      nx = grid%nx_total
      ny = grid%ny_total
      ng = grid%nghost

      allocate (h_layer(nx, ny, NZ))
      allocate (hu_layer(nx, ny, NZ))
      allocate (hv_layer(nx, ny, NZ))
      allocate (hTr_layer(nx, ny, NZ))
      allocate (mass_flux_x(nx, ny))
      allocate (mass_flux_y(nx, ny))

      ! ---- uniform thickness everywhere (incl. ghosts) ----
      h_layer = H_VAL

      ! ---- checkerboard tracer: conc = 10 or 30, large gradient ----
      do k = 1, NZ
         do j = 1, ny
            do i = 1, nx
               if (mod(i + j, 2) == 0) then
                  conc_ij = 30.0_wp
               else
                  conc_ij = 10.0_wp
               end if
               hTr_layer(i, j, k) = H_VAL*conc_ij
            end do
         end do
      end do

      ! ---- uniform velocity: hu = U_FACE * H_VAL per layer ----
      ! A uniform-velocity uniform-thickness field makes the per-face
      ! CWC split exact: each layer carries exactly (1/nz) of the column
      ! mass flux, frac_active = 1/nz everywhere.
      hu_layer = U_FACE*H_VAL
      hv_layer = U_FACE*H_VAL

      ! ---- barotropic mass flux = U * H_total ----
      ! For uniform-thickness uniform-velocity: column velocity = U_FACE,
      ! total depth = NZ * H_VAL.
      mass_flux_x = U_FACE*real(NZ, wp)*H_VAL
      mass_flux_y = U_FACE*real(NZ, wp)*H_VAL

      ! ---- periodic ghost fill: wrap interior edges into ghost rows ----
      ! x-direction: ghost columns [1..ng] <- mirror of interior columns
      ! [nx-2*ng+1..nx-ng]; ghost columns [nx-ng+1..nx] <- [ng+1..2*ng].
      do k = 1, NZ
         do j = 1, ny
            do i = 1, ng
               ! left ghost (i=1..ng) <- right interior (nx-ng-ng+i..nx-ng+i)
               h_layer(i, j, k) = h_layer(nx - 2*ng + i, j, k)
               hTr_layer(i, j, k) = hTr_layer(nx - 2*ng + i, j, k)
               hu_layer(i, j, k) = hu_layer(nx - 2*ng + i, j, k)
               hv_layer(i, j, k) = hv_layer(nx - 2*ng + i, j, k)
               ! right ghost (nx-ng+i) <- left interior (ng+i)
               h_layer(nx - ng + i, j, k) = h_layer(ng + i, j, k)
               hTr_layer(nx - ng + i, j, k) = hTr_layer(ng + i, j, k)
               hu_layer(nx - ng + i, j, k) = hu_layer(ng + i, j, k)
               hv_layer(nx - ng + i, j, k) = hv_layer(ng + i, j, k)
            end do
         end do
      end do
      ! y-direction
      do k = 1, NZ
         do j = 1, ng
            do i = 1, nx
               h_layer(i, j, k) = h_layer(i, ny - 2*ng + j, k)
               hTr_layer(i, j, k) = hTr_layer(i, ny - 2*ng + j, k)
               hu_layer(i, j, k) = hu_layer(i, ny - 2*ng + j, k)
               hv_layer(i, j, k) = hv_layer(i, ny - 2*ng + j, k)
               h_layer(i, ny - ng + j, k) = h_layer(i, ng + j, k)
               hTr_layer(i, ny - ng + j, k) = hTr_layer(i, ng + j, k)
               hu_layer(i, ny - ng + j, k) = hu_layer(i, ng + j, k)
               hv_layer(i, ny - ng + j, k) = hv_layer(i, ng + j, k)
            end do
         end do
      end do
      ! mass flux ghost fill
      do j = 1, ng
         do i = 1, nx
            mass_flux_x(i, j) = mass_flux_x(i, ny - 2*ng + j)
            mass_flux_y(i, j) = mass_flux_y(i, ny - 2*ng + j)
            mass_flux_x(i, ny - ng + j) = mass_flux_x(i, ng + j)
            mass_flux_y(i, ny - ng + j) = mass_flux_y(i, ng + j)
         end do
      end do
      do j = 1, ny
         do i = 1, ng
            mass_flux_x(i, j) = mass_flux_x(nx - 2*ng + i, j)
            mass_flux_y(i, j) = mass_flux_y(nx - 2*ng + i, j)
            mass_flux_x(nx - ng + i, j) = mass_flux_x(ng + i, j)
            mass_flux_y(nx - ng + i, j) = mass_flux_y(ng + i, j)
         end do
      end do

      ! ---- total interior tracer mass before the step ----
      m0 = 0.0_wp
      do k = 1, NZ
         do j = ng + 1, ny - ng
            do i = ng + 1, nx - ng
               m0 = m0 + hTr_layer(i, j, k)
            end do
         end do
      end do

      ! ---- one advection step ----
      call tracer_scratch_ensure(scratch, nx, ny, NZ, &
                                 want_ppm=.false., want_upwind=.true., want_hdiff=.false.)
      call ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                            mass_flux_x, mass_flux_y, &
                            hTr_layer, scratch%adv_hTr_snap, &
                            grid, NZ, DT_VAL)

      ! ---- total interior tracer mass after the step ----
      m1 = 0.0_wp
      do k = 1, NZ
         do j = ng + 1, ny - ng
            do i = ng + 1, nx - ng
               m1 = m1 + hTr_layer(i, j, k)
            end do
         end do
      end do

      ! ---- conservation check ----
      ! On a periodic domain first-order upwind is exactly conservative
      ! to floating-point round-off.  The racy kernel leaks ~1e-7 relative
      ! on a sequential backend; the fixed kernel passes at 1e-9 relative.
      call check(error, abs(m1 - m0) < TOL_CONS*abs(m0), &
                 "ml_advect_tracer: interior tracer mass not conserved on periodic domain")
      if (allocated(error)) return

      ! ---- release scratch so subsequent tests start clean ----
      call tracer_scratch_destroy(scratch)

   end subroutine test_advect_tracer_periodic_conservation

end module test_ml_tracer_advect_race
