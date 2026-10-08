!! Parallel array fill helpers for stdpar (GPU/multicore).
module rki_array_utils
   !! `fill_gpu(arr, val)` replaces `allocate(arr, source=val)` (serial host
   !! zero-fill) with a `do concurrent` kernel. Generic over rank 1D/2D/3D;
   !! uses `lbound`/`ubound` so any lower bound works (e.g. `z_ref(0:nz,...)`).
   use rki_constants, only: wp
   implicit none
   private

   public :: fill_gpu

   interface fill_gpu
      module procedure fill_gpu_1d
      module procedure fill_gpu_2d
      module procedure fill_gpu_3d
   end interface fill_gpu

contains

   pure subroutine fill_gpu_1d(arr, val)
      ! assumed-shape-ok: generic fill utility; must accept any rank-1 slice
      real(wp), intent(inout) :: arr(:)
      real(wp), intent(in)    :: val
      integer :: i, lo, hi
      lo = lbound(arr, 1)
      hi = ubound(arr, 1)
      do concurrent(i=lo:hi)
         arr(i) = val
      end do
   end subroutine fill_gpu_1d

   pure subroutine fill_gpu_2d(arr, val)
      ! assumed-shape-ok: generic fill utility; must accept any rank-2 slice
      real(wp), intent(inout) :: arr(:, :)
      real(wp), intent(in)    :: val
      integer :: i, j, lo1, hi1, lo2, hi2
      lo1 = lbound(arr, 1)
      hi1 = ubound(arr, 1)
      lo2 = lbound(arr, 2)
      hi2 = ubound(arr, 2)
      do concurrent(j=lo2:hi2, i=lo1:hi1)
         arr(i, j) = val
      end do
   end subroutine fill_gpu_2d

   pure subroutine fill_gpu_3d(arr, val)
      ! assumed-shape-ok: generic fill utility; must accept any rank-3 slice
      real(wp), intent(inout) :: arr(:, :, :)
      real(wp), intent(in)    :: val
      integer :: i, j, k, lo1, hi1, lo2, hi2, lo3, hi3
      lo1 = lbound(arr, 1)
      hi1 = ubound(arr, 1)
      lo2 = lbound(arr, 2)
      hi2 = ubound(arr, 2)
      lo3 = lbound(arr, 3)
      hi3 = ubound(arr, 3)
      do concurrent(k=lo3:hi3, j=lo2:hi2, i=lo1:hi1)
         arr(i, j, k) = val
      end do
   end subroutine fill_gpu_3d

end module rki_array_utils
