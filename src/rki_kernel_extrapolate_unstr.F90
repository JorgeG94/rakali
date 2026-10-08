!! Second-order extrapolation kernel for the unstructured solver
module rki_kernel_extrapolate_unstr
   !! Reconstruct cell-centroid values to edge midpoints with a per-cell
   !! linear gradient and slope limiter, in three passes (centroid prep,
   !! per-cell extrapolation, centroid restore). Cell loops are
   !! `do concurrent`: each iteration writes only to its own cell's
   !! edge-value arrays, and the pointwise helpers are all `pure`.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   implicit none
   private

   public :: extrapolate_second_order_edge_unstr

   ! Slope-limiter / wet-dry blending constants. Defined locally so the
   ! `*_impl` helpers below inline in NVHPC device codegen (same TU).
   real(wp), parameter :: LIMITER_EPSILON = 1.0e-100_wp
   real(wp), parameter :: HFACTOR_A = 0.3_wp
   real(wp), parameter :: HFACTOR_B = 0.1_wp
   real(wp), parameter :: HFACTOR_C = 1.0_wp/(HFACTOR_A - HFACTOR_B)
   real(wp), parameter :: HFACTOR_D = 1.0_wp - HFACTOR_C*HFACTOR_A
   real(wp), parameter :: HFACTOR_HC_FLOOR = 1.0e-6_wp

contains

   subroutine extrapolate_second_order_edge_unstr(state, mesh)
      !! Reconstruct edge-midpoint values from centroid values.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      call centroid_pass(state, mesh)
      call per_cell_extrapolation(state, mesh)
      call centroid_restore(state, mesh)
   end subroutine extrapolate_second_order_edge_unstr

   pure subroutine centroid_pass(state, mesh)
      !! Compute heights, zero momenta on dry cells, prep velocity-form.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k
      real(wp) :: dk, xmom_out, ymom_out, inv_dk
      logical :: is_dry, do_extrap

      do concurrent(k=1:mesh%ncells_total)
         dk = max(state%barotropic%stage_cv(k) - state%barotropic%bed_cv(k), 0.0_wp)
         state%barotropic%height_cv(k) = dk

         is_dry = (dk <= state%minimum_allowed_height)
         do_extrap = state%extrapolate_velocity_second_order .and. .not. is_dry

         if (is_dry) then
            xmom_out = 0.0_wp
            ymom_out = 0.0_wp
         else
            xmom_out = state%barotropic%xmom_cv(k)
            ymom_out = state%barotropic%ymom_cv(k)
         end if

         if (do_extrap) then
            inv_dk = 1.0_wp/dk
            state%barotropic%x_centroid_work(k) = xmom_out
            state%barotropic%y_centroid_work(k) = ymom_out
         else
            inv_dk = 1.0_wp
            state%barotropic%x_centroid_work(k) = 0.0_wp
            state%barotropic%y_centroid_work(k) = 0.0_wp
         end if

         state%barotropic%xmom_cv(k) = xmom_out*inv_dk
         state%barotropic%ymom_cv(k) = ymom_out*inv_dk
      end do
   end subroutine centroid_pass

   pure subroutine per_cell_extrapolation(state, mesh)
      !! Compute edge-midpoint values for stage, height, xmom, ymom and
      !! recover bed_ev from the well-balanced identity.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k

      do concurrent(k=1:mesh%ncells_total)
         call extrapolate_one_cell(state, mesh, k)
      end do
   end subroutine per_cell_extrapolation

   pure subroutine extrapolate_one_cell(state, mesh, k)
      !! Per-cell body of the extrapolation pass. Pure: never writes to
      !! any cell other than `k`.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: k

      real(wp) :: dxv0, dxv1, dxv2, dyv0, dyv1, dyv2
      real(wp) :: stage_ev_local(3), xmom_ev_local(3), ymom_ev_local(3)
      real(wp) :: height_ev_local(3)
      integer :: i

      ! Edge-midpoint offsets from this cell's centroid.
      dxv0 = mesh%edge_mid_x(1, k) - mesh%centroid_x(k)
      dxv1 = mesh%edge_mid_x(2, k) - mesh%centroid_x(k)
      dxv2 = mesh%edge_mid_x(3, k) - mesh%centroid_x(k)
      dyv0 = mesh%edge_mid_y(1, k) - mesh%centroid_y(k)
      dyv1 = mesh%edge_mid_y(2, k) - mesh%centroid_y(k)
      dyv2 = mesh%edge_mid_y(3, k) - mesh%centroid_y(k)

      select case (mesh%n_boundaries(k))
      case (3)
         ! Isolated triangle -- constant (zeroth-order) edge values.
         call set_constant_edge_values_impl(state%barotropic%stage_cv(k), stage_ev_local)
         call set_constant_edge_values_impl(state%barotropic%xmom_cv(k), xmom_ev_local)
         call set_constant_edge_values_impl(state%barotropic%ymom_cv(k), ymom_ev_local)
         call set_constant_edge_values_impl(state%barotropic%height_cv(k), height_ev_local)

      case (2)
         call extrapolate_one_neighbour_branch( &
            state, mesh, k, &
            dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
            stage_ev_local, xmom_ev_local, ymom_ev_local, height_ev_local)

      case default
         ! 0 or 1 boundary edges -- full 3-neighbour reconstruction.
         call extrapolate_three_neighbour_branch( &
            state, mesh, k, &
            dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
            stage_ev_local, xmom_ev_local, ymom_ev_local, height_ev_local)
      end select

      ! Velocity-form recovery: if on, xmom_ev/ymom_ev hold velocities;
      ! multiply by reconstructed height to recover momenta.
      if (state%extrapolate_velocity_second_order) then
         do i = 1, 3
            xmom_ev_local(i) = xmom_ev_local(i)*height_ev_local(i)
            ymom_ev_local(i) = ymom_ev_local(i)*height_ev_local(i)
         end do
      end if

      ! Commit edge values for this cell.
      do i = 1, 3
         state%barotropic%stage_ev(i, k) = stage_ev_local(i)
         state%barotropic%xmom_ev(i, k) = xmom_ev_local(i)
         state%barotropic%ymom_ev(i, k) = ymom_ev_local(i)
         state%barotropic%height_ev(i, k) = height_ev_local(i)
         ! Recover bed from the well-balanced identity, preserving
         ! stage_ev = bed_ev + height_ev exactly (lake-at-rest).
         state%barotropic%bed_ev(i, k) = stage_ev_local(i) - height_ev_local(i)
      end do
   end subroutine extrapolate_one_cell

   pure subroutine extrapolate_three_neighbour_branch( &
      state, mesh, k, &
      dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
      stage_ev_local, xmom_ev_local, ymom_ev_local, height_ev_local)
      !! Full 3-neighbour reconstruction via the closed-form 2x2
      !! gradient. The limiter β is blended wet->dry via `hfactor` so
      !! the slope dials down as cell/neighbours dry out.
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: k
      real(wp), intent(in) :: dxv0, dxv1, dxv2, dyv0, dyv1, dyv2
      real(wp), intent(out) :: stage_ev_local(3)
      real(wp), intent(out) :: xmom_ev_local(3)
      real(wp), intent(out) :: ymom_ev_local(3)
      real(wp), intent(out) :: height_ev_local(3)

      integer :: k0, k1, k2
      real(wp) :: x0, y0, x1, y1, x2, y2
      real(wp) :: dx1, dy1, dx2, dy2, area2, inv_area2
      real(wp) :: hfactor
      real(wp) :: beta_stage, beta_xmom, beta_ymom

      k0 = mesh%surrogate_neighbours(1, k)
      k1 = mesh%surrogate_neighbours(2, k)
      k2 = mesh%surrogate_neighbours(3, k)

      x0 = mesh%centroid_x(k0)
      y0 = mesh%centroid_y(k0)
      x1 = mesh%centroid_x(k1)
      y1 = mesh%centroid_y(k1)
      x2 = mesh%centroid_x(k2)
      y2 = mesh%centroid_y(k2)

      dx1 = x1 - x0
      dy1 = y1 - y0
      dx2 = x2 - x0
      dy2 = y2 - y0
      area2 = dy2*dx1 - dy1*dx2
      ! Collinear/coincident neighbours -> area2 vanishes; fall back to
      ! constant edge values.
      if (abs(area2) <= 1.0e-30_wp) then
         call set_constant_edge_values_impl(state%barotropic%stage_cv(k), stage_ev_local)
         call set_constant_edge_values_impl(state%barotropic%xmom_cv(k), xmom_ev_local)
         call set_constant_edge_values_impl(state%barotropic%ymom_cv(k), ymom_ev_local)
         call set_constant_edge_values_impl(state%barotropic%height_cv(k), height_ev_local)
         return
      end if
      inv_area2 = 1.0_wp/area2

      hfactor = compute_hfactor_impl(state%barotropic%height_cv(k), &
                                     state%barotropic%height_cv(k0), &
                                     state%barotropic%height_cv(k1), &
                                     state%barotropic%height_cv(k2), &
                                     state%minimum_allowed_height)

      beta_stage = state%beta_w_dry + (state%beta_w - state%beta_w_dry)*hfactor
      beta_xmom = state%beta_uh_dry + (state%beta_uh - state%beta_uh_dry)*hfactor
      beta_ymom = state%beta_vh_dry + (state%beta_vh - state%beta_vh_dry)*hfactor

      call calc_edge_values_three_neighbours_impl( &
        state%barotropic%stage_cv(k), state%barotropic%stage_cv(k0), state%barotropic%stage_cv(k1), state%barotropic%stage_cv(k2), &
         dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
         dx1, dx2, dy1, dy2, inv_area2, &
         beta_stage, stage_ev_local)

      call calc_edge_values_three_neighbours_impl( &
    state%barotropic%height_cv(k), state%barotropic%height_cv(k0), state%barotropic%height_cv(k1), state%barotropic%height_cv(k2), &
         dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
         dx1, dx2, dy1, dy2, inv_area2, &
         beta_stage, height_ev_local)

      call calc_edge_values_three_neighbours_impl( &
         state%barotropic%xmom_cv(k), state%barotropic%xmom_cv(k0), state%barotropic%xmom_cv(k1), state%barotropic%xmom_cv(k2), &
         dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
         dx1, dx2, dy1, dy2, inv_area2, &
         beta_xmom, xmom_ev_local)

      call calc_edge_values_three_neighbours_impl( &
         state%barotropic%ymom_cv(k), state%barotropic%ymom_cv(k0), state%barotropic%ymom_cv(k1), state%barotropic%ymom_cv(k2), &
         dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
         dx1, dx2, dy1, dy2, inv_area2, &
         beta_ymom, ymom_ev_local)
   end subroutine extrapolate_three_neighbour_branch

   pure subroutine extrapolate_one_neighbour_branch( &
      state, mesh, k, &
      dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
      stage_ev_local, xmom_ev_local, ymom_ev_local, height_ev_local)
      !! Two-boundary fallback: gradient taken toward the single internal
      !! neighbour, displacement normalised by squared distance:
      !! `dq * (dx, dy)/|dx,dy|^2 . (xv, yv)`.
      type(unstr_state_t), intent(in) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: k
      real(wp), intent(in) :: dxv0, dxv1, dxv2, dyv0, dyv1, dyv2
      real(wp), intent(out) :: stage_ev_local(3)
      real(wp), intent(out) :: xmom_ev_local(3)
      real(wp), intent(out) :: ymom_ev_local(3)
      real(wp), intent(out) :: height_ev_local(3)

      integer :: i, kn
      real(wp) :: dx, dy, dist2
      real(wp) :: grad_dx2, grad_dy2
      real(wp) :: dq1, qmin, qmax
      real(wp) :: dqv(3)

      ! Find the single internal neighbour.
      kn = k
      do i = 1, 3
         if (mesh%surrogate_neighbours(i, k) /= k) then
            kn = mesh%surrogate_neighbours(i, k)
            exit
         end if
      end do

      if (kn == k) then
         ! Should not happen for n_boundaries == 2, but be defensive.
         call set_constant_edge_values_impl(state%barotropic%stage_cv(k), stage_ev_local)
         call set_constant_edge_values_impl(state%barotropic%xmom_cv(k), xmom_ev_local)
         call set_constant_edge_values_impl(state%barotropic%ymom_cv(k), ymom_ev_local)
         call set_constant_edge_values_impl(state%barotropic%height_cv(k), height_ev_local)
         return
      end if

      dx = mesh%centroid_x(kn) - mesh%centroid_x(k)
      dy = mesh%centroid_y(kn) - mesh%centroid_y(k)
      dist2 = dx*dx + dy*dy
      if (dist2 > 0.0_wp) then
         grad_dx2 = dx/dist2
         grad_dy2 = dy/dist2
      else
         grad_dx2 = 0.0_wp
         grad_dy2 = 0.0_wp
      end if

      ! Stage
      dq1 = state%barotropic%stage_cv(kn) - state%barotropic%stage_cv(k)
      call calc_dqv_one_neighbour_impl(dq1, grad_dx2, grad_dy2, &
                                       dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv)
      call qrange_one_diff_impl(dq1, qmin, qmax)
      call limit_gradient_impl(dqv, qmin, qmax, state%beta_w)
      stage_ev_local(1) = state%barotropic%stage_cv(k) + dqv(1)
      stage_ev_local(2) = state%barotropic%stage_cv(k) + dqv(2)
      stage_ev_local(3) = state%barotropic%stage_cv(k) + dqv(3)

      ! Height
      dq1 = state%barotropic%height_cv(kn) - state%barotropic%height_cv(k)
      call calc_dqv_one_neighbour_impl(dq1, grad_dx2, grad_dy2, &
                                       dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv)
      call qrange_one_diff_impl(dq1, qmin, qmax)
      call limit_gradient_impl(dqv, qmin, qmax, state%beta_w)
      height_ev_local(1) = state%barotropic%height_cv(k) + dqv(1)
      height_ev_local(2) = state%barotropic%height_cv(k) + dqv(2)
      height_ev_local(3) = state%barotropic%height_cv(k) + dqv(3)

      ! x-momentum
      dq1 = state%barotropic%xmom_cv(kn) - state%barotropic%xmom_cv(k)
      call calc_dqv_one_neighbour_impl(dq1, grad_dx2, grad_dy2, &
                                       dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv)
      call qrange_one_diff_impl(dq1, qmin, qmax)
      call limit_gradient_impl(dqv, qmin, qmax, state%beta_uh)
      xmom_ev_local(1) = state%barotropic%xmom_cv(k) + dqv(1)
      xmom_ev_local(2) = state%barotropic%xmom_cv(k) + dqv(2)
      xmom_ev_local(3) = state%barotropic%xmom_cv(k) + dqv(3)

      ! y-momentum
      dq1 = state%barotropic%ymom_cv(kn) - state%barotropic%ymom_cv(k)
      call calc_dqv_one_neighbour_impl(dq1, grad_dx2, grad_dy2, &
                                       dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, dqv)
      call qrange_one_diff_impl(dq1, qmin, qmax)
      call limit_gradient_impl(dqv, qmin, qmax, state%beta_vh)
      ymom_ev_local(1) = state%barotropic%ymom_cv(k) + dqv(1)
      ymom_ev_local(2) = state%barotropic%ymom_cv(k) + dqv(2)
      ymom_ev_local(3) = state%barotropic%ymom_cv(k) + dqv(3)
   end subroutine extrapolate_one_neighbour_branch

   pure subroutine centroid_restore(state, mesh)
      !! If velocity-form extrapolation was on, restore the original
      !! centroid momenta from the work arrays.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k

      if (.not. state%extrapolate_velocity_second_order) return

      do concurrent(k=1:mesh%ncells_total)
         state%barotropic%xmom_cv(k) = state%barotropic%x_centroid_work(k)
         state%barotropic%ymom_cv(k) = state%barotropic%y_centroid_work(k)
      end do
   end subroutine centroid_restore

   ! Inlined pointwise helpers (`*_impl`): self-contained pure helpers
   ! in the same TU as `extrapolate_one_cell` so NVHPC device codegen
   ! inlines them.

   pure subroutine qrange_one_diff_impl(dq1, qmin, qmax)
      !$acc routine seq
      real(wp), intent(in) :: dq1
      real(wp), intent(out) :: qmin, qmax

      if (dq1 >= 0.0_wp) then
         qmin = 0.0_wp
         qmax = dq1
      else
         qmin = dq1
         qmax = 0.0_wp
      end if
   end subroutine qrange_one_diff_impl

   pure subroutine qrange_three_diff_impl(dq0, dq1, dq2, qmin, qmax)
      !$acc routine seq
      real(wp), intent(in) :: dq0, dq1, dq2
      real(wp), intent(out) :: qmin, qmax

      qmax = max(dq0, dq0 + dq1, dq0 + dq2, 0.0_wp)
      qmin = min(dq0, dq0 + dq1, dq0 + dq2, 0.0_wp)
   end subroutine qrange_three_diff_impl

   pure subroutine limit_gradient_impl(dqv, qmin, qmax, beta)
      !$acc routine seq
      real(wp), intent(inout) :: dqv(3)
      real(wp), intent(in) :: qmin, qmax, beta

      real(wp) :: r, phi, d
      integer :: i

      r = 1000.0_wp
      do i = 1, 3
         d = dqv(i)
         if (d < -LIMITER_EPSILON) then
            r = min(r, qmin/d)
         else if (d > LIMITER_EPSILON) then
            r = min(r, qmax/d)
         end if
      end do

      phi = min(r*beta, 1.0_wp)
      dqv(1) = dqv(1)*phi
      dqv(2) = dqv(2)*phi
      dqv(3) = dqv(3)*phi
   end subroutine limit_gradient_impl

   pure subroutine calc_edge_values_three_neighbours_impl( &
      cv_k, cv_k0, cv_k1, cv_k2, &
      dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
      dx1, dx2, dy1, dy2, inv_area2, &
      beta, edge_values)
      !$acc routine seq
      real(wp), intent(in) :: cv_k, cv_k0, cv_k1, cv_k2
      real(wp), intent(in) :: dxv0, dxv1, dxv2
      real(wp), intent(in) :: dyv0, dyv1, dyv2
      real(wp), intent(in) :: dx1, dx2, dy1, dy2
      real(wp), intent(in) :: inv_area2
      real(wp), intent(in) :: beta
      real(wp), intent(out) :: edge_values(3)

      real(wp) :: dq0, dq1, dq2
      real(wp) :: a, b
      real(wp) :: dqv(3)
      real(wp) :: qmin, qmax

      dq0 = cv_k0 - cv_k
      dq1 = cv_k1 - cv_k0
      dq2 = cv_k2 - cv_k0

      a = (dy2*dq1 - dy1*dq2)*inv_area2
      b = (dx1*dq2 - dx2*dq1)*inv_area2

      dqv(1) = a*dxv0 + b*dyv0
      dqv(2) = a*dxv1 + b*dyv1
      dqv(3) = a*dxv2 + b*dyv2

      call qrange_three_diff_impl(dq0, dq1, dq2, qmin, qmax)
      call limit_gradient_impl(dqv, qmin, qmax, beta)

      edge_values(1) = cv_k + dqv(1)
      edge_values(2) = cv_k + dqv(2)
      edge_values(3) = cv_k + dqv(3)
   end subroutine calc_edge_values_three_neighbours_impl

   pure subroutine calc_dqv_one_neighbour_impl( &
      dq1, dx2, dy2, &
      dxv0, dxv1, dxv2, dyv0, dyv1, dyv2, &
      dqv)
      !$acc routine seq
      real(wp), intent(in) :: dq1
      real(wp), intent(in) :: dx2, dy2
      real(wp), intent(in) :: dxv0, dxv1, dxv2
      real(wp), intent(in) :: dyv0, dyv1, dyv2
      real(wp), intent(out) :: dqv(3)

      real(wp) :: a, b

      a = dq1*dx2
      b = dq1*dy2

      dqv(1) = a*dxv0 + b*dyv0
      dqv(2) = a*dxv1 + b*dyv1
      dqv(3) = a*dxv2 + b*dyv2
   end subroutine calc_dqv_one_neighbour_impl

   pure subroutine set_constant_edge_values_impl(cv_k, edge_values)
      !$acc routine seq
      real(wp), intent(in) :: cv_k
      real(wp), intent(out) :: edge_values(3)

      edge_values(1) = cv_k
      edge_values(2) = cv_k
      edge_values(3) = cv_k
   end subroutine set_constant_edge_values_impl

   pure function compute_hfactor_impl(hc, h0, h1, h2, min_h) result(hfactor)
      !$acc routine seq
      real(wp), intent(in) :: hc, h0, h1, h2, min_h
      real(wp) :: hfactor

      real(wp) :: hmin, hmax, tmp1, tmp2

      hmin = min(hc, h0, h1, h2)
      hmax = max(hc, h0, h1, h2)

      tmp1 = HFACTOR_C*max(hmin, 0.0_wp)/max(hc, HFACTOR_HC_FLOOR) + HFACTOR_D
      tmp2 = HFACTOR_C*max(hc, 0.0_wp)/max(hmax, HFACTOR_HC_FLOOR) + HFACTOR_D

      hfactor = max(0.0_wp, min(tmp1, tmp2, 1.0_wp))

      hfactor = min(1.2_wp*max(hmin - min_h, 0.0_wp) &
                    /(max(hmin, 0.0_wp) + min_h), hfactor)
   end function compute_hfactor_impl

end module rki_kernel_extrapolate_unstr
