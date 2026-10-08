!! Central-upwind flux kernel for the unstructured solver
module rki_kernel_flux_unstr
   !! Net flux divergence at every owned cell via the
   !! Kurganov-Noelle-Petrova central-upwind scheme on triangles.
   !!
   !! Cell-centric: each triangle loops its three edges, computes the
   !! central-upwind flux from both sides, and accumulates ONLY into its
   !! own update arrays. Each interior edge is visited twice (once per
   !! neighbour) -- no scatter, no race, no atomics -- which is what
   !! makes the kernel safe for `do concurrent` on GPU. Boundary edges
   !! (`mesh%neighbours(i, k) < 0`) dispatch per-edge BC state; group 0
   !! falls through to BC_WALL (normal-momentum mirror). The
   !! well-balanced gravity bed-slope source is added per cell, then
   !! accumulators are normalised by area to flux divergence (per unit
   !! area), ready to be multiplied by `dt` in the update step.
   use rki_constants, only: wp, GRAVITY
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_state_unstr, only: unstr_state_t
   use rki_boundary_unstr, only: BC_WALL, BC_OPEN, BC_INFLOW, BC_OUTFLOW, &
                                 BC_TIDAL, BC_NESTED, BC_CLAMPED, BC_DISCHARGE, BC_CHAPMAN
   implicit none
   private

   public :: compute_fluxes_central_unstr

   ! Pointwise central-upwind constants. Defined locally so the `*_impl`
   ! helpers below inline in NVHPC device codegen (same TU). The
   ! LOW_FROUDE_* values duplicate the `rki_state_unstr` params -- keep
   ! in sync.
   real(wp), parameter :: FLUX_EPSILON = 1.0e-15_wp
   real(wp), parameter :: FLUX_DENOM_FLOOR = 1.0e-100_wp
   integer, parameter :: LOW_FROUDE_OFF = 0
   integer, parameter :: LOW_FROUDE_1 = 1
   integer, parameter :: LOW_FROUDE_2 = 2

contains

   subroutine compute_fluxes_central_unstr(state, mesh)
      !! Top-level driver: zero `(stage_eu, xmom_eu, ymom_eu)` and
      !! recompute them from the current edge values produced by the
      !! extrapolation kernel.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: k

      !$omp target teams distribute parallel do
      do k = 1, mesh%ncells_owned
         call flux_one_cell(state, mesh, k)
      end do
      !$omp end target teams distribute parallel do
   end subroutine compute_fluxes_central_unstr

   pure subroutine flux_one_cell(state, mesh, k)
      !! Per-cell body of the flux pass. Pure: never writes to any cell
      !! other than `k`.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: k

      integer :: i, neighbour, m
      integer :: bc_group, bc_kind
      real(wp) :: bc_eta_ext, bc_un_ext, bc_h_ext, bc_Qn_ext
      real(wp) :: q_left(3), q_right(3)
      real(wp) :: edgeflux(3)
      real(wp) :: hc, zc, hle, hre, zl, zr, h_left, h_right, z_half
      real(wp) :: n1, n2, length
      real(wp) :: max_speed_local, pressure_flux, speed_max_last
      real(wp) :: pressuregrad_work
      real(wp) :: stage_eu_acc, xmom_eu_acc, ymom_eu_acc
      real(wp) :: inv_area

      ! Centroid values for this cell -- needed by the bed-slope source.
      hc = state%barotropic%height_cv(k)
      zc = state%barotropic%bed_cv(k)

      ! Loop-private accumulators; avoid touching state%*_eu in the inner
      ! loop so the writes coalesce.
      stage_eu_acc = 0.0_wp
      xmom_eu_acc = 0.0_wp
      ymom_eu_acc = 0.0_wp
      speed_max_last = 0.0_wp

      do i = 1, 3
         ! Left state (this cell's edge values from the extrapolation pass).
         q_left(1) = state%barotropic%stage_ev(i, k)
         q_left(2) = state%barotropic%xmom_ev(i, k)
         q_left(3) = state%barotropic%ymom_ev(i, k)
         zl = state%barotropic%bed_ev(i, k)
         hle = state%barotropic%height_ev(i, k)

         ! Edge geometry.
         length = mesh%edge_len(i, k)
         n1 = mesh%edge_n1(i, k)
         n2 = mesh%edge_n2(i, k)

         neighbour = mesh%neighbours(i, k)

         if (neighbour < 0) then
            ! Boundary edge -- dispatch the edge's BC group. Group <= 0
            ! falls through to BC_WALL (inline wall mirror).
            bc_group = mesh%bc_group(i, k)
            if (bc_group <= 0) then
               bc_kind = BC_WALL
               bc_eta_ext = 0.0_wp
               bc_un_ext = 0.0_wp
               bc_h_ext = 0.0_wp
               bc_Qn_ext = 0.0_wp
            else
               bc_kind = state%bc_type(bc_group)
               bc_eta_ext = state%bc_eta_ext(bc_group)
               bc_un_ext = state%bc_un_ext(bc_group)
               bc_h_ext = state%bc_h_ext(bc_group)
               bc_Qn_ext = state%bc_Qn_ext(bc_group)
            end if
            if (bc_kind == BC_CHAPMAN) then
               ! Orlanski radiation: advect the persistent eta_old(i, k)
               ! along the estimated outward phase speed, then dispatch
               ! through the BC_TIDAL/Flather algebra with that ghost
               ! target.
               call compute_chapman_eta_target(state, mesh, k, i, &
                                               hle, bc_eta_ext)
               bc_kind = BC_TIDAL
            end if
            call compute_boundary_state_impl(q_left, n1, n2, zl, hle, bc_kind, &
                                             bc_eta_ext, bc_un_ext, bc_h_ext, &
                                             bc_Qn_ext, q_right, zr, hre)
         else
            ! Interior edge -- read the matching edge of the neighbour.
            m = mesh%neighbour_edges(i, k)
            q_right(1) = state%barotropic%stage_ev(m, neighbour)
            q_right(2) = state%barotropic%xmom_ev(m, neighbour)
            q_right(3) = state%barotropic%ymom_ev(m, neighbour)
            zr = state%barotropic%bed_ev(m, neighbour)
            hre = state%barotropic%height_ev(m, neighbour)
         end if

         ! KP draining: shared edge bed is the higher of the two sides;
         ! effective heights clipped at zero from below.
         z_half = max(zl, zr)
         h_left = max(hle + zl - z_half, 0.0_wp)
         h_right = max(hre + zr - z_half, 0.0_wp)

         if (h_left == 0.0_wp .and. h_right == 0.0_wp) then
            ! Both sides dry -- no flux, no pressure contribution.
            edgeflux(1) = 0.0_wp
            edgeflux(2) = 0.0_wp
            edgeflux(3) = 0.0_wp
            max_speed_local = 0.0_wp
            pressure_flux = 0.0_wp
         else
            call cu_flux_pointwise_impl(q_left, q_right, &
                                        h_left, h_right, &
                                        hle, hre, &
                                        n1, n2, &
                                        state%low_froude_mode, &
                                        edgeflux, max_speed_local, pressure_flux)
         end if

         ! Edge length * negate: accumulate -div(F) so the update reads
         ! `Q^{n+1} = Q^n + dt * stage_eu`.
         edgeflux(1) = -length*edgeflux(1)
         edgeflux(2) = -length*edgeflux(2)
         edgeflux(3) = -length*edgeflux(3)

         ! Stash the barotropic mass flux (net inflow into cell `k`
         ! across edge `i`, m^3/s, positive = gaining) for the tracer
         ! advection kernel to reuse (consistent-with-continuity).
         state%barotropic%mass_flux_e(i, k) = edgeflux(1)

         stage_eu_acc = stage_eu_acc + edgeflux(1)
         xmom_eu_acc = xmom_eu_acc + edgeflux(2)
         ymom_eu_acc = ymom_eu_acc + edgeflux(3)

         speed_max_last = max(speed_max_last, max_speed_local)

         ! Well-balanced bed-slope source: discrete hydrostatic-pressure
         ! work plus centroid-to-edge bed-slope term, plus `pressure_flux`
         ! from the central-upwind helper; telescopes to zero for lake at
         ! rest.
         pressuregrad_work = length*(-GRAVITY*0.5_wp* &
                                     (h_left*h_left - hle*hle &
                                      - (hle + hc)*(zl - zc)) &
                                     + pressure_flux)
         xmom_eu_acc = xmom_eu_acc - n1*pressuregrad_work
         ymom_eu_acc = ymom_eu_acc - n2*pressuregrad_work
      end do

      ! Normalise by cell area so the accumulators hold flux divergence
      ! per unit area, ready for the explicit update step.
      if (mesh%area(k) > 0.0_wp) then
         inv_area = 1.0_wp/mesh%area(k)
      else
         inv_area = 0.0_wp
      end if
      state%barotropic%stage_eu(k) = stage_eu_acc*inv_area
      state%barotropic%xmom_eu(k) = xmom_eu_acc*inv_area
      state%barotropic%ymom_eu(k) = ymom_eu_acc*inv_area

      state%barotropic%max_speed(k) = speed_max_last
   end subroutine flux_one_cell

   pure subroutine compute_chapman_eta_target(state, mesh, k, i, hle, eta_target)
      !$omp declare target
      !! Orlanski-style radiating ghost stage at boundary edge `i` of
      !! cell `k`. Advects the persistent `chapman_eta_old(i, k)` toward
      !! the interior centroid stage along the phase speed, writes back
      !! the new `eta_old`, and returns `eta_target` as the BC_TIDAL
      !! ghost eta. Per-edge writes to `chapman_eta_old(i, k)` are
      !! race-free: each (i, k) is touched only by cell `k`.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: k, i
      real(wp), intent(in) :: hle
      real(wp), intent(out) :: eta_target

      real(wp) :: eta_int, eta_at_edge, eta_old_val
      real(wp) :: deta_dt, deta_dn, c_phase, d_perp, edge_len, dt_local
      real(wp), parameter :: H_FLOOR = 1.0e-6_wp
      real(wp), parameter :: GRAD_FLOOR = 1.0e-10_wp

      eta_int = state%barotropic%stage_cv(k)
      eta_at_edge = state%barotropic%stage_ev(i, k)
      eta_old_val = state%barotropic%chapman_eta_old(i, k)
      dt_local = state%dt
      edge_len = mesh%edge_len(i, k)

      if (edge_len > 0.0_wp .and. mesh%area(k) > 0.0_wp) then
         d_perp = 2.0_wp*mesh%area(k)/(3.0_wp*edge_len)
      else
         d_perp = 0.0_wp
      end if

      if (dt_local > 0.0_wp .and. d_perp > 0.0_wp) then
         deta_dt = (eta_int - eta_old_val)/dt_local
         deta_dn = (eta_int - eta_at_edge)/d_perp  ! along inward normal
         if (abs(deta_dn) > GRAD_FLOOR) then
            c_phase = -deta_dt/deta_dn
         else
            c_phase = 0.0_wp
         end if
         c_phase = max(0.0_wp, min(c_phase, sqrt(GRAVITY*max(hle, H_FLOOR))))
         eta_target = eta_old_val - c_phase*dt_local/d_perp* &
                      (eta_old_val - eta_int)
      else
         ! First step or degenerate triangle -- seed with interior eta.
         eta_target = eta_int
      end if

      state%barotropic%chapman_eta_old(i, k) = eta_int
   end subroutine compute_chapman_eta_target

   ! Inlined pointwise helpers (`*_impl`): self-contained pure helpers
   ! in the same TU as `flux_one_cell` so NVHPC device codegen inlines
   ! them (cross-module `!$acc routine seq` would be a real call).

   pure subroutine rotate_to_normal_impl(q, n1, n2)
      !$omp declare target
      real(wp), intent(inout) :: q(3)
      real(wp), intent(in) :: n1, n2
      real(wp) :: q1, q2

      q1 = q(2)
      q2 = q(3)
      q(2) = n1*q1 + n2*q2
      q(3) = -n2*q1 + n1*q2
   end subroutine rotate_to_normal_impl

   pure subroutine compute_velocity_terms_impl(h, h_edge, uh_raw, vh_raw, &
                                               u, uh, v, vh)
      !$omp declare target
      real(wp), intent(in) :: h, h_edge, uh_raw, vh_raw
      real(wp), intent(out) :: u, uh, v, vh
      real(wp) :: inv_h_edge

      if (h_edge > 0.0_wp) then
         inv_h_edge = 1.0_wp/h_edge
         u = uh_raw*inv_h_edge
         uh = h*u
         v = vh_raw*inv_h_edge
         vh = h*inv_h_edge*vh_raw
      else
         u = 0.0_wp
         uh = 0.0_wp
         v = 0.0_wp
         vh = 0.0_wp
      end if
   end subroutine compute_velocity_terms_impl

   pure function compute_s_max_impl(u_left, u_right, c_left, c_right) result(s)
      !$omp declare target
      real(wp), intent(in) :: u_left, u_right, c_left, c_right
      real(wp) :: s

      s = max(u_left + c_left, u_right + c_right)
      if (s < 0.0_wp) s = 0.0_wp
   end function compute_s_max_impl

   pure function compute_s_min_impl(u_left, u_right, c_left, c_right) result(s)
      !$omp declare target
      real(wp), intent(in) :: u_left, u_right, c_left, c_right
      real(wp) :: s

      s = min(u_left - c_left, u_right - c_right)
      if (s > 0.0_wp) s = 0.0_wp
   end function compute_s_min_impl

   pure function compute_local_froude_impl(low_froude_mode, &
                                           u_left, u_right, &
                                           v_left, v_right, &
                                           c_left, c_right) result(fr)
      !$omp declare target
      integer, intent(in) :: low_froude_mode
      real(wp), intent(in) :: u_left, u_right, v_left, v_right, c_left, c_right
      real(wp) :: fr
      real(wp) :: numerator, denominator, ratio

      select case (low_froude_mode)
      case (LOW_FROUDE_1)
         numerator = u_right*u_right + u_left*u_left &
                     + v_right*v_right + v_left*v_left
         denominator = c_left*c_left + c_right*c_right + 1.0e-10_wp
         ratio = numerator/denominator
         fr = sqrt(max(0.001_wp, min(1.0_wp, ratio)))
      case (LOW_FROUDE_2)
         numerator = u_right*u_right + u_left*u_left &
                     + v_right*v_right + v_left*v_left
         denominator = c_left*c_left + c_right*c_right + 1.0e-10_wp
         ratio = sqrt(numerator/denominator)
         fr = sqrt(min(1.0_wp, 0.01_wp + max(ratio - 0.01_wp, 0.0_wp)))
      case default
         fr = 1.0_wp
      end select
   end function compute_local_froude_impl

   pure subroutine cu_flux_pointwise_impl( &
      q_left, q_right, &
      h_left, h_right, &
      hle, hre, &
      n1, n2, &
      low_froude_mode, &
      edgeflux, max_speed, pressure_flux)
      !$omp declare target
      real(wp), intent(in) :: q_left(3), q_right(3)
      real(wp), intent(in) :: h_left, h_right
      real(wp), intent(in) :: hle, hre
      real(wp), intent(in) :: n1, n2
      integer, intent(in) :: low_froude_mode
      real(wp), intent(out) :: edgeflux(3)
      real(wp), intent(out) :: max_speed
      real(wp), intent(out) :: pressure_flux

      real(wp) :: q_left_rot(3), q_right_rot(3)
      real(wp) :: u_left, uh_left, v_left, vh_left
      real(wp) :: u_right, uh_right, v_right, vh_right
      real(wp) :: c_left, c_right
      real(wp) :: s_max, s_min, denom, inv_denom, s_max_s_min
      real(wp) :: flux_left(3), flux_right(3)
      real(wp) :: local_fr

      q_left_rot(:) = q_left(:)
      q_right_rot(:) = q_right(:)
      call rotate_to_normal_impl(q_left_rot, n1, n2)
      call rotate_to_normal_impl(q_right_rot, n1, n2)

      call compute_velocity_terms_impl(h_left, hle, q_left_rot(2), q_left_rot(3), &
                                       u_left, uh_left, v_left, vh_left)
      call compute_velocity_terms_impl(h_right, hre, q_right_rot(2), q_right_rot(3), &
                                       u_right, uh_right, v_right, vh_right)

      c_left = sqrt(GRAVITY*h_left)
      c_right = sqrt(GRAVITY*h_right)
      s_max = compute_s_max_impl(u_left, u_right, c_left, c_right)
      s_min = compute_s_min_impl(u_left, u_right, c_left, c_right)

      local_fr = compute_local_froude_impl(low_froude_mode, &
                                           u_left, u_right, &
                                           v_left, v_right, &
                                           c_left, c_right)

      flux_left(1) = u_left*h_left
      flux_left(2) = u_left*uh_left
      flux_left(3) = u_left*vh_left

      flux_right(1) = u_right*h_right
      flux_right(2) = u_right*uh_right
      flux_right(3) = u_right*vh_right

      denom = s_max - s_min

      if (denom < FLUX_EPSILON) then
         edgeflux(1) = 0.0_wp
         edgeflux(2) = 0.0_wp
         edgeflux(3) = 0.0_wp
         max_speed = 0.0_wp
         pressure_flux = 0.5_wp*GRAVITY*0.5_wp*(h_left*h_left + h_right*h_right)
      else
         max_speed = max(s_max, -s_min)
         inv_denom = 1.0_wp/max(denom, FLUX_DENOM_FLOOR)
         s_max_s_min = s_max*s_min

         edgeflux(1) = (s_max*flux_left(1) - s_min*flux_right(1) &
                        + s_max_s_min*(h_right - h_left))*inv_denom

         edgeflux(2) = (s_max*flux_left(2) - s_min*flux_right(2) &
                        + local_fr*s_max_s_min*(uh_right - uh_left))*inv_denom

         edgeflux(3) = (s_max*flux_left(3) - s_min*flux_right(3) &
                        + local_fr*s_max_s_min*(vh_right - vh_left))*inv_denom

         pressure_flux = 0.5_wp*GRAVITY*(s_max*h_left*h_left &
                                         - s_min*h_right*h_right)*inv_denom

         call rotate_to_normal_impl(edgeflux, n1, -n2)
      end if
   end subroutine cu_flux_pointwise_impl

   pure subroutine compute_boundary_state_impl( &
      q_left, n1, n2, zl, hle, bc_type, &
      eta_ext, un_ext, h_ext, Qn_ext, &
      q_right, zr, hre)
      !$omp declare target
      real(wp), intent(in) :: q_left(3)
      real(wp), intent(in) :: n1, n2
      real(wp), intent(in) :: zl, hle
      integer, intent(in) :: bc_type
      real(wp), intent(in) :: eta_ext, un_ext, h_ext, Qn_ext
      real(wp), intent(out) :: q_right(3)
      real(wp), intent(out) :: zr
      real(wp), intent(out) :: hre

      real(wp) :: normal_mom
      real(wp) :: h_guard, n_mom, u_n_int, u_n_ghost, c_speed

      select case (bc_type)
      case (BC_OPEN)
         q_right(1) = q_left(1)
         q_right(2) = q_left(2)
         q_right(3) = q_left(3)
         zr = zl
         hre = hle

      case (BC_CLAMPED)
         q_right(1) = eta_ext
         q_right(2) = q_left(2)
         q_right(3) = q_left(3)
         zr = zl
         hre = max(eta_ext - zl, 0.0_wp)

      case (BC_TIDAL, BC_NESTED)
         ! BC_NESTED: smoke-level — Flather radiation driven by per-group
         ! prescribed elevation, identical to BC_TIDAL; parent-model
         ! coupling (live eta from a coarser grid) is deferred.
         zr = zl
         hre = max(eta_ext - zl, 0.0_wp)
         q_right(1) = eta_ext
         h_guard = max(hle, 1.0e-6_wp)
         if (hle > 1.0e-6_wp) then
            n_mom = q_left(2)*n1 + q_left(3)*n2
            u_n_int = n_mom/h_guard
            c_speed = sqrt(GRAVITY/h_guard)
            u_n_ghost = u_n_int + c_speed*(q_left(1) - eta_ext)
            q_right(2) = hre*u_n_ghost*n1 + (q_left(2) - n_mom*n1)
            q_right(3) = hre*u_n_ghost*n2 + (q_left(3) - n_mom*n2)
         else
            q_right(2) = q_left(2)
            q_right(3) = q_left(3)
         end if

      case (BC_INFLOW)
         q_right(1) = q_left(1)
         q_right(2) = hle*un_ext*n1
         q_right(3) = hle*un_ext*n2
         zr = zl
         hre = hle

      case (BC_DISCHARGE)
         q_right(1) = q_left(1)
         q_right(2) = Qn_ext*n1
         q_right(3) = Qn_ext*n2
         zr = zl
         hre = hle

      case (BC_OUTFLOW)
         zr = zl
         hre = h_ext
         q_right(1) = zr + h_ext
         q_right(2) = q_left(2)
         q_right(3) = q_left(3)

      case default
         ! BC_WALL (and any unknown type as a safe fallback).
         q_right(1) = q_left(1)
         normal_mom = q_left(2)*n1 + q_left(3)*n2
         q_right(2) = q_left(2) - 2.0_wp*normal_mom*n1
         q_right(3) = q_left(3) - 2.0_wp*normal_mom*n2
         zr = zl
         hre = hle
      end select
   end subroutine compute_boundary_state_impl

end module rki_kernel_flux_unstr
