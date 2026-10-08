!! Boundary condition application for ghost cells
module rki_boundary
   !! Fills ghost cells per boundary type: wall (mirror), open/tidal
   !! (Flather radiation), inflow/discharge (prescribed velocity/flux),
   !! clamped (Dirichlet eta), sponge (falls through to wall), Chapman
   !! (Orlanski radiation).
   use rki_constants, only: wp, GRAVITY
   use rki_config, only: MAX_TIDAL_CONSTITUENTS
   use rki_boundary_types, only: bc_params_t, BC_WALL, BC_OPEN, BC_TIDAL, &
                                 BC_NESTED, BC_INFLOW, BC_DISCHARGE, &
                                 BC_CLAMPED, BC_SPONGE, BC_CHAPMAN, &
                                 bc_type_from_string, TWO_PI
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: apply_boundary_conditions
   public :: apply_nh_boundary_conditions
   public :: bc_params_t
   public :: BC_WALL, BC_OPEN, BC_TIDAL, BC_NESTED
   public :: BC_INFLOW, BC_DISCHARGE, BC_CLAMPED, BC_SPONGE, BC_CHAPMAN
   public :: bc_type_from_string

   integer, parameter :: N_EDGES = 4
      !! Number of domain edges (west, east, south, north).

   ! Per-edge geometry descriptor for the unified boundary dispatch.
   type :: edge_t
      integer  :: bc_type      !! BC tag for this edge (bc%west / east / south / north)
      integer  :: n_tan        !! tangential extent (ny on W/E, nx on S/N)
      integer  :: ig0, jg0     !! first ghost cell (k=1, t=1)
      integer  :: di_n, dj_n   !! outward step per ghost layer k
      integer  :: di_t, dj_t   !! step along the edge per tangential index t
      real(wp) :: nrm_sign     !! +1 on high edges (E/N), -1 on low edges (W/S)
      logical  :: normal_is_u  !! normal momentum is hu (W/E) vs hv (S/N)
   end type edge_t

contains

   subroutine apply_boundary_conditions(h, hu, hv, b, grid, bc)
      !! Fill ghost cells based on boundary condition type
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: b(grid%nx_total, grid%ny_total)
         !! INOUT: wall/sponge edges mirror the bed into their ghost cells
         !! (see `bc_wall_edge`).  Open/tidal/inflow edges leave it untouched
         !! — there the domain genuinely continues past the boundary.
      type(bc_params_t), intent(inout) :: bc
         !! Boundary condition parameters

      integer :: k, j, i, n
      integer :: nx, ny, nghost
      integer :: i_first, i_last, j_first, j_last
      type(edge_t) :: edges(N_EDGES)
      real(wp) :: tidal_eta, h_int, eta_int, speed
      real(wp) :: inflow_vel, disch_flux, clamp_eta
      real(wp) :: c_phase, deta_dt, deta_dx, eta_ghost_new, eta_old_val, bc_dt, bc_dx, bc_dy
      integer :: n_tc
      real(wp) :: tc_amp(MAX_TIDAL_CONSTITUENTS)
      real(wp) :: tc_phase(MAX_TIDAL_CONSTITUENTS)
      real(wp) :: tc_omega(MAX_TIDAL_CONSTITUENTS)
      real(wp) :: bc_time

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      i_first = nghost + 1
      i_last = nx - nghost
      j_first = nghost + 1
      j_last = ny - nghost

      ! Per-edge geometry (west / east / south / north). W/E loop over j,
      ! normal momentum hu; S/N over i, normal momentum hv.
      edges(1) = edge_t(bc_type=bc%west, n_tan=ny, ig0=i_first - 1, jg0=1, &
                        di_n=-1, dj_n=0, di_t=0, dj_t=1, nrm_sign=-1.0_wp, normal_is_u=.true.)
      edges(2) = edge_t(bc_type=bc%east, n_tan=ny, ig0=i_last + 1, jg0=1, &
                        di_n=1, dj_n=0, di_t=0, dj_t=1, nrm_sign=1.0_wp, normal_is_u=.true.)
      edges(3) = edge_t(bc_type=bc%south, n_tan=nx, ig0=1, jg0=j_first - 1, &
                        di_n=0, dj_n=-1, di_t=1, dj_t=0, nrm_sign=-1.0_wp, normal_is_u=.false.)
      edges(4) = edge_t(bc_type=bc%north, n_tan=nx, ig0=1, jg0=j_last + 1, &
                        di_n=0, dj_n=1, di_t=1, dj_t=0, nrm_sign=1.0_wp, normal_is_u=.false.)

      n_tc = bc%n_tidal_constituents
      bc_time = bc%time
      tc_amp = bc%tidal_amp
      tc_phase = bc%tidal_phase
      tc_omega = bc%tidal_omega

      tidal_eta = 0.0_wp
      if (n_tc > 0) then
         ! eta(t) = sum A_n cos(omega_n*t - phi_n)
         do n = 1, n_tc
            tidal_eta = tidal_eta + tc_amp(n)*cos(tc_omega(n)*bc_time - tc_phase(n))
         end do
      else if (bc%tidal_period > 0.0_wp) then
         ! Single-constituent fallback
         tidal_eta = bc%tidal_amplitude*sin(TWO_PI*bc_time/bc%tidal_period)
      end if

      inflow_vel = bc%inflow_velocity
      disch_flux = bc%discharge_flux
      clamp_eta = bc%clamped_elevation
      bc_dt = bc%dt
      bc_dx = bc%dx
      bc_dy = bc%dy

      ! ---- West boundary ----
      select case (bc%west)
      case (BC_WALL, BC_SPONGE)
         call bc_wall_edge(h, hu, hv, b, grid, edges(1))

      case (BC_OPEN)
         call bc_flather_edge(h, hu, hv, b, grid, edges(1), 0.0_wp, .false.)

      case (BC_TIDAL)
         call bc_flather_edge(h, hu, hv, b, grid, edges(1), tidal_eta, .true.)

      case (BC_INFLOW)
         call bc_inflow_edge(h, hu, hv, b, grid, edges(1), inflow_vel)

      case (BC_DISCHARGE)
         call bc_disch_edge(h, hu, hv, b, grid, edges(1), disch_flux)

      case (BC_CLAMPED)
         call bc_clamped_edge(h, hu, hv, b, grid, edges(1), clamp_eta)

      case (BC_CHAPMAN)
         call bc_chapman_edge(h, hu, hv, b, grid, edges(1), bc%eta_old_west, bc_dt, bc_dx, bc_dy)

      case default
         continue
      end select

      ! ---- East boundary ----
      select case (bc%east)
      case (BC_WALL, BC_SPONGE)
         call bc_wall_edge(h, hu, hv, b, grid, edges(2))

      case (BC_OPEN)
         call bc_flather_edge(h, hu, hv, b, grid, edges(2), 0.0_wp, .false.)

      case (BC_TIDAL)
         call bc_flather_edge(h, hu, hv, b, grid, edges(2), tidal_eta, .true.)

      case (BC_INFLOW)
         call bc_inflow_edge(h, hu, hv, b, grid, edges(2), inflow_vel)

      case (BC_DISCHARGE)
         call bc_disch_edge(h, hu, hv, b, grid, edges(2), disch_flux)

      case (BC_CLAMPED)
         call bc_clamped_edge(h, hu, hv, b, grid, edges(2), clamp_eta)

      case (BC_CHAPMAN)
         call bc_chapman_edge(h, hu, hv, b, grid, edges(2), bc%eta_old_east, bc_dt, bc_dx, bc_dy)

      case default
         continue
      end select

      ! ---- South boundary ----
      select case (bc%south)
      case (BC_WALL, BC_SPONGE)
         call bc_wall_edge(h, hu, hv, b, grid, edges(3))

      case (BC_OPEN)
         call bc_flather_edge(h, hu, hv, b, grid, edges(3), 0.0_wp, .false.)

      case (BC_TIDAL)
         call bc_flather_edge(h, hu, hv, b, grid, edges(3), tidal_eta, .true.)

      case (BC_INFLOW)
         call bc_inflow_edge(h, hu, hv, b, grid, edges(3), inflow_vel)

      case (BC_DISCHARGE)
         call bc_disch_edge(h, hu, hv, b, grid, edges(3), disch_flux)

      case (BC_CLAMPED)
         call bc_clamped_edge(h, hu, hv, b, grid, edges(3), clamp_eta)

      case (BC_CHAPMAN)
         call bc_chapman_edge(h, hu, hv, b, grid, edges(3), bc%eta_old_south, bc_dt, bc_dx, bc_dy)

      case default
         continue
      end select

      ! ---- North boundary ----
      select case (bc%north)
      case (BC_WALL, BC_SPONGE)
         call bc_wall_edge(h, hu, hv, b, grid, edges(4))

      case (BC_OPEN)
         call bc_flather_edge(h, hu, hv, b, grid, edges(4), 0.0_wp, .false.)

      case (BC_TIDAL)
         call bc_flather_edge(h, hu, hv, b, grid, edges(4), tidal_eta, .true.)

      case (BC_INFLOW)
         call bc_inflow_edge(h, hu, hv, b, grid, edges(4), inflow_vel)

      case (BC_DISCHARGE)
         call bc_disch_edge(h, hu, hv, b, grid, edges(4), disch_flux)

      case (BC_CLAMPED)
         call bc_clamped_edge(h, hu, hv, b, grid, edges(4), clamp_eta)

      case (BC_CHAPMAN)
         call bc_chapman_edge(h, hu, hv, b, grid, edges(4), bc%eta_old_north, bc_dt, bc_dx, bc_dy)

      case default
         continue
      end select

   end subroutine apply_boundary_conditions

   pure subroutine bc_wall_edge(h, hu, hv, b, grid, ed)
      !! Wall / sponge ghost fill for one edge: mirror the BED and the
      !! surface, negate the NORMAL momentum, copy the tangential.
      !! Geometry comes from `ed`; its scalars are hoisted to locals so the
      !! do concurrent carries no derived-type dereference.
      !!
      !! ## Why the bed is mirrored too
      !!
      !! A wall IS a reflection: the ghost column is a fictitious mirror of
      !! the interior, so its bed must equal its mirror's bed.  Anything else
      !! describes a domain that CONTINUES past the wall with different
      !! topography — which is not a wall.
      !!
      !! Formula bathymetry setters fill the ghost rows by evaluating the
      !! formula there (they must, or the EOS sees h_layer = 0 at ghosts and
      !! falls back to rho_0), so without this line `b(ghost)` is the
      !! CONTINUATION of the topography rather than its reflection, and the
      !! reconstruction across the wall face is asymmetric.
      !!
      !! At rest that is harmless — the Audusse hydrostatic reconstruction
      !! still returns exactly zero flux, which is why the C-property held to
      !! 1e-12 over irregular bathymetry.  But WITH FLOW the asymmetry gives a
      !! non-zero HLL mass flux through a closed wall.  Measured on a closed
      !! box, stratified, sinusoidal bed, 400 steps: relative mass change
      !! 1.364e-6 with the formula bed in the ghosts vs -2.1e-14 with it
      !! mirrored, and the solver's own `mass_out` accumulator attributed
      !! 100 % of it to boundary flux (budget residual 2e-14).
      !!
      !! Mirroring here rather than once at setup means it also survives a bed
      !! that EVOLVES — which is exactly what the Exner update in the
      !! morphodynamics plan will do.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: b(grid%nx_total, grid%ny_total)
      type(edge_t), intent(in) :: ed

      integer :: t, k, ig, jg, im, jm, nghost
      integer :: ig0, jg0, din, djn, dit, djt, n_tan
      logical :: nu

      nghost = grid%nghost
      ig0 = ed%ig0
      jg0 = ed%jg0
      din = ed%di_n
      djn = ed%dj_n
      dit = ed%di_t
      djt = ed%dj_t
      n_tan = ed%n_tan
      nu = ed%normal_is_u

      do concurrent(t=1:n_tan, k=1:nghost) local(ig, jg, im, jm)
         ig = ig0 + (k - 1)*din + (t - 1)*dit       ! ghost cell
         jg = jg0 + (k - 1)*djn + (t - 1)*djt
         im = ig - (2*k - 1)*din                    ! mirror interior cell
         jm = jg - (2*k - 1)*djn
         ! Mirror the bed FIRST (see the note above), then the surface.  With
         ! b symmetric the eta-mirror below reduces to h(ig,jg) = h(im,jm);
         ! the general form is kept so the routine stays correct if a caller
         ! ever supplies an already-mirrored bed.
         b(ig, jg) = b(im, jm)
         h(ig, jg) = max(0.0_wp, h(im, jm) + b(im, jm) - b(ig, jg))
         if (nu) then
            hu(ig, jg) = -hu(im, jm)
            hv(ig, jg) = hv(im, jm)
         else
            hv(ig, jg) = -hv(im, jm)
            hu(ig, jg) = hu(im, jm)
         end if
      end do
   end subroutine bc_wall_edge

   pure subroutine bc_flather_edge(h, hu, hv, b, grid, ed, eta_ext, prescribe_h)
      !! Flather radiation for one edge.  Shared by OPEN (eta_ext = 0,
      !! prescribe_h = .false. -> zero-gradient ghost depth) and TIDAL
      !! (eta_ext = tidal_eta, prescribe_h = .true. -> ghost depth from
      !! the prescribed surface).  The NORMAL momentum radiates the
      !! surface anomaly outward (sign from ed%nrm_sign); the tangential
      !! momentum is copied from the first interior cell.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      type(edge_t), intent(in) :: ed
      real(wp), intent(in) :: eta_ext
      logical, intent(in) :: prescribe_h

      integer :: t, k, ig, jg, ii, ji, nghost
      integer :: ig0, jg0, din, djn, dit, djt, n_tan
      real(wp) :: sgn, h_int, eta_int, speed, hmax, h_ghost
      logical :: nu

      nghost = grid%nghost
      ig0 = ed%ig0
      jg0 = ed%jg0
      din = ed%di_n
      djn = ed%dj_n
      dit = ed%di_t
      djt = ed%dj_t
      n_tan = ed%n_tan
      nu = ed%normal_is_u
      sgn = ed%nrm_sign

      do concurrent(t=1:n_tan, k=1:nghost) &
         local(ig, jg, ii, ji, h_int, eta_int, speed, hmax, h_ghost)
         ig = ig0 + (k - 1)*din + (t - 1)*dit
         jg = jg0 + (k - 1)*djn + (t - 1)*djt
         ii = ig0 + (t - 1)*dit - din       ! first interior cell at this tangential index
         ji = jg0 + (t - 1)*djt - djn
         h_int = h(ii, ji)
         eta_int = h_int + b(ii, ji)
         hmax = max(h_int, 1.0e-6_wp)
         speed = sqrt(GRAVITY*hmax)
         h_ghost = merge(max(0.0_wp, eta_ext - b(ig, jg)), h_int, prescribe_h)
         h(ig, jg) = h_ghost
         if (nu) then
            hu(ig, jg) = h_ghost*(hu(ii, ji)/hmax + sgn*speed*(eta_int - eta_ext)/hmax)
            hv(ig, jg) = hv(ii, ji)
         else
            hv(ig, jg) = h_ghost*(hv(ii, ji)/hmax + sgn*speed*(eta_int - eta_ext)/hmax)
            hu(ig, jg) = hu(ii, ji)
         end if
      end do
   end subroutine bc_flather_edge

   pure subroutine bc_inflow_edge(h, hu, hv, b, grid, ed, vel)
      !! Inflow: zero-gradient ghost depth; the NORMAL momentum is the
      !! prescribed inward transport h_ghost*vel (inward = -ed%nrm_sign);
      !! tangential momentum copied from the first interior cell.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      type(edge_t), intent(in) :: ed
      real(wp), intent(in) :: vel

      integer :: t, k, ig, jg, ii, ji, nghost
      integer :: ig0, jg0, din, djn, dit, djt, n_tan
      real(wp) :: in_sgn, h_int
      logical :: nu

      nghost = grid%nghost
      ig0 = ed%ig0
      jg0 = ed%jg0
      din = ed%di_n
      djn = ed%dj_n
      dit = ed%di_t
      djt = ed%dj_t
      n_tan = ed%n_tan
      nu = ed%normal_is_u
      in_sgn = -ed%nrm_sign

      do concurrent(t=1:n_tan, k=1:nghost) local(ig, jg, ii, ji, h_int)
         ig = ig0 + (k - 1)*din + (t - 1)*dit
         jg = jg0 + (k - 1)*djn + (t - 1)*djt
         ii = ig0 + (t - 1)*dit - din
         ji = jg0 + (t - 1)*djt - djn
         h_int = h(ii, ji)
         h(ig, jg) = h_int
         if (nu) then
            hu(ig, jg) = in_sgn*h_int*vel
            hv(ig, jg) = hv(ii, ji)
         else
            hv(ig, jg) = in_sgn*h_int*vel
            hu(ig, jg) = hu(ii, ji)
         end if
      end do
   end subroutine bc_inflow_edge

   pure subroutine bc_disch_edge(h, hu, hv, b, grid, ed, flux)
      !! Discharge: zero-gradient ghost depth; the NORMAL momentum (the
      !! depth-integrated transport) is set directly to the prescribed
      !! inward flux (inward = -ed%nrm_sign); tangential momentum copied.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      type(edge_t), intent(in) :: ed
      real(wp), intent(in) :: flux

      integer :: t, k, ig, jg, ii, ji, nghost
      integer :: ig0, jg0, din, djn, dit, djt, n_tan
      real(wp) :: in_sgn
      logical :: nu

      nghost = grid%nghost
      ig0 = ed%ig0
      jg0 = ed%jg0
      din = ed%di_n
      djn = ed%dj_n
      dit = ed%di_t
      djt = ed%dj_t
      n_tan = ed%n_tan
      nu = ed%normal_is_u
      in_sgn = -ed%nrm_sign

      do concurrent(t=1:n_tan, k=1:nghost) local(ig, jg, ii, ji)
         ig = ig0 + (k - 1)*din + (t - 1)*dit
         jg = jg0 + (k - 1)*djn + (t - 1)*djt
         ii = ig0 + (t - 1)*dit - din
         ji = jg0 + (t - 1)*djt - djn
         h(ig, jg) = h(ii, ji)
         if (nu) then
            hu(ig, jg) = in_sgn*flux
            hv(ig, jg) = hv(ii, ji)
         else
            hv(ig, jg) = in_sgn*flux
            hu(ig, jg) = hu(ii, ji)
         end if
      end do
   end subroutine bc_disch_edge

   pure subroutine bc_clamped_edge(h, hu, hv, b, grid, ed, eta)
      !! Clamped (Dirichlet surface): prescribe the ghost depth from the
      !! clamped elevation; copy BOTH momenta from the first interior
      !! cell (no normal/tangential distinction).
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      type(edge_t), intent(in) :: ed
      real(wp), intent(in) :: eta

      integer :: t, k, ig, jg, ii, ji, nghost
      integer :: ig0, jg0, din, djn, dit, djt, n_tan

      nghost = grid%nghost
      ig0 = ed%ig0
      jg0 = ed%jg0
      din = ed%di_n
      djn = ed%dj_n
      dit = ed%di_t
      djt = ed%dj_t
      n_tan = ed%n_tan

      do concurrent(t=1:n_tan, k=1:nghost) local(ig, jg, ii, ji)
         ig = ig0 + (k - 1)*din + (t - 1)*dit
         jg = jg0 + (k - 1)*djn + (t - 1)*djt
         ii = ig0 + (t - 1)*dit - din
         ji = jg0 + (t - 1)*djt - djn
         h(ig, jg) = max(0.0_wp, eta - b(ig, jg))
         hu(ig, jg) = hu(ii, ji)
         hv(ig, jg) = hv(ii, ji)
      end do
   end subroutine bc_clamped_edge

   pure subroutine bc_chapman_edge(h, hu, hv, b, grid, ed, eta_old, dt, dx, dy)
      !! Chapman / Orlanski radiation for one edge.  The ghost surface
      !! advects at a clamped, non-negative phase speed estimated from the
      !! interior eta tendency and along-normal gradient; the momentum is
      !! the Flather radiation about that surface (sign ed%nrm_sign).
      !! `eta_old` is this edge's persistent surface (length n_tan),
      !! refreshed to the interior eta on the k==1 layer (per-t writes are
      !! disjoint).  `dd = dx` on W/E, `dy` on S/N.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(inout) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      type(edge_t), intent(in) :: ed
      ! assumed-shape-ok: eta_old is a 1D edge-count-sized array; its length
      ! varies per boundary edge and is not a grid-level nx/ny dimension.
      real(wp), intent(inout) :: eta_old(:)
      real(wp), intent(in) :: dt, dx, dy

      integer :: t, k, ig, jg, ii, ji, i2, j2, nx, ny, nghost
      integer :: ig0, jg0, din, djn, dit, djt, n_tan
      real(wp) :: sgn, dd, h_int, eta_int, eta_old_val, speed, hmax
      real(wp) :: deta_dt, deta_dx, c_phase, eta_ghost_new, h_ghost
      logical :: nu

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      ig0 = ed%ig0
      jg0 = ed%jg0
      din = ed%di_n
      djn = ed%dj_n
      dit = ed%di_t
      djt = ed%dj_t
      n_tan = ed%n_tan
      nu = ed%normal_is_u
      sgn = ed%nrm_sign
      dd = merge(dx, dy, nu)

      do concurrent(t=1:n_tan, k=1:nghost) &
         local(ig, jg, ii, ji, i2, j2, h_int, eta_int, eta_old_val, speed, &
               hmax, deta_dt, deta_dx, c_phase, eta_ghost_new, h_ghost)
         ig = ig0 + (k - 1)*din + (t - 1)*dit
         jg = jg0 + (k - 1)*djn + (t - 1)*djt
         ii = ig0 + (t - 1)*dit - din          ! first interior cell
         ji = jg0 + (t - 1)*djt - djn
         i2 = ii - din                          ! second interior (along normal)
         j2 = ji - djn
         h_int = h(ii, ji)
         eta_int = h_int + b(ii, ji)
         eta_old_val = eta_old(t)
         hmax = max(h_int, 1.0e-6_wp)
         speed = sqrt(GRAVITY*hmax)
         if (dt > 0.0_wp) then
            deta_dt = (eta_int - eta_old_val)/dt
            if (i2 >= 1 .and. i2 <= nx .and. j2 >= 1 .and. j2 <= ny) then
               deta_dx = (h(i2, j2) + b(i2, j2) - eta_int)/dd
            else
               deta_dx = 0.0_wp
            end if
            if (abs(deta_dx) > 1.0e-10_wp) then
               c_phase = -deta_dt/deta_dx
            else
               c_phase = 0.0_wp
            end if
            c_phase = max(0.0_wp, min(c_phase, speed))
            eta_ghost_new = eta_old_val - c_phase*dt/dd*(eta_old_val - eta_int)
         else
            eta_ghost_new = eta_int
         end if
         if (k == 1) eta_old(t) = eta_int
         h_ghost = max(0.0_wp, eta_ghost_new - b(ig, jg))
         h(ig, jg) = h_ghost
         if (nu) then
            hu(ig, jg) = h_ghost*(hu(ii, ji)/hmax + sgn*speed*(eta_int - eta_ghost_new)/hmax)
            hv(ig, jg) = hv(ii, ji)
         else
            hv(ig, jg) = h_ghost*(hv(ii, ji)/hmax + sgn*speed*(eta_int - eta_ghost_new)/hmax)
            hu(ig, jg) = hu(ii, ji)
         end if
      end do
   end subroutine bc_chapman_edge

   pure subroutine apply_nh_boundary_conditions(q_nh, u_layer, v_layer, w, h, b, &
                                                grid, nz, bc)
      !! Fill ghost cells for non-hydrostatic variables
      !!
      !! q_nh: zero-gradient (Neumann) at all boundaries
      !! u_layer, v_layer: mirror/negate at walls, extrapolate at open
      !! w: zero-gradient laterally
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: q_nh(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: u_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: v_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: b(grid%nx_total, grid%ny_total)
      type(bc_params_t), intent(in) :: bc

      integer :: i, j, k, g, nx, ny, nghost
      integer :: i_first, i_last, j_first, j_last
      integer :: bc_w, bc_e, bc_s, bc_n

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      i_first = nghost + 1
      i_last = nx - nghost
      j_first = nghost + 1
      j_last = ny - nghost

      ! Extract bc types into local scalars for GPU accessibility
      bc_w = bc%west
      bc_e = bc%east
      bc_s = bc%south
      bc_n = bc%north

      ! ---- West boundary ----
      do concurrent(k=1:nz, j=1:ny, g=1:nghost)
         select case (bc_w)
         case (BC_WALL)
            q_nh(i_first - g, j, k) = q_nh(i_first + g - 1, j, k)
            u_layer(i_first - g, j, k) = -u_layer(i_first + g - 1, j, k)
            v_layer(i_first - g, j, k) = v_layer(i_first + g - 1, j, k)
         case default
            q_nh(i_first - g, j, k) = q_nh(i_first + g - 1, j, k)
            u_layer(i_first - g, j, k) = u_layer(i_first + g - 1, j, k)
            v_layer(i_first - g, j, k) = v_layer(i_first + g - 1, j, k)
         end select
      end do
      do concurrent(k=1:nz + 1, j=1:ny, g=1:nghost)
         w(i_first - g, j, k) = w(i_first + g - 1, j, k)
      end do

      ! ---- East boundary ----
      do concurrent(k=1:nz, j=1:ny, g=1:nghost)
         select case (bc_e)
         case (BC_WALL)
            q_nh(i_last + g, j, k) = q_nh(i_last - g + 1, j, k)
            u_layer(i_last + g, j, k) = -u_layer(i_last - g + 1, j, k)
            v_layer(i_last + g, j, k) = v_layer(i_last - g + 1, j, k)
         case default
            q_nh(i_last + g, j, k) = q_nh(i_last - g + 1, j, k)
            u_layer(i_last + g, j, k) = u_layer(i_last - g + 1, j, k)
            v_layer(i_last + g, j, k) = v_layer(i_last - g + 1, j, k)
         end select
      end do
      do concurrent(k=1:nz + 1, j=1:ny, g=1:nghost)
         w(i_last + g, j, k) = w(i_last - g + 1, j, k)
      end do

      ! ---- South boundary ----
      do concurrent(k=1:nz, i=1:nx, g=1:nghost)
         select case (bc_s)
         case (BC_WALL)
            q_nh(i, j_first - g, k) = q_nh(i, j_first + g - 1, k)
            u_layer(i, j_first - g, k) = u_layer(i, j_first + g - 1, k)
            v_layer(i, j_first - g, k) = -v_layer(i, j_first + g - 1, k)
         case default
            q_nh(i, j_first - g, k) = q_nh(i, j_first + g - 1, k)
            u_layer(i, j_first - g, k) = u_layer(i, j_first + g - 1, k)
            v_layer(i, j_first - g, k) = v_layer(i, j_first + g - 1, k)
         end select
      end do
      do concurrent(k=1:nz + 1, i=1:nx, g=1:nghost)
         w(i, j_first - g, k) = w(i, j_first + g - 1, k)
      end do

      ! ---- North boundary ----
      do concurrent(k=1:nz, i=1:nx, g=1:nghost)
         select case (bc_n)
         case (BC_WALL)
            q_nh(i, j_last + g, k) = q_nh(i, j_last - g + 1, k)
            u_layer(i, j_last + g, k) = u_layer(i, j_last - g + 1, k)
            v_layer(i, j_last + g, k) = -v_layer(i, j_last - g + 1, k)
         case default
            q_nh(i, j_last + g, k) = q_nh(i, j_last - g + 1, k)
            u_layer(i, j_last + g, k) = u_layer(i, j_last - g + 1, k)
            v_layer(i, j_last + g, k) = v_layer(i, j_last - g + 1, k)
         end select
      end do
      do concurrent(k=1:nz + 1, i=1:nx, g=1:nghost)
         w(i, j_last + g, k) = w(i, j_last - g + 1, k)
      end do

   end subroutine apply_nh_boundary_conditions

end module rki_boundary
