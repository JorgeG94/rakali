!! Per-layer scalar-tracer transport for the structured multilayer solver.
module rki_ml_tracers
   !! Horizontal advection (PPM and first-order upwind), vertical
   !! exchange and diffusion, distribution under the barotropic step,
   !! horizontal Laplacian diffusion, surface fluxes, RK2 plumbing,
   !! ghost-cell boundary application, and bounds clamping.
   !!
   !! Both prognostic tracers (salinity hS, temperature hT) share these
   !! generic kernels — callers pick which array to pass.
   !!
   !! The advection / PPM / horizontal-diffusion kernels are `pure`: the
   !! per-kernel scratch (PPM reconstruction snapshots, Jacobi neighbour
   !! snapshots) is caller-owned through `tracer_scratch_t`, allocated
   !! once up front by the solver (`tracer_scratch_ensure`) and released
   !! by `tracer_scratch_destroy`.  Writing an `intent(out)` scratch dummy
   !! is pure-legal, so the kernels keep purity while still avoiding the
   !! per-step allocate/deallocate that would stall the GPU pipeline.
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            RHO_WATER, LATENT_HEAT_VAPORIZATION
#else
   use rki_constants, only: NZ_STACK_MAX, wp, GRAVITY, DRY_TOLERANCE, THIN_LAYER_THRESHOLD, &
                            RHO_WATER, LATENT_HEAT_VAPORIZATION
#endif
   use rki_grid, only: hgrid_t
   use rki_boundary_types, only: bc_params_t, BC_INFLOW, BC_TIDAL, BC_CLAMPED
   use rki_ml_common, only: ML_MIN_DEPTH, VANISHING_LAYER_TOL
   use rki_remap_column, only: remap_column_ppm
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: ml_advect_tracer
   public :: ml_advect_tracer_ppm
   public :: ml_distribute_tracer
   public :: ml_vertical_exchange_tracer
   public :: ml_init_tracer
   public :: ml_init_tracer_stratified
   public :: ml_seed_tracer
   public :: ml_save_tracer
   public :: ml_rk2_average_tracer
   public :: ml_apply_tracer_bcs
   public :: ml_clamp_tracer
   public :: ml_horizontal_diffusion_tracer
   public :: ml_apply_surface_tracer_fluxes
   public :: ml_vertical_diffusion_tracer
   public :: ml_settling_tracer
   public :: ml_erode_bed_sediment

   ! Caller-owned per-kernel scratch (keeps the kernels `pure`).
   public :: tracer_scratch_t
   public :: tracer_scratch_ensure
   public :: tracer_scratch_destroy

   type :: tracer_scratch_t
      !! Persistent per-kernel scratch owned by the caller (the solver, or
      !! a unit test).  Allocated once up front to (nx,ny,nz), device-
      !! attached, and reused for every kernel call so the `pure` kernels
      !! never allocate.  Each group is guarded by a `want_*` flag in
      !! `tracer_scratch_ensure`, so a PPM-only run does not allocate the
      !! upwind snapshot and vice-versa.

      ! ml_advect_tracer (first-order upwind) — Jacobi neighbour snapshot.
      real(wp), allocatable :: adv_hTr_snap(:, :, :)

      ! ml_advect_tracer_ppm — reconstruction snapshots + limited edges.
      real(wp), allocatable :: ppm_tr_snap(:, :, :)
      real(wp), allocatable :: ppm_ql_x(:, :, :)
      real(wp), allocatable :: ppm_qr_x(:, :, :)
      real(wp), allocatable :: ppm_ql_y(:, :, :)
      real(wp), allocatable :: ppm_qr_y(:, :, :)

      ! ml_advect_tracer_weno (WENO5/7/9 reconstruction) — concentration snapshot.
      real(wp), allocatable :: weno_tr_snap(:, :, :)

      ! ml_horizontal_diffusion_tracer — Jacobi pre-step snapshot.
      real(wp), allocatable :: hdiff_hTr_pre(:, :, :)
   end type tracer_scratch_t

contains

   subroutine ml_advect_tracer(h_layer, hu_layer, hv_layer, &
                                    mass_flux_x, mass_flux_y, &
                                    hTr_layer, adv_hTr_snap, &
                                    grid, nz, dt)
      !! Per-layer horizontal advection of a scalar tracer hTr = h*Tr
      !! using first-order upwind reconstruction.
      !!
      !! Per-face per-layer mass flux uses centred per-layer velocity
      !! with a vanishing-layer GUARD that excludes layers where either
      !! side has h_layer(k) <= VANISHING_LAYER_TOL from the per-layer
      !! split, redistributing the residual across the active layers
      !! via a per-face symmetric frac:
      !!
      !!   ! For each face (one cell-side shown, the other is symmetric):
      !!   if h_self(k) > tol .and. h_nbr(k) > tol then
      !!      v_face(k)         = 0.5 * (hu_self(k) + hu_nbr(k))
      !!      h_face_active(k)  = 0.5 * (h_self(k)  + h_nbr(k))
      !!   else
      !!      v_face(k)         = 0
      !!      h_face_active(k)  = 0     ! vanishing layer carries no flux
      !!   end if
      !!   H_active   = sum_k h_face_active(k)
      !!   corr       = mass_flux_face - sum_k v_face(k)
      !!   frac(k)    = h_face_active(k) / H_active   ! 0 if vanishing
      !!   f_layer(k) = v_face(k) + frac(k) * corr
      !!
      !! `v_face`, `h_face_active`, and `frac_face_active` are all
      !! computed from `0.5*(self + nbr)` quantities, so they are
      !! IDENTICAL when computed from cell A's right-face perspective
      !! versus cell A+1's left-face perspective → per-layer face mass
      !! conservation holds.  `sum_k(f_layer) = sum_k v_face + corr =
      !! mass_flux_face` (CWC).
      !!
      !! For sigma / zsigma / zstar-lite where every layer is non-
      !! vanishing on both sides, `frac_active(k) = h_face(k)/H_face =
      !! 1/nz` for uniform sigma, and the formula degenerates to the
      !! original centred CWC scheme — zero risk to existing tests.
      !!
      !! For ZSTAR_FULL deep-shallow interfaces, vanishing bed layers
      !! in the shallow cell are excluded from the per-layer split, so
      !! the column flux is redistributed across the (active) non-
      !! vanishing layers instead of being centred-averaged into the
      !! tiny vanishing receivers.  No bloat of S in vanishing layers,
      !! no clamp clipping, no salt-mass leak.
      !!
      !! History: the original formulation used `v_centred = 0.5 *
      !! (hu_self + hu_nbr)` plus a frac_self-weighted correction.
      !! Under sharp h_layer mismatches between cells (ZSTAR_FULL deep-
      !! shallow faces) the centred average over-allocated flux to
      !! vanishing layers, producing per-layer S drift that the
      !! [S_min, S_max] clamp then destroyed (~4 % salt-mass leak on
      !! the LE 4 h tidal validation, 2026-05-11).
      !!
      !! Jacobi snapshot: all neighbour and own-cell tracer-value reads
      !! (tr_here, tr_up for each of the four upwind faces) are taken
      !! from `adv_hTr_snap`, a caller-owned scratch dummy that this
      !! kernel fills with `hTr_layer` (frozen at entry via a full-array
      !! `do concurrent`) before the main kernel.
      !! Without this, the read-neighbour / write-own pattern inside
      !! `do concurrent` is a race (UB): on sequential backends (gfortran
      !! CPU, nvfortran serial) later iterations see already-updated
      !! neighbour values, so the two cell-visits of a shared face
      !! disagree and ~1e-7 tracer mass leaks per step.  The snapshot
      !! gives Jacobi semantics: CPU-serial and GPU-parallel launches
      !! apply the same frozen input and produce bit-identical,
      !! conservative results.  Mirrors the identical fix in
      !! `ml_horizontal_diffusion_tracer` (hdiff_hTr_pre).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: mass_flux_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: mass_flux_y(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: adv_hTr_snap(grid%nx_total, grid%ny_total, nz)
         !! Caller-owned Jacobi snapshot scratch (filled then read here).
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: h_val, tr_here, h_col
      real(wp) :: f_layer, tr_up
      real(wp) :: flux_x_left, flux_x_right, flux_y_bot, flux_y_top
      real(wp) :: inv_dx, inv_dy
      real(wp) :: corr_xL, corr_xR, corr_yB, corr_yT
      real(wp) :: H_active_xL, H_active_xR, H_active_yB, H_active_yT
      real(wp) :: inv_H_xL, inv_H_xR, inv_H_yB, inv_H_yT
      real(wp) :: v_xL(NZ_STACK_MAX), v_xR(NZ_STACK_MAX)
      real(wp) :: v_yB(NZ_STACK_MAX), v_yT(NZ_STACK_MAX)
      real(wp) :: hf_xL(NZ_STACK_MAX), hf_xR(NZ_STACK_MAX)
      real(wp) :: hf_yB(NZ_STACK_MAX), hf_yT(NZ_STACK_MAX)
      real(wp) :: h_iL, h_iR, h_jB, h_jT

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      inv_dx = 1.0_wp/dx
      inv_dy = 1.0_wp/dy

      ! Jacobi snapshot: freeze hTr_layer before the main kernel so that
      ! all neighbour reads come from the pre-step values.  See docstring.
      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         adv_hTr_snap(i, j, k) = hTr_layer(i, j, k)
      end do
      end do
      end do
      !$omp end parallel do

      !$omp parallel do collapse(2) private(k, kk, h_val, tr_here, h_col, f_layer, tr_up, flux_x_left, flux_x_right, flux_y_bot, flux_y_top, corr_xL, corr_xR, corr_yB, corr_yT, H_active_xL, H_active_xR, H_active_yB, H_active_yT, inv_H_xL, inv_H_xR, inv_H_yB, inv_H_yT, v_xL, v_xR, v_yB, v_yT, hf_xL, hf_xR, hf_yB, hf_yT, h_iL, h_iR, h_jB, h_jT)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then

            ! Vanishing-layer-gated centred CWC.  For each face, layers
            ! where either side has h_layer(k) <= VANISHING_LAYER_TOL
            ! contribute zero per-layer flux; the residual is
            ! redistributed across active layers via a per-face
            ! symmetric frac (h_face_active(k) / sum h_face_active).
            ! See kernel docstring.
            corr_xL = mass_flux_x(i - 1, j)
            corr_xR = mass_flux_x(i, j)
            corr_yB = mass_flux_y(i, j - 1)
            corr_yT = mass_flux_y(i, j)
            H_active_xL = 0.0_wp
            H_active_xR = 0.0_wp
            H_active_yB = 0.0_wp
            H_active_yT = 0.0_wp
            do k = 1, nz
               h_val = h_layer(i, j, k)
               h_iL = h_layer(i - 1, j, k)
               h_iR = h_layer(i + 1, j, k)
               h_jB = h_layer(i, j - 1, k)
               h_jT = h_layer(i, j + 1, k)

               if (h_val > VANISHING_LAYER_TOL .and. h_iL > VANISHING_LAYER_TOL) then
                  v_xL(k) = 0.5_wp*(hu_layer(i - 1, j, k) + hu_layer(i, j, k))
                  hf_xL(k) = 0.5_wp*(h_iL + h_val)
               else
                  v_xL(k) = 0.0_wp
                  hf_xL(k) = 0.0_wp
               end if
               if (h_val > VANISHING_LAYER_TOL .and. h_iR > VANISHING_LAYER_TOL) then
                  v_xR(k) = 0.5_wp*(hu_layer(i, j, k) + hu_layer(i + 1, j, k))
                  hf_xR(k) = 0.5_wp*(h_val + h_iR)
               else
                  v_xR(k) = 0.0_wp
                  hf_xR(k) = 0.0_wp
               end if
               if (h_val > VANISHING_LAYER_TOL .and. h_jB > VANISHING_LAYER_TOL) then
                  v_yB(k) = 0.5_wp*(hv_layer(i, j - 1, k) + hv_layer(i, j, k))
                  hf_yB(k) = 0.5_wp*(h_jB + h_val)
               else
                  v_yB(k) = 0.0_wp
                  hf_yB(k) = 0.0_wp
               end if
               if (h_val > VANISHING_LAYER_TOL .and. h_jT > VANISHING_LAYER_TOL) then
                  v_yT(k) = 0.5_wp*(hv_layer(i, j, k) + hv_layer(i, j + 1, k))
                  hf_yT(k) = 0.5_wp*(h_val + h_jT)
               else
                  v_yT(k) = 0.0_wp
                  hf_yT(k) = 0.0_wp
               end if

               corr_xL = corr_xL - v_xL(k)
               corr_xR = corr_xR - v_xR(k)
               corr_yB = corr_yB - v_yB(k)
               corr_yT = corr_yT - v_yT(k)
               H_active_xL = H_active_xL + hf_xL(k)
               H_active_xR = H_active_xR + hf_xR(k)
               H_active_yB = H_active_yB + hf_yB(k)
               H_active_yT = H_active_yT + hf_yT(k)
            end do

            if (H_active_xL > 0.0_wp) then
               inv_H_xL = 1.0_wp/H_active_xL
            else
               inv_H_xL = 0.0_wp
            end if
            if (H_active_xR > 0.0_wp) then
               inv_H_xR = 1.0_wp/H_active_xR
            else
               inv_H_xR = 0.0_wp
            end if
            if (H_active_yB > 0.0_wp) then
               inv_H_yB = 1.0_wp/H_active_yB
            else
               inv_H_yB = 0.0_wp
            end if
            if (H_active_yT > 0.0_wp) then
               inv_H_yT = 1.0_wp/H_active_yT
            else
               inv_H_yT = 0.0_wp
            end if

            do k = 1, nz
               h_val = h_layer(i, j, k)
               if (h_val > DRY_TOLERANCE) then
                  tr_here = adv_hTr_snap(i, j, k)/h_val
               else
                  tr_here = 0.0_wp
               end if

               ! Left face (i-1/2)
               f_layer = v_xL(k) + hf_xL(k)*inv_H_xL*corr_xL
               if (f_layer > 0.0_wp) then
                  if (h_layer(i - 1, j, k) > DRY_TOLERANCE) then
                     tr_up = adv_hTr_snap(i - 1, j, k)/h_layer(i - 1, j, k)
                  else
                     tr_up = 0.0_wp
                  end if
               else
                  tr_up = tr_here
               end if
               flux_x_left = f_layer*tr_up

               ! Right face (i+1/2)
               f_layer = v_xR(k) + hf_xR(k)*inv_H_xR*corr_xR
               if (f_layer > 0.0_wp) then
                  tr_up = tr_here
               else
                  if (h_layer(i + 1, j, k) > DRY_TOLERANCE) then
                     tr_up = adv_hTr_snap(i + 1, j, k)/h_layer(i + 1, j, k)
                  else
                     tr_up = 0.0_wp
                  end if
               end if
               flux_x_right = f_layer*tr_up

               ! Bottom face (j-1/2)
               f_layer = v_yB(k) + hf_yB(k)*inv_H_yB*corr_yB
               if (f_layer > 0.0_wp) then
                  if (h_layer(i, j - 1, k) > DRY_TOLERANCE) then
                     tr_up = adv_hTr_snap(i, j - 1, k)/h_layer(i, j - 1, k)
                  else
                     tr_up = 0.0_wp
                  end if
               else
                  tr_up = tr_here
               end if
               flux_y_bot = f_layer*tr_up

               ! Top face (j+1/2)
               f_layer = v_yT(k) + hf_yT(k)*inv_H_yT*corr_yT
               if (f_layer > 0.0_wp) then
                  tr_up = tr_here
               else
                  if (h_layer(i, j + 1, k) > DRY_TOLERANCE) then
                     tr_up = adv_hTr_snap(i, j + 1, k)/h_layer(i, j + 1, k)
                  else
                     tr_up = 0.0_wp
                  end if
               end if
               flux_y_top = f_layer*tr_up

               hTr_layer(i, j, k) = hTr_layer(i, j, k) &
                                    - dt*((flux_x_right - flux_x_left)*inv_dx &
                                          + (flux_y_top - flux_y_bot)*inv_dy)
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine ml_advect_tracer

   subroutine ml_advect_tracer_ppm(h_layer, hu_layer, hv_layer, &
                                        mass_flux_x, mass_flux_y, &
                                        hTr_layer, &
                                        ppm_tr_snap, ppm_ql_x, ppm_qr_x, &
                                        ppm_ql_y, ppm_qr_y, &
                                        grid, nz, dt)
      !! Per-layer horizontal advection of a scalar tracer using the
      !! Piecewise Parabolic Method (Colella & Woodward 1984).
      !!
      !! Two-phase unsplit approach:
      !!   Phase 1 — PPM reconstruction: 4th-order edge interpolation with
      !!             Colella-Woodward monotonicity limiting.
      !!   Phase 2 — Flux + update: time-averaged parabolic integration
      !!             over the CFL-swept region at each face.
      !!
      !! Falls back to piecewise-constant (first-order upwind) at ghost
      !! cells where the 5-point stencil extends outside the domain.
      !!
      !! Per-face mass flux uses the same vanishing-layer-gated centred
      !! CWC split as `ml_advect_tracer` — see that kernel's docstring
      !! for the algebra and the ZSTAR_FULL deep-shallow motivation.
      !!
      !! `h_face(k) = 0.5 * (h_layer(here) + h_layer(nbr))` is symmetric
      !! across the face: the PPM CFL `sigma = u_eff*dt/dx` (with
      !! `u_eff = f_layer / h_face`) drives the parabolic-profile
      !! integration, and an asymmetric `h_face` would make the two
      !! cells at a face see different `sigma` and integrate the
      !! upwind profile to different `tr_avg_up` values, leaking mass
      !! for any non-uniform tracer profile (the nz=1 closed-basin
      !! mass conservation test at 1e-12 catches this).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: mass_flux_x(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: mass_flux_y(grid%nx_total, grid%ny_total)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: ppm_tr_snap(grid%nx_total, grid%ny_total, nz)
         !! Caller-owned concentration snapshot scratch.
      real(wp), intent(out) :: ppm_ql_x(grid%nx_total, grid%ny_total, nz)
         !! Caller-owned limited left/right edge scratch (x/y directions).
      real(wp), intent(out) :: ppm_qr_x(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: ppm_ql_y(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: ppm_qr_y(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: inv_dx, inv_dy, h_col, h_val

      ! PPM reconstruction variables
      real(wp) :: q_im2, q_im1, q_i, q_ip1, q_ip2
      real(wp) :: q_L, q_R, dq, q6

      ! Flux computation variables
      real(wp) :: f_layer, u_eff, h_face, sigma
      real(wp) :: flux_left, flux_right, flux_bot, flux_top
      real(wp) :: q_L_up, q_R_up, dq_up, q6_up, tr_avg_up
      real(wp) :: corr_xL, corr_xR, corr_yB, corr_yT
      real(wp) :: H_active_xL, H_active_xR, H_active_yB, H_active_yT
      real(wp) :: inv_H_xL, inv_H_xR, inv_H_yB, inv_H_yT
      real(wp) :: v_xL(NZ_STACK_MAX), v_xR(NZ_STACK_MAX)
      real(wp) :: v_yB(NZ_STACK_MAX), v_yT(NZ_STACK_MAX)
      real(wp) :: hf_xL(NZ_STACK_MAX), hf_xR(NZ_STACK_MAX)
      real(wp) :: hf_yB(NZ_STACK_MAX), hf_yT(NZ_STACK_MAX)
      real(wp) :: h_iL, h_iR, h_jB, h_jT

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      inv_dx = 1.0_wp/dx
      inv_dy = 1.0_wp/dy

      ! Short-name aliases onto the caller-owned scratch dummies.
      associate (tr_snap => ppm_tr_snap, &
                 ql_x => ppm_ql_x, qr_x => ppm_qr_x, &
                 ql_y => ppm_ql_y, qr_y => ppm_qr_y)

         ! ---- Phase 0: concentration snapshot + fallback initialisation ----
         ! Ghost cells keep ql = qr = cell average (piecewise constant),
         ! so boundary faces degrade gracefully to first-order upwind.
         !$omp parallel do collapse(3)
         do j = 1, ny
         do i = 1, nx
         do k = 1, nz
            if (h_layer(i, j, k) > DRY_TOLERANCE) then
               tr_snap(i, j, k) = hTr_layer(i, j, k)/h_layer(i, j, k)
            else
               tr_snap(i, j, k) = 0.0_wp
            end if
            ql_x(i, j, k) = tr_snap(i, j, k)
            qr_x(i, j, k) = tr_snap(i, j, k)
            ql_y(i, j, k) = tr_snap(i, j, k)
            qr_y(i, j, k) = tr_snap(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do

         ! ---- Phase 1: PPM reconstruction for interior cells ----
         ! Stencil: i-2..i+2, fits within nghost=2.
         !$omp parallel do collapse(2) private(k, kk, q_im2, q_im1, q_i, q_ip1, q_ip2, q_L, q_R, dq, q6, h_col)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            h_col = 0.0_wp
            do kk = 1, nz
               h_col = h_col + h_layer(i, j, kk)
            end do
            if (h_col > ML_MIN_DEPTH) then
               do k = 1, nz

                  ! --- x-direction ---
                  q_im2 = tr_snap(i - 2, j, k)
                  q_im1 = tr_snap(i - 1, j, k)
                  q_i = tr_snap(i, j, k)
                  q_ip1 = tr_snap(i + 1, j, k)
                  q_ip2 = tr_snap(i + 2, j, k)

                  ! 4th-order interpolation to cell edges
                  q_L = (7.0_wp/12.0_wp)*(q_im1 + q_i) &
                        - (1.0_wp/12.0_wp)*(q_im2 + q_ip1)
                  q_R = (7.0_wp/12.0_wp)*(q_i + q_ip1) &
                        - (1.0_wp/12.0_wp)*(q_im1 + q_ip2)

                  ! Colella-Woodward monotonicity limiting
                  if ((q_R - q_i)*(q_i - q_L) <= 0.0_wp) then
                     ! Local extremum: flatten to constant
                     q_L = q_i
                     q_R = q_i
                  else
                     dq = q_R - q_L
                     q6 = 6.0_wp*(q_i - 0.5_wp*(q_L + q_R))
                     if (dq*q6 > dq*dq) then
                        q_L = 3.0_wp*q_i - 2.0_wp*q_R
                     end if
                     if (dq*q6 < -(dq*dq)) then
                        q_R = 3.0_wp*q_i - 2.0_wp*q_L
                     end if
                  end if

                  ql_x(i, j, k) = q_L
                  qr_x(i, j, k) = q_R

                  ! --- y-direction ---
                  q_im2 = tr_snap(i, j - 2, k)
                  q_im1 = tr_snap(i, j - 1, k)
                  q_i = tr_snap(i, j, k)
                  q_ip1 = tr_snap(i, j + 1, k)
                  q_ip2 = tr_snap(i, j + 2, k)

                  q_L = (7.0_wp/12.0_wp)*(q_im1 + q_i) &
                        - (1.0_wp/12.0_wp)*(q_im2 + q_ip1)
                  q_R = (7.0_wp/12.0_wp)*(q_i + q_ip1) &
                        - (1.0_wp/12.0_wp)*(q_im1 + q_ip2)

                  if ((q_R - q_i)*(q_i - q_L) <= 0.0_wp) then
                     q_L = q_i
                     q_R = q_i
                  else
                     dq = q_R - q_L
                     q6 = 6.0_wp*(q_i - 0.5_wp*(q_L + q_R))
                     if (dq*q6 > dq*dq) then
                        q_L = 3.0_wp*q_i - 2.0_wp*q_R
                     end if
                     if (dq*q6 < -(dq*dq)) then
                        q_R = 3.0_wp*q_i - 2.0_wp*q_L
                     end if
                  end if

                  ql_y(i, j, k) = q_L
                  qr_y(i, j, k) = q_R
               end do
            end if
         end do
         end do
         !$omp end parallel do

         ! ---- Phase 2: flux computation and conservative update ----
         !$omp parallel do collapse(2) private(k, kk, h_val, f_layer, u_eff, h_face, sigma, q_L_up, q_R_up, dq_up, q6_up, tr_avg_up, flux_left, flux_right, flux_bot, flux_top, h_col, corr_xL, corr_xR, corr_yB, corr_yT, H_active_xL, H_active_xR, H_active_yB, H_active_yT, inv_H_xL, inv_H_xR, inv_H_yB, inv_H_yT, v_xL, v_xR, v_yB, v_yT, hf_xL, hf_xR, hf_yB, hf_yT, h_iL, h_iR, h_jB, h_jT)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            h_col = 0.0_wp
            do kk = 1, nz
               h_col = h_col + h_layer(i, j, kk)
            end do
            if (h_col > ML_MIN_DEPTH) then

               ! Vanishing-layer-gated centred CWC (matches
               ! ml_advect_tracer above).  Layers with either side
               ! h_layer(k) <= VANISHING_LAYER_TOL contribute zero per-
               ! layer flux; residual goes to active layers via per-face
               ! symmetric frac (h_face_active(k) / sum h_face_active).
               corr_xL = mass_flux_x(i - 1, j)
               corr_xR = mass_flux_x(i, j)
               corr_yB = mass_flux_y(i, j - 1)
               corr_yT = mass_flux_y(i, j)
               H_active_xL = 0.0_wp
               H_active_xR = 0.0_wp
               H_active_yB = 0.0_wp
               H_active_yT = 0.0_wp
               do k = 1, nz
                  h_val = h_layer(i, j, k)
                  h_iL = h_layer(i - 1, j, k)
                  h_iR = h_layer(i + 1, j, k)
                  h_jB = h_layer(i, j - 1, k)
                  h_jT = h_layer(i, j + 1, k)

                  if (h_val > VANISHING_LAYER_TOL .and. h_iL > VANISHING_LAYER_TOL) then
                     v_xL(k) = 0.5_wp*(hu_layer(i - 1, j, k) + hu_layer(i, j, k))
                     hf_xL(k) = 0.5_wp*(h_iL + h_val)
                  else
                     v_xL(k) = 0.0_wp
                     hf_xL(k) = 0.0_wp
                  end if
                  if (h_val > VANISHING_LAYER_TOL .and. h_iR > VANISHING_LAYER_TOL) then
                     v_xR(k) = 0.5_wp*(hu_layer(i, j, k) + hu_layer(i + 1, j, k))
                     hf_xR(k) = 0.5_wp*(h_val + h_iR)
                  else
                     v_xR(k) = 0.0_wp
                     hf_xR(k) = 0.0_wp
                  end if
                  if (h_val > VANISHING_LAYER_TOL .and. h_jB > VANISHING_LAYER_TOL) then
                     v_yB(k) = 0.5_wp*(hv_layer(i, j - 1, k) + hv_layer(i, j, k))
                     hf_yB(k) = 0.5_wp*(h_jB + h_val)
                  else
                     v_yB(k) = 0.0_wp
                     hf_yB(k) = 0.0_wp
                  end if
                  if (h_val > VANISHING_LAYER_TOL .and. h_jT > VANISHING_LAYER_TOL) then
                     v_yT(k) = 0.5_wp*(hv_layer(i, j, k) + hv_layer(i, j + 1, k))
                     hf_yT(k) = 0.5_wp*(h_val + h_jT)
                  else
                     v_yT(k) = 0.0_wp
                     hf_yT(k) = 0.0_wp
                  end if

                  corr_xL = corr_xL - v_xL(k)
                  corr_xR = corr_xR - v_xR(k)
                  corr_yB = corr_yB - v_yB(k)
                  corr_yT = corr_yT - v_yT(k)
                  H_active_xL = H_active_xL + hf_xL(k)
                  H_active_xR = H_active_xR + hf_xR(k)
                  H_active_yB = H_active_yB + hf_yB(k)
                  H_active_yT = H_active_yT + hf_yT(k)
               end do

               if (H_active_xL > 0.0_wp) then
                  inv_H_xL = 1.0_wp/H_active_xL
               else
                  inv_H_xL = 0.0_wp
               end if
               if (H_active_xR > 0.0_wp) then
                  inv_H_xR = 1.0_wp/H_active_xR
               else
                  inv_H_xR = 0.0_wp
               end if
               if (H_active_yB > 0.0_wp) then
                  inv_H_yB = 1.0_wp/H_active_yB
               else
                  inv_H_yB = 0.0_wp
               end if
               if (H_active_yT > 0.0_wp) then
                  inv_H_yT = 1.0_wp/H_active_yT
               else
                  inv_H_yT = 0.0_wp
               end if

               do k = 1, nz

                  ! ---- x-direction fluxes ----

                  ! Left face (i-1/2)
                  f_layer = v_xL(k) + hf_xL(k)*inv_H_xL*corr_xL
                  h_face = 0.5_wp*(h_layer(i - 1, j, k) + h_layer(i, j, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dx, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     ! Upwind is cell i-1: sweep from right edge of i-1
                     q_L_up = ql_x(i - 1, j, k)
                     q_R_up = qr_x(i - 1, j, k)
                     tr_avg_up = tr_snap(i - 1, j, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_left = f_layer*(q_R_up &
                                          - 0.5_wp*sigma*(dq_up - q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else if (f_layer < 0.0_wp) then
                     ! Upwind is cell i: sweep from left edge of i
                     q_L_up = ql_x(i, j, k)
                     q_R_up = qr_x(i, j, k)
                     tr_avg_up = tr_snap(i, j, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_left = f_layer*(q_L_up &
                                          + 0.5_wp*sigma*(dq_up + q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else
                     flux_left = 0.0_wp
                  end if

                  ! Right face (i+1/2)
                  f_layer = v_xR(k) + hf_xR(k)*inv_H_xR*corr_xR
                  h_face = 0.5_wp*(h_layer(i, j, k) + h_layer(i + 1, j, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dx, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     q_L_up = ql_x(i, j, k)
                     q_R_up = qr_x(i, j, k)
                     tr_avg_up = tr_snap(i, j, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_right = f_layer*(q_R_up &
                                           - 0.5_wp*sigma*(dq_up - q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else if (f_layer < 0.0_wp) then
                     q_L_up = ql_x(i + 1, j, k)
                     q_R_up = qr_x(i + 1, j, k)
                     tr_avg_up = tr_snap(i + 1, j, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_right = f_layer*(q_L_up &
                                           + 0.5_wp*sigma*(dq_up + q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else
                     flux_right = 0.0_wp
                  end if

                  ! ---- y-direction fluxes ----

                  ! Bottom face (j-1/2)
                  f_layer = v_yB(k) + hf_yB(k)*inv_H_yB*corr_yB
                  h_face = 0.5_wp*(h_layer(i, j - 1, k) + h_layer(i, j, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dy, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     q_L_up = ql_y(i, j - 1, k)
                     q_R_up = qr_y(i, j - 1, k)
                     tr_avg_up = tr_snap(i, j - 1, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_bot = f_layer*(q_R_up &
                                         - 0.5_wp*sigma*(dq_up - q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else if (f_layer < 0.0_wp) then
                     q_L_up = ql_y(i, j, k)
                     q_R_up = qr_y(i, j, k)
                     tr_avg_up = tr_snap(i, j, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_bot = f_layer*(q_L_up &
                                         + 0.5_wp*sigma*(dq_up + q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else
                     flux_bot = 0.0_wp
                  end if

                  ! Top face (j+1/2)
                  f_layer = v_yT(k) + hf_yT(k)*inv_H_yT*corr_yT
                  h_face = 0.5_wp*(h_layer(i, j, k) + h_layer(i, j + 1, k))
                  if (h_face > DRY_TOLERANCE) then
                     u_eff = f_layer/h_face
                     sigma = min(abs(u_eff)*dt*inv_dy, 1.0_wp)
                  else
                     u_eff = 0.0_wp
                     sigma = 0.0_wp
                  end if

                  if (f_layer > 0.0_wp) then
                     q_L_up = ql_y(i, j, k)
                     q_R_up = qr_y(i, j, k)
                     tr_avg_up = tr_snap(i, j, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_top = f_layer*(q_R_up &
                                         - 0.5_wp*sigma*(dq_up - q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else if (f_layer < 0.0_wp) then
                     q_L_up = ql_y(i, j + 1, k)
                     q_R_up = qr_y(i, j + 1, k)
                     tr_avg_up = tr_snap(i, j + 1, k)
                     dq_up = q_R_up - q_L_up
                     q6_up = 6.0_wp*(tr_avg_up - 0.5_wp*(q_L_up + q_R_up))
                     flux_top = f_layer*(q_L_up &
                                         + 0.5_wp*sigma*(dq_up + q6_up*(1.0_wp - (2.0_wp/3.0_wp)*sigma)))
                  else
                     flux_top = 0.0_wp
                  end if

                  ! Conservative update
                  hTr_layer(i, j, k) = hTr_layer(i, j, k) &
                                       - dt*((flux_right - flux_left)*inv_dx &
                                             + (flux_top - flux_bot)*inv_dy)
               end do
            end if
         end do
         end do
         !$omp end parallel do

      end associate

   end subroutine ml_advect_tracer_ppm

   subroutine tracer_scratch_ensure(scratch, nx, ny, nz, &
                                    want_ppm, want_upwind, want_hdiff, want_weno)
      !! Allocate (or resize) the caller-owned tracer-kernel scratch and
      !! attach it to the GPU.  Idempotent: a no-op once each requested
      !! group exists at the requested size; resizes if (nx,ny,nz) changes
      !! (e.g. a unit test re-running at a different grid size).  Each
      !! group is guarded so only the kernels actually in use allocate
      !! (PPM-only run skips the upwind snapshot and vice-versa).
      !!
      !! Call once up front (solver_enter_data), never per step — the
      !! kernels themselves are `pure` and must not allocate.
      type(tracer_scratch_t), intent(inout) :: scratch
      integer, intent(in) :: nx, ny, nz
      logical, intent(in) :: want_ppm, want_upwind, want_hdiff
      logical, intent(in), optional :: want_weno
         !! WENO5/7/9 reconstruction snapshot (`tracer_recon` /= ppm).
         !! Optional so pre-WENO callers stay source-compatible.

      if (want_upwind) then
         call scratch_array_ensure(scratch%adv_hTr_snap, nx, ny, nz)
      end if

      if (want_ppm) then
         call scratch_array_ensure(scratch%ppm_tr_snap, nx, ny, nz)
         call scratch_array_ensure(scratch%ppm_ql_x, nx, ny, nz)
         call scratch_array_ensure(scratch%ppm_qr_x, nx, ny, nz)
         call scratch_array_ensure(scratch%ppm_ql_y, nx, ny, nz)
         call scratch_array_ensure(scratch%ppm_qr_y, nx, ny, nz)
      end if

      if (present(want_weno)) then
         if (want_weno) then
            call scratch_array_ensure(scratch%weno_tr_snap, nx, ny, nz)
         end if
      end if

      if (want_hdiff) then
         call scratch_array_ensure(scratch%hdiff_hTr_pre, nx, ny, nz)
      end if
   end subroutine tracer_scratch_ensure

   subroutine tracer_scratch_destroy(scratch)
      !! Detach and deallocate every allocated scratch group.  Called by
      !! the solver at simulation end (mirror of `tracer_scratch_ensure`).
      type(tracer_scratch_t), intent(inout) :: scratch

      call scratch_array_destroy(scratch%adv_hTr_snap)
      call scratch_array_destroy(scratch%ppm_tr_snap)
      call scratch_array_destroy(scratch%ppm_ql_x)
      call scratch_array_destroy(scratch%ppm_qr_x)
      call scratch_array_destroy(scratch%ppm_ql_y)
      call scratch_array_destroy(scratch%ppm_qr_y)
      call scratch_array_destroy(scratch%weno_tr_snap)
      call scratch_array_destroy(scratch%hdiff_hTr_pre)
   end subroutine tracer_scratch_destroy

   subroutine scratch_array_ensure(arr, nx, ny, nz)
      !! Allocate `arr` to (nx,ny,nz) + attach to the device, resizing if
      !! it already exists at a different shape.  The array is passed as a
      !! plain dummy so the `!$acc` directives carry no derived-type
      !! dereference (the leaf form the OpenMP-portability rule requires).
      real(wp), allocatable, intent(inout) :: arr(:, :, :)
      integer, intent(in) :: nx, ny, nz

      if (allocated(arr)) then
         if (size(arr, 1) == nx .and. size(arr, 2) == ny .and. &
             size(arr, 3) == nz) return
         ! [acc->omp-cpu] dropped device directive (host no-op): exit data
         deallocate (arr)
      end if
      allocate (arr(nx, ny, nz))
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
   end subroutine scratch_array_ensure

   subroutine scratch_array_destroy(arr)
      !! Detach + deallocate a single scratch array (no-op if unallocated).
      real(wp), allocatable, intent(inout) :: arr(:, :, :)

      if (.not. allocated(arr)) return
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      deallocate (arr)
   end subroutine scratch_array_destroy

   subroutine ml_distribute_tracer(h_layer, hTr_layer, dsig, h_new, &
                                        grid, nz)
      !! Update hTr through the barotropic step: preserve tracer Tr,
      !! set h_k = dsig_k * H_new, hTr_k = h_k_new * Tr_k_old.
      !!
      !! Call BEFORE ml_distribute_barotropic (which overwrites h_layer).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
         !! Old layer depths (before barotropic update)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
         !! Tracer h*Tr (updated using old Tr and new h)
      real(wp), intent(in) :: dsig(nz)
      real(wp), intent(in) :: h_new(grid%nx_total, grid%ny_total)
         !! Total depth AFTER barotropic update

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: tr_k, h_k_old, h_k_new

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(2) private(k, tr_k, h_k_old, h_k_new)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         do k = 1, nz
            h_k_old = h_layer(i, j, k)
            h_k_new = dsig(k)*h_new(i, j)
            if (h_k_old > DRY_TOLERANCE) then
               tr_k = hTr_layer(i, j, k)/h_k_old
            else
               tr_k = 0.0_wp
            end if
            hTr_layer(i, j, k) = h_k_new*tr_k
         end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_distribute_tracer

   subroutine ml_vertical_exchange_tracer(h_layer, hTr_layer, w, &
                                               grid, nz, dt)
      !! Vertical tracer exchange using upwind advection at interfaces.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: w(grid%nx_total, grid%ny_total, nz + 1)
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: w_face, tr_up, flux_tr, h_val, h_col

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (nz < 2) return

      !$omp parallel do collapse(2) private(k, kk, w_face, tr_up, flux_tr, h_val, h_col)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then
            do k = 1, nz - 1
               w_face = w(i, j, k + 1)

               if (w_face > 0.0_wp) then
                  h_val = h_layer(i, j, k)
                  if (h_val > DRY_TOLERANCE) then
                     tr_up = hTr_layer(i, j, k)/h_val
                  else
                     tr_up = 0.0_wp
                  end if
               else
                  h_val = h_layer(i, j, k + 1)
                  if (h_val > DRY_TOLERANCE) then
                     tr_up = hTr_layer(i, j, k + 1)/h_val
                  else
                     tr_up = 0.0_wp
                  end if
               end if

               flux_tr = w_face*tr_up*dt
               ! Limit: don't extract more than 50% of donor layer tracer
               if (w_face > 0.0_wp) then
                  if (abs(flux_tr) > 0.5_wp*abs(hTr_layer(i, j, k))) then
                     flux_tr = sign(0.5_wp*abs(hTr_layer(i, j, k)), flux_tr)
                  end if
               else
                  if (abs(flux_tr) > 0.5_wp*abs(hTr_layer(i, j, k + 1))) then
                     flux_tr = sign(0.5_wp*abs(hTr_layer(i, j, k + 1)), flux_tr)
                  end if
               end if

               hTr_layer(i, j, k) = hTr_layer(i, j, k) - flux_tr
               hTr_layer(i, j, k + 1) = hTr_layer(i, j, k + 1) + flux_tr
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine ml_vertical_exchange_tracer

   subroutine ml_settling_tracer(h_layer, hTr_layer, bed, w_s, grid, nz, dt)
      !! Constant-velocity downward settling of one tracer (e.g. suspended
      !! sediment), conservative + positivity-preserving.  Bottom-up
      !! convention: k=1 is the bed layer, k=nz the surface.  The whole
      !! suspended profile sinks at `w_s`; whatever crosses the bed deposits
      !! into the 2D `bed` reservoir, so total `sum_k hTr + bed` is conserved.
      !!
      !! Implemented as a conservative PPM remap (Colella & Woodward 1984) on
      !! an *augmented* column.  Settling by `s = w_s*dt` is a downward
      !! advection of the concentration profile by `s`; expressed as a remap,
      !! the new content of fixed layer k is the old reconstruction integrated
      !! over [z_{k-1}+s, z_k+s].  To keep that a same-total (hence exactly
      !! conservative) remap we add two cells: a clean cap (conc 0, thickness
      !! `s`) above the surface to supply the zero-inflow surface BC, and a
      !! deposition slab (thickness `s`) below the bed to catch settled-out
      !! mass.  `remap_column_ppm` then does the parabolic reconstruction +
      !! Colella-Woodward limiting; target cell 1 is the deposition, target
      !! cell k+1 is the new layer k.
      !!
      !! Unlike the old donor-cell flux this is 2nd-order (less numerical
      !! diffusion) and has no one-layer-per-step limit: `s` may span several
      !! layers (or the whole column) in a single call and still conserve and
      !! stay non-negative.  nz=1 reduces to the same linear sink as before;
      !! dry / vanishing layers (zero thickness) carry zero mass.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: bed(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: w_s, dt

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: settle, h_col, mass_old, mass_new
      real(wp) :: dz_src(NZ_STACK_MAX), dz_tgt(NZ_STACK_MAX)
      real(wp) :: c_src(NZ_STACK_MAX), c_tgt(NZ_STACK_MAX)

      if (w_s <= 0.0_wp) return

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      settle = w_s*dt

      !$omp parallel do collapse(2) private(k, h_col, mass_old, mass_new, dz_src, dz_tgt, c_src, c_tgt)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         ! Skip dry columns — nothing to settle.
         h_col = 0.0_wp
         do k = 1, nz
            h_col = h_col + h_layer(i, j, k)
         end do
         if (h_col > DRY_TOLERANCE) then
            ! Source: real layers (bed..surface) then a clean cap of
            ! thickness `settle` (concentration 0).
            mass_old = 0.0_wp
            do k = 1, nz
               dz_src(k) = h_layer(i, j, k)
               if (h_layer(i, j, k) > 0.0_wp) then
                  c_src(k) = hTr_layer(i, j, k)/h_layer(i, j, k)
               else
                  c_src(k) = 0.0_wp
               end if
               mass_old = mass_old + hTr_layer(i, j, k)
            end do
            dz_src(nz + 1) = settle
            c_src(nz + 1) = 0.0_wp

            ! Target: deposition slab (thickness `settle`) at the bed, then the
            ! real layer thicknesses shifted up one index.  Same total height
            ! as the source, so remapping onto it advects the profile down by
            ! `settle` exactly conservatively.
            dz_tgt(1) = settle
            do k = 1, nz
               dz_tgt(k + 1) = h_layer(i, j, k)
            end do

            call remap_column_ppm(nz + 1, dz_src, dz_tgt, c_src, c_tgt)

            ! Target cell k+1 is the new content of real layer k; target
            ! cell 1 is the deposition.  Bank the column-mass deficit into the
            ! bed so total (suspended + bed) is conserved to round-off.
            mass_new = 0.0_wp
            do k = 1, nz
               hTr_layer(i, j, k) = c_tgt(k + 1)*h_layer(i, j, k)
               mass_new = mass_new + hTr_layer(i, j, k)
            end do
            bed(i, j) = bed(i, j) + (mass_old - mass_new)
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine ml_settling_tracer

   subroutine ml_erode_bed_sediment(h_layer, hu_layer, hv_layer, hTr_sed, bed, &
                                         manning_n, rho_0, tau_ce, e_rate, porosity, &
                                         grid, nz, dt)
      !! Ariathurai-Arulanandan (1978) bed erosion of suspended sediment.
      !! Where the bed shear stress exceeds the critical erosion stress
      !! `tau_ce`, lift bed sediment into the bottom (lowest active) layer:
      !!
      !!   E = e_rate * (1 - porosity) * (tau_b/tau_ce - 1)   [kg/m^2/s], +ve only
      !!
      !! integrated over `dt` and capped at the available bed mass for
      !! positivity.  Conservative against the bed: the suspended bottom
      !! layer gains exactly what `bed` loses, so total (suspended + bed) is
      !! preserved.  The bed stress is the *same* one the momentum drag
      !! applies (Manning-consistent, mirrors `ml_apply_bottom_friction`):
      !!
      !!   cf = g * n^2 / H^(1/3),   tau_b = rho_0 * cf * |u_b|^2
      !!
      !! so erosion responds to the resolved bed shear without a separate
      !! sediment-only drag law.  Bottom-up convention: k=1 is the bed; the
      !! lowest non-vanishing layer receives the eroded mass (handles
      !! ZSTAR_FULL parked bottom layers like the friction kernel).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hTr_sed(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: bed(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: manning_n, rho_0, tau_ce, e_rate, porosity, dt

      integer :: i, j, kk, k_bot, nx, ny, nghost
      real(wp) :: h_k, h_total, u_b, v_b, speed_b, cf, tau_b, dm

      if (e_rate <= 0.0_wp .or. tau_ce <= 0.0_wp .or. manning_n <= 0.0_wp) return

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(2) private(kk, k_bot, h_k, h_total, u_b, v_b, speed_b, cf, tau_b, dm)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         ! Lowest active layer (mirrors ml_apply_bottom_friction).
         k_bot = 1
         do kk = 1, nz
            if (h_layer(i, j, kk) > VANISHING_LAYER_TOL) then
               k_bot = kk
               exit
            end if
         end do

         h_k = h_layer(i, j, k_bot)
         if (h_k < THIN_LAYER_THRESHOLD) cycle

         h_total = 0.0_wp
         do kk = 1, nz
            h_total = h_total + h_layer(i, j, kk)
         end do
         if (h_total < THIN_LAYER_THRESHOLD) cycle

         u_b = hu_layer(i, j, k_bot)/h_k
         v_b = hv_layer(i, j, k_bot)/h_k
         speed_b = sqrt(u_b*u_b + v_b*v_b)

         cf = GRAVITY*manning_n*manning_n/(h_total**(1.0_wp/3.0_wp))
         tau_b = rho_0*cf*speed_b*speed_b
         if (tau_b <= tau_ce) cycle

         ! Eroded mass this step, capped at what the bed holds.
         dm = e_rate*(1.0_wp - porosity)*(tau_b/tau_ce - 1.0_wp)*dt
         dm = min(dm, bed(i, j))
         bed(i, j) = bed(i, j) - dm
         hTr_sed(i, j, k_bot) = hTr_sed(i, j, k_bot) + dm
      end do
      end do
      !$omp end parallel do

   end subroutine ml_erode_bed_sediment

   subroutine ml_clamp_tracer(hTr_layer, h_layer, tr_min, tr_max, &
                                   grid, nz)
      !! Enforce physical tracer bounds tr_min <= Tr <= tr_max
      !! conservatively within each column.
      !!
      !! Two-pass column-conservative limiter (MOM6/NEMO style):
      !!
      !!   Pass 1 — clip violators.  For any layer with `S = hTr/h`
      !!            outside [tr_min, tr_max], reset hTr to the bound and
      !!            accumulate the rejected `excess_above` (above tr_max,
      !!            destroyed) and `deficit_below` (below tr_min, created)
      !!            in PSU·m.  Track the per-cell `cap_to_add` (room
      !!            still available to push concentration up toward
      !!            tr_max in non-saturated layers) and `cap_to_remove`
      !!            (room to pull concentration down toward tr_min).
      !!
      !!   Pass 2 — redistribute the net imbalance `net = excess_above
      !!            - deficit_below` into the column's spare capacity.
      !!            If net > 0 (more destroyed than created): push it
      !!            back into layers with cap_to_add, proportionally to
      !!            each layer's room-to-fill.  If net < 0: pull it
      !!            from layers with cap_to_remove, proportionally to
      !!            each layer's room-to-empty.  When the column-mean
      !!            S is itself outside [tr_min, tr_max] (rare for
      !!            normal seawater) the redistribution saturates and
      !!            the residual imbalance is unavoidably destroyed /
      !!            created — this is the irreducible cost of holding
      !!            a hard bound.
      !!
      !! Per-cell tracer total is preserved exactly when the column-
      !! mean concentration is in bounds, even though individual layers
      !! were violating.  This is the right behaviour for the rescaled-
      !! S inflation produced by `ml_distribute_with_remap_full` under
      !! VCOORD_ZSTAR_FULL: the per-layer S after the conservative
      !! remap can sit above tr_max as a numerical artefact of the
      !! `dz_old_col = remap_dz_old * scale_old` rescaling, while the
      !! physical column-mean is still well inside bounds; the previous
      !! hard clip would have destroyed the very tracer the conservative
      !! remap had just preserved.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: tr_min, tr_max

      integer :: i, j, k, nx, ny, nghost
      real(wp) :: tr_k, h_k, room_k
      real(wp) :: excess_above, deficit_below, cap_to_add, cap_to_remove
      real(wp) :: net, frac
      real(wp) :: col_hTr, col_h
      real(wp), parameter :: TINY_CAP = 1.0e-30_wp

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(2) private(k, tr_k, h_k, room_k, excess_above, deficit_below, cap_to_add, cap_to_remove, net, frac, col_hTr, col_h)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost

         ! Pass 0: column mean.  The two-pass clip+redistribute below is
         ! EXACTLY conservative if and only if the column-mean concentration
         ! lies within [tr_min, tr_max].  Proof for the upper bound:
         !   mean <= tr_max
         !     <=> sum_k (tr_max - tr_k) h_k >= 0
         !     <=> -excess_above + cap_unsat + (cap_below + deficit_below) >= 0
         !     <=> excess_above - deficit_below <= cap_unsat + cap_below
         !     <=> net <= cap_to_add
         ! so `frac = min(1, net/cap_to_add)` never saturates and every unit
         ! removed in pass 1 is returned in pass 2.  Symmetrically for tr_min.
         !
         ! When the mean is OUTSIDE the bounds no clamp can both hold the
         ! bound and conserve, and the old code silently destroyed the
         ! un-absorbable residual `net - cap_to_add`.  With the default
         ! S_max sitting at 40 PSU that fired continuously at sharp fronts
         ! and leaked ~9e-6 of the domain salt per 400 steps.  Conservation
         ! wins: skip the column instead.  Leaving the value visibly out of
         ! range also surfaces the upstream problem rather than hiding it.
         col_hTr = 0.0_wp
         col_h = 0.0_wp
         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k > DRY_TOLERANCE) then
               col_hTr = col_hTr + hTr_layer(i, j, k)
               col_h = col_h + h_k
            else
               hTr_layer(i, j, k) = 0.0_wp
            end if
         end do
         if (col_h > DRY_TOLERANCE) then
            tr_k = col_hTr/col_h
            if (tr_k > tr_max .or. tr_k < tr_min) cycle
         else
            cycle
         end if

         ! Pass 1: hard clip violators, accumulate per-cell excess/deficit
         ! and the spare capacity remaining in non-saturated layers.
         excess_above = 0.0_wp
         deficit_below = 0.0_wp
         cap_to_add = 0.0_wp
         cap_to_remove = 0.0_wp
         do k = 1, nz
            h_k = h_layer(i, j, k)
            if (h_k > DRY_TOLERANCE) then
               tr_k = hTr_layer(i, j, k)/h_k
               if (tr_k > tr_max) then
                  excess_above = excess_above + (tr_k - tr_max)*h_k
                  hTr_layer(i, j, k) = tr_max*h_k
                  ! Saturated above: no room to add, all room to remove.
                  cap_to_remove = cap_to_remove + (tr_max - tr_min)*h_k
               else if (tr_k < tr_min) then
                  deficit_below = deficit_below + (tr_min - tr_k)*h_k
                  hTr_layer(i, j, k) = tr_min*h_k
                  ! Saturated below: all room to add, none to remove.
                  cap_to_add = cap_to_add + (tr_max - tr_min)*h_k
               else
                  cap_to_add = cap_to_add + (tr_max - tr_k)*h_k
                  cap_to_remove = cap_to_remove + (tr_k - tr_min)*h_k
               end if
            else
               hTr_layer(i, j, k) = 0.0_wp
            end if
         end do

         ! Pass 2: redistribute net imbalance into available capacity.
         ! Conservative within the cell when capacity is sufficient;
         ! falls back to leaving the residual outside bounds when the
         ! column-mean itself is out of range.
         net = excess_above - deficit_below
         if (net > 0.0_wp .and. cap_to_add > TINY_CAP) then
            frac = min(1.0_wp, net/cap_to_add)
            do k = 1, nz
               h_k = h_layer(i, j, k)
               if (h_k > DRY_TOLERANCE) then
                  tr_k = hTr_layer(i, j, k)/h_k
                  room_k = (tr_max - tr_k)*h_k
                  if (room_k > 0.0_wp) then
                     hTr_layer(i, j, k) = hTr_layer(i, j, k) + frac*room_k
                  end if
               end if
            end do
         else if (net < 0.0_wp .and. cap_to_remove > TINY_CAP) then
            frac = min(1.0_wp, (-net)/cap_to_remove)
            do k = 1, nz
               h_k = h_layer(i, j, k)
               if (h_k > DRY_TOLERANCE) then
                  tr_k = hTr_layer(i, j, k)/h_k
                  room_k = (tr_k - tr_min)*h_k
                  if (room_k > 0.0_wp) then
                     hTr_layer(i, j, k) = hTr_layer(i, j, k) - frac*room_k
                  end if
               end if
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine ml_clamp_tracer

   subroutine ml_horizontal_diffusion_tracer(h_layer, hTr_layer, kappa, &
                                                  hdiff_hTr_pre, &
                                                  grid, nz, dt)
      !! Horizontal Laplacian diffusion of a scalar tracer in flux form,
      !! conservative per-face under non-uniform layer thicknesses.
      !!
      !! Per face the diffusive flux is
      !!   flux_face_x = h_face * kappa_eff * (Tr_nbr - Tr_self) / dx
      !! where `h_face = 0.5 * (h_self + h_nbr)` is symmetric across the
      !! face.  Applied with opposite signs to the two adjacent cells
      !! through the divergence form
      !!   d(hTr)/dt = (flux_x_right - flux_x_left)/dx + (flux_y_top - flux_y_bot)/dy
      !! so what cell A gains from face A-B equals what cell B loses
      !! through the same face; mass conservation holds even when the
      !! per-cell layer thicknesses differ across the face.
      !!
      !! History: the previous formulation computed `dtr = kappa·dt·∇²Tr`
      !! and updated `hTr_new = h_self · (Tr_self + dtr)`.  Each cell
      !! used its OWN h on the RHS, so adjacent cells with different
      !! `h_layer(k)` (e.g., zsigma / zstar at bathymetric transitions)
      !! saw different mass changes through the same shared face,
      !! leaking total mass.  Currently dormant in production runs
      !! (`hdiff_kappa = 0` in the LE / Sydney Harbour configs) but
      !! would bite the moment any user enables it.
      !!
      !! Skips dry columns (sum(h_layer) < ML_MIN_DEPTH) and treats dry
      !! neighbours as zero-flux walls (face flux = 0), preserving the
      !! symmetric per-face conservation property.  Includes a
      !! diffusive stability limiter: kappa_eff <= 0.5 / (dt·(1/dx² + 1/dy²)).
      !!
      !! Jacobi vs Gauss-Seidel: snapshot `hTr_layer` into a local
      !! scratch on entry and read neighbour concentrations from the
      !! frozen copy, so CPU-sequential and GPU-parallel builds give
      !! bit-identical results.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hdiff_hTr_pre(grid%nx_total, grid%ny_total, nz)
         !! Caller-owned Jacobi pre-step snapshot scratch.
      real(wp), intent(in) :: kappa, dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: dx, dy
      real(wp) :: tr_c, tr_ip, tr_im, tr_jp, tr_jm, h_col, h_self
      real(wp) :: kappa_eff, inv_dx, inv_dy, inv_dx2, inv_dy2
      real(wp) :: flux_x_right, flux_x_left, flux_y_top, flux_y_bot

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost
      dx = grid%dx
      dy = grid%dy

      if (kappa <= 0.0_wp) return

      ! Stability limiter: kappa * dt * (1/dx^2 + 1/dy^2) <= 0.5
      inv_dx = 1.0_wp/dx
      inv_dy = 1.0_wp/dy
      inv_dx2 = inv_dx*inv_dx
      inv_dy2 = inv_dy*inv_dy
      kappa_eff = min(kappa, 0.5_wp/(dt*(inv_dx2 + inv_dy2) + 1.0e-30_wp))

      associate (hTr_pre => hdiff_hTr_pre)

         !$omp parallel do collapse(3)
         do k = 1, nz
         do j = 1, ny
         do i = 1, nx
            hTr_pre(i, j, k) = hTr_layer(i, j, k)
         end do
         end do
         end do
         !$omp end parallel do

         !$omp parallel do collapse(2) private(k, kk, tr_c, tr_ip, tr_im, tr_jp, tr_jm, h_col, h_self, flux_x_right, flux_x_left, flux_y_top, flux_y_bot)
         do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost

            h_col = 0.0_wp
            do kk = 1, nz
               h_col = h_col + h_layer(i, j, kk)
            end do
            if (h_col < ML_MIN_DEPTH) cycle

            do k = 1, nz
               h_self = h_layer(i, j, k)
               if (h_self < THIN_LAYER_THRESHOLD) cycle

               tr_c = hTr_pre(i, j, k)/h_self

               ! Per-face flux: h_face = 0.5*(h_self + h_nbr) is symmetric,
               ! face flux = h_face*kappa*(tr_nbr - tr_self)/dx.  Dry
               ! neighbours close the face (zero-flux wall).
               if (h_layer(i + 1, j, k) > THIN_LAYER_THRESHOLD) then
                  tr_ip = hTr_pre(i + 1, j, k)/h_layer(i + 1, j, k)
                  flux_x_right = 0.5_wp*(h_self + h_layer(i + 1, j, k)) &
                                 *kappa_eff*(tr_ip - tr_c)*inv_dx
               else
                  flux_x_right = 0.0_wp
               end if
               if (h_layer(i - 1, j, k) > THIN_LAYER_THRESHOLD) then
                  tr_im = hTr_pre(i - 1, j, k)/h_layer(i - 1, j, k)
                  flux_x_left = 0.5_wp*(h_layer(i - 1, j, k) + h_self) &
                                *kappa_eff*(tr_c - tr_im)*inv_dx
               else
                  flux_x_left = 0.0_wp
               end if
               if (h_layer(i, j + 1, k) > THIN_LAYER_THRESHOLD) then
                  tr_jp = hTr_pre(i, j + 1, k)/h_layer(i, j + 1, k)
                  flux_y_top = 0.5_wp*(h_self + h_layer(i, j + 1, k)) &
                               *kappa_eff*(tr_jp - tr_c)*inv_dy
               else
                  flux_y_top = 0.0_wp
               end if
               if (h_layer(i, j - 1, k) > THIN_LAYER_THRESHOLD) then
                  tr_jm = hTr_pre(i, j - 1, k)/h_layer(i, j - 1, k)
                  flux_y_bot = 0.5_wp*(h_layer(i, j - 1, k) + h_self) &
                               *kappa_eff*(tr_c - tr_jm)*inv_dy
               else
                  flux_y_bot = 0.0_wp
               end if

               hTr_layer(i, j, k) = hTr_layer(i, j, k) &
                                    + dt*((flux_x_right - flux_x_left)*inv_dx &
                                          + (flux_y_top - flux_y_bot)*inv_dy)
            end do
         end do
         end do
         !$omp end parallel do

      end associate

   end subroutine ml_horizontal_diffusion_tracer

   subroutine ml_init_tracer(h_layer, hTr_layer, tr_init, &
                                  grid, nz)
      !! Initialise hTr = h_k * tr_init for all layers
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: tr_init

      integer :: i, j, k, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         hTr_layer(i, j, k) = h_layer(i, j, k)*tr_init
      end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_init_tracer

   subroutine ml_init_tracer_stratified(h_layer, hTr_layer, tr_bottom, &
                                             tr_surface, grid, nz)
      !! Seed a tracer with a linear vertical profile: `tr_bottom` at k=1
      !! (bed) to `tr_surface` at k=nz (surface), per the bottom-up layer
      !! convention.  The profile is computed inline from the two scalars
      !! (no host-side per-layer array) so the seed runs correctly on the
      !! device after `enter_data`.  Caller guarantees nz > 1.  Mirrors the
      !! unstructured `ml_init_tracer_stratified_unstr`.
      !!   hTr_k = h_layer_k * [tr_bottom + dtr*(k - 1)]
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: tr_bottom, tr_surface

      integer :: i, j, k, nx, ny
      real(wp) :: dtr

      nx = grid%nx_total
      ny = grid%ny_total
      dtr = (tr_surface - tr_bottom)/real(nz - 1, wp)

      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         hTr_layer(i, j, k) = h_layer(i, j, k) &
                              *(tr_bottom + dtr*real(k - 1, wp))
      end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_init_tracer_stratified

   subroutine ml_seed_tracer(h_layer, hTr_layer, tr_uniform, &
                                  tr_bottom, tr_surface, grid, nz)
      !! Seed one tracer's `hTr` field, dispatching between the uniform and
      !! the linear-stratified profile from the tracer's IC scalars.  When
      !! both `tr_bottom` and `tr_surface` are set (non-zero) and nz > 1 a
      !! linear profile is laid down from `tr_bottom` at the bed (k=1) to
      !! `tr_surface` at the surface (k=nz); otherwise every layer is filled
      !! with `tr_uniform`.  This is the single stratified-vs-uniform gate
      !! shared by all IC seed sites (rki_init, rki_api), so the polarity
      !! convention lives in exactly one place.  Mirrors `ml_seed_tracer_unstr`.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: tr_uniform, tr_bottom, tr_surface

      if (tr_bottom /= 0.0_wp .and. tr_surface /= 0.0_wp .and. nz > 1) then
         call ml_init_tracer_stratified(h_layer, hTr_layer, tr_bottom, &
                                        tr_surface, grid, nz)
      else
         call ml_init_tracer(h_layer, hTr_layer, tr_uniform, grid, nz)
      end if

   end subroutine ml_seed_tracer

   subroutine ml_save_tracer(hTr_layer, hTr_layer0, grid, nz)
      !! Save tracer field for RK2 averaging
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: hTr_layer0(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         hTr_layer0(i, j, k) = hTr_layer(i, j, k)
      end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_save_tracer

   subroutine ml_rk2_average_tracer(hTr_layer, hTr_layer0, &
                                         grid, nz)
      !! SSP-RK2 averaging for a tracer: hTr = 0.5*(hTr0 + hTr).
      !! Clamps the result to >= 0 inline; additional bounds (upper and/or
      !! tracer-specific lower) are enforced by ml_clamp_tracer if called
      !! after this kernel.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hTr_layer0(grid%nx_total, grid%ny_total, nz)

      integer :: i, j, k, nx, ny, nghost

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      !$omp parallel do collapse(3)
      do k = 1, nz
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         hTr_layer(i, j, k) = 0.5_wp*(hTr_layer0(i, j, k) + hTr_layer(i, j, k))
         if (hTr_layer(i, j, k) < 0.0_wp) hTr_layer(i, j, k) = 0.0_wp
      end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_rk2_average_tracer

   pure function edge_prescribes(edge_type) result(prescribes)
      !! An edge accepts a prescribed inflow tracer value when it is an
      !! inflow (BC_INFLOW) or a Dirichlet-eta OPEN edge (BC_TIDAL /
      !! BC_CLAMPED) — the edges through which water can enter the domain
      !! carrying a specified concentration.  The prescription still only
      !! takes effect when tr_inflow >= 0 (checked at the call site), and the
      !! prescribed ghost is only READ on inflow (upwind), so a reversing
      !! tidal edge is inflow-value on flood and zero-gradient on ebb.
      integer, intent(in) :: edge_type
      logical :: prescribes
      prescribes = (edge_type == BC_INFLOW .or. edge_type == BC_TIDAL &
                    .or. edge_type == BC_CLAMPED)
   end function edge_prescribes

   subroutine ml_apply_tracer_bcs(hTr_layer, h_layer, tr_inflow, &
                                       grid, nz, bc)
      !! Set a scalar tracer in ghost cells.
      !!
      !! Default: zero-gradient extrapolation (copy hTr from nearest interior
      !! cell, scaled by local layer thickness).  For OPEN boundaries with a
      !! prescribed tracer value (tr_inflow >= 0) — inflow (BC_INFLOW) and the
      !! Dirichlet-eta open edges (BC_TIDAL / BC_CLAMPED) — ghost cells are set
      !! to h_layer * tr_inflow so that INFLOWING water carries the prescribed
      !! value rather than the interior one.  The ghost value is only READ by
      !! the upwind advection when the face flux is INTO the domain (flood on
      !! a tidal edge); on OUTFLOW (ebb) the interior cell is upwind, so the
      !! prescribed ghost is ignored and the edge is zero-gradient — the
      !! reversing tidal edge gets the correct inflow-value / zero-gradient-
      !! outflow behaviour automatically.  With tr_inflow < 0 (the default,
      !! unset) every edge is zero-gradient — bit-identical to before.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: tr_inflow
      type(bc_params_t), intent(in) :: bc

      integer :: i, j, k, g, i_int, nx, ny, nghost
      logical :: prescribe_west, prescribe_east, prescribe_south, prescribe_north

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      prescribe_west = (edge_prescribes(bc%west) .and. tr_inflow >= 0.0_wp)
      prescribe_east = (edge_prescribes(bc%east) .and. tr_inflow >= 0.0_wp)
      prescribe_south = (edge_prescribes(bc%south) .and. tr_inflow >= 0.0_wp)
      prescribe_north = (edge_prescribes(bc%north) .and. tr_inflow >= 0.0_wp)

      ! West
      do g = 1, nghost
         i_int = nghost + 1
         if (prescribe_west) then
            !$omp parallel do collapse(2)
            do j = 1, ny
            do k = 1, nz
               hTr_layer(g, j, k) = h_layer(g, j, k)*tr_inflow
            end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do collapse(2)
            do j = 1, ny
            do k = 1, nz
               hTr_layer(g, j, k) = hTr_layer(i_int, j, k) &
                                    *h_layer(g, j, k) &
                                    /max(h_layer(i_int, j, k), DRY_TOLERANCE)
            end do
            end do
            !$omp end parallel do
         end if
      end do
      ! East
      do g = nx - nghost + 1, nx
         i_int = nx - nghost
         if (prescribe_east) then
            !$omp parallel do collapse(2)
            do j = 1, ny
            do k = 1, nz
               hTr_layer(g, j, k) = h_layer(g, j, k)*tr_inflow
            end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do collapse(2)
            do j = 1, ny
            do k = 1, nz
               hTr_layer(g, j, k) = hTr_layer(i_int, j, k) &
                                    *h_layer(g, j, k) &
                                    /max(h_layer(i_int, j, k), DRY_TOLERANCE)
            end do
            end do
            !$omp end parallel do
         end if
      end do
      ! South
      do g = 1, nghost
         i_int = nghost + 1
         if (prescribe_south) then
            !$omp parallel do collapse(2)
            do i = 1, nx
            do k = 1, nz
               hTr_layer(i, g, k) = h_layer(i, g, k)*tr_inflow
            end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do collapse(2)
            do i = 1, nx
            do k = 1, nz
               hTr_layer(i, g, k) = hTr_layer(i, i_int, k) &
                                    *h_layer(i, g, k) &
                                    /max(h_layer(i, i_int, k), DRY_TOLERANCE)
            end do
            end do
            !$omp end parallel do
         end if
      end do
      ! North
      do g = ny - nghost + 1, ny
         i_int = ny - nghost
         if (prescribe_north) then
            !$omp parallel do collapse(2)
            do i = 1, nx
            do k = 1, nz
               hTr_layer(i, g, k) = h_layer(i, g, k)*tr_inflow
            end do
            end do
            !$omp end parallel do
         else
            !$omp parallel do collapse(2)
            do i = 1, nx
            do k = 1, nz
               hTr_layer(i, g, k) = hTr_layer(i, i_int, k) &
                                    *h_layer(i, g, k) &
                                    /max(h_layer(i, i_int, k), DRY_TOLERANCE)
            end do
            end do
            !$omp end parallel do
         end if
      end do

   end subroutine ml_apply_tracer_bcs

   subroutine ml_vertical_diffusion_tracer(h_layer, hu_layer, hv_layer, &
                                                hTr_layer, rho_layer, &
                                                rho_0, nu_max, kappa_bg, ri_alpha, &
                                                grid, nz, dt)
      !! Ri-dependent vertical diffusion of a scalar tracer (Pacanowski-Philander 1981).
      !! Uses stronger Ri suppression for tracers: kappa_t = nu_max / (1 + alpha*Ri)^3 + kappa_bg.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0
         !! Reference density (kg/m^3)
      real(wp), intent(in) :: nu_max
         !! Maximum turbulent viscosity (m^2/s)
      real(wp), intent(in) :: kappa_bg
         !! Background tracer diffusivity (m^2/s)
      real(wp), intent(in) :: ri_alpha
         !! Richardson number scaling coefficient
      real(wp), intent(in) :: dt

      integer :: i, j, k, kk, nx, ny, nghost
      real(wp) :: h_avg, h_k, h_k1, u_k, v_k, u_k1, v_k1, tr_k, tr_k1, h_col
      real(wp) :: du, dv, dtr, drho, Ri, denom, kappa_t, flux_tr

      real(wp), parameter :: EPS = 1.0e-20_wp

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      if (nz < 2) return

      !$omp parallel do collapse(2) private(k, kk, h_k, h_k1, h_avg, u_k, v_k, u_k1, v_k1, tr_k, tr_k1, du, dv, dtr, drho, Ri, denom, kappa_t, flux_tr, h_col)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then
            do k = 1, nz - 1
               h_k = h_layer(i, j, k)
               h_k1 = h_layer(i, j, k + 1)
               h_avg = 0.5_wp*(h_k + h_k1)

               ! VCOORD_ZSTAR_FULL: vanishing layers on either side make
               ! kappa_t*dtr/h_avg pathological.  No-flux through such
               ! interfaces.
               if (h_avg > VANISHING_LAYER_TOL .and. &
                   h_k > VANISHING_LAYER_TOL .and. &
                   h_k1 > VANISHING_LAYER_TOL) then
                  if (h_k > DRY_TOLERANCE) then
                     u_k = hu_layer(i, j, k)/h_k
                     v_k = hv_layer(i, j, k)/h_k
                     tr_k = hTr_layer(i, j, k)/h_k
                  else
                     u_k = 0.0_wp
                     v_k = 0.0_wp
                     tr_k = 0.0_wp
                  end if
                  if (h_k1 > DRY_TOLERANCE) then
                     u_k1 = hu_layer(i, j, k + 1)/h_k1
                     v_k1 = hv_layer(i, j, k + 1)/h_k1
                     tr_k1 = hTr_layer(i, j, k + 1)/h_k1
                  else
                     u_k1 = 0.0_wp
                     v_k1 = 0.0_wp
                     tr_k1 = 0.0_wp
                  end if

                  du = u_k1 - u_k
                  dv = v_k1 - v_k
                  drho = rho_layer(i, j, k) - rho_layer(i, j, k + 1)
                  Ri = GRAVITY*(drho/rho_0)*h_avg/(du*du + dv*dv + EPS)
                  if (Ri < 0.0_wp) Ri = 0.0_wp

                  denom = 1.0_wp + ri_alpha*Ri
                  kappa_t = nu_max/(denom*denom*denom) + kappa_bg

                  ! Limit diffusive flux: stability requires kappa_t*dt/h_avg^2 < 0.5
                  if (kappa_t*dt > 0.5_wp*h_avg*h_avg) then
                     kappa_t = 0.5_wp*h_avg*h_avg/dt
                  end if

                  dtr = tr_k1 - tr_k
                  flux_tr = kappa_t*dtr/h_avg*dt

                  ! Conservative + dimensionally correct: hTr has units [tr]*m,
                  ! flux_tr = kappa*dtr/h_avg*dt has units [tr]*m. Same quantity
                  ! moves across the interface so both layers see the same mass
                  ! change (previously multiplied each side by its own h, which
                  ! leaks mass under non-uniform layers, e.g. VCOORD_ZSTAR_FULL).
                  hTr_layer(i, j, k) = hTr_layer(i, j, k) + flux_tr
                  hTr_layer(i, j, k + 1) = hTr_layer(i, j, k + 1) - flux_tr
               end if
            end do
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine ml_vertical_diffusion_tracer

   subroutine ml_apply_surface_tracer_fluxes(hS_layer, hT_layer, h_layer, &
                                                  hu_layer, hv_layer, &
                                                  heat_flux, evap_minus_precip, &
                                                  evap_mass_sink, evap_latent_cooling, &
                                                  rho_0, cp, grid, nz, &
                                                  dt)
      !! Apply prescribed surface fluxes to the top layer (k=nz).
      !!
      !! Two modes, selected by `evap_mass_sink`:
      !!
      !! **Default (evap_mass_sink = .false.) — virtual salt flux:**
      !!   hT(nz) += dt * heat_flux / (rho_0 * cp)    ! thermal
      !!   hS(nz) += dt * (E - P) * S_top             ! virtual salt flux
      !! Neither term changes h_layer — real volume forcing is the job of
      !! the existing precipitation field (applied to the barotropic h in
      !! apply_friction).
      !!
      !! **Mass-sink mode (evap_mass_sink = .true.) — real water removal:**
      !! net evaporation dh = E*dt is removed COLUMN-PROPORTIONALLY — every
      !! layer's h, hu, hv, hT is scaled by fac = 1 - dh/h_col while hS is
      !! held.  Consequences (all exact, any tracer profile):
      !!   - water mass:  column loses exactly dh;
      !!   - salt:        evaporation removes WATER not salt — hS unchanged,
      !!                  so S = hS/h rises by 1/fac;
      !!   - heat:        the evaporated water carries its own temperature
      !!                  away (hT scales with h), so T = hT/h is invariant —
      !!                  no spurious warming.  This is the SENSIBLE heat of
      !!                  the removed liquid mass;
      !!   - latent heat: with evap_latent_cooling, evaporation ALSO withdraws
      !!                  the latent heat of vaporization of the removed mass,
      !!                  Q_lat = rho_water*E_eff*L_v [W/m²], applied as an
      !!                  additional surface cooling to hT(nz).  hT(nz) is
      !!                  fac-scaled by the sensible sink FIRST; the latent
      !!                  term is the additive increment applied AFTER that
      !!                  scaling (it is NOT the full hT(nz) update):
      !!                    hT(nz) += -(rho_water*L_v/(rho_0*cp))*dh_eff  [°C·m].
      !!                  E_eff = dh_eff/dt is the effective (clamped) rate, so
      !!                  a dried-out column withdraws no latent heat.  This is
      !!                  DISTINCT from the sensible term above (phase-change
      !!                  energy vs. enthalpy of the departing liquid) — no
      !!                  double-count;
      !!   - momentum:    the water leaves at its own velocity (hu, hv scale
      !!                  with h), so u, v are invariant — no spurious
      !!                  acceleration.
      !! The proportional form deliberately preserves the CURRENT layer
      !! layout (sigma fractions are untouched), so the vcoord distribute /
      !! remap machinery sees a consistent column and no spurious
      !! layer-reshape transport arises.  This mirrors how the barotropic
      !! precipitation source dilutes the column as a whole; a
      !! surface-INTENSIFIED evaporative salinity structure (sink localised
      !! to k=nz with w-coupled tracer transport) is a documented follow-up.
      !!
      !! Precipitation (evap_minus_precip < 0) is NOT added in mass-sink
      !! mode — real freshwater input remains the job of the barotropic
      !! precipitation field.
      !!
      !! WET/DRY: the removal is clamped to the water available above
      !! ML_MIN_DEPTH (dh_eff = min(dh, h_col - ML_MIN_DEPTH)); a column at
      !! or below ML_MIN_DEPTH stops evaporating.  The clamped (unsatisfied)
      !! fraction of the demand is dropped honestly — the water simply was
      !! not there to evaporate; no compensating source is invented.
      !!
      !! Layer-index convention: the structured solver indexes layers
      !! bottom-up (k=1 is the bottom layer, k=nz is the surface layer
      !! just below the free surface).  This matches the BPG, rho_ref
      !! builder, and NH init.  Surface forcings therefore land at k=nz.
      !!
      !! Heat flux is in W/m^2, positive into the ocean.  Evaporation
      !! minus precipitation is in m/s, positive for net evaporation.
      !!
      !! Skips columns shallower than ML_MIN_DEPTH and cells with a
      !! thin top layer (h_top < THIN_LAYER_THRESHOLD).  The clamp
      !! kernel that follows enforces physical bounds.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(inout) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(inout) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: heat_flux(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: evap_minus_precip(grid%nx_total, grid%ny_total)
      logical, intent(in) :: evap_mass_sink
      logical, intent(in) :: evap_latent_cooling
      real(wp), intent(in) :: rho_0, cp, dt

      integer :: i, j, kk, nx, ny, nghost
      real(wp) :: h_top, kin_heat_flux, S_top, inv_rho0cp, h_col
      real(wp) :: dh_req, dh_eff, fac, latent_gain

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      inv_rho0cp = 1.0_wp/(rho_0*cp)
      ! Kinematic latent-heat withdrawal per metre of evaporated water:
      ! d(hT) = -latent_gain*dh_eff  [°C·m].  Zero unless latent cooling is on.
      latent_gain = 0.0_wp
      if (evap_latent_cooling) latent_gain = RHO_WATER*LATENT_HEAT_VAPORIZATION*inv_rho0cp

      !$omp parallel do collapse(2) private(kk, h_top, kin_heat_flux, S_top, h_col, dh_req, dh_eff, fac)
      do j = nghost + 1, ny - nghost
      do i = nghost + 1, nx - nghost
         h_col = 0.0_wp
         do kk = 1, nz
            h_col = h_col + h_layer(i, j, kk)
         end do
         if (h_col > ML_MIN_DEPTH) then
            h_top = h_layer(i, j, nz)
            if (h_top > THIN_LAYER_THRESHOLD) then
               ! Thermal (unchanged in both modes): dt * Q_net / (rho_0 * cp)
               kin_heat_flux = heat_flux(i, j)*inv_rho0cp
               hT_layer(i, j, nz) = hT_layer(i, j, nz) + dt*kin_heat_flux

               if (evap_mass_sink) then
                  ! Column-proportional water-mass sink (see docstring).
                  ! Salt mass hS is held; h, hu, hv, hT scale by fac so S
                  ! concentrates while T, u, v stay invariant.
                  dh_req = dt*max(evap_minus_precip(i, j), 0.0_wp)
                  dh_eff = min(dh_req, max(h_col - ML_MIN_DEPTH, 0.0_wp))
                  if (dh_eff > 0.0_wp) then
                     fac = 1.0_wp - dh_eff/h_col
                     do kk = 1, nz
                        h_layer(i, j, kk) = h_layer(i, j, kk)*fac
                        hu_layer(i, j, kk) = hu_layer(i, j, kk)*fac
                        hv_layer(i, j, kk) = hv_layer(i, j, kk)*fac
                        hT_layer(i, j, kk) = hT_layer(i, j, kk)*fac
                        ! hS_layer held -> concentration rises naturally.
                     end do
                     ! Latent-heat cooling (additive, top layer only; not
                     ! scaled by fac so the withdrawal equals Q_lat exactly).
                     ! latent_gain == 0 when the knob is off -> bit-identical.
                     hT_layer(i, j, nz) = hT_layer(i, j, nz) - latent_gain*dh_eff
                  end if
               else
                  ! Legacy virtual salt flux (bit-identical default).
                  S_top = hS_layer(i, j, nz)/h_top
                  hS_layer(i, j, nz) = hS_layer(i, j, nz) &
                                       + dt*evap_minus_precip(i, j)*S_top
               end if
            end if
         end if
      end do
      end do
      !$omp end parallel do

   end subroutine ml_apply_surface_tracer_fluxes

end module rki_ml_tracers
