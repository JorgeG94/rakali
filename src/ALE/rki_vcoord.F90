!! Vertical coordinate type for hybrid z-sigma and ALE support
module rki_vcoord
   !! Defines `vcoord_t`, a lightweight config type for the vertical coordinate.
   !! Type-bound procedures dispatch via `select case` on an integer enum (no
   !! runtime polymorphism) so the type is GPU-safe.
   !!
   !! Implemented (`VCOORD_*`):
   !!   - SIGMA       — pure terrain-following (no remap).
   !!   - ZSIGMA      — smoothstep blend sigma (shallow) → fixed z-levels (deep);
   !!                   reduces sigma PGE on steep bathymetry. Conservative remap.
   !!   - ZSTAR       — z*-lite SSH-tracking: a global `z_ref` stretched per column
   !!                   by H/z_ref(nz) so sum(dz)=H, relative spacing preserved as
   !!                   η changes (no sigma distortion at large SSH).
   !!   - ZSTAR_FULL  — per-column z_ref from local bathymetry (vanishing layers).
   !!   - ZSTAR_SIGMA — smoothstep blend of sigma (shallow) and z*-lite (deep).
#ifdef LFORTRAN_PASSING
   use rki_constants, only: wp, &
                            nz_stack_required, nz_stack_is_sufficient, &
                            VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, &
                            VCOORD_ZSTAR, VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED, VCOORD_RHO, VCOORD_HYCOM, &
                            REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM
#else
   use rki_constants, only: NZ_STACK_MAX, wp, &
                            nz_stack_required, nz_stack_is_sufficient, &
                            VCOORD_LAGRANGIAN, VCOORD_EULERIAN_Z, &
                            VCOORD_SIGMA, VCOORD_ZSIGMA, &
                            VCOORD_ZSTAR, VCOORD_ZSTAR_FULL, VCOORD_ZSTAR_SIGMA, &
                            VCOORD_Z_FIXED, VCOORD_RHO, VCOORD_HYCOM, &
                            REMAP_PCM, REMAP_PLM, REMAP_PPM, REMAP_PPM_H4, REMAP_PQM
#endif
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

#ifdef LFORTRAN_PASSING
   integer, parameter :: NZ_STACK_MAX = 64
      !! LFortran 0.64 workaround: module-local copy of the rki_constants value
      !! (an imported parameter used as an explicit-shape dummy bound inside a
      !! PURE call becomes an impure getter under LFortran). Keep in sync (=64).
#endif

   public :: vcoord_t
   public :: vcoord_target_dz_column
   public :: vcoord_target_dz_column_zstar_full
   public :: zstar_full_build_column
   public :: parse_vcoord_type
   public :: parse_remap_method
   public :: parse_stretching_mode

   integer, parameter, public :: STRETCH_UNIFORM = 0
   integer, parameter, public :: STRETCH_LOG = 1

   type :: vcoord_t
      integer :: coord_type = VCOORD_SIGMA
         !! Vertical coordinate type (VCOORD_SIGMA, VCOORD_ZSIGMA, etc.)
      integer :: remap_method = REMAP_PPM
         !! Remap reconstruction order (REMAP_PCM/PLM/PPM/PPM_H4/PQM). Default PPM.
         !! PQM falls back to PPM for nz < 5 (see `remap_column_pqm`).
      integer :: nz = 0
         !! Number of vertical layers

      real(wp), allocatable :: dsig_target(:)
         !! (nz) Target sigma-like layer fractions, sum = 1.0

      real(wp), allocatable :: z_ref(:)
         !! (nz+1) Reference interface depths for z-levels (m, positive down).
         !! Only allocated for coord types that use z-levels.

      ! z-sigma hybrid parameters
      real(wp) :: depth_transition = 100.0_wp
         !! Depth (m) where sigma-to-z blending begins
      real(wp) :: blend_width = 50.0_wp
         !! Width of the blending zone (m)

      ! Full MOM6 z* parameters (VCOORD_ZSTAR_FULL)
      real(wp) :: zstar_h_surf_target = 0.0_wp
         !! Target physical thickness of the surface layer (m).
         !! 0 => auto: uniform per-column (falls back to lite behaviour
         !! per column).  Set > 0 to anchor the surface layer at a
         !! fixed thickness regardless of total depth H.
      real(wp) :: zstar_h_min = 1.0e-4_wp
         !! Vanishing-layer floor (m).  Layers that would land below the
         !! local bed get clipped to this thickness rather than going to zero.
      integer :: zstar_stretching = 1
         !! Surface-concentration stretching: 1=log, 0=uniform
      integer :: zstar_n_surf = 0
         !! Number of "fine" near-surface layers using stretching.
         !! 0 => auto (use max(1, nz/3))
   contains
      procedure, non_overridable :: init => vcoord_init
      procedure, non_overridable :: needs_remap => vcoord_needs_remap
      procedure, non_overridable :: enter_data => vcoord_enter_data
      procedure, non_overridable :: exit_data => vcoord_exit_data
      procedure, non_overridable :: cleanup => vcoord_cleanup
   end type vcoord_t

contains

   subroutine vcoord_init(self, nz, coord_type, remap_method)
      !! Initialise a vertical coordinate definition. Allocates + populates
      !! `dsig_target` for the given type. Unknown types: error stop.
      class(vcoord_t), intent(inout) :: self
      integer, intent(in) :: nz
         !! Number of vertical layers
      integer, intent(in) :: coord_type
         !! Coordinate type constant (VCOORD_SIGMA, etc.)
      integer, intent(in) :: remap_method
         !! Remapping method constant (REMAP_PLM, etc.)

      self%nz = nz
      self%coord_type = coord_type
      self%remap_method = remap_method

      ! (Re)allocate target distribution
      if (allocated(self%dsig_target)) deallocate (self%dsig_target)
      allocate (self%dsig_target(nz))

      ! (Re)allocate z_ref if needed
      if (allocated(self%z_ref)) deallocate (self%z_ref)

      select case (coord_type)
      case (VCOORD_SIGMA)
         self%dsig_target = 1.0_wp/real(nz, wp)

      case (VCOORD_ZSIGMA)
         ! z-sigma hybrid: dsig_target is used for the sigma component,
         ! z_ref defines the reference z-level interfaces.
         self%dsig_target = 1.0_wp/real(nz, wp)
         allocate (self%z_ref(0:nz))
         ! Reference z-levels: uniform spacing over the blending depth.
         ! Columns shallower than depth_transition use pure sigma.
         ! Columns deeper blend toward these fixed horizontal levels.
         block
            real(wp) :: h_ref
            integer :: kk
            h_ref = self%depth_transition + self%blend_width
            self%z_ref(0) = 0.0_wp
            do kk = 1, nz
               self%z_ref(kk) = real(kk, wp)*h_ref/real(nz, wp)
            end do
         end block

      case (VCOORD_ZSTAR_FULL)
         ! Per-column z_ref built later by zstar_full_build_column from local
         ! bathymetry; this uniform template keeps allocated(z_ref) checks valid.
         self%dsig_target = 1.0_wp/real(nz, wp)
         allocate (self%z_ref(0:nz))
         block
            real(wp) :: h_ref
            integer :: kk
            h_ref = self%depth_transition + self%blend_width
            self%z_ref(0) = 0.0_wp
            do kk = 1, nz
               self%z_ref(kk) = real(kk, wp)*h_ref/real(nz, wp)
            end do
         end block

      case (VCOORD_ZSTAR_SIGMA)
         ! z*/sigma hybrid: same template as VCOORD_ZSTAR / VCOORD_ZSIGMA.
         ! The blend is applied per-column in `vcoord_target_dz_column`
         ! based on the local H vs depth_transition / blend_width.
         self%dsig_target = 1.0_wp/real(nz, wp)
         allocate (self%z_ref(0:nz))
         block
            real(wp) :: h_ref
            integer :: kk
            h_ref = self%depth_transition + self%blend_width
            self%z_ref(0) = 0.0_wp
            do kk = 1, nz
               self%z_ref(kk) = real(kk, wp)*h_ref/real(nz, wp)
            end do
         end block

      case (VCOORD_ZSTAR)
         ! z*-lite: uniform reference spacing over depth_transition+blend_width;
         ! stretched per column by H/z_ref(nz) in vcoord_target_dz_column.
         self%dsig_target = 1.0_wp/real(nz, wp)
         allocate (self%z_ref(0:nz))
         block
            real(wp) :: h_ref
            integer :: kk
            h_ref = self%depth_transition + self%blend_width
            self%z_ref(0) = 0.0_wp
            do kk = 1, nz
               self%z_ref(kk) = real(kk, wp)*h_ref/real(nz, wp)
            end do
         end block

      case default
         call logger%error("vcoord_init: unknown coord_type "//to_string(coord_type)// &
                           ". Valid options are 'sigma', 'zsigma', and 'zstar'.")
         error stop "vcoord_init: unknown coord_type"
      end select

      ! Defence in depth.  `validate_config` refuses this outright, so a
      ! namelist-driven run never reaches here; this catches the FFI /
      ! two-phase `rki_create` paths that build a vcoord without going
      ! through config validation.  Threshold is `nz + 1` — see
      ! `nz_stack_required`.
      if (.not. nz_stack_is_sufficient(nz)) then
         call logger%warning("nz = "//to_string(nz)//" needs NZ_STACK_MAX >= "// &
                             to_string(nz_stack_required(nz))//" but this binary was "// &
                             "compiled with NZ_STACK_MAX = "//to_string(NZ_STACK_MAX)// &
                             ". Per-column stack kernels (BPG, remap) will overrun "// &
                             "thread-local storage and silently produce wrong answers. "// &
                             "Rebuild with -DRAKALI_NZ_STACK_MAX="// &
                             to_string(nz_stack_required(nz))//" (or larger).")
      end if
   end subroutine vcoord_init

   pure logical function vcoord_needs_remap(self) result(needs_remap)
      !! Returns `.true.` if this coordinate type requires conservative
      !! vertical remapping after the barotropic step.  Pure sigma does
      !! not — layers are simply rescaled by dsig(k) * H.
      class(vcoord_t), intent(in) :: self
      needs_remap = (self%coord_type /= VCOORD_SIGMA)
   end function vcoord_needs_remap

   subroutine vcoord_enter_data(self)
      !! Map read-only coordinate arrays to GPU.
      class(vcoord_t), intent(inout) :: self
      if (allocated(self%dsig_target)) then
         ! [acc->omp-cpu] dropped device directive (host no-op): enter data
      end if
      if (allocated(self%z_ref)) then
         ! [acc->omp-cpu] dropped device directive (host no-op): enter data
      end if
   end subroutine vcoord_enter_data

   subroutine vcoord_exit_data(self)
      !! Unmap coordinate arrays from GPU.
      class(vcoord_t), intent(inout) :: self
      if (allocated(self%dsig_target)) then
         ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      end if
      if (allocated(self%z_ref)) then
         ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      end if
   end subroutine vcoord_exit_data

   subroutine vcoord_cleanup(self)
      !! Deallocate all arrays.  Safe to call on uninitialised instances.
      class(vcoord_t), intent(inout) :: self
      if (allocated(self%dsig_target)) deallocate (self%dsig_target)
      if (allocated(self%z_ref)) deallocate (self%z_ref)
      self%nz = 0
   end subroutine vcoord_cleanup

   ! ---- Per-column target thickness computation ----

   pure subroutine vcoord_target_dz_column(coord_type, nz, H, dsig, z_ref, &
                                           depth_transition, blend_width, dz)
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      !! Compute target layer thicknesses for a single water column.
      !! Pure, called from `do concurrent` (one thread per column).
      !!
      !! Output `dz` is bottom-up (ROMS): dz(1) = BOTTOM, dz(nz) = SURFACE;
      !! matches the solvers' `h_layer` indexing (consumers write
      !! h_layer(k,...) = dz(k) with no reversal). sum(dz) = H in all cases.
      !!   SIGMA:       dz(k) = dsig(k)*H (terrain-following).
      !!   ZSIGMA:      smooth blend sigma (shallow) → fixed z-levels (deep).
      !!   ZSTAR:       z*-lite, dz(k) = (z_ref(nz-k+1)-z_ref(nz-k))*H/z_ref(nz);
      !!                interfaces stay at fixed relative position as η changes.
      !!   ZSTAR_SIGMA: smoothstep blend of sigma (shallow) and z*-lite (deep);
      !!                each branch sums to H so no surface trim needed.
      integer, intent(in) :: coord_type
      integer, intent(in) :: nz
      real(wp), intent(in) :: H
         !! Total water depth at this column (m)
      real(wp), intent(in) :: dsig(nz)
         !! Reference sigma fractions (sum = 1), ROMS-ordered: dsig(1) bottom,
         !! dsig(nz) surface. Currently uniform 1/nz.
      real(wp), intent(in) :: z_ref(0:nz)
         !! Reference z-level interface depths (m, positive down).
         !! `z_ref(0) = 0` is the surface; `z_ref(nz)` is the deepest
         !! reference interface.  Used for VCOORD_ZSIGMA and VCOORD_ZSTAR.
      real(wp), intent(in) :: depth_transition
         !! Depth (m) below which blending begins
      real(wp), intent(in) :: blend_width
         !! Width of the blending zone (m)
      real(wp), intent(out) :: dz(nz)
         !! Output target layer thicknesses, ROMS-ordered (sum = H,
         !! `dz(1)` = bottom, `dz(nz)` = surface)

      real(wp) :: alpha, x, z_top_k, z_bot_k, dz_z, dz_sum, deficit
      integer :: k

      select case (coord_type)

      case (VCOORD_SIGMA)
         ! Pure terrain-following.  Uniform dsig means orientation is
         ! immaterial — the output is the same in either direction.
         do k = 1, nz
            dz(k) = dsig(k)*H
         end do

      case (VCOORD_ZSTAR)
         ! z*-lite: stretch the global reference pattern uniformly so
         ! sum(dz) = H (H = h_bed + η). ROMS-ordered: dz(1) deepest, dz(nz)
         ! surface. Reduces to uniform sigma for uniform z_ref.
         if (z_ref(nz) > 0.0_wp) then
            do k = 1, nz
               dz_z = z_ref(nz - k + 1) - z_ref(nz - k)
               dz(k) = dz_z*H/z_ref(nz)
            end do
         else
            ! Degenerate z_ref — fall back to uniform sigma
            do k = 1, nz
               dz(k) = dsig(k)*H
            end do
         end if

      case (VCOORD_ZSIGMA)
         ! Smooth blend sigma↔z-levels: alpha=0 shallow (H<=depth_transition,
         ! pure sigma) → alpha=1 deep (H>=depth_transition+blend_width, z-levels).
         ! ROMS order: layer k spans z_ref(nz-k)..z_ref(nz-k+1).
         if (H <= depth_transition) then
            ! Pure sigma — uniform dsig, orientation immaterial
            do k = 1, nz
               dz(k) = dsig(k)*H
            end do
         else
            ! Compute blending factor
            if (blend_width > 0.0_wp) then
               x = (H - depth_transition)/blend_width
               x = max(0.0_wp, min(1.0_wp, x))
               alpha = x*x*(3.0_wp - 2.0_wp*x)  ! smoothstep
            else
               alpha = 1.0_wp
            end if

            ! Compute z-level thicknesses (clip to column depth)
            dz_sum = 0.0_wp
            do k = 1, nz
               z_top_k = min(z_ref(nz - k), H)        ! shallower interface
               z_bot_k = min(z_ref(nz - k + 1), H)    ! deeper interface
               dz_z = max(z_bot_k - z_top_k, 0.0_wp)
               ! Blend: (1-alpha)*sigma + alpha*z-level
               dz(k) = (1.0_wp - alpha)*dsig(k)*H + alpha*dz_z
               dz_sum = dz_sum + dz(k)
            end do

            ! When H is shallower than the deepest z_ref, layers near the
            ! bottom (k=1 under ROMS) get clipped to zero.  Put the
            ! deficit back into the bottom layer so sum(dz) = H.
            deficit = H - dz_sum
            if (abs(deficit) > 0.0_wp) then
               dz(1) = dz(1) + deficit
            end if
         end if

      case (VCOORD_ZSTAR_SIGMA)
         ! z*/sigma smoothstep blend: alpha=0 shallow (pure sigma) → alpha=1 deep
         ! (pure z*-lite). Both branches sum to H so the blend does, no trim.
         if (H <= depth_transition .or. z_ref(nz) <= 0.0_wp) then
            ! Pure sigma — including the degenerate-z_ref fallback so
            ! mass conservation holds without an extra dz_sum fix-up.
            do k = 1, nz
               dz(k) = dsig(k)*H
            end do
         else
            if (blend_width > 0.0_wp) then
               x = (H - depth_transition)/blend_width
               x = max(0.0_wp, min(1.0_wp, x))
               alpha = x*x*(3.0_wp - 2.0_wp*x)  ! smoothstep
            else
               alpha = 1.0_wp
            end if
            ! z*-lite contribution: dz_z = (z_ref(nz-k+1)-z_ref(nz-k)) * H/z_ref(nz)
            do k = 1, nz
               dz_z = (z_ref(nz - k + 1) - z_ref(nz - k))*H/z_ref(nz)
               dz(k) = (1.0_wp - alpha)*dsig(k)*H + alpha*dz_z
            end do
         end if

      case default
         ! Fallback to sigma (uniform, orientation immaterial)
         do k = 1, nz
            dz(k) = dsig(k)*H
         end do

      end select
   end subroutine vcoord_target_dz_column

   ! ---- Full MOM6 z* per-column reference builder ----

   pure subroutine zstar_full_build_column(nz, h_bed, h_surf_target, &
                                           n_surf, stretching, &
                                           z_ref_col)
      !! Build a per-column reference z-level pattern for VCOORD_ZSTAR_FULL.
      !! Output `z_ref_col(0:nz)` monotonically increasing (positive-down),
      !! z_ref_col(0)=0 surface, z_ref_col(nz)=h_bed. Built top-down; the caller
      !! handles ROMS ordering. Top `n_surf_use` layers use the stretching mode
      !! (log/uniform) toward h_surf_target; below that uniform to the bed;
      !! degenerate (h_bed≤0, nz≤0) ⇒ uniform.
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_bed
         !! Local bed depth (m, positive down).  Must be > 0 in normal use.
      real(wp), intent(in) :: h_surf_target
         !! Target thickness of the surface layer (m).  ≤ 0 means "auto"
         !! and falls back to uniform spacing.
      integer, intent(in) :: n_surf
         !! Number of fine near-surface layers.  ≤ 0 means auto.
      integer, intent(in) :: stretching
         !! STRETCH_LOG or STRETCH_UNIFORM
      real(wp), intent(out) :: z_ref_col(0:nz)

      integer :: k, n_surf_use, n_coarse
      real(wp) :: h_fine, h_coarse, dz_uniform, r, base, w

      if (nz <= 0) return

      ! Degenerate / dry column: output zeros (never negative interfaces, which
      ! would propagate as negative dz_new and crash the solver via NaN CFL).
      if (h_bed <= 0.0_wp) then
         do k = 0, nz
            z_ref_col(k) = 0.0_wp
         end do
         return
      end if

      ! Auto n_surf: 1/3 of layers, at least 1, at most nz-1.
      if (n_surf <= 0) then
         n_surf_use = max(1, nz/3)
      else
         n_surf_use = max(1, min(nz - 1, n_surf))
      end if

      ! If no surface target requested, just uniform.
      if (h_surf_target <= 0.0_wp .or. nz == 1) then
         z_ref_col(0) = 0.0_wp
         do k = 1, nz
            z_ref_col(k) = h_bed*real(k, wp)/real(nz, wp)
         end do
         return
      end if

      n_coarse = nz - n_surf_use
      z_ref_col(0) = 0.0_wp

      select case (stretching)
      case (STRETCH_LOG)
         ! Geometric fine zone: layer k thickness = h_surf_target*r^(k-1). Pick r
         ! (bisection, monotonic) so the bottom fine layer matches the coarse
         ! thickness and fine+coarse sums = h_bed.
         if (n_surf_use == 1 .or. n_coarse == 0) then
            ! Single fine layer or no coarse — fall through to UNIFORM behaviour
            r = 1.0_wp
         else
            block
               real(wp) :: r_lo, r_hi, r_mid, f_mid, target_ratio
               integer :: it
               target_ratio = h_bed/h_surf_target
               r_lo = 1.000001_wp
               r_hi = 10.0_wp
               do it = 1, 60
                  r_mid = 0.5_wp*(r_lo + r_hi)
                  ! f(r) = (r^n - 1)/(r - 1) + r^(n-1) * n_coarse
                  f_mid = (r_mid**n_surf_use - 1.0_wp)/(r_mid - 1.0_wp) &
                          + r_mid**(n_surf_use - 1)*real(n_coarse, wp)
                  if (f_mid > target_ratio) then
                     r_hi = r_mid
                  else
                     r_lo = r_mid
                  end if
                  if (r_hi - r_lo < 1.0e-9_wp) exit
               end do
               r = 0.5_wp*(r_lo + r_hi)
            end block
         end if
         w = 1.0_wp
         base = 0.0_wp
         do k = 1, n_surf_use
            base = base + h_surf_target*w
            z_ref_col(k) = base
            w = w*r
         end do
         ! Coarse zone: continue with dz = h_surf_target * r^(n_surf-1)
         dz_uniform = h_surf_target*r**(n_surf_use - 1)
         do k = n_surf_use + 1, nz
            z_ref_col(k) = z_ref_col(k - 1) + dz_uniform
         end do
      case default  ! STRETCH_UNIFORM
         ! Fine zone: n_surf_use layers each of thickness h_surf_target.
         ! Coarse zone: uniform fill of the remainder.
         h_fine = h_surf_target*real(n_surf_use, wp)
         if (h_fine >= h_bed) then
            ! Column too shallow to honour h_surf_target for all fine layers
            ! — fall back to uniform spacing across the whole column.
            do k = 1, nz
               z_ref_col(k) = h_bed*real(k, wp)/real(nz, wp)
            end do
            return
         end if
         h_coarse = h_bed - h_fine
         do k = 1, n_surf_use
            z_ref_col(k) = h_surf_target*real(k, wp)
         end do
         if (n_coarse > 0) then
            dz_uniform = h_coarse/real(n_coarse, wp)
            do k = n_surf_use + 1, nz
               z_ref_col(k) = h_fine + dz_uniform*real(k - n_surf_use, wp)
            end do
         end if
      end select

      ! Lock the bed interface exactly to h_bed (in case of round-off).
      z_ref_col(nz) = h_bed
   end subroutine zstar_full_build_column

   ! ---- Full MOM6 z* per-step dz builder (called per column on GPU) ----

   pure subroutine vcoord_target_dz_column_zstar_full(nz, H, z_ref_col, &
                                                      h_min, dz)
      ! [acc->omp-cpu] dropped device directive (host no-op): routine seq
      !! Compute target layer thicknesses for VCOORD_ZSTAR_FULL.
      !! Inputs: z_ref_col(0:nz) local reference (top-down, 0=surface, nz=h_bed);
      !! H current total depth (m) = h_bed + η; h_min vanishing-layer floor (m).
      !! Output: dz(1:nz) ROMS-ordered (dz(1) bottom, dz(nz) surface),
      !! sum(dz) = H exactly, vanishing rows set to h_min.
      !! Surface layer absorbs η; if H < h_bed the deepest layers clip to h_min
      !! and the surface is trimmed to keep sum = H.
      integer, intent(in) :: nz
      real(wp), intent(in) :: H
      real(wp), intent(in) :: z_ref_col(0:nz)
      real(wp), intent(in) :: h_min
      real(wp), intent(out) :: dz(nz)

      real(wp) :: dz_top_kt, z_upper, z_lower, H_eff, h_bed_ref, eta
      real(wp) :: sum_dz, deficit
      integer :: kt, k

      h_bed_ref = z_ref_col(nz)
      eta = H - h_bed_ref
      H_eff = max(H, 0.0_wp)

      if (eta >= 0.0_wp) then
         ! Column at or above reference depth.  Subsurface layers keep
         ! their reference thicknesses; the surface layer absorbs the
         ! SSH offset.  sum(dz) = h_bed_ref + eta = H exactly.
         do kt = 1, nz
            k = nz - kt + 1
            dz_top_kt = max(z_ref_col(kt) - z_ref_col(kt - 1), 0.0_wp)
            if (kt == 1) then
               dz(k) = dz_top_kt + eta     ! ROMS surface gets +η
            else
               dz(k) = dz_top_kt
            end if
         end do
      else
         ! Column shallower than reference (H < h_bed_ref). Walk top-down:
         ! layer above bed keeps full thickness; straddling layer gets
         ! H_eff - z_upper; below-bed layer is vanishing (h_min). Then trim the
         ! surface so sum = H exactly (remap mass conservation needs this).
         do kt = 1, nz
            k = nz - kt + 1
            z_upper = z_ref_col(kt - 1)
            z_lower = z_ref_col(kt)
            if (z_lower <= H_eff) then
               dz(k) = z_lower - z_upper
            else if (z_upper < H_eff) then
               dz(k) = H_eff - z_upper
            else
               dz(k) = h_min
            end if
         end do
         ! Re-balance via the surface layer.
         sum_dz = 0.0_wp
         do k = 1, nz
            sum_dz = sum_dz + dz(k)
         end do
         deficit = sum_dz - H_eff
         if (deficit > 0.0_wp) then
            if (dz(nz) - deficit >= h_min) then
               dz(nz) = dz(nz) - deficit
            else
               ! Too thin for the floor — set surface to h_min, let downstream
               ! physics guards (DRY_TOLERANCE) handle it; sum = H not enforced.
               dz(nz) = h_min
            end if
         end if
      end if
   end subroutine vcoord_target_dz_column_zstar_full

   ! ---- String-to-enum parsers for namelist config ----

   pure integer function parse_vcoord_type(str, default_code) result(coord_type)
      !! Convert a namelist string to a `VCOORD_*` constant (shared by coastal
      !! and ocean backends). Unrecognised ⇒ `default_code` if given, else
      !! `VCOORD_SIGMA`. Ocean callers pass `default_code = VCOORD_EULERIAN_Z`.
      character(len=*), intent(in) :: str
      integer, intent(in), optional :: default_code
      integer :: fallback
      fallback = VCOORD_SIGMA
      if (present(default_code)) fallback = default_code
      select case (trim(adjustl(str)))
      case ("lagrangian", "LAGRANGIAN", "isopycnal", "ISOPYCNAL")
         coord_type = VCOORD_LAGRANGIAN
      case ("eulerian_z", "EULERIAN_Z", "z", "Z")
         coord_type = VCOORD_EULERIAN_Z
      case ("sigma", "SIGMA")
         coord_type = VCOORD_SIGMA
      case ("zsigma", "ZSIGMA", "z-sigma", "z_sigma")
         coord_type = VCOORD_ZSIGMA
      case ("zstar", "ZSTAR", "z-star", "z_star", "zstar_lite", "ZSTAR_LITE")
         coord_type = VCOORD_ZSTAR
      case ("zstar_full", "ZSTAR_FULL", "z-star-full", "z_star_full", "zstarfull")
         coord_type = VCOORD_ZSTAR_FULL
      case ("zstar_sigma", "ZSTAR_SIGMA", "z-star-sigma", "z_star_sigma", "zstarsigma")
         coord_type = VCOORD_ZSTAR_SIGMA
      case ("z_fixed", "Z_FIXED", "z_levels", "Z_LEVELS", "gprime", "GPRIME")
         coord_type = VCOORD_Z_FIXED
      case ("rho", "RHO", "isopycnic", "ISOPYCNIC", "rho_target", "RHO_TARGET")
         coord_type = VCOORD_RHO
      case ("hycom", "HYCOM", "hybrid", "HYBRID")
         coord_type = VCOORD_HYCOM
      case default
         coord_type = fallback
      end select
   end function parse_vcoord_type

   pure integer function parse_stretching_mode(str) result(mode)
      !! Convert a namelist string to a STRETCH_* constant.
      character(len=*), intent(in) :: str
      select case (trim(adjustl(str)))
      case ("log")
         mode = STRETCH_LOG
      case ("uniform")
         mode = STRETCH_UNIFORM
      case default
         mode = STRETCH_LOG
      end select
   end function parse_stretching_mode

   pure integer function parse_remap_method(str) result(method)
      !! Convert a namelist string to a REMAP_* constant.
      character(len=*), intent(in) :: str
      select case (trim(adjustl(str)))
      case ("pcm")
         method = REMAP_PCM
      case ("plm")
         method = REMAP_PLM
      case ("ppm")
         method = REMAP_PPM
      case ("ppm_h4")
         method = REMAP_PPM_H4
      case ("pqm")
         method = REMAP_PQM
      case default
         method = REMAP_PLM
      end select
   end function parse_remap_method

end module rki_vcoord
