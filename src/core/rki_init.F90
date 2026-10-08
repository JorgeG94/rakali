!! Initial condition routines and diagnostics
module rki_init
   !! Provides initial condition setup and diagnostic functions.
   use rki_constants, only: wp, GRAVITY, DRY_TOLERANCE
   use rki_config, only: config_t
   use rki_state, only: state_t
   use rki_ml_common, only: ml_set_min_depth
   use rki_ml_eos, only: ml_equation_of_state, ml_equation_of_state_shared, &
                         ml_compute_reference_density
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   implicit none
   private

   public :: init_dam_break
   public :: init_gaussian_hump
   public :: init_castro_dam_break
   public :: init_still_water
   public :: init_from_config
   public :: init_multilayer_from_barotropic
   public :: refresh_multilayer_derived
   public :: restore_multilayer_after_restart
   public :: compute_total_mass
   public :: compute_layer_sums
   public :: diagnose_bathymetry

contains

   subroutine init_dam_break(s)
      !! Initialise a 1D dam break (left half h=2m, right half h=1m, no
      !! velocity). Bathymetry preserved if already loaded, else zeroed.
      !! Public only for the unit-test suite.
      type(state_t), intent(inout) :: s

      integer :: i, j, mid_x
      real(wp), parameter :: H_LEFT = 2.0_wp
      real(wp), parameter :: H_RIGHT = 1.0_wp
      logical :: has_bathymetry

      mid_x = s%grid%nghost + s%grid%nx_phys/2

      ! Check if bathymetry was already loaded (non-zero values present)
      has_bathymetry = any(s%barotropic%b /= 0.0_wp)

      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            if (.not. has_bathymetry) s%barotropic%b(i, j) = 0.0_wp
            s%barotropic%hu(i, j) = 0.0_wp
            s%barotropic%hv(i, j) = 0.0_wp
            if (i <= mid_x) then
               s%barotropic%h(i, j) = H_LEFT
            else
               s%barotropic%h(i, j) = H_RIGHT
            end if
         end do
      end do

   end subroutine init_dam_break

   subroutine init_gaussian_hump(s, h0, amplitude, sigma, x_off, y_off, u_shear)
      !! Initialise a Gaussian hump perturbation. When u_shear /= 0 and NH
      !! is enabled, the top layer gets +u_shear and the bottom -u_shear
      !! (vertical shear → non-zero vorticity). Bathymetry preserved if
      !! already loaded. Public only for the unit-test suite.
      type(state_t), intent(inout) :: s
      real(wp), intent(in) :: h0
         !! Background water depth (m)
      real(wp), intent(in) :: amplitude
         !! Hump amplitude (m)
      real(wp), intent(in) :: sigma
         !! Gaussian width (m)
      real(wp), intent(in) :: x_off
         !! Hump centre x-offset from domain centre (m)
      real(wp), intent(in) :: y_off
         !! Hump centre y-offset from domain centre (m)
      real(wp), intent(in) :: u_shear
         !! Initial vertical shear: top +val, bottom -val (m/s)

      integer :: i, j, k
      real(wp) :: x, y, x_c, y_c, r2, u_bar
      real(wp) :: layer_sign
      logical :: has_bathymetry

      x_c = real(s%grid%nx_phys, wp)*s%grid%dx*0.5_wp + x_off
      y_c = real(s%grid%ny_phys, wp)*s%grid%dy*0.5_wp + y_off

      ! Check if bathymetry was already loaded (non-zero values present)
      has_bathymetry = any(s%barotropic%b /= 0.0_wp)

      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            x = (real(i - s%grid%nghost, wp) - 0.5_wp)*s%grid%dx
            y = (real(j - s%grid%nghost, wp) - 0.5_wp)*s%grid%dy
            r2 = (x - x_c)**2 + (y - y_c)**2

            if (.not. has_bathymetry) s%barotropic%b(i, j) = 0.0_wp
            s%barotropic%h(i, j) = h0 + amplitude*exp(-r2/(2.0_wp*sigma**2))
            s%barotropic%hu(i, j) = 0.0_wp
            s%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

      ! NH layer velocities with optional vertical shear
      if (s%use_nonhydrostatic) then
         do k = 1, s%nh%nz
            ! Linear shear: top layer positive, bottom negative
            ! layer_sign goes from -(nz-1)/nz to +(nz-1)/nz
            layer_sign = (2.0_wp*real(k, wp) - real(s%nh%nz + 1, wp))/real(s%nh%nz, wp)
            do j = 1, s%grid%ny_total
               do i = 1, s%grid%nx_total
                  s%nh%u_layer(i, j, k) = u_shear*layer_sign
                  s%nh%v_layer(i, j, k) = 0.0_wp
               end do
            end do
         end do
         ! Set depth-averaged momentum consistent with layer velocities
         ! hu = h * sum(dsig(k) * u_layer(k))
         do j = 1, s%grid%ny_total
            do i = 1, s%grid%nx_total
               u_bar = 0.0_wp
               do k = 1, s%nh%nz
                  u_bar = u_bar + s%nh%dsig(k)*s%nh%u_layer(i, j, k)
               end do
               s%barotropic%hu(i, j) = s%barotropic%h(i, j)*u_bar
            end do
         end do
      end if

   end subroutine init_gaussian_hump

   subroutine init_castro_dam_break(s)
      !! Initialise the circular dam break from Castro et al. (2011):
      !! domain [-5,5]^2, Gaussian bump bathymetry, circular water column
      !! (r <= 0.6: eta = 3, else eta = 1). Public only for the unit-test suite.
      type(state_t), intent(inout) :: s

      integer :: i, j
      real(wp) :: x, y, r2, h_depth, x0, y0

      ! Domain origin: cell centres span [-5, 5] in both directions
      x0 = -0.5_wp*real(s%grid%nx_phys, wp)*s%grid%dx
      y0 = -0.5_wp*real(s%grid%ny_phys, wp)*s%grid%dy

      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            ! Cell centre coordinates
            x = x0 + (real(i - s%grid%nghost, wp) - 0.5_wp)*s%grid%dx
            y = y0 + (real(j - s%grid%nghost, wp) - 0.5_wp)*s%grid%dy

            ! Castro's depth function: H = 1 - 0.4 * exp(-x^2 - y^2)
            h_depth = 1.0_wp - 0.4_wp*exp(-x*x - y*y)

            ! Tidal bathymetry: B = -H (positive up convention)
            s%barotropic%b(i, j) = -h_depth

            ! Radial distance squared from origin
            r2 = x*x + y*y

            ! Initial water depth: h = eta - B
            if (r2 <= 0.36_wp) then
               ! Inside circle (r <= 0.6): free surface eta = 3.0
               s%barotropic%h(i, j) = 3.0_wp - s%barotropic%b(i, j)
            else
               ! Outside circle: free surface eta = 1.0
               s%barotropic%h(i, j) = 1.0_wp - s%barotropic%b(i, j)
            end if

            ! Zero initial momentum
            s%barotropic%hu(i, j) = 0.0_wp
            s%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

   end subroutine init_castro_dam_break

   subroutine init_still_water(s, eta0)
      !! Initialise a flat free surface at elevation eta0:
      !! h = max(eta0 - b, 0) everywhere, zero velocity. Requires
      !! bathymetry already loaded into s%barotropic%b. Public only for tests.
      type(state_t), intent(inout) :: s
      real(wp), intent(in) :: eta0
         !! Initial free-surface elevation (m, same datum as bathymetry)

      integer :: i, j

      do j = 1, s%grid%ny_total
         do i = 1, s%grid%nx_total
            s%barotropic%h(i, j) = max(eta0 - s%barotropic%b(i, j), 0.0_wp)
            s%barotropic%hu(i, j) = 0.0_wp
            s%barotropic%hv(i, j) = 0.0_wp
         end do
      end do

   end subroutine init_still_water

   subroutine init_from_config(s, cfg)
      !! Select and apply initial condition based on config
      type(state_t), intent(inout) :: s
      type(config_t), intent(in) :: cfg

      select case (trim(cfg%init_condition))
      case ("still_water")
         call init_still_water(s, cfg%eta_initial)
      case ("gaussian_hump")
         call init_gaussian_hump(s, cfg%h0, cfg%hump_amplitude, cfg%hump_sigma, &
                                 cfg%hump_x_offset, cfg%hump_y_offset, cfg%init_u_shear)
      case ("castro_dam_break")
         call init_castro_dam_break(s)
      case default
         call init_dam_break(s)
      end select

      ! Initialise multilayer arrays from barotropic state
      if (s%use_multilayer) call init_multilayer_from_barotropic(s, cfg)

   end subroutine init_from_config

   subroutine init_multilayer_from_barotropic(s, cfg)
      !! Cold-start multilayer initialisation from the barotropic state:
      !! lay out the layers, seed the tracers from their IC scalars, then
      !! derive density and the reference-density profile.
      type(state_t), intent(inout) :: s
      type(config_t), intent(in) :: cfg

      call ml_set_min_depth(cfg%ml_min_depth)
      ! Lay out the initial layer distribution from the barotropic state
      ! (builds the VCOORD_ZSTAR_FULL reference z-levels first if needed).
      call s%multilayer%init_layers(s%barotropic, s%grid)
      ! Seed salinity + temperature from their IC scalars (each dispatching
      ! uniform-vs-stratified; S_init_bottom > S_init_surface = salty bed).
      call s%multilayer%seed_tracers(s%grid)
      call refresh_multilayer_derived(s)
      call logger%info("Multilayer initialised: "//to_string(s%multilayer%nz_ml)// &
                       " layers, S = "//to_string(s%multilayer%initial_salinity)// &
                       " PSU, T = "//to_string(s%multilayer%initial_temperature)// &
                       " degC, ml_min_depth = "//to_string(cfg%ml_min_depth)// &
                       " m, BPG = "//trim(cfg%bpg_method))
      if (cfg%ml_min_depth < 0.1_wp*s%multilayer%nz_ml) then
         call logger%warning("ml_min_depth = "//to_string(cfg%ml_min_depth)// &
                             " m may be too small for "//to_string(s%multilayer%nz_ml)// &
                             " layers (min layer = "// &
                             to_string(cfg%ml_min_depth/s%multilayer%nz_ml)//" m). "// &
                             "Consider increasing to "// &
                             to_string(0.1_wp*s%multilayer%nz_ml)//" m if instabilities appear.")
      end if
   end subroutine init_multilayer_from_barotropic

   subroutine refresh_multilayer_derived(s)
      !! Recompute the multilayer fields that are derived from the layer
      !! thicknesses and tracers: the equation of state (`rho_layer`) and
      !! the reference-density profile (`rho_ref`) that the baroclinic
      !! pressure-gradient kernels subtract.
      type(state_t), intent(inout) :: s

      if (s%multilayer%eos_use_shared) then
         call ml_equation_of_state_shared(s%multilayer%h_layer, &
                                          s%multilayer%tracers(s%multilayer%idx_salinity)%hTr, &
                                          s%multilayer%tracers(s%multilayer%idx_temperature)%hTr, &
                                          s%multilayer%rho_layer, s%multilayer%eos, &
                                          s%grid, s%multilayer%nz_ml)
      else
         call ml_equation_of_state(s%multilayer%h_layer, &
                                   s%multilayer%tracers(s%multilayer%idx_salinity)%hTr, &
                                   s%multilayer%tracers(s%multilayer%idx_temperature)%hTr, &
                                   s%multilayer%rho_layer, &
                                   s%multilayer%rho_0, s%multilayer%beta_S, s%multilayer%S_ref, &
                                   s%multilayer%alpha_T, s%multilayer%T_ref, &
                                   s%grid, s%multilayer%nz_ml)
      end if
      call ml_compute_reference_density(s%multilayer%h_layer, s%multilayer%rho_layer, &
                                        s%multilayer%rho_ref, s%multilayer%rho_ref_dz, s%multilayer%nz_ref, &
                                        s%multilayer%rho_0, &
                                        s%grid, s%multilayer%nz_ml)
   end subroutine refresh_multilayer_derived

   subroutine restore_multilayer_after_restart(s)
      !! A restart file restores the prognostic multilayer fields (layer
      !! thicknesses, momenta, tracers) but not the derived ones.  Rebuild
      !! the z*-full reference levels and the density / reference-density
      !! profile so a warm start sees the same state a continuous run does.
      type(state_t), intent(inout) :: s

      call s%multilayer%build_vcoord_reference(s%barotropic, s%grid)
      call refresh_multilayer_derived(s)
   end subroutine restore_multilayer_after_restart

   function compute_total_mass(s) result(mass)
      !! Sum h * dx * dy over interior cells via OpenACC reduction (only
      !! the scalar crosses the PCIe bus). Runs on host or device; NaN /
      !! dry-cell garbage is filtered.
      type(state_t), intent(inout) :: s
      real(wp) :: mass

      integer :: i, j, ng, i_lo, i_hi, j_lo, j_hi
      real(wp) :: cell_area, h_ij, acc

      ng = s%grid%nghost
      i_lo = ng + 1
      i_hi = s%grid%nx_total - ng
      j_lo = ng + 1
      j_hi = s%grid%ny_total - ng
      cell_area = s%grid%dx*s%grid%dy
      acc = 0.0_wp

      ! present_or_copyin makes this runnable in production (s mapped via
      ! solver_enter_data) and in standalone unit tests (nothing mapped):
      ! OpenMP target is strict and demands explicit map() when the data
      ! isn't already on device, unlike OpenACC's implicit copy.
      !$omp parallel do collapse(2) private(h_ij)
      do j = j_lo, j_hi
      do i = i_lo, i_hi
         h_ij = s%barotropic%h(i, j)
         if (h_ij > 0.0_wp .and. h_ij < 1.0e6_wp) then
            acc = acc + h_ij*cell_area
         end if
      end do
      end do
      !$omp end parallel do

      mass = acc
   end function compute_total_mass

   subroutine compute_layer_sums(s, sum_h, sum_hS, sum_hT)
      !! Sum h_layer, hS_layer, hT_layer over interior wet cells and layers,
      !! reducing on device (three scalars cross the PCIe bus). Caller forms
      !! sum_hS/sum_h and sum_hT/sum_h for volume-weighted mean S and T.
      !!
      !! Outer shim: dereferences the tracer registry on host and hands flat
      !! 3D arrays to the inner impl — accessing
      !! `s%multilayer%tracers(it)%hTr(i,j,k)` inside the kernel triggers
      !! CUDA_ERROR_ILLEGAL_ADDRESS on NVHPC (array-of-derived-types
      !! descriptor not reachable device-side).
      type(state_t), intent(inout) :: s
      real(wp), intent(out) :: sum_h, sum_hS, sum_hT

      sum_h = 0.0_wp
      sum_hS = 0.0_wp
      sum_hT = 0.0_wp
      if (.not. s%use_multilayer) return

      call compute_layer_sums_impl(s%multilayer%h_layer, &
                                   s%multilayer%tracers(s%multilayer%idx_salinity)%hTr, &
                                   s%multilayer%tracers(s%multilayer%idx_temperature)%hTr, &
                                   s%grid, s%multilayer%nz_ml, &
                                   sum_h, sum_hS, sum_hT)
   end subroutine compute_layer_sums

   subroutine compute_layer_sums_impl(h_layer, hS_layer, hT_layer, grid, nz, &
                                      sum_h, sum_hS, sum_hT)
      !! Inner reduction kernel for compute_layer_sums; takes flat 3D arrays
      !! so the device kernel avoids the tracer-registry indirection.
      use rki_grid, only: hgrid_t
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: sum_h, sum_hS, sum_hT

      integer :: i, j, k, ng, i_lo, i_hi, j_lo, j_hi
      real(wp) :: h_ijk, acc_h, acc_S, acc_T

      ng = grid%nghost
      i_lo = ng + 1
      i_hi = grid%nx_total - ng
      j_lo = ng + 1
      j_hi = grid%ny_total - ng

      acc_h = 0.0_wp
      acc_S = 0.0_wp
      acc_T = 0.0_wp

      ! Explicit present_or_copyin (as in compute_total_mass) so the kernel
      ! works in both production and unit-test contexts.
      !$omp parallel do collapse(3) private(h_ijk)
      do k = 1, nz
      do j = j_lo, j_hi
      do i = i_lo, i_hi
         h_ijk = h_layer(i, j, k)
         ! Skip NaN / dry-cell garbage with the same guard as mass.
         if (h_ijk > 0.0_wp .and. h_ijk < 1.0e6_wp) then
            acc_h = acc_h + h_ijk
            acc_S = acc_S + hS_layer(i, j, k)
            acc_T = acc_T + hT_layer(i, j, k)
         end if
      end do
      end do
      end do
      !$omp end parallel do

      sum_h = acc_h
      sum_hS = acc_S
      sum_hT = acc_T
   end subroutine compute_layer_sums_impl

   subroutine diagnose_bathymetry(s, cfg)
      !! Report bathymetry diagnostics and warn about potential CFL issues:
      !! extreme depths, steep wet-to-wet gradients, shallow-next-to-deep
      !! cells (velocity blowup risk), and the estimated CFL timestep.
      type(state_t), intent(in) :: s
      type(config_t), intent(in) :: cfg

      integer :: i, j, ng
      real(wp) :: b_min, b_max, h_max, dt_est, dx_min
      real(wp) :: grad, max_grad_i, max_grad_j
      real(wp) :: h_here, h_nbr, depth_ratio, worst_ratio
      integer :: n_wet, n_steep, n_shallow_deep
      integer :: steep_i, steep_j, ratio_i, ratio_j
      real(wp), parameter :: STEEP_THRESHOLD = 10.0_wp   !! m/cell
      real(wp), parameter :: DEPTH_RATIO_WARN = 20.0_wp  !! dimensionless

      ng = s%grid%nghost
      b_min = huge(1.0_wp)
      b_max = -huge(1.0_wp)
      h_max = 0.0_wp
      n_wet = 0
      n_steep = 0
      n_shallow_deep = 0
      max_grad_i = 0.0_wp
      max_grad_j = 0.0_wp
      steep_i = 0
      steep_j = 0
      worst_ratio = 0.0_wp
      ratio_i = 0
      ratio_j = 0

      do j = ng + 1, s%grid%ny_total - ng
         do i = ng + 1, s%grid%nx_total - ng
            ! Basic stats
            if (s%barotropic%b(i, j) < b_min) b_min = s%barotropic%b(i, j)
            if (s%barotropic%b(i, j) > b_max) b_max = s%barotropic%b(i, j)
            if (s%barotropic%h(i, j) > h_max) h_max = s%barotropic%h(i, j)
            if (s%barotropic%h(i, j) > 0.0_wp) n_wet = n_wet + 1

            ! Wet-to-wet gradient in i-direction
            if (s%barotropic%h(i, j) > 0.0_wp .and. i < s%grid%nx_total - ng) then
               if (s%barotropic%h(i + 1, j) > 0.0_wp) then
                  grad = abs(s%barotropic%b(i + 1, j) - s%barotropic%b(i, j))
                  if (grad > max_grad_i) then
                     max_grad_i = grad
                     steep_i = i - ng
                     steep_j = j - ng
                  end if
                  if (grad > STEEP_THRESHOLD) n_steep = n_steep + 1
               end if
            end if

            ! Wet-to-wet gradient in j-direction
            if (s%barotropic%h(i, j) > 0.0_wp .and. j < s%grid%ny_total - ng) then
               if (s%barotropic%h(i, j + 1) > 0.0_wp) then
                  grad = abs(s%barotropic%b(i, j + 1) - s%barotropic%b(i, j))
                  if (grad > max_grad_j) then
                     max_grad_j = grad
                  end if
                  if (grad > STEEP_THRESHOLD) n_steep = n_steep + 1
               end if
            end if

            ! Shallow cell next to deep cell (velocity blowup risk)
            h_here = s%barotropic%h(i, j)
            if (h_here > 0.01_wp .and. h_here < 2.0_wp) then
               ! Check all 4 neighbours
               if (i > ng + 1) then
                  h_nbr = s%barotropic%h(i - 1, j)
                  if (h_nbr > 0.01_wp) then
                     depth_ratio = h_nbr/h_here
                     if (depth_ratio > DEPTH_RATIO_WARN) then
                        n_shallow_deep = n_shallow_deep + 1
                        if (depth_ratio > worst_ratio) then
                           worst_ratio = depth_ratio
                           ratio_i = i - ng
                           ratio_j = j - ng
                        end if
                     end if
                  end if
               end if
               if (i < s%grid%nx_total - ng) then
                  h_nbr = s%barotropic%h(i + 1, j)
                  if (h_nbr > 0.01_wp) then
                     depth_ratio = h_nbr/h_here
                     if (depth_ratio > DEPTH_RATIO_WARN) then
                        n_shallow_deep = n_shallow_deep + 1
                        if (depth_ratio > worst_ratio) then
                           worst_ratio = depth_ratio
                           ratio_i = i - ng
                           ratio_j = j - ng
                        end if
                     end if
                  end if
               end if
            end if
         end do
      end do

      ! Report
      call logger%info("  Bathymetry diagnostics:")
      call logger%info("    Elevation: ["//to_string(b_min)//", "// &
                       to_string(b_max)//"] m")
      call logger%info("    Wet cells: "//to_string(n_wet)//" / "// &
                       to_string((s%grid%nx_total - 2*ng)*(s%grid%ny_total - 2*ng)))
      call logger%info("    Max water depth: "//to_string(h_max)//" m")

      ! Estimated CFL timestep
      dx_min = min(cfg%dx, cfg%dy)
      if (h_max > 0.0_wp) then
         dt_est = cfg%cfl*dx_min/sqrt(GRAVITY*h_max)
         call logger%info("    Est. CFL dt (from depth alone): "// &
                          to_string(dt_est)//" s")
      end if

      ! Gradient warnings
      if (max_grad_i > STEEP_THRESHOLD .or. max_grad_j > STEEP_THRESHOLD) then
         call logger%warning("    Steep wet-to-wet bathymetry gradient detected!")
         call logger%warning("      Max i-gradient: "//to_string(max_grad_i)// &
                             " m/cell at ("//to_string(steep_i)//", "// &
                             to_string(steep_j)//")")
         call logger%warning("      Max j-gradient: "//to_string(max_grad_j)//" m/cell")
         call logger%warning("      Total cells > "//to_string(STEEP_THRESHOLD)// &
                             " m/cell: "//to_string(n_steep))
      else
         call logger%info("    Max gradient: "//to_string(max(max_grad_i, max_grad_j))// &
                          " m/cell (OK)")
      end if

      ! Shallow-deep warnings
      if (n_shallow_deep > 0) then
         call logger%warning("    Shallow cells adjacent to deep cells: "// &
                             to_string(n_shallow_deep))
         call logger%warning("      Worst depth ratio: "//to_string(worst_ratio)// &
                             ":1 at ("//to_string(ratio_i)//", "// &
                             to_string(ratio_j)//")")
         call logger%warning("      This can cause extreme velocities and tiny timesteps")
      end if

   end subroutine diagnose_bathymetry

end module rki_init
