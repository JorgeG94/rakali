!! Equation of state and reference-density profile for the structured
!! multilayer solver.
module rki_ml_eos
   !! Linear two-tracer EOS `rho = rho_0 + beta_S*(S-S_ref) - alpha_T*(T-T_ref)`
   !! and the horizontally-uniform `rho_ref(z)` profile used by the
   !! Jacobian baroclinic pressure gradient.
   use rki_constants, only: wp, DRY_TOLERANCE
   use rki_grid, only: hgrid_t
   use rki_ml_common, only: ML_MIN_DEPTH
   use rki_eos, only: eos_t, eos_compute_arrays, &
                      EOS_VARIANT_LINEAR, EOS_VARIANT_WRIGHT_97, &
                      EOS_VARIANT_ROQUET_SPV
   implicit none
   private

   public :: ml_equation_of_state
   public :: ml_equation_of_state_shared
   public :: ml_parse_eos_form
   public :: ml_compute_reference_density

contains

   subroutine ml_parse_eos_form(name, use_shared, variant)
      !! Translate the coastal `&tracer_nml eos_form` knob into the
      !! (use_shared, EOS_VARIANT_*) pair held on the multilayer state.
      !!
      !! `"linear"` (the default) keeps the legacy coastal linear kernel
      !! (`ml_equation_of_state` / `ml_equation_of_state_unstr`) — the
      !! bit-identity-by-construction path.  `"wright"` / `"roquet_spv"`
      !! opt in to the shared nonlinear EOS (`rki_eos`) via
      !! `ml_equation_of_state_shared{,_unstr}`.  The shared
      !! `EOS_VARIANT_LINEAR` is deliberately NOT reachable from this
      !! knob: it is a different linear formulation from the coastal
      !! kernel and exposing both invites silent near-miss configs.
      !! Unknown strings fail loud (configure time).
      character(len=*), intent(in) :: name
      logical, intent(out) :: use_shared
      integer, intent(out) :: variant

      select case (trim(adjustl(name)))
      case ("linear", "")
         use_shared = .false.
         variant = EOS_VARIANT_LINEAR
      case ("wright", "wright_97", "wright97")
         use_shared = .true.
         variant = EOS_VARIANT_WRIGHT_97
      case ("roquet_spv", "roquet")
         use_shared = .true.
         variant = EOS_VARIANT_ROQUET_SPV
      case default
         error stop "ml_parse_eos_form: unknown eos_form (use linear | "// &
            "wright | roquet_spv)"
      end select
   end subroutine ml_parse_eos_form

   subroutine ml_equation_of_state_shared(h_layer, hS_layer, hT_layer, &
                                          rho_layer, eos, grid, nz)
      !! Opt-in coastal entry to the SHARED (ocean) equation of state —
      !! selected by `&tracer_nml eos_form = "wright" | "roquet_spv"`.
      !! Thin forwarding shim onto the regime-agnostic
      !! `eos_compute_arrays` variant dispatch; evaluation is at
      !! `eos%p_ref` (default 0 Pa ⇒ surface/potential density — the
      !! right choice for the coastal depth range) with the shared
      !! `H_VANISHED` vanishing-layer fallback to `eos%rho0` (note: the
      !! legacy coastal kernel gates on `DRY_TOLERANCE` instead).  The
      !! `eos_t` handle is a flat POD — its scalars are read on
      !! the host here and passed by value into the device kernels, so
      !! no device mapping of the handle is required.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      type(eos_t), intent(in) :: eos
         !! Shared EOS handle (variant + scalar coefficients).

      call eos_compute_arrays(eos, h_layer, hS_layer, hT_layer, &
                              rho_layer, grid%nx_total, grid%ny_total, nz)
   end subroutine ml_equation_of_state_shared

   subroutine ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                        rho_0, beta_S, S_ref, &
                                        alpha_T, T_ref, &
                                        grid, nz)
      !! Linear two-tracer equation of state:
      !!   rho_k = rho_0 + beta_S * (S_k - S_ref) - alpha_T * (T_k - T_ref)
      !! where S_k = hS_k / h_k and T_k = hT_k / h_k.
      !!
      !! Setting alpha_T = 0 makes temperature a passive tracer (no density
      !! feedback).  Setting S_ref = 0 reproduces the one-tracer legacy form
      !! rho = rho_0 + beta_S*S.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hS_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hT_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_0
         !! Reference density (kg/m^3)
      real(wp), intent(in) :: beta_S
         !! Haline contraction coefficient (kg/m^3 per PSU)
      real(wp), intent(in) :: S_ref
         !! Reference salinity (PSU)
      real(wp), intent(in) :: alpha_T
         !! Thermal expansion coefficient (kg/m^3 per degC)
      real(wp), intent(in) :: T_ref
         !! Reference temperature (degC)

      integer :: i, j, k, nx, ny

      real(wp) :: S_k, T_k, inv_h

      nx = grid%nx_total
      ny = grid%ny_total

      ! Compute rho for ALL cells including ghost cells, because the
      ! baroclinic pressure gradient uses central differences that access
      ! ghost cell rho values.  Ghost cells already have valid h_layer
      ! (from ml_apply_layer_bcs) and hS/hT_layer (from ml_apply_tracer_bcs).
      !$omp parallel do collapse(3) private(S_k, T_k, inv_h)
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         if (h_layer(i, j, k) > DRY_TOLERANCE) then
            inv_h = 1.0_wp/h_layer(i, j, k)
            S_k = hS_layer(i, j, k)*inv_h
            T_k = hT_layer(i, j, k)*inv_h
         else
            S_k = S_ref
            T_k = T_ref
         end if
         rho_layer(i, j, k) = rho_0 + beta_S*(S_k - S_ref) - alpha_T*(T_k - T_ref)
      end do
      end do
      end do
      !$omp end parallel do

   end subroutine ml_equation_of_state

   subroutine ml_compute_reference_density(h_layer, rho_layer, &
                                           rho_ref, rho_ref_dz, nz_ref, &
                                           rho_0, &
                                           grid, nz)
      !! Build the horizontally-uniform reference-density profile
      !! `rho_ref(z)` from the current `rho_layer`.  Sampled on a uniform
      !! z-grid `0..rho_ref_dz*nz_ref`, area-weighted across interior cells.
      !!
      !! Subtracted by both BPG kernels so the integrand becomes the
      !! density anomaly `rho - rho_ref(z)` -- 10-100x smaller than
      !! `rho - rho_0` in stratified water (Phase B2).
      !!
      !! Computed on the host: copies layer fields back from the device,
      !! builds the profile, and pushes `rho_ref` to the device.  Single
      !! one-time D2H+H2D, negligible at init.
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz, nz_ref
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: rho_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(out) :: rho_ref(0:nz_ref)
      real(wp), intent(out) :: rho_ref_dz
      real(wp), intent(in) :: rho_0

      integer :: i, j, k, iz, ncells, nx, ny, nghost
      real(wp) :: H_max, H_col, dz, z, z_top, z_bot, rho_at_z
      real(wp), allocatable :: sum_rho(:)

      nx = grid%nx_total
      ny = grid%ny_total
      nghost = grid%nghost

      ! [acc->omp-cpu] dropped device directive (host no-op): update

      H_max = 0.0_wp
      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            H_col = 0.0_wp
            do k = 1, nz
               H_col = H_col + h_layer(i, j, k)
            end do
            if (H_col > H_max) H_max = H_col
         end do
      end do

      if (H_max <= 0.0_wp) then
         rho_ref = rho_0
         rho_ref_dz = 0.0_wp
         ! [acc->omp-cpu] dropped device directive (host no-op): update
         return
      end if

      H_max = 1.05_wp*H_max
      dz = H_max/real(nz_ref, wp)
      rho_ref_dz = dz

      allocate (sum_rho(0:nz_ref), source=0.0_wp)
      ncells = 0

      do j = nghost + 1, ny - nghost
         do i = nghost + 1, nx - nghost
            H_col = 0.0_wp
            do k = 1, nz
               H_col = H_col + h_layer(i, j, k)
            end do
            if (H_col < ML_MIN_DEPTH) cycle

            ncells = ncells + 1

            do iz = 0, nz_ref
               z = real(iz, wp)*dz
               z_top = 0.0_wp
               rho_at_z = rho_layer(i, j, nz)
               if (z >= H_col) then
                  rho_at_z = rho_layer(i, j, 1)
               else
                  do k = nz, 1, -1
                     z_bot = z_top + h_layer(i, j, k)
                     if (z >= z_top .and. z <= z_bot) then
                        rho_at_z = rho_layer(i, j, k)
                        exit
                     end if
                     z_top = z_bot
                  end do
               end if
               sum_rho(iz) = sum_rho(iz) + rho_at_z
            end do
         end do
      end do

      if (ncells > 0) then
         do iz = 0, nz_ref
            rho_ref(iz) = sum_rho(iz)/real(ncells, wp)
         end do
      else
         rho_ref = rho_0
      end if

      deallocate (sum_rho)

      ! [acc->omp-cpu] dropped device directive (host no-op): update
   end subroutine ml_compute_reference_density

end module rki_ml_eos
