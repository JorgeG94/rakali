!! Surface forcing fields for the structured solver
module rki_surface_forcing
   !! Container for the spatially-varying surface forcing fields applied
   !! by the barotropic and multilayer pipelines: wind stress, atmospheric
   !! pressure, precipitation, net heat flux, and the virtual freshwater
   !! flux (E - P) used as a virtual salt flux on the top layer.
   !!
   !! Distinct from `forcing_t` in `rki_forcing` (the NetCDF reader that
   !! *fills* these); this is the on-state container kernels read each step.
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_mem_report, only: arr_bytes
   implicit none
   private

   public :: surface_forcing_t

   type :: surface_forcing_t
      !! Surface flux fields applied by the barotropic friction kernel
      !! and the multilayer surface-flux kernel.

      ! Spatially-varying fields (nx_total, ny_total)
      real(wp), allocatable :: wind_x(:, :)
         !! Surface wind stress in x-direction (Pa)
      real(wp), allocatable :: wind_y(:, :)
         !! Surface wind stress in y-direction (Pa)
      real(wp), allocatable :: atm_pressure(:, :)
         !! Atmospheric pressure (Pa)
      real(wp), allocatable :: precipitation(:, :)
         !! Precipitation rate (m/s), adds volume to barotropic h
      real(wp), allocatable :: heat_flux(:, :)
         !! Net surface heat flux (W/m^2), positive into ocean.
         !! Applied to top layer hT only; does not change h_layer.
      real(wp), allocatable :: evap_minus_precip(:, :)
         !! Virtual surface freshwater flux (m/s), positive for net
         !! evaporation.  Applied as a virtual salt flux on the top
         !! layer: hS(1) += dt * (E-P) * S_top.  Does not change h_layer.

      ! Uniform fallback wind stress (Pa) — initialises wind_x/wind_y at
      ! startup when no forcing is supplied; unused once those are written.
      real(wp) :: wind_stress_x = 0.0_wp
      real(wp) :: wind_stress_y = 0.0_wp

      logical :: active = .false.
         !! True when forcing arrays contain non-zero data
   contains
      procedure, non_overridable :: init => surface_forcing_init
      procedure, non_overridable :: destroy => surface_forcing_destroy
      procedure, non_overridable :: bytes => surface_forcing_bytes
   end type surface_forcing_t

contains

   subroutine surface_forcing_init(this, grid)
      !! Allocate the surface-forcing arrays sized to the grid, zero-filled.
      !! Uniform fallbacks (wind_stress_x/y) are set by the caller from config.
      class(surface_forcing_t), intent(inout) :: this
      type(hgrid_t), intent(in)    :: grid

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      allocate (this%wind_x(nx, ny), source=0.0_wp)
      allocate (this%wind_y(nx, ny), source=0.0_wp)
      allocate (this%atm_pressure(nx, ny), source=0.0_wp)
      allocate (this%precipitation(nx, ny), source=0.0_wp)
      allocate (this%heat_flux(nx, ny), source=0.0_wp)
      allocate (this%evap_minus_precip(nx, ny), source=0.0_wp)
   end subroutine surface_forcing_init

   subroutine surface_forcing_destroy(this)
      class(surface_forcing_t), intent(inout) :: this

      if (allocated(this%wind_x)) deallocate (this%wind_x)
      if (allocated(this%wind_y)) deallocate (this%wind_y)
      if (allocated(this%atm_pressure)) deallocate (this%atm_pressure)
      if (allocated(this%precipitation)) deallocate (this%precipitation)
      if (allocated(this%heat_flux)) deallocate (this%heat_flux)
      if (allocated(this%evap_minus_precip)) deallocate (this%evap_minus_precip)
      this%active = .false.
   end subroutine surface_forcing_destroy

   pure function surface_forcing_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the surface-forcing slot.
      class(surface_forcing_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%wind_x) + arr_bytes(this%wind_y) &
               + arr_bytes(this%atm_pressure) + arr_bytes(this%precipitation) &
               + arr_bytes(this%heat_flux) + arr_bytes(this%evap_minus_precip)
   end function surface_forcing_bytes

end module rki_surface_forcing
