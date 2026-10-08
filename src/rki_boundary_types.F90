!! Boundary condition types, constants, and utility functions
module rki_boundary_types
   !! Backend-agnostic boundary-condition definitions: `bc_params_t`,
   !! the BC type constants, and helper functions.
   use rki_constants, only: wp, TWO_PI
   use rki_config, only: config_t, MAX_TIDAL_CONSTITUENTS
   implicit none
   private

   public :: bc_params_t
   public :: BC_WALL, BC_OPEN, BC_TIDAL, BC_NESTED
   public :: BC_INFLOW, BC_DISCHARGE, BC_CLAMPED, BC_SPONGE, BC_CHAPMAN
   public :: bc_type_from_string
   public :: bc_init_from_config
   public :: bc_init_chapman
   public :: TWO_PI

   integer, parameter :: BC_WALL = 1
      !! Reflective wall boundary
   integer, parameter :: BC_OPEN = 2
      !! Open (Flather radiation) boundary
   integer, parameter :: BC_TIDAL = 3
      !! Tidal forcing (prescribed elevation) boundary
   integer, parameter :: BC_NESTED = 4
      !! Nested boundary (ghost cells filled by nesting module)
   integer, parameter :: BC_INFLOW = 5
      !! Prescribed normal velocity boundary
   integer, parameter :: BC_DISCHARGE = 6
      !! Prescribed volume flux boundary
   integer, parameter :: BC_CLAMPED = 7
      !! Clamped (Dirichlet) elevation boundary
   integer, parameter :: BC_SPONGE = 8
      !! Sponge relaxation zone (ghost cells use wall BC)
   integer, parameter :: BC_CHAPMAN = 9
      !! Chapman (Orlanski radiation) boundary

   type :: bc_params_t
      !! Boundary condition parameters for all four edges
      integer :: west = BC_WALL
         !! West boundary type
      integer :: east = BC_WALL
         !! East boundary type
      integer :: south = BC_WALL
         !! South boundary type
      integer :: north = BC_WALL
         !! North boundary type
      real(wp) :: tidal_amplitude = 0.0_wp
         !! Tidal elevation amplitude (m) — single-constituent legacy
      real(wp) :: tidal_period = 44712.0_wp
         !! Tidal period (s), default M2 = 12.42 h — single-constituent legacy
      real(wp) :: time = 0.0_wp
         !! Current simulation time (s)

      ! Multi-constituent tidal forcing
      integer :: n_tidal_constituents = 0
         !! Number of active tidal constituents (0 = use legacy single)
      real(wp) :: tidal_amp(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent amplitudes (m)
      real(wp) :: tidal_phase(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent phases (radians)
      real(wp) :: tidal_omega(MAX_TIDAL_CONSTITUENTS) = 0.0_wp
         !! Constituent angular frequencies (rad/s)
      real(wp) :: inflow_velocity = 0.0_wp
         !! Prescribed normal velocity for BC_INFLOW (m/s)
      real(wp) :: inflow_salinity = -1.0_wp
         !! Prescribed salinity for BC_INFLOW (PSU); <0 = zero-gradient
      real(wp) :: inflow_temperature = -999.0_wp
         !! Prescribed temperature for BC_INFLOW (degC); <0 = zero-gradient
      real(wp) :: discharge_flux = 0.0_wp
         !! Prescribed volume flux for BC_DISCHARGE (m^2/s)
      real(wp) :: clamped_elevation = 0.0_wp
         !! Prescribed surface elevation for BC_CLAMPED (m)
      integer :: sponge_width = 0
         !! Sponge layer width in cells
      real(wp) :: sponge_strength = 0.0_wp
         !! Sponge relaxation rate (1/s)
      real(wp) :: sponge_h_ref = 0.0_wp
         !! Sponge reference depth (m)
      real(wp), allocatable :: eta_old_west(:)
         !! Chapman: previous-step eta at west boundary (ny)
      real(wp), allocatable :: eta_old_east(:)
         !! Chapman: previous-step eta at east boundary (ny)
      real(wp), allocatable :: eta_old_south(:)
         !! Chapman: previous-step eta at south boundary (nx)
      real(wp), allocatable :: eta_old_north(:)
         !! Chapman: previous-step eta at north boundary (nx)
      real(wp) :: dx = 1.0_wp
         !! Cell size in x-direction (m), needed for Chapman
      real(wp) :: dy = 1.0_wp
         !! Cell size in y-direction (m), needed for Chapman
      real(wp) :: dt = 0.0_wp
         !! Current timestep (s), needed for Chapman
   end type bc_params_t

contains

   pure function bc_type_from_string(name) result(bc_type)
      !! Convert boundary type string to integer constant
      character(len=*), intent(in) :: name
      integer :: bc_type

      select case (trim(name))
      case ("wall")
         bc_type = BC_WALL
      case ("open")
         bc_type = BC_OPEN
      case ("tidal")
         bc_type = BC_TIDAL
      case ("nested")
         bc_type = BC_NESTED
      case ("inflow")
         bc_type = BC_INFLOW
      case ("discharge")
         bc_type = BC_DISCHARGE
      case ("clamped")
         bc_type = BC_CLAMPED
      case ("sponge")
         bc_type = BC_SPONGE
      case ("chapman")
         bc_type = BC_CHAPMAN
      case default
         bc_type = BC_WALL
      end select

   end function bc_type_from_string

   subroutine bc_init_from_config(bc, cfg)
      !! Initialise boundary condition parameters from config
      type(bc_params_t), intent(inout) :: bc
      type(config_t), intent(in) :: cfg

      bc%west = bc_type_from_string(cfg%bc_west)
      bc%east = bc_type_from_string(cfg%bc_east)
      bc%south = bc_type_from_string(cfg%bc_south)
      bc%north = bc_type_from_string(cfg%bc_north)
      bc%tidal_amplitude = cfg%tidal_amplitude
      bc%tidal_period = cfg%tidal_period
      bc%n_tidal_constituents = cfg%n_tidal_constituents
      bc%tidal_amp = cfg%tidal_amp
      bc%tidal_phase = cfg%tidal_phase
      bc%tidal_omega = cfg%tidal_omega
      bc%inflow_velocity = cfg%inflow_velocity
      bc%inflow_salinity = cfg%inflow_salinity
      bc%inflow_temperature = cfg%inflow_temperature
      bc%discharge_flux = cfg%discharge_flux
      bc%clamped_elevation = cfg%clamped_elevation
      bc%sponge_width = cfg%sponge_width
      bc%sponge_strength = cfg%sponge_strength
      bc%sponge_h_ref = cfg%sponge_h_ref
      bc%dx = cfg%dx
      bc%dy = cfg%dy

   end subroutine bc_init_from_config

   subroutine bc_init_chapman(bc, h, b, nx_total, ny_total, nghost)
      !! Allocate Chapman eta_old arrays for all four edges (so they are
      !! available for unconditional GPU mapping in `solver_enter_data`)
      !! and seed the active edges from current eta. Inactive edges keep
      !! their zero source value -- the kernel never reads them when the
      !! corresponding `bc%<side>` is not `BC_CHAPMAN`, so the zeros are
      !! harmless.
      type(bc_params_t), intent(inout) :: bc
      real(wp), intent(in) :: h(:, :)
         !! Water depth array (nx_total, ny_total)
      real(wp), intent(in) :: b(:, :)
         !! Bottom elevation array (nx_total, ny_total)
      integer, intent(in) :: nx_total
         !! Total cells in x including ghost cells
      integer, intent(in) :: ny_total
         !! Total cells in y including ghost cells
      integer, intent(in) :: nghost
         !! Ghost cell width

      integer :: ii, jj

      if (.not. allocated(bc%eta_old_west)) then
         allocate (bc%eta_old_west(ny_total), source=0.0_wp)
      end if
      if (.not. allocated(bc%eta_old_east)) then
         allocate (bc%eta_old_east(ny_total), source=0.0_wp)
      end if
      if (.not. allocated(bc%eta_old_south)) then
         allocate (bc%eta_old_south(nx_total), source=0.0_wp)
      end if
      if (.not. allocated(bc%eta_old_north)) then
         allocate (bc%eta_old_north(nx_total), source=0.0_wp)
      end if

      if (bc%west == BC_CHAPMAN .or. bc%east == BC_CHAPMAN) then
         do jj = 1, ny_total
            bc%eta_old_west(jj) = h(nghost + 1, jj) + b(nghost + 1, jj)
            bc%eta_old_east(jj) = h(nx_total - nghost, jj) + b(nx_total - nghost, jj)
         end do
      end if
      if (bc%south == BC_CHAPMAN .or. bc%north == BC_CHAPMAN) then
         do ii = 1, nx_total
            bc%eta_old_south(ii) = h(ii, nghost + 1) + b(ii, nghost + 1)
            bc%eta_old_north(ii) = h(ii, ny_total - nghost) + b(ii, ny_total - nghost)
         end do
      end if

   end subroutine bc_init_chapman

end module rki_boundary_types
