!! Multilayer tracer descriptor
module rki_tracer
   !! Per-tracer state and config for the structured multilayer solver.
   !!
   !! Holds the prognostic field `hTr` (= h * Tr), the RK2 save buffer
   !! `hTr0`, and all scalar config the kernels consume (initial value,
   !! BC inflow value, physical clamp bounds, vertical/horizontal
   !! background diffusivities, and a linear-EOS contribution
   !! `eos_coeff * (Tr - eos_ref)`).
   !!
   !! Lets the solver iterate `do it = 1, size(state%tracers)` and dispatch
   !! the generic `ml_*_tracer` kernels per element.  Special-shape physics
   !! (surface heat/salt flux, sediment) indexes by name instead.
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_mem_report, only: arr_bytes
   implicit none
   private

   public :: tracer_t
   public :: TRACER_BUDGET_NONE, TRACER_BUDGET_HEAT, TRACER_BUDGET_SALT

   integer, parameter :: TRACER_BUDGET_NONE = 0
      !! No conservation-budget contributor slot (passive tracers, ideal
      !! age, pseudo-salt).
   integer, parameter :: TRACER_BUDGET_HEAT = 1
      !! Tracer feeds the OCEAN path's `heat_budget_*` accumulators.
   integer, parameter :: TRACER_BUDGET_SALT = 2
      !! Tracer feeds the OCEAN path's `salt_budget_*` accumulators.

   type :: tracer_t
      ! Identity (used for output/config lookup).  All four strings are
      ! written verbatim onto the NetCDF variable; output loops over the
      ! registry, so a new tracer needs no changes to `rki_output`.
      character(len=32) :: name = ""
      character(len=64) :: long_name = ""
      character(len=16) :: units = ""
      character(len=64) :: standard_name = ""
         !! CF-1.11 standard_name (e.g. "sea_water_salinity").  Empty
         !! string => no standard_name attribute is written.

      ! Prognostic field and RK2 save (nx_total, ny_total, nz_ml)
      real(wp), allocatable :: hTr(:, :, :)
      real(wp), allocatable :: hTr0(:, :, :)

      ! Scalar config
      real(wp) :: tr_init = 0.0_wp
         !! Uniform initial concentration (set by init from cfg)
      real(wp) :: tr_inflow = 0.0_wp
         !! BC inflow concentration for BC_INFLOW / BC_DISCHARGE
      real(wp) :: tr_min = -huge(1.0_wp)
         !! Physical lower clamp bound
      real(wp) :: tr_max = huge(1.0_wp)
         !! Physical upper clamp bound
      real(wp) :: kappa_bg = 0.0_wp
         !! Background vertical diffusivity (m^2/s)
      real(wp) :: hdiff_kappa = 0.0_wp
         !! Horizontal Laplacian diffusivity (m^2/s)

      ! Linear EOS contribution: rho_anomaly += eos_coeff * (Tr - eos_ref)
      ! Salinity: eos_coeff = +beta_S, eos_ref = S_ref
      ! Temperature: eos_coeff = -alpha_T, eos_ref = T_ref
      ! Passive scalar: eos_coeff = 0
      real(wp) :: eos_coeff = 0.0_wp
      real(wp) :: eos_ref = 0.0_wp

      integer :: budget_id = TRACER_BUDGET_NONE
         !! Which conservation-budget contributor slot this tracer's
         !! kernels fill on the OCEAN path: NONE (no budget — passive
         !! tracers, ideal age, pseudo-salt), HEAT (`heat_budget_*`) or
         !! SALT (`salt_budget_*`).  Read HOST-side only, at the seven
         !! budget dispatch blocks; never dereferenced inside a device
         !! loop.  Unused on the coastal paths (no budget slots there) —
         !! left at NONE.

      ! Pipeline opt-outs (skip work for diagnostic/passive tracers)
      logical :: do_horizontal_advection = .true.
      logical :: do_vertical_exchange = .true.
      logical :: do_vertical_diffusion = .true.
      logical :: do_horizontal_diffusion = .true.
      logical :: do_clamp = .true.
   contains
      procedure, non_overridable :: init => tracer_init
      procedure, non_overridable :: destroy => tracer_destroy
      procedure, non_overridable :: bytes => tracer_bytes
   end type tracer_t

contains

   subroutine tracer_init(this, grid, nz)
      !! Allocate hTr / hTr0 at the multilayer grid size, zero-filled.
      !! Caller is responsible for populating the prognostic field
      !! (e.g. `hTr = h_layer * tr_init`) after layer thicknesses are set.
      class(tracer_t), intent(inout) :: this
      type(hgrid_t), intent(in)    :: grid
      integer, intent(in)    :: nz

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      allocate (this%hTr(nx, ny, nz), source=0.0_wp)
      allocate (this%hTr0(nx, ny, nz), source=0.0_wp)
   end subroutine tracer_init

   subroutine tracer_destroy(this)
      class(tracer_t), intent(inout) :: this

      if (allocated(this%hTr)) deallocate (this%hTr)
      if (allocated(this%hTr0)) deallocate (this%hTr0)
   end subroutine tracer_destroy

   pure function tracer_bytes(this) result(nbytes)
      !! Counted allocatable footprint of one tracer (hTr + RK2 save).
      class(tracer_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%hTr) + arr_bytes(this%hTr0)
   end function tracer_bytes

end module rki_tracer
