!! Barotropic (depth-integrated) solution state for the structured solver.
module rki_barotropic_state
   !! Depth-integrated conserved variables, bathymetry, flux-divergence
   !! and per-face mass-flux scratch, RK2 save buffers, and the scalar
   !! physics params the barotropic kernels read (Manning `n`, Coriolis
   !! `f`).  Lives at `state%barotropic`; multilayer code reads `h, hu,
   !! hv` and `mass_flux_x/y` from here (produced by the barotropic step,
   !! consumed by layer-distribute + CWC tracer advection).
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_mem_report, only: arr_bytes
   implicit none
   private

   public :: barotropic_state_t

   type :: barotropic_state_t
      !! Depth-integrated conserved variables and the workspace they
      !! produce / consume each step.

      ! Conserved variables -- SoA layout: (nx_total, ny_total)
      real(wp), allocatable :: h(:, :)
         !! Water depth (m)
      real(wp), allocatable :: hu(:, :)
         !! x-momentum, h * u (m^2/s)
      real(wp), allocatable :: hv(:, :)
         !! y-momentum, h * v (m^2/s)

      ! Bathymetry (static after initialisation)
      real(wp), allocatable :: b(:, :)
         !! Bottom elevation (m, positive up)

      ! Flux divergence (workspace, recomputed each step)
      real(wp), allocatable :: flux_h(:, :)
         !! Net flux divergence for h
      real(wp), allocatable :: flux_hu(:, :)
         !! Net flux divergence for hu
      real(wp), allocatable :: flux_hv(:, :)
         !! Net flux divergence for hv

      ! Per-face HLL mass flux (workspace).  mass_flux_x(i, j) is the
      ! HLL h-flux at the right face of cell (i, j) -- i.e. between
      ! cells (i, j) and (i+1, j).  mass_flux_y(i, j) is the top-face
      ! flux between cells (i, j) and (i, j+1).  Reused by the tracer
      ! advection kernels to enforce consistency-with-continuity (CWC).
      real(wp), allocatable :: mass_flux_x(:, :)
      real(wp), allocatable :: mass_flux_y(:, :)

      ! RK2 storage (saves U^n for the averaging step)
      real(wp), allocatable :: h0(:, :)
         !! Saved h at start of RK step
      real(wp), allocatable :: hu0(:, :)
         !! Saved hu at start of RK step
      real(wp), allocatable :: hv0(:, :)
         !! Saved hv at start of RK step

      ! Physics scalars used by the barotropic friction / Coriolis
      ! kernels.  Multilayer-side per-layer bottom-friction (when
      ! enabled) reads `manning_n` from here too.
      real(wp) :: manning_n = 0.0_wp
         !! Manning roughness coefficient (s/m^{1/3})
      logical :: use_linear_bottom_drag = .false.
         !! Use linear (Rayleigh) bottom drag instead of quadratic Manning.
      real(wp) :: bottom_drag_r = 0.0_wp
         !! Linear bottom-drag rate (1/s); used when use_linear_bottom_drag.
      logical :: wind_surface_layer = .false.
         !! Apply wind stress to the top layer (surface-concentrated) instead
         !! of the depth-uniform barotropic mode.  Multilayer only.
      real(wp) :: coriolis_f = 0.0_wp
         !! Coriolis parameter f (1/s), typically 2*Omega*sin(lat)

      ! Conservation-budget accumulators (host scalars; not device-mapped).
      real(wp) :: mass_out = 0.0_wp
         !! Cumulative net volume that has left the domain through its
         !! boundaries since t=0 (m³, positive = outflow), accumulated per
         !! RK stage from the scheme's own flux divergence so the console
         !! conservation `Error` closes to round-off even with open BCs.
         !! See `rki_console_stats` / the driver's budget wiring.
      real(wp) :: mass_src = 0.0_wp
         !! Cumulative volume ADDED to the domain by surface sources since
         !! t=0 (m³, positive = added), accumulated per RK stage from
         !! precipitation (`h += dt·P`), the evaporative mass sink
         !! (`evap_mass_sink=.true.`; negative contribution), and the sponge
         !! h-relaxation Δh.  Feeds `budget%mass_src` in `rki_console_stats`
         !! so the residual `(total−ref) + out − src` stays ~round-off on
         !! forced runs instead of growing monotonically with the forcing.
      logical :: mass_out_tracked = .false.
         !! Set `.true.` once the solver has accumulated `mass_out`.
         !! Covered paths: 2D-barotropic SSP-RK2 fused path + coastal
         !! multilayer nz>1.  Not set by: nz=1 ML fused fast path (safe
         !! raw-drift fallback — see solver note), semi-implicit, generic-RK.
      logical :: mass_src_tracked = .false.
         !! Set `.true.` once the solver has accumulated `mass_src` for the
         !! current run.  Same coverage as `mass_out_tracked`: 2D-barotropic
         !! fused path + coastal multilayer nz>1.  The closed budget
         !! (`budget%mass_active`) is enabled in the console only when BOTH
         !! flags are set.

      ! Scheme flags
      logical :: use_hllc = .false.
         !! Use HLLC (with contact wave for transverse momentum) instead
         !! of HLL for the barotropic Riemann flux.  Structured-only —
         !! the unstructured solver always uses KNP central-upwind, which
         !! resolves contact discontinuities natively.  Default `.false.`.
   contains
      procedure, non_overridable :: init => barotropic_state_init
      procedure, non_overridable :: destroy => barotropic_state_destroy
      procedure, non_overridable :: bytes => barotropic_state_bytes
   end type barotropic_state_t

contains

   subroutine barotropic_state_init(this, grid)
      !! Allocate all arrays at the grid size, zero-filled.  Scalars
      !! (`manning_n`, `coriolis_f`) are populated by `state_init_from_config`
      !! from the namelist after this routine returns.
      class(barotropic_state_t), intent(inout) :: this
      type(hgrid_t), intent(in)    :: grid

      integer :: nx, ny

      nx = grid%nx_total
      ny = grid%ny_total

      allocate (this%h(nx, ny), source=0.0_wp)
      allocate (this%hu(nx, ny), source=0.0_wp)
      allocate (this%hv(nx, ny), source=0.0_wp)
      allocate (this%b(nx, ny), source=0.0_wp)
      allocate (this%flux_h(nx, ny), source=0.0_wp)
      allocate (this%flux_hu(nx, ny), source=0.0_wp)
      allocate (this%flux_hv(nx, ny), source=0.0_wp)
      allocate (this%mass_flux_x(nx, ny), source=0.0_wp)
      allocate (this%mass_flux_y(nx, ny), source=0.0_wp)
      allocate (this%h0(nx, ny), source=0.0_wp)
      allocate (this%hu0(nx, ny), source=0.0_wp)
      allocate (this%hv0(nx, ny), source=0.0_wp)
   end subroutine barotropic_state_init

   subroutine barotropic_state_destroy(this)
      class(barotropic_state_t), intent(inout) :: this

      if (allocated(this%h)) deallocate (this%h)
      if (allocated(this%hu)) deallocate (this%hu)
      if (allocated(this%hv)) deallocate (this%hv)
      if (allocated(this%b)) deallocate (this%b)
      if (allocated(this%flux_h)) deallocate (this%flux_h)
      if (allocated(this%flux_hu)) deallocate (this%flux_hu)
      if (allocated(this%flux_hv)) deallocate (this%flux_hv)
      if (allocated(this%mass_flux_x)) deallocate (this%mass_flux_x)
      if (allocated(this%mass_flux_y)) deallocate (this%mass_flux_y)
      if (allocated(this%h0)) deallocate (this%h0)
      if (allocated(this%hu0)) deallocate (this%hu0)
      if (allocated(this%hv0)) deallocate (this%hv0)
      this%manning_n = 0.0_wp
      this%use_linear_bottom_drag = .false.
      this%bottom_drag_r = 0.0_wp
      this%wind_surface_layer = .false.
      this%coriolis_f = 0.0_wp
   end subroutine barotropic_state_destroy

   pure function barotropic_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the barotropic slot (0 for any
      !! unallocated array).  One `arr_bytes` term per allocatable — add a
      !! term here whenever a new array joins the type.
      class(barotropic_state_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%h) + arr_bytes(this%hu) + arr_bytes(this%hv) &
               + arr_bytes(this%b) + arr_bytes(this%flux_h) + arr_bytes(this%flux_hu) &
               + arr_bytes(this%flux_hv) + arr_bytes(this%mass_flux_x) &
               + arr_bytes(this%mass_flux_y) + arr_bytes(this%h0) &
               + arr_bytes(this%hu0) + arr_bytes(this%hv0)
   end function barotropic_state_bytes

end module rki_barotropic_state
