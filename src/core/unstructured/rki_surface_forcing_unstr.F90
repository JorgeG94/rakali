!! Surface forcing fields for the unstructured solver
module rki_surface_forcing_unstr
   !! Container for the spatially-varying surface forcing fields applied
   !! by the unstructured multilayer pipeline.  Currently only the two
   !! tracer-side forcings are supported:
   !!
   !!   `heat_flux`         net surface heat flux (W/m^2, positive into
   !!                       ocean) applied to the top layer's hT field
   !!                       via `hT(1,c) += dt * Q / (rho_0 * cp)`.
   !!
   !!   `evap_minus_precip` virtual surface freshwater flux (m/s, positive
   !!                       for net evaporation) applied as a virtual
   !!                       salt flux on the top layer:
   !!                       `hS(1,c) += dt * (E-P) * S_top`.
   !!                       Does not change `h_layer` -- pure tracer
   !!                       forcing in the ROMS / MOM convention.
   !!
   !!   `wind_x` / `wind_y`   per-cell surface wind STRESS (Pa) added to the
   !!                       barotropic momentum in the fused update kernel:
   !!                       `hu += dt * wind_x / RHO_WATER`, on top of the
   !!                       spatially-uniform `barotropic%wind_tau_*`.
   !!                       Pushed per owned cell via `rki_unstr_set_wind_field`.
   !!
   !! Atmospheric pressure / precipitation will be added when those ports
   !! land; this container will gain them to mirror `surface_forcing_t`.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   implicit none
   private

   public :: unstr_surface_forcing_t

   type :: unstr_surface_forcing_t
      !! Surface flux fields applied by the unstructured multilayer
      !! surface-flux kernel.  Sized to `mesh%ncells_total` so halo
      !! cells carry whatever the user pushed in (forcing is read-only
      !! from the kernel's perspective and only owned cells are written
      !! into tracer fields, so halo content is irrelevant for
      !! correctness -- but the allocation matches the rest of the
      !! state for consistency).

      real(wp), allocatable :: heat_flux(:)
         !! Net surface heat flux (W/m^2), positive into ocean.
         !! Applied to top layer hT only; does not change h_layer.
      real(wp), allocatable :: evap_minus_precip(:)
         !! Virtual surface freshwater flux (m/s), positive for net
         !! evaporation.  Applied as a virtual salt flux on the top
         !! layer: hS(1,c) += dt * (E-P) * S_top.
      real(wp), allocatable :: wind_x(:)
         !! Per-cell surface wind stress, x-component (Pa).  Added to the
         !! barotropic momentum (on top of the uniform barotropic%wind_tau_x)
         !! by the fused update kernel.  Zero => no per-cell contribution.
      real(wp), allocatable :: wind_y(:)
         !! Per-cell surface wind stress, y-component (Pa).

      logical :: wind_active = .false.
         !! True once a non-zero per-cell wind field has been pushed via
         !! `rki_unstr_set_wind_field`.  Gates the per-cell wind read in the
         !! fused kernel so the no-wind hot path skips it entirely.
      logical :: active = .false.
         !! True when forcing arrays contain non-zero data.  Currently
         !! advisory only -- the kernel runs unconditionally because
         !! zero forcing is a bit-exact no-op.
   contains
      procedure, non_overridable :: init => unstr_surface_forcing_init
      procedure, non_overridable :: destroy => unstr_surface_forcing_destroy
      procedure, non_overridable :: enter_data => unstr_surface_forcing_enter_data
      procedure, non_overridable :: exit_data => unstr_surface_forcing_exit_data
   end type unstr_surface_forcing_t

contains

   subroutine unstr_surface_forcing_init(this, mesh)
      !! Allocate the surface-forcing arrays sized to the mesh,
      !! zero-initialised.  Always safe to call -- destroys first.
      class(unstr_surface_forcing_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in)    :: mesh

      integer :: n

      call this%destroy()

      n = mesh%ncells_total

      allocate (this%heat_flux(n), source=0.0_wp)
      allocate (this%evap_minus_precip(n), source=0.0_wp)
      allocate (this%wind_x(n), source=0.0_wp)
      allocate (this%wind_y(n), source=0.0_wp)
   end subroutine unstr_surface_forcing_init

   subroutine unstr_surface_forcing_destroy(this)
      class(unstr_surface_forcing_t), intent(inout) :: this

      if (allocated(this%heat_flux)) deallocate (this%heat_flux)
      if (allocated(this%evap_minus_precip)) deallocate (this%evap_minus_precip)
      if (allocated(this%wind_x)) deallocate (this%wind_x)
      if (allocated(this%wind_y)) deallocate (this%wind_y)
      this%active = .false.
      this%wind_active = .false.
   end subroutine unstr_surface_forcing_destroy

   subroutine unstr_surface_forcing_enter_data(this)
      !! Map the forcing allocatables onto the device.  The parent
      !! `unstr_state_t%forcing` value-type member rides over with the
      !! parent `copyin(state)`; only the inner allocatable arrays need
      !! explicit attaches.  Both fields are read-only in the kernel
      !! and updated from host via `set_*` setters, so `copyin` is the
      !! right initial direction.
      class(unstr_surface_forcing_t), intent(inout) :: this

      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
   end subroutine unstr_surface_forcing_enter_data

   subroutine unstr_surface_forcing_exit_data(this)
      !! Counterpart to `enter_data`.  No copy-out: forcing is set from
      !! host, never produced on device, so the device copies have no
      !! information the host doesn't already have.
      class(unstr_surface_forcing_t), intent(inout) :: this

      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
   end subroutine unstr_surface_forcing_exit_data

end module rki_surface_forcing_unstr
