!! Barotropic (depth-integrated) solution state for the unstructured solver.
module rki_barotropic_state_unstr
   !! Container for the centroid conserved variables, the edge
   !! reconstruction scratch, the per-step accumulators, RK2 saves,
   !! Chapman persistent state, the per-edge mass flux consumed by
   !! tracer CWC, and the small set of barotropic physics scalars
   !! (Manning's `n`, uniform wind stress).
   !!
   !! Lives at `state%barotropic`.  Multilayer code reads `mass_flux_e`
   !! from here (CWC tracer advection); the layer-distribute kernels
   !! also read `stage_cv` / `xmom_cv` / `ymom_cv` / `bed_cv` to
   !! compute new layer states.  Mirrors structured `barotropic_state_t`.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   implicit none
   private

   public :: unstr_barotropic_state_t

   integer, parameter :: EDGES_PER_TRI = 3

   type :: unstr_barotropic_state_t
      !! Centroid + edge conserved variables and the per-step workspace.

      ! ----- centroid conserved variables -----
      real(wp), allocatable :: stage_cv(:)
         !! Stage w = h + b at the centroid (m).
      real(wp), allocatable :: xmom_cv(:)
         !! x-momentum h*u at the centroid (m^2/s).
      real(wp), allocatable :: ymom_cv(:)
         !! y-momentum h*v at the centroid (m^2/s).

      ! ----- bathymetry and derived height -----
      real(wp), allocatable :: bed_cv(:)
         !! Bed elevation at the centroid (m). Static after init.
      real(wp), allocatable :: height_cv(:)
         !! Water depth h = max(stage - bed, 0) (m). Recomputed each step.

      ! ----- reconstructed values at edge midpoints (3, ncells_total) -----
      real(wp), allocatable :: stage_ev(:, :)
         !! Stage at edge midpoints, after second-order extrapolation.
      real(wp), allocatable :: xmom_ev(:, :)
         !! x-momentum at edge midpoints.
      real(wp), allocatable :: ymom_ev(:, :)
         !! y-momentum at edge midpoints.
      real(wp), allocatable :: bed_ev(:, :)
         !! Bed at edge midpoints (= stage_ev - height_ev for well-balancing).
      real(wp), allocatable :: height_ev(:, :)
         !! Water depth at edge midpoints.

      real(wp), allocatable :: chapman_eta_old(:, :)
         !! Persistent previous-step eta for `BC_CHAPMAN`. Shape
         !! `(EDGES_PER_TRI, ncells_total)`. Allocated unconditionally
         !! so the array is GPU-mappable regardless of whether Chapman
         !! is active on a given mesh.

      ! ----- explicit update accumulators -----
      real(wp), allocatable :: stage_eu(:)
         !! d(stage)/dt accumulator (cell-centric, no scatter).
      real(wp), allocatable :: xmom_eu(:)
         !! d(xmom)/dt accumulator.
      real(wp), allocatable :: ymom_eu(:)
         !! d(ymom)/dt accumulator.

      ! ----- semi-implicit friction update -----
      real(wp), allocatable :: xmom_siu(:)
         !! Semi-implicit increment for x-momentum (Manning friction).
      real(wp), allocatable :: ymom_siu(:)
         !! Semi-implicit increment for y-momentum.

      ! ----- RK2 backup of U^n at the start of each step -----
      real(wp), allocatable :: stage_backup(:)
      real(wp), allocatable :: xmom_backup(:)
      real(wp), allocatable :: ymom_backup(:)

      ! ----- velocity-form extrapolation scratch -----
      real(wp), allocatable :: x_centroid_work(:)
         !! Stash for original xmom during velocity-form reconstruction.
      real(wp), allocatable :: y_centroid_work(:)
         !! Stash for original ymom during velocity-form reconstruction.

      ! ----- per-cell max wave speed for CFL -----
      real(wp), allocatable :: max_speed(:)
         !! Maximum wave speed seen on any edge of this cell during the
         !! last flux pass.

      ! ----- per-edge barotropic mass flux -----
      real(wp), allocatable :: mass_flux_e(:, :)
         !! Net mass inflow to this cell across edge i, written by the
         !! central-upwind flux kernel and reused by
         !! `ml_advect_tracer_unstr` for CWC-consistent advection.
         !! Shape `(3, ncells_total)`, units m^3/s, positive INTO the cell.

      ! ----- scalar physics knobs -----
      real(wp) :: manning_n = 0.0_wp
         !! Manning roughness coefficient (s/m^{1/3}). Uniform.
      real(wp) :: wind_tau_x = 0.0_wp
         !! Spatially-uniform wind stress in x (Pa).
      real(wp) :: wind_tau_y = 0.0_wp
         !! Spatially-uniform wind stress in y (Pa).
      real(wp) :: coriolis_f = 0.0_wp
         !! Coriolis parameter f (1/s), typically 2*Omega*sin(lat). Uniform.
         !! Consumed per-layer by `ml_apply_coriolis_unstr` (multilayer path).
   contains
      procedure, non_overridable :: init => unstr_barotropic_init
      procedure, non_overridable :: destroy => unstr_barotropic_destroy
      procedure, non_overridable :: enter_data => unstr_barotropic_enter_data
      procedure, non_overridable :: exit_data => unstr_barotropic_exit_data
   end type unstr_barotropic_state_t

contains

   subroutine unstr_barotropic_init(this, mesh)
      !! Allocate every array at the cell counts in `mesh`. All arrays
      !! are zero-initialised.  Calling on an already-initialised
      !! sub-object is safe -- it destroys first.
      class(unstr_barotropic_state_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in)    :: mesh

      integer :: n

      call this%destroy()

      n = mesh%ncells_total

      allocate (this%stage_cv(n), source=0.0_wp)
      allocate (this%xmom_cv(n), source=0.0_wp)
      allocate (this%ymom_cv(n), source=0.0_wp)
      allocate (this%bed_cv(n), source=0.0_wp)
      allocate (this%height_cv(n), source=0.0_wp)

      allocate (this%stage_ev(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%xmom_ev(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%ymom_ev(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%bed_ev(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%height_ev(EDGES_PER_TRI, n), source=0.0_wp)
      allocate (this%chapman_eta_old(EDGES_PER_TRI, n), source=0.0_wp)

      allocate (this%stage_eu(n), source=0.0_wp)
      allocate (this%xmom_eu(n), source=0.0_wp)
      allocate (this%ymom_eu(n), source=0.0_wp)

      allocate (this%xmom_siu(n), source=0.0_wp)
      allocate (this%ymom_siu(n), source=0.0_wp)

      allocate (this%stage_backup(n), source=0.0_wp)
      allocate (this%xmom_backup(n), source=0.0_wp)
      allocate (this%ymom_backup(n), source=0.0_wp)

      allocate (this%x_centroid_work(n), source=0.0_wp)
      allocate (this%y_centroid_work(n), source=0.0_wp)

      allocate (this%max_speed(n), source=0.0_wp)

      allocate (this%mass_flux_e(EDGES_PER_TRI, n), source=0.0_wp)
   end subroutine unstr_barotropic_init

   subroutine unstr_barotropic_destroy(this)
      class(unstr_barotropic_state_t), intent(inout) :: this

      if (allocated(this%stage_cv)) deallocate (this%stage_cv)
      if (allocated(this%xmom_cv)) deallocate (this%xmom_cv)
      if (allocated(this%ymom_cv)) deallocate (this%ymom_cv)
      if (allocated(this%bed_cv)) deallocate (this%bed_cv)
      if (allocated(this%height_cv)) deallocate (this%height_cv)

      if (allocated(this%stage_ev)) deallocate (this%stage_ev)
      if (allocated(this%xmom_ev)) deallocate (this%xmom_ev)
      if (allocated(this%ymom_ev)) deallocate (this%ymom_ev)
      if (allocated(this%bed_ev)) deallocate (this%bed_ev)
      if (allocated(this%height_ev)) deallocate (this%height_ev)
      if (allocated(this%chapman_eta_old)) deallocate (this%chapman_eta_old)

      if (allocated(this%stage_eu)) deallocate (this%stage_eu)
      if (allocated(this%xmom_eu)) deallocate (this%xmom_eu)
      if (allocated(this%ymom_eu)) deallocate (this%ymom_eu)

      if (allocated(this%xmom_siu)) deallocate (this%xmom_siu)
      if (allocated(this%ymom_siu)) deallocate (this%ymom_siu)

      if (allocated(this%stage_backup)) deallocate (this%stage_backup)
      if (allocated(this%xmom_backup)) deallocate (this%xmom_backup)
      if (allocated(this%ymom_backup)) deallocate (this%ymom_backup)

      if (allocated(this%x_centroid_work)) deallocate (this%x_centroid_work)
      if (allocated(this%y_centroid_work)) deallocate (this%y_centroid_work)

      if (allocated(this%max_speed)) deallocate (this%max_speed)
      if (allocated(this%mass_flux_e)) deallocate (this%mass_flux_e)
   end subroutine unstr_barotropic_destroy

   subroutine unstr_barotropic_enter_data(this)
      !! Map the barotropic allocatables onto the device (the parent
      !! value-type member rides over with `copyin(state)`; only the inner
      !! arrays need explicit attaches).  Read+write fields use `copyin`,
      !! pure-workspace fields use `create`.
      class(unstr_barotropic_state_t), intent(inout) :: this

      ! Conserved + bathymetry + height: read first, then written.  Backup
      ! arrays are also copyin so test setups can pre-fill them.
      !$acc enter data copyin(this%stage_cv, this%xmom_cv, this%ymom_cv, &
      !$acc&                  this%bed_cv, this%height_cv, &
      !$acc&                  this%stage_eu, this%xmom_eu, this%ymom_eu, &
      !$acc&                  this%stage_backup, this%xmom_backup, this%ymom_backup, &
      !$acc&                  this%chapman_eta_old)

      ! Pure scratch / written-before-read: create-only.
      !$acc enter data create(this%x_centroid_work, this%y_centroid_work, &
      !$acc&                  this%stage_ev, this%xmom_ev, this%ymom_ev, &
      !$acc&                  this%bed_ev, this%height_ev, &
      !$acc&                  this%xmom_siu, this%ymom_siu, &
      !$acc&                  this%mass_flux_e)
      ! `max_speed` feeds the CFL reduction; copyin so a defined value
      ! is present if the timestep kernel runs before the flux kernel
      ! (production paths overwrite it every flux pass).
      !$acc enter data copyin(this%max_speed)
   end subroutine unstr_barotropic_enter_data

   subroutine unstr_barotropic_exit_data(this)
      !! Reverse of `enter_data` — copy out host-visible fields, drop
      !! pure scratch.  Components first, parent last (handled by the
      !! caller `state_unstr_exit_data`).
      class(unstr_barotropic_state_t), intent(inout) :: this

      !$acc exit data copyout(this%stage_cv, this%xmom_cv, this%ymom_cv, &
      !$acc&                  this%height_cv, &
      !$acc&                  this%stage_ev, this%xmom_ev, this%ymom_ev, &
      !$acc&                  this%bed_ev, this%height_ev, &
      !$acc&                  this%stage_eu, this%xmom_eu, this%ymom_eu, &
      !$acc&                  this%stage_backup, this%xmom_backup, this%ymom_backup, &
      !$acc&                  this%max_speed)
      !$acc exit data delete(this%chapman_eta_old)

      !$acc exit data delete(this%bed_cv, &
      !$acc&                 this%x_centroid_work, this%y_centroid_work, &
      !$acc&                 this%xmom_siu, this%ymom_siu, &
      !$acc&                 this%mass_flux_e)
   end subroutine unstr_barotropic_exit_data

end module rki_barotropic_state_unstr
