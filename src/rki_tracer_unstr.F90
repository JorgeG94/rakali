!! Multilayer tracer descriptor (unstructured grid).
module rki_tracer_unstr
   !! Per-tracer state and config for the unstructured multilayer
   !! solver.  Mirrors `rki_tracer` (structured) field-for-field; only
   !! the array axis order differs.  Cell-centric layout puts `nz`
   !! first to match the rest of the unstructured multilayer fields
   !! (`h_layer(nz, ncells_total)` etc.).
   !!
   !! Holds the prognostic field `hTr` (= h * Tr), the RK2 save buffer
   !! `hTr0`, and all scalar config the kernels consume.  Designed so
   !! the solver can iterate `do it = 1, size(state%tracers)` and
   !! dispatch the generic `ml_*_tracer_unstr` kernels per element.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_boundary_unstr, only: MAX_BC_GROUPS
   implicit none
   private

   public :: unstr_tracer_t

   type :: unstr_tracer_t
      ! Identity (used for output/config lookup).  Mirrors structured
      ! `tracer_t`: output loops over the registry and writes these
      ! verbatim onto the NetCDF variable.
      character(len=32) :: name = ""
      character(len=64) :: long_name = ""
      character(len=16) :: units = ""
      character(len=64) :: standard_name = ""
         !! CF-1.11 standard_name (empty string => no attribute written).

      ! Prognostic field and RK2 save (nz, ncells_total)
      real(wp), allocatable :: hTr(:, :)
      real(wp), allocatable :: hTr0(:, :)

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
      real(wp) :: bc_ext(MAX_BC_GROUPS) = 0.0_wp
         !! Per-BC-group inflow concentration consumed by the advection
         !! kernel (`ml_advect_tracer{,_muscl}_unstr`).  For salinity /
         !! temperature the solver syncs this from `bc_S_ext` / `bc_T_ext`
         !! each step; a registered passive tracer gets its scalar
         !! `tr_inflow` broadcast across all groups at registration.  Read
         !! host-side (copied into the kernel-local before the device
         !! loop), so it needs no separate device mapping.

      ! Linear EOS contribution: rho_anomaly += eos_coeff * (Tr - eos_ref)
      ! Salinity:    eos_coeff = +beta_S,  eos_ref = S_ref
      ! Temperature: eos_coeff = -alpha_T, eos_ref = T_ref
      ! Passive scalar: eos_coeff = 0
      real(wp) :: eos_coeff = 0.0_wp
      real(wp) :: eos_ref = 0.0_wp

      ! Pipeline opt-outs (skip work for diagnostic/passive tracers)
      logical :: do_horizontal_advection = .true.
      logical :: do_vertical_exchange = .true.
      logical :: do_vertical_diffusion = .true.
      logical :: do_horizontal_diffusion = .true.
      logical :: do_clamp = .true.
   contains
      procedure, non_overridable :: init => unstr_tracer_init
      procedure, non_overridable :: destroy => unstr_tracer_destroy
   end type unstr_tracer_t

contains

   subroutine unstr_tracer_init(this, mesh, nz)
      !! Allocate hTr / hTr0 at the multilayer mesh size, zero-filled.
      !! Caller is responsible for populating the prognostic field
      !! (e.g. `hTr = h_layer * tr_init`) after layer thicknesses are set.
      class(unstr_tracer_t), intent(inout) :: this
      type(mesh_unstr_t), intent(in)    :: mesh
      integer, intent(in)    :: nz

      allocate (this%hTr(nz, mesh%ncells_total), source=0.0_wp)
      allocate (this%hTr0(nz, mesh%ncells_total), source=0.0_wp)
   end subroutine unstr_tracer_init

   subroutine unstr_tracer_destroy(this)
      class(unstr_tracer_t), intent(inout) :: this

      if (allocated(this%hTr)) deallocate (this%hTr)
      if (allocated(this%hTr0)) deallocate (this%hTr0)
   end subroutine unstr_tracer_destroy

end module rki_tracer_unstr
