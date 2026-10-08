!! Non-hydrostatic (multi-layer pressure correction) solution state.
module rki_nonhydrostatic_state
   !! Per-layer NH prognostics (`u_layer`, `v_layer`, `w`, `q_nh`) + RK2
   !! saves, CG Poisson workspaces, barotropic pre-HLL momentum saves,
   !! mode-splitting time-average + barotropic-tendency arrays, and the
   !! NH/CG/split scalar config.  Lives at `state%nh`; allocated by
   !! `state_init` only when `use_nonhydrostatic = .true.` (else scalar
   !! defaults).  Mode-splitting arrays further gated on `split_ratio /= 1`.
   use, intrinsic :: iso_fortran_env, only: int64
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_mem_report, only: arr_bytes
   implicit none
   private

   public :: nonhydrostatic_state_t

   type :: nonhydrostatic_state_t
      !! Per-layer NH prognostics + CG workspaces + mode-split buffers.

      ! Layer count
      integer :: nz = 1
         !! Number of NH sigma layers

      ! Sigma vertical coordinate
      real(wp), allocatable :: dsig(:)
         !! Sigma layer thicknesses (nz), sum = 1

      ! Per-layer prognostic velocities + pressure (nx_total, ny_total, nz)
      real(wp), allocatable :: u_layer(:, :, :)
         !! Per-layer x-velocity
      real(wp), allocatable :: v_layer(:, :, :)
         !! Per-layer y-velocity
      real(wp), allocatable :: w(:, :, :)
         !! Vertical velocity at interfaces (nx_total, ny_total, nz+1)
      real(wp), allocatable :: q_nh(:, :, :)
         !! Non-hydrostatic pressure

      ! RK2 saves
      real(wp), allocatable :: u_layer0(:, :, :)
      real(wp), allocatable :: v_layer0(:, :, :)
      real(wp), allocatable :: w0(:, :, :)

      ! CG Poisson solver workspaces (nx_total, ny_total, nz)
      real(wp), allocatable :: cg_r(:, :, :)
         !! Residual
      real(wp), allocatable :: cg_p(:, :, :)
         !! Search direction
      real(wp), allocatable :: cg_Ap(:, :, :)
         !! Operator application
      real(wp), allocatable :: cg_z(:, :, :)
         !! Preconditioned residual (PCG only)

      ! Workspace: barotropic momentum before HLL update
      real(wp), allocatable :: hu_pre(:, :)
      real(wp), allocatable :: hv_pre(:, :)

      ! ----- CG / NH scalar config -----
      real(wp) :: nh_tolerance = 1.0e-6_wp
         !! CG convergence tolerance
      integer :: nh_max_iterations = 100
         !! Maximum CG iterations
      logical :: nh_preconditioner = .false.
         !! Enable vertical-line preconditioner for CG solver
      integer :: nh_solve_interval = 1
         !! Solve NH pressure every N steps (1 = every step)

      ! ----- Mode-splitting config + arrays -----
      ! Allocated only when `split_ratio /= 1`.  N fast barotropic
      ! substeps per slow baroclinic step, ROMS-style time averaging.
      integer :: split_ratio = 1
         !! Barotropic-baroclinic split: 0=auto, 1=unsplit, N>1=fixed
      logical :: split_averaging = .true.
         !! ROMS-style time averaging over barotropic substeps
      real(wp), allocatable :: h_avg(:, :)
         !! Time-averaged barotropic depth
      real(wp), allocatable :: hu_avg(:, :)
         !! Time-averaged barotropic x-momentum
      real(wp), allocatable :: hv_avg(:, :)
         !! Time-averaged barotropic y-momentum
      real(wp), allocatable :: h_bt0(:, :)
         !! Barotropic RK2 saved depth
      real(wp), allocatable :: hu_bt0(:, :)
         !! Barotropic RK2 saved x-momentum
      real(wp), allocatable :: hv_bt0(:, :)
         !! Barotropic RK2 saved y-momentum
      real(wp), allocatable :: baro_tend_x(:, :)
         !! NH barotropic tendency in x
      real(wp), allocatable :: baro_tend_y(:, :)
         !! NH barotropic tendency in y

   contains
      procedure, non_overridable :: init => nh_state_init
      procedure, non_overridable :: destroy => nh_state_destroy
      procedure, non_overridable :: bytes => nh_state_bytes
   end type nonhydrostatic_state_t

contains

   subroutine nh_state_init(this, grid)
      !! Allocate all NH arrays at the grid size and the configured
      !! layer count.  Mode-splitting arrays are allocated only when
      !! `split_ratio /= 1`.  PCG preconditioner workspace `cg_z` is
      !! always allocated; the OpenACC `enter data` directive is
      !! conditional on `nh_preconditioner` so unused PCG state isn't
      !! mapped to the device.
      !!
      !! Caller must have set `this%nz`, the CG/NH scalar config, and
      !! `this%split_ratio` / `this%split_averaging` before calling.
      class(nonhydrostatic_state_t), intent(inout) :: this
      type(hgrid_t), intent(in)    :: grid

      integer :: nx, ny, nz

      nx = grid%nx_total
      ny = grid%ny_total
      nz = this%nz

      allocate (this%dsig(nz))
      this%dsig = 1.0_wp/real(nz, wp)

      allocate (this%u_layer(nx, ny, nz), source=0.0_wp)
      allocate (this%v_layer(nx, ny, nz), source=0.0_wp)
      allocate (this%w(nx, ny, nz + 1), source=0.0_wp)
      allocate (this%q_nh(nx, ny, nz), source=0.0_wp)

      allocate (this%u_layer0(nx, ny, nz), source=0.0_wp)
      allocate (this%v_layer0(nx, ny, nz), source=0.0_wp)
      allocate (this%w0(nx, ny, nz + 1), source=0.0_wp)

      allocate (this%cg_r(nx, ny, nz), source=0.0_wp)
      allocate (this%cg_p(nx, ny, nz), source=0.0_wp)
      allocate (this%cg_Ap(nx, ny, nz), source=0.0_wp)
      allocate (this%cg_z(nx, ny, nz), source=0.0_wp)

      allocate (this%hu_pre(nx, ny), source=0.0_wp)
      allocate (this%hv_pre(nx, ny), source=0.0_wp)

      ! Mode-splitting arrays only when split is on
      if (this%split_ratio /= 1) then
         allocate (this%h_avg(nx, ny), source=0.0_wp)
         allocate (this%hu_avg(nx, ny), source=0.0_wp)
         allocate (this%hv_avg(nx, ny), source=0.0_wp)
         allocate (this%h_bt0(nx, ny), source=0.0_wp)
         allocate (this%hu_bt0(nx, ny), source=0.0_wp)
         allocate (this%hv_bt0(nx, ny), source=0.0_wp)
         allocate (this%baro_tend_x(nx, ny), source=0.0_wp)
         allocate (this%baro_tend_y(nx, ny), source=0.0_wp)
      end if
   end subroutine nh_state_init

   subroutine nh_state_destroy(this)
      class(nonhydrostatic_state_t), intent(inout) :: this

      if (allocated(this%dsig)) deallocate (this%dsig)
      if (allocated(this%u_layer)) deallocate (this%u_layer)
      if (allocated(this%v_layer)) deallocate (this%v_layer)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%q_nh)) deallocate (this%q_nh)
      if (allocated(this%u_layer0)) deallocate (this%u_layer0)
      if (allocated(this%v_layer0)) deallocate (this%v_layer0)
      if (allocated(this%w0)) deallocate (this%w0)
      if (allocated(this%cg_r)) deallocate (this%cg_r)
      if (allocated(this%cg_p)) deallocate (this%cg_p)
      if (allocated(this%cg_Ap)) deallocate (this%cg_Ap)
      if (allocated(this%cg_z)) deallocate (this%cg_z)
      if (allocated(this%hu_pre)) deallocate (this%hu_pre)
      if (allocated(this%hv_pre)) deallocate (this%hv_pre)
      if (allocated(this%h_avg)) deallocate (this%h_avg)
      if (allocated(this%hu_avg)) deallocate (this%hu_avg)
      if (allocated(this%hv_avg)) deallocate (this%hv_avg)
      if (allocated(this%h_bt0)) deallocate (this%h_bt0)
      if (allocated(this%hu_bt0)) deallocate (this%hu_bt0)
      if (allocated(this%hv_bt0)) deallocate (this%hv_bt0)
      if (allocated(this%baro_tend_x)) deallocate (this%baro_tend_x)
      if (allocated(this%baro_tend_y)) deallocate (this%baro_tend_y)
   end subroutine nh_state_destroy

   pure function nh_state_bytes(this) result(nbytes)
      !! Counted allocatable footprint of the non-hydrostatic slot: layer
      !! velocities + NH pressure + RK2 saves + CG Poisson workspace + the
      !! gated mode-splitting buffers (0 when split_ratio == 1).
      class(nonhydrostatic_state_t), intent(in) :: this
      integer(int64) :: nbytes
      nbytes = arr_bytes(this%dsig) &
               + arr_bytes(this%u_layer) + arr_bytes(this%v_layer) &
               + arr_bytes(this%w) + arr_bytes(this%q_nh) &
               + arr_bytes(this%u_layer0) + arr_bytes(this%v_layer0) &
               + arr_bytes(this%w0) &
               + arr_bytes(this%cg_r) + arr_bytes(this%cg_p) &
               + arr_bytes(this%cg_Ap) + arr_bytes(this%cg_z) &
               + arr_bytes(this%hu_pre) + arr_bytes(this%hv_pre) &
               + arr_bytes(this%h_avg) + arr_bytes(this%hu_avg) &
               + arr_bytes(this%hv_avg) &
               + arr_bytes(this%h_bt0) + arr_bytes(this%hu_bt0) &
               + arr_bytes(this%hv_bt0) &
               + arr_bytes(this%baro_tend_x) + arr_bytes(this%baro_tend_y)
   end function nh_state_bytes

end module rki_nonhydrostatic_state
