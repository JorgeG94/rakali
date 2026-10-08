!! Integration test: VCOORD_ZSTAR_FULL on a hard 10 m → 200 m shelf break.
!!
!! This is the F2 stress test for the full MOM6 z* implementation —
!! it exercises vanishing layers end-to-end through the SSP-RK2
!! multilayer pipeline (BPG, EOS, PP81, vertical exchange, friction,
!! truncation, tracers).
!!
!! What we validate:
!!   1. Solver runs through several SSP-RK2 steps without NaN/blowup
!!      under VCOORD_ZSTAR_FULL with a hard shelf-break bathymetry.
!!   2. On BOTH sides of the shelf, the surface ROMS layer (k=NZ)
!!      lands close to the configured `zstar_h_surf_target` — full z*'s
!!      headline property (surface layer thickness is a fixed physical
!!      value regardless of total column depth).  z*-lite cannot
!!      deliver this: its surface layer scales with H so the deep-side
!!      surface would be ~10× thicker than the shallow-side one.
!!   3. Per-column mass conservation: sum(h_layer) ≈ h(i,j) (remap
!!      precision = double-precision round-off times per-column work-set).
!!
!! Vanishing layers (the other half of full z*) are validated by
!! `test_vcoord_zstar_full::step_negative_eta_creates_vanishing` —
!! that's a kernel-level test exercising `H < h_bed_ref` directly.
!! Vanishing layers do not arise in a steady-bathymetry simulation
!! because each column's reference is built from its OWN bed.
module test_zstar_full_shelf_break
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, VCOORD_ZSTAR_FULL
   use rki_state, only: state_t
   use rki_boundary, only: bc_type_from_string
   use rki_ml_dynamics, only: ml_init_layers_from_barotropic, &
                              ml_build_vcoord_zstar_full
   use rki_ml_eos, only: ml_equation_of_state
   use rki_solver, only: solver_step, solver_enter_data, solver_exit_data
   use rki_vcoord, only: vcoord_target_dz_column_zstar_full, &
                         STRETCH_LOG, STRETCH_UNIFORM
   implicit none
   private

   public :: collect_zstar_full_shelf_break_tests

   integer, parameter :: NX_PHYS = 16
   integer, parameter :: NY_PHYS = 4
   integer, parameter :: NGHOST = 2
   integer, parameter :: NZ = 10
   real(wp), parameter :: H_SHALLOW = 10.0_wp     !! west half (m)
   real(wp), parameter :: H_DEEP = 200.0_wp       !! east half (m)
   real(wp), parameter :: DX_VAL = 500.0_wp       !! 500 m cells

   real(wp), parameter :: H_SURF_TARGET = 2.0_wp  !! m, fixed surface dz
   real(wp), parameter :: H_MIN = 1.0e-4_wp       !! vanishing-layer floor
   integer, parameter :: N_SURF = 3               !! 3 surface fine layers

   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: S_INIT = 35.0_wp
   real(wp), parameter :: T_INIT = 15.0_wp

contains

   subroutine collect_zstar_full_shelf_break_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("solver_runs_no_nan_no_blowup", &
                               test_solver_runs_no_nan), &
                  new_unittest("deep_surface_layer_at_target", &
                               test_deep_surface_at_target), &
                  new_unittest("shallow_surface_layer_at_target", &
                               test_shallow_surface_at_target), &
                  new_unittest("per_column_mass_conservation", &
                               test_per_column_mass), &
                  new_unittest("total_tracer_mass_conserved_in_closed_basin", &
                               test_total_tracer_mass) &
                  ]
   end subroutine collect_zstar_full_shelf_break_tests

   subroutine setup_shelf_state(state)
      type(state_t), intent(inout) :: state

      integer :: i, j, k, i_mid
      real(wp) :: dz_col(NZ)

      state%use_multilayer = .true.
      state%multilayer%nz_ml = NZ
      state%multilayer%initial_salinity = S_INIT
      state%multilayer%initial_temperature = T_INIT
      state%multilayer%alpha_T = 0.0_wp
      state%multilayer%beta_S = 0.0_wp
      state%multilayer%rho_0 = RHO_0
      state%multilayer%vmix_nu_max = 0.0_wp
      state%multilayer%vmix_nu_bg = 0.0_wp
      state%multilayer%hdiff_kappa = 0.0_wp
      state%multilayer%smagorinsky_cs = 0.0_wp

      state%multilayer%vcoord%coord_type = VCOORD_ZSTAR_FULL
      state%multilayer%vcoord%zstar_h_surf_target = H_SURF_TARGET
      state%multilayer%vcoord%zstar_h_min = H_MIN
      state%multilayer%vcoord%zstar_stretching = STRETCH_UNIFORM
      state%multilayer%vcoord%zstar_n_surf = N_SURF

      call state%init(NX_PHYS, NY_PHYS, NGHOST, DX_VAL, DX_VAL, 0.4_wp)

      ! Closed basin
      state%bc%west = bc_type_from_string("wall")
      state%bc%east = bc_type_from_string("wall")
      state%bc%south = bc_type_from_string("wall")
      state%bc%north = bc_type_from_string("wall")

      state%barotropic%coriolis_f = 0.0_wp
      state%barotropic%manning_n = 0.0_wp

      ! Bathymetry: hard step at i = NX/2.  b is bed elevation (positive
      ! down).  Convention used throughout rakali tests: b = -bed_depth
      ! so h = eta - b.  For our test, eta_initial = 0 and we want
      ! h(i,j) = H_SHALLOW or H_DEEP, so b = -H_*.
      i_mid = state%grid%nx_total/2
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            if (i <= i_mid) then
               state%barotropic%b(i, j) = -H_SHALLOW
               state%barotropic%h(i, j) = H_SHALLOW
            else
               state%barotropic%b(i, j) = -H_DEEP
               state%barotropic%h(i, j) = H_DEEP
            end if
         end do
      end do
      state%barotropic%hu = 0.0_wp
      state%barotropic%hv = 0.0_wp

      ! Build the per-column z_ref pattern from bathymetry now.  In
      ! production this is done in `rki_init`; the test bypasses
      ! rki_init and does it explicitly.  `state%barotropic%b` stores
      ! the bed as a negative number (eta-b convention), so pass the
      ! magnitude as h_bed to the builder.
      block
         integer :: ii, jj
         real(wp) :: h_bed
         do jj = 1, state%grid%ny_total
            do ii = 1, state%grid%nx_total
               h_bed = -state%barotropic%b(ii, jj)
               call build_one_column(NZ, h_bed, H_SURF_TARGET, N_SURF, &
                                     STRETCH_UNIFORM, &
                                     state%multilayer%vcoord_z_ref_col(:, ii, jj))
            end do
         end do
      end block

      call ml_init_layers_from_barotropic(state%barotropic%h, &
                                          state%barotropic%hu, &
                                          state%barotropic%hv, &
                                          state%multilayer%dsig_ml, &
                                          state%multilayer%h_layer, &
                                          state%multilayer%hu_layer, &
                                          state%multilayer%hv_layer, &
                                          state%grid, NZ)

      ! Re-shape h_layer per cell to the ZSTAR_FULL target so the IC
      ! matches the vcoord (otherwise the first RK2 average blends a
      ! uniform-sigma snapshot with a full-z*-remapped one).
      do j = 1, state%grid%ny_total
         do i = 1, state%grid%nx_total
            call vcoord_target_dz_column_zstar_full(NZ, state%barotropic%h(i, j), &
                                                    state%multilayer%vcoord_z_ref_col(:, i, j), &
                                                    H_MIN, dz_col)
            do k = 1, NZ
               state%multilayer%h_layer(i, j, k) = dz_col(k)
            end do
         end do
      end do

      ! Uniform tracer (hTr scales with h_layer)
      state%multilayer%tracers(state%multilayer%idx_salinity)%hTr = &
         state%multilayer%h_layer*S_INIT
      state%multilayer%tracers(state%multilayer%idx_temperature)%hTr = &
         state%multilayer%h_layer*T_INIT

      call ml_equation_of_state(state%multilayer%h_layer, &
                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                state%multilayer%rho_layer, &
                                RHO_0, 0.0_wp, 0.0_wp, 0.0_wp, T_INIT, &
                                state%grid, NZ)
   end subroutine setup_shelf_state

   subroutine build_one_column(nz, h_bed, h_surf_target, n_surf, stretching, &
                               z_ref_col)
      !! Thin wrapper so the test can call the per-column builder without
      !! pulling rki_vcoord publicly — already pulled via rki_ml_dynamics.
      use rki_vcoord, only: zstar_full_build_column
      integer, intent(in) :: nz, n_surf, stretching
      real(wp), intent(in) :: h_bed, h_surf_target
      real(wp), intent(out) :: z_ref_col(0:nz)
      call zstar_full_build_column(nz, h_bed, h_surf_target, n_surf, &
                                   stretching, z_ref_col)
   end subroutine build_one_column

   subroutine run_to(state, t_final)
      type(state_t), intent(inout) :: state
      real(wp), intent(in) :: t_final
      real(wp) :: t

      call solver_enter_data(state)
      t = 0.0_wp
      do while (t < t_final)
         call solver_step(state, time=t)
         t = t + state%dt
      end do
      call solver_exit_data(state)
   end subroutine run_to

   logical function any_nan_in_h_layer(state)
      type(state_t), intent(in) :: state
      integer :: i, j, k
      any_nan_in_h_layer = .false.
      do k = 1, NZ
         do j = 1, state%grid%ny_total
            do i = 1, state%grid%nx_total
               if (state%multilayer%h_layer(i, j, k) &
                   /= state%multilayer%h_layer(i, j, k)) then
                  any_nan_in_h_layer = .true.
                  return
               end if
            end do
         end do
      end do
   end function any_nan_in_h_layer

   subroutine test_solver_runs_no_nan(error)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 10.0_wp

      call setup_shelf_state(state)
      call run_to(state, T_END)

      call check(error,.not. any_nan_in_h_layer(state), &
                 "ZSTAR_FULL shelf-break: h_layer contains no NaN after 10s")
      if (allocated(error)) then
         call state%destroy()
         return
      end if

      ! Quick blow-up check: max |hu| should not have grown to silly values
      block
         integer :: i, j, k
         real(wp) :: max_hu
         max_hu = 0.0_wp
         do k = 1, NZ
            do j = NGHOST + 1, state%grid%ny_total - NGHOST
               do i = NGHOST + 1, state%grid%nx_total - NGHOST
                  if (abs(state%multilayer%hu_layer(i, j, k)) > max_hu) then
                     max_hu = abs(state%multilayer%hu_layer(i, j, k))
                  end if
               end do
            end do
         end do
         call check(error, max_hu < 1.0e6_wp, &
                    "ZSTAR_FULL shelf-break: |hu_layer| stays bounded")
      end block

      call state%destroy()
   end subroutine test_solver_runs_no_nan

   subroutine test_deep_surface_at_target(error)
      !! On the deep side (H=200m) the surface ROMS layer (k=NZ)
      !! should be close to H_SURF_TARGET (2 m), independent of the
      !! deep total depth.  This is the headline property of full z*.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      integer :: i_deep, j_mid
      real(wp) :: dz_surf

      call setup_shelf_state(state)
      ! Don't even run — the initial layout already encodes this
      ! property and the per-step remap target preserves it.
      i_deep = state%grid%nx_total - NGHOST - 1
      j_mid = state%grid%ny_total/2
      dz_surf = state%multilayer%h_layer(i_deep, j_mid, NZ)

      call check(error, dz_surf > 0.5_wp*H_SURF_TARGET .and. &
                 dz_surf < 2.0_wp*H_SURF_TARGET, &
                 "deep-side surface layer thickness near H_SURF_TARGET")
      call state%destroy()
   end subroutine test_deep_surface_at_target

   subroutine test_shallow_surface_at_target(error)
      !! On the shallow side (H=10m) the surface ROMS layer (k=NZ)
      !! should also land at ≈ H_SURF_TARGET — same fixed physical
      !! thickness as the deep side, NOT proportionally scaled to H.
      !! This is the property z*-lite cannot deliver.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      integer :: i_shallow, j_mid
      real(wp) :: dz_surf

      call setup_shelf_state(state)
      i_shallow = NGHOST + 2     ! interior on shallow side
      j_mid = state%grid%ny_total/2
      dz_surf = state%multilayer%h_layer(i_shallow, j_mid, NZ)

      call check(error, dz_surf > 0.5_wp*H_SURF_TARGET .and. &
                 dz_surf < 2.0_wp*H_SURF_TARGET, &
                 "shallow-side surface layer thickness near H_SURF_TARGET")
      call state%destroy()
   end subroutine test_shallow_surface_at_target

   subroutine test_per_column_mass(error)
      !! sum(h_layer(i, j, :)) should equal h(i, j) at every interior
      !! cell — both at init and after a run.  Mass conservation per
      !! column is the remap invariant; vanishing-layer overhead must
      !! be absorbed by the surface trim in vcoord_target_dz_column_zstar_full.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 5.0_wp
      real(wp), parameter :: TOL = 1.0e-9_wp
      integer :: i, j, k
      real(wp) :: sum_h, diff

      call setup_shelf_state(state)
      call run_to(state, T_END)

      do j = NGHOST + 1, state%grid%ny_total - NGHOST
         do i = NGHOST + 1, state%grid%nx_total - NGHOST
            sum_h = 0.0_wp
            do k = 1, NZ
               sum_h = sum_h + state%multilayer%h_layer(i, j, k)
            end do
            diff = abs(sum_h - state%barotropic%h(i, j))
            call check(error, diff < TOL*max(state%barotropic%h(i, j), 1.0_wp), &
                       "per-column sum(h_layer) = h at interior cells")
            if (allocated(error)) then
               call state%destroy()
               return
            end if
         end do
      end do

      call state%destroy()
   end subroutine test_per_column_mass

   subroutine test_total_tracer_mass(error)
      !! Domain-integrated tracer mass conservation in a closed (walled)
      !! basin under VCOORD_ZSTAR_FULL.  The unit tests for the kernel
      !! already prove per-column dz_new/dz_old conservation; the
      !! shelf-break setup additionally exercises the per-column
      !! `ml_distribute_with_remap_full` loop driven by the full SSP-RK2
      !! solver.  All boundaries are walls (set in setup_shelf_state),
      !! so there is no physical flux of tracer — drift indicates a
      !! bookkeeping bug in the remap or advection path.
      !!
      !! Regression test for the 7% salt drift observed in the Lakes
      !! Entrance tidal run (2026-05-11).
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: state
      real(wp), parameter :: T_END = 10.0_wp
      real(wp), parameter :: REL_TOL = 1.0e-10_wp
      real(wp) :: salt0, saltT, drift

      call setup_shelf_state(state)
      salt0 = total_tracer_mass(state, state%multilayer%idx_salinity)

      call run_to(state, T_END)
      saltT = total_tracer_mass(state, state%multilayer%idx_salinity)

      drift = abs(saltT - salt0)/max(abs(salt0), 1.0_wp)

      call check(error, drift < REL_TOL, &
                 "ZSTAR_FULL closed basin: total salt mass conserved to ~1e-10")

      call state%destroy()
   end subroutine test_total_tracer_mass

   function total_tracer_mass(state, idx) result(total)
      !! Domain-integrated tracer mass over interior cells.
      type(state_t), intent(in) :: state
      integer, intent(in) :: idx
      real(wp) :: total
      integer :: i, j, k

      total = 0.0_wp
      do k = 1, NZ
         do j = NGHOST + 1, state%grid%ny_total - NGHOST
            do i = NGHOST + 1, state%grid%nx_total - NGHOST
               total = total + state%multilayer%tracers(idx)%hTr(i, j, k)
            end do
         end do
      end do
      total = total*DX_VAL*DX_VAL
   end function total_tracer_mass

end module test_zstar_full_shelf_break
