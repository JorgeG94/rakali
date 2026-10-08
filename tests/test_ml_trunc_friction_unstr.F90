!! Unit tests for the unstructured per-layer velocity truncation and
!! per-layer bottom-friction twins (parity with structured
!! `ml_truncate_velocity` / `ml_apply_bottom_friction`).
module test_ml_trunc_friction_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, GRAVITY, THIN_LAYER_THRESHOLD
   use rki_mesh_unstr, only: mesh_unstr_t, &
                             mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy, &
                             mesh_unstr_enter_data, &
                             mesh_unstr_exit_data
   use rki_state_unstr, only: unstr_state_t, &
                              state_unstr_init, &
                              state_unstr_destroy, &
                              state_unstr_enter_data, &
                              state_unstr_exit_data
   use rki_ml_truncate_unstr, only: ml_truncate_velocity_unstr
   use rki_ml_friction_unstr, only: ml_apply_bottom_friction_unstr
   implicit none
   private

   public :: collect_ml_trunc_friction_unstr_tests

   integer, parameter :: NZ = 4
   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_ml_trunc_friction_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("trunc_clips_fast_layer", test_trunc_clips), &
                  new_unittest("trunc_leaves_slow_layer", test_trunc_slow_untouched), &
                  new_unittest("trunc_skips_thin_layer", test_trunc_thin_skipped), &
                  new_unittest("trunc_zero_cfl_no_op", test_trunc_zero_cfl_noop), &
                  new_unittest("friction_zero_manning_no_op", test_fric_zero_manning), &
                  new_unittest("friction_hand_computed_bottom", test_fric_hand_computed), &
                  new_unittest("friction_only_bottom_layer", test_fric_only_bottom), &
                  new_unittest("friction_monotonic_decay", test_fric_monotonic) &
                  ]
   end subroutine collect_ml_trunc_friction_unstr_tests

   ! ----------------------------------------------------------------------
   ! Helpers
   ! ----------------------------------------------------------------------

   subroutine build_two_triangle_square(mesh)
      !! Two right triangles filling the unit square; radius = inscribed
      !! circle of a right triangle with legs 1,1 (r = area/s).
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      tv(:, 2) = [1, 3, 4]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_two_triangle_square

   subroutine run_trunc(state, mesh, dt, cfl_trunc)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt, cfl_trunc

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_truncate_velocity_unstr(state, mesh, dt, cfl_trunc)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_trunc

   subroutine run_friction(state, mesh, dt)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh
      real(wp), intent(in) :: dt

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call ml_apply_bottom_friction_unstr(state, dt)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_friction

   ! ----------------------------------------------------------------------
   ! Truncation tests
   ! ----------------------------------------------------------------------

   subroutine test_trunc_clips(error)
      !! A layer velocity above vel_max = 0.9*cfl*radius/dt is clipped to
      !! exactly that speed (magnitude), sign preserved.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: dt, cfl, vel_max, h_k, u_fast
      integer :: cell, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      dt = 0.5_wp
      cfl = 0.8_wp
      h_k = 2.0_wp
      ! Fast velocity well above the cap for any triangle radius.
      u_fast = 1.0e6_wp

      do cell = 1, mesh%ncells_total
         do k = 1, NZ
            state%multilayer%h_layer(k, cell) = h_k
            state%multilayer%hu_layer(k, cell) = u_fast*h_k
            state%multilayer%hv_layer(k, cell) = -u_fast*h_k
         end do
      end do

      call run_trunc(state, mesh, dt, cfl)

      do cell = 1, mesh%ncells_owned
         vel_max = 0.9_wp*cfl*mesh%radius(cell)/dt
         do k = 1, NZ
            call check(error, &
                       abs(state%multilayer%hu_layer(k, cell) - vel_max*h_k) < TOL*max(vel_max*h_k, 1.0_wp), &
                       "hu clipped to +vel_max*h")
            if (allocated(error)) return
            call check(error, &
                       abs(state%multilayer%hv_layer(k, cell) + vel_max*h_k) < TOL*max(vel_max*h_k, 1.0_wp), &
                       "hv clipped to -vel_max*h")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_trunc_clips

   subroutine test_trunc_slow_untouched(error)
      !! A velocity below the cap is left bit-for-bit unchanged.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: hu0, hv0
      integer :: cell, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      hu0 = 1.0e-6_wp
      hv0 = -2.0e-6_wp
      do cell = 1, mesh%ncells_total
         do k = 1, NZ
            state%multilayer%h_layer(k, cell) = 1.0_wp
            state%multilayer%hu_layer(k, cell) = hu0
            state%multilayer%hv_layer(k, cell) = hv0
         end do
      end do

      call run_trunc(state, mesh, 0.5_wp, 0.8_wp)

      do cell = 1, mesh%ncells_owned
         do k = 1, NZ
            call check(error, abs(state%multilayer%hu_layer(k, cell) - hu0) < TOL, "hu unchanged")
            if (allocated(error)) return
            call check(error, abs(state%multilayer%hv_layer(k, cell) - hv0) < TOL, "hv unchanged")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_trunc_slow_untouched

   subroutine test_trunc_thin_skipped(error)
      !! A layer thinner than THIN_LAYER_THRESHOLD is not divided/clipped.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: hu0
      integer :: cell, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      hu0 = 1.0e6_wp
      do cell = 1, mesh%ncells_total
         do k = 1, NZ
            state%multilayer%h_layer(k, cell) = 0.5_wp*THIN_LAYER_THRESHOLD
            state%multilayer%hu_layer(k, cell) = hu0
            state%multilayer%hv_layer(k, cell) = 0.0_wp
         end do
      end do

      call run_trunc(state, mesh, 0.5_wp, 0.8_wp)

      do cell = 1, mesh%ncells_owned
         do k = 1, NZ
            call check(error, abs(state%multilayer%hu_layer(k, cell) - hu0) < TOL, "thin layer untouched")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_trunc_thin_skipped

   subroutine test_trunc_zero_cfl_noop(error)
      !! cfl_trunc <= 0 is a no-op (matches the gate).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: hu0
      integer :: cell, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      hu0 = 1.0e6_wp
      do cell = 1, mesh%ncells_total
         do k = 1, NZ
            state%multilayer%h_layer(k, cell) = 1.0_wp
            state%multilayer%hu_layer(k, cell) = hu0
            state%multilayer%hv_layer(k, cell) = 0.0_wp
         end do
      end do

      call run_trunc(state, mesh, 0.5_wp, 0.0_wp)

      do cell = 1, mesh%ncells_owned
         do k = 1, NZ
            call check(error, abs(state%multilayer%hu_layer(k, cell) - hu0) < TOL, "zero cfl no-op")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_trunc_zero_cfl_noop

   ! ----------------------------------------------------------------------
   ! Friction tests
   ! ----------------------------------------------------------------------

   subroutine test_fric_zero_manning(error)
      !! Default manning_n = 0 -> per-layer friction is a no-op.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: hu0
      integer :: cell, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      hu0 = 0.4_wp
      do cell = 1, mesh%ncells_total
         do k = 1, NZ
            state%multilayer%h_layer(k, cell) = 1.0_wp
            state%multilayer%hu_layer(k, cell) = hu0
            state%multilayer%hv_layer(k, cell) = 0.0_wp
         end do
      end do
      ! manning_n stays default 0.

      call run_friction(state, mesh, 0.5_wp)

      do cell = 1, mesh%ncells_owned
         do k = 1, NZ
            call check(error, abs(state%multilayer%hu_layer(k, cell) - hu0) < TOL, "zero manning no-op")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_fric_zero_manning

   subroutine test_fric_hand_computed(error)
      !! Bottom layer (k=1) drag = hu / (1 + dt*C_d*|u|/h_k),
      !! C_d = g*n^2 / H^{1/3}. Verify against a hand computation.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: dt, n, h_k, hu0, u_b, h_total, cd, denom, expected
      integer :: cell

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      dt = 10.0_wp
      n = 0.03_wp
      h_k = 1.5_wp             ! bottom-layer thickness
      hu0 = 0.6_wp             ! bottom-layer momentum
      h_total = NZ*h_k         ! uniform layers

      state%barotropic%manning_n = n
      do cell = 1, mesh%ncells_total
         state%multilayer%h_layer(:, cell) = h_k
         state%multilayer%hu_layer(:, cell) = hu0
         state%multilayer%hv_layer(:, cell) = 0.0_wp
      end do

      call run_friction(state, mesh, dt)

      u_b = hu0/h_k
      cd = GRAVITY*n*n/(h_total**(1.0_wp/3.0_wp))
      denom = 1.0_wp + dt*cd*abs(u_b)/h_k
      expected = hu0/denom

      do cell = 1, mesh%ncells_owned
         call check(error, abs(state%multilayer%hu_layer(1, cell) - expected) < 1.0e-12_wp, &
                    "bottom-layer drag matches hand computation")
         if (allocated(error)) return
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_fric_hand_computed

   subroutine test_fric_only_bottom(error)
      !! Friction acts on the bottom layer only; layers k=2..nz untouched.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: hu0
      integer :: cell, k

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      hu0 = 0.6_wp
      state%barotropic%manning_n = 0.03_wp
      do cell = 1, mesh%ncells_total
         state%multilayer%h_layer(:, cell) = 1.5_wp
         state%multilayer%hu_layer(:, cell) = hu0
         state%multilayer%hv_layer(:, cell) = 0.0_wp
      end do

      call run_friction(state, mesh, 10.0_wp)

      do cell = 1, mesh%ncells_owned
         call check(error, state%multilayer%hu_layer(1, cell) < hu0 - TOL, "bottom layer reduced")
         if (allocated(error)) return
         do k = 2, NZ
            call check(error, abs(state%multilayer%hu_layer(k, cell) - hu0) < TOL, "upper layers untouched")
            if (allocated(error)) return
         end do
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_fric_only_bottom

   subroutine test_fric_monotonic(error)
      !! Repeated application monotonically shrinks the bottom-layer speed
      !! and never flips its sign.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: prev
      integer :: it, cell

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh, nz_ml=NZ)

      state%barotropic%manning_n = 0.03_wp
      do cell = 1, mesh%ncells_total
         state%multilayer%h_layer(:, cell) = 1.5_wp
         state%multilayer%hu_layer(:, cell) = 0.6_wp
         state%multilayer%hv_layer(:, cell) = 0.0_wp
      end do

      prev = 0.6_wp
      do it = 1, 20
         call run_friction(state, mesh, 10.0_wp)
         call check(error, state%multilayer%hu_layer(1, 1) < prev + TOL, "monotonic non-increasing")
         if (allocated(error)) return
         call check(error, state%multilayer%hu_layer(1, 1) > 0.0_wp, "sign preserved")
         if (allocated(error)) return
         prev = state%multilayer%hu_layer(1, 1)
      end do

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_fric_monotonic

end module test_ml_trunc_friction_unstr
