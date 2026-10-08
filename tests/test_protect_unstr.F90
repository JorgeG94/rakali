!! Unit tests for the protect (negative-depth guard) kernel
module test_protect_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, DESING_EPS
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
   use rki_kernel_protect_unstr, only: protect_unstr
   implicit none
   private

   public :: collect_protect_tests

   real(wp), parameter :: TOL = 1.0e-12_wp

contains

   subroutine collect_protect_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("wet_cell_unchanged", test_wet_cell_unchanged), &
                  new_unittest("shallow_cell_desingularizes_momentum", test_shallow_zeros_momentum), &
                  new_unittest("dry_cell_zeros_momentum", test_dry_zeros_momentum), &
                  new_unittest("negative_depth_clipped_to_bed", test_negative_clipped), &
                  new_unittest("mass_error_starts_zero", test_mass_error_zero_initially), &
                  new_unittest("mass_error_accumulates", test_mass_error_accumulates), &
                  new_unittest("all_wet_no_mass_error", test_all_wet_no_mass_error), &
                  new_unittest("multiple_calls_accumulate", test_multiple_calls) &
                  ]
   end subroutine collect_protect_tests

   ! ----------------------------------------------------------------------
   ! Helpers
   ! ----------------------------------------------------------------------

   subroutine build_two_triangle_square(mesh)
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp) :: vx(4), vy(4)
      integer :: tv(3, 2)

      vx = [0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
      vy = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp]
      tv(:, 1) = [1, 2, 3]
      tv(:, 2) = [1, 3, 4]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)
   end subroutine build_two_triangle_square

   subroutine run_protect(state, mesh)
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(inout) :: mesh

      call mesh_unstr_enter_data(mesh)
      call state_unstr_enter_data(state)
      call protect_unstr(state, mesh)
      call state_unstr_exit_data(state)
      call mesh_unstr_exit_data(mesh)
   end subroutine run_protect

   ! ----------------------------------------------------------------------
   ! Tests
   ! ----------------------------------------------------------------------

   subroutine test_wet_cell_unchanged(error)
      !! Comfortably wet cells (depth >> min_h) must be left untouched.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Both cells well above min_h.
      state%barotropic%bed_cv = [0.0_wp, 0.0_wp]
      state%barotropic%stage_cv = [2.0_wp, 1.5_wp]
      state%barotropic%xmom_cv = [0.5_wp, -0.3_wp]
      state%barotropic%ymom_cv = [-0.2_wp, 0.4_wp]

      call run_protect(state, mesh)

      call check(error, abs(state%barotropic%stage_cv(1) - 2.0_wp) < TOL, "stage(1) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(2) - 1.5_wp) < TOL, "stage(2) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - 0.5_wp) < TOL, "xmom(1) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(2) - 0.4_wp) < TOL, "ymom(2) unchanged")
      if (allocated(error)) return
      call check(error, abs(state%mass_error) < TOL, "no mass error")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_wet_cell_unchanged

   subroutine test_shallow_zeros_momentum(error)
      !! Cell with DRY_TOLERANCE <= depth < eps: stage stays where it is,
      !! momentum gets the Kurganov & Petrova (2007) desingularized value
      !!     uh := h * sqrt(2) h uh / sqrt(h^4 + eps^4),
      !! eps = DESING_EPS — NOT the old hard zero, which froze thin water
      !! at wet/dry fronts.  Fixture momenta are small enough that the
      !! FROUDE_BACKSTOP does not bind, isolating the KP07 formula.
      !! mass_error is NOT incremented.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: h, eps, scale

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Cell 1: shallow (depth = 1e-4, below eps).  Cell 2: well wet.
      h = 1.0e-4_wp
      state%barotropic%bed_cv = [0.0_wp, 0.0_wp]
      state%barotropic%stage_cv = [h, 1.0_wp]
      state%barotropic%xmom_cv = [2.0e-4_wp, 0.5_wp]
      state%barotropic%ymom_cv = [-1.0e-4_wp, -0.3_wp]

      call run_protect(state, mesh)

      eps = DESING_EPS
      scale = sqrt(2.0_wp)*h**2/sqrt(h**4 + eps**4)

      ! Cell 1: stage unchanged (still positive), momentum desingularized.
      call check(error, abs(state%barotropic%stage_cv(1) - h) < TOL, &
                 "shallow cell stage unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1) - scale*2.0e-4_wp) < TOL, &
                 "shallow cell xmom must equal the KP07 desingularized value")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1) + scale*1.0e-4_wp) < TOL, &
                 "shallow cell ymom must equal the KP07 desingularized value")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1)) < 2.0e-4_wp, &
                 "desingularization must reduce shallow-cell momentum")
      if (allocated(error)) return

      ! Cell 2 (depth >= eps): untouched exactly.
      call check(error, abs(state%barotropic%stage_cv(2) - 1.0_wp) < TOL, "wet cell unchanged")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(2) - 0.5_wp) < TOL, "wet xmom unchanged")
      if (allocated(error)) return

      ! No negative depths -> no mass error.
      call check(error, abs(state%mass_error) < TOL, "shallow but >=0 -> no mass error")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_shallow_zeros_momentum

   subroutine test_dry_zeros_momentum(error)
      !! Cell with depth below DRY_TOLERANCE: momentum is zeroed exactly
      !! (a genuinely dry cell has no meaningful velocity), stage stays.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = [0.0_wp, 0.0_wp]
      state%barotropic%stage_cv = [1.0e-7_wp, 1.0_wp]
      state%barotropic%xmom_cv = [0.5_wp, 0.5_wp]
      state%barotropic%ymom_cv = [-0.3_wp, -0.3_wp]

      call run_protect(state, mesh)

      call check(error, state%barotropic%xmom_cv(1) == 0.0_wp, &
                 "dry cell xmom zeroed exactly")
      if (allocated(error)) return
      call check(error, state%barotropic%ymom_cv(1) == 0.0_wp, &
                 "dry cell ymom zeroed exactly")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%stage_cv(1) - 1.0e-7_wp) < TOL, &
                 "dry cell stage (mass) preserved")
      if (allocated(error)) return
      call check(error, abs(state%mass_error) < TOL, "dry but >=0 -> no mass error")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_dry_zeros_momentum

   subroutine test_negative_clipped(error)
      !! Cell with negative depth: stage gets clipped UP to bed,
      !! momentum gets zeroed (because depth < DRY_TOLERANCE), and the
      !! removed unphysical mass is added to state%mass_error.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: expected_error

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Cell 1 has stage = 0.4, bed = 0.5 -> depth = -0.1.
      ! Cell 2 wet, untouched.
      state%barotropic%bed_cv = [0.5_wp, 0.0_wp]
      state%barotropic%stage_cv = [0.4_wp, 1.0_wp]
      state%barotropic%xmom_cv = [0.5_wp, 0.5_wp]
      state%barotropic%ymom_cv = [-0.3_wp, -0.3_wp]

      call run_protect(state, mesh)

      call check(error, abs(state%barotropic%stage_cv(1) - 0.5_wp) < TOL, &
                 "negative-depth cell clipped to bed")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%xmom_cv(1)) < TOL, "xmom zeroed")
      if (allocated(error)) return
      call check(error, abs(state%barotropic%ymom_cv(1)) < TOL, "ymom zeroed")
      if (allocated(error)) return

      ! Cell 2 unchanged.
      call check(error, abs(state%barotropic%stage_cv(2) - 1.0_wp) < TOL, "wet stage unchanged")
      if (allocated(error)) return

      ! Mass error: (-h) * area = 0.1 * 0.5 (right triangle area = 1/2)
      expected_error = 0.1_wp*0.5_wp
      call check(error, abs(state%mass_error - expected_error) < 1.0e-13_wp, &
                 "mass_error = (-h)*area = 0.1*0.5 = 0.05")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_negative_clipped

   subroutine test_mass_error_zero_initially(error)
      !! Right after state_unstr_init, mass_error must be 0.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      call check(error, abs(state%mass_error) < TOL, "mass_error starts at 0")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_mass_error_zero_initially

   subroutine test_mass_error_accumulates(error)
      !! Both cells with negative depths -> mass_error is the sum.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: expected

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      ! Both cells negative.
      state%barotropic%bed_cv = [0.5_wp, 1.0_wp]
      state%barotropic%stage_cv = [0.3_wp, 0.7_wp]
      state%barotropic%xmom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]

      call run_protect(state, mesh)

      ! Cell 1: -h = 0.2, area = 0.5 -> 0.1
      ! Cell 2: -h = 0.3, area = 0.5 -> 0.15
      ! Total = 0.25
      expected = 0.2_wp*0.5_wp + 0.3_wp*0.5_wp
      call check(error, abs(state%mass_error - expected) < 1.0e-13_wp, &
                 "mass_error sums over both cells")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_mass_error_accumulates

   subroutine test_all_wet_no_mass_error(error)
      !! Lake at rest with deep water -> mass_error should remain 0
      !! after a protect pass (nothing to clip).
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = 0.0_wp
      state%barotropic%stage_cv = 5.0_wp
      state%barotropic%xmom_cv = 0.0_wp
      state%barotropic%ymom_cv = 0.0_wp

      call run_protect(state, mesh)

      call check(error, abs(state%mass_error) < TOL, "all-wet -> no mass error")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_all_wet_no_mass_error

   subroutine test_multiple_calls(error)
      !! Idempotence: protect is a projection onto the "h >= 0" manifold.
      !! After the first call cell 1 has stage = bed exactly (depth = 0),
      !! which is NOT strictly negative, so a second call must not clip
      !! again. Everything the first call produced — the `mass_error`
      !! tally AND the full cell state — must be unchanged by the
      !! second invocation, to machine precision.
      type(error_type), allocatable, intent(out) :: error
      type(mesh_unstr_t) :: mesh
      type(unstr_state_t) :: state
      real(wp) :: error_after_first
      real(wp) :: stage_after_first(2), xmom_after_first(2), ymom_after_first(2)

      call build_two_triangle_square(mesh)
      call state_unstr_init(state, mesh)

      state%barotropic%bed_cv = [0.5_wp, 0.0_wp]
      state%barotropic%stage_cv = [0.3_wp, 1.0_wp]
      state%barotropic%xmom_cv = [0.0_wp, 0.0_wp]
      state%barotropic%ymom_cv = [0.0_wp, 0.0_wp]

      call run_protect(state, mesh)
      error_after_first = state%mass_error
      stage_after_first = state%barotropic%stage_cv
      xmom_after_first = state%barotropic%xmom_cv
      ymom_after_first = state%barotropic%ymom_cv

      call check(error, error_after_first > 0.0_wp, "first call should record error")
      if (allocated(error)) return

      ! Second call on the already-clipped state must be a strict no-op.
      call run_protect(state, mesh)

      call check(error, abs(state%mass_error - error_after_first) < TOL, &
                 "second protect on already-clipped state must not grow mass_error")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%stage_cv - stage_after_first)) < TOL, &
                 "second protect must leave stage_cv bit-identical")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%xmom_cv - xmom_after_first)) < TOL, &
                 "second protect must leave xmom_cv bit-identical")
      if (allocated(error)) return
      call check(error, maxval(abs(state%barotropic%ymom_cv - ymom_after_first)) < TOL, &
                 "second protect must leave ymom_cv bit-identical")

      call state_unstr_destroy(state)
      call mesh_unstr_destroy(mesh)
   end subroutine test_multiple_calls

end module test_protect_unstr
