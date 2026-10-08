!! Unit tests for the state type initialisation and destruction
module test_state
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_state, only: state_t
   use rki_solver, only: solver_enter_data, solver_exit_data, solver_step
   implicit none
   private

   public :: collect_state_tests

contains

   subroutine collect_state_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("init_allocates_arrays", test_init_allocates), &
                  new_unittest("init_sets_dimensions", test_init_dimensions), &
                  new_unittest("init_zeroes_arrays", test_init_zeroes), &
                  new_unittest("destroy_deallocates", test_destroy_deallocates), &
                  new_unittest("double_init_safe", test_double_init), &
                  new_unittest("device_cycle_reentry", test_device_cycle_reentry) &
                  ]
   end subroutine collect_state_tests

   subroutine test_init_allocates(error)
      !! All arrays should be allocated after init
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 8, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      call check(error, allocated(s%barotropic%h), "h should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%hu), "hu should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%hv), "hv should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%b), "b should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%flux_h), "flux_h should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%flux_hu), "flux_hu should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%flux_hv), "flux_hv should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%h0), "h0 should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%hu0), "hu0 should be allocated")
      if (allocated(error)) return
      call check(error, allocated(s%barotropic%hv0), "hv0 should be allocated")

      call s%destroy()
   end subroutine test_init_allocates

   subroutine test_init_dimensions(error)
      !! Grid metadata should be set correctly
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 8, 2, 5.0_wp, 3.0_wp, 0.4_wp)

      call check(error, s%grid%nx_phys == 10, "nx_phys should be 10")
      if (allocated(error)) return
      call check(error, s%grid%ny_phys == 8, "ny_phys should be 8")
      if (allocated(error)) return
      call check(error, s%grid%nx_total == 14, "nx_total should be 10 + 2*2 = 14")
      if (allocated(error)) return
      call check(error, s%grid%ny_total == 12, "ny_total should be 8 + 2*2 = 12")
      if (allocated(error)) return
      call check(error, s%grid%nghost == 2, "nghost should be 2")
      if (allocated(error)) return
      call check(error, abs(s%grid%dx - 5.0_wp) < 1.0e-15_wp, "dx should be 5.0")
      if (allocated(error)) return
      call check(error, abs(s%grid%dy - 3.0_wp) < 1.0e-15_wp, "dy should be 3.0")
      if (allocated(error)) return
      call check(error, abs(s%cfl - 0.4_wp) < 1.0e-15_wp, "cfl should be 0.4")
      if (allocated(error)) return

      ! Check array shapes match nx_total x ny_total
      call check(error, size(s%barotropic%h, 1) == 14, "h dim 1 should be nx_total")
      if (allocated(error)) return
      call check(error, size(s%barotropic%h, 2) == 12, "h dim 2 should be ny_total")

      call s%destroy()
   end subroutine test_init_dimensions

   subroutine test_init_zeroes(error)
      !! All arrays should be zero-initialised
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)

      call check(error, maxval(abs(s%barotropic%h)) < 1.0e-15_wp, "h should be zero")
      if (allocated(error)) return
      call check(error, maxval(abs(s%barotropic%hu)) < 1.0e-15_wp, "hu should be zero")
      if (allocated(error)) return
      call check(error, maxval(abs(s%barotropic%hv)) < 1.0e-15_wp, "hv should be zero")
      if (allocated(error)) return
      call check(error, maxval(abs(s%barotropic%b)) < 1.0e-15_wp, "b should be zero")

      call s%destroy()
   end subroutine test_init_zeroes

   subroutine test_destroy_deallocates(error)
      !! All arrays should be deallocated after destroy
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call s%destroy()

      call check(error,.not. allocated(s%barotropic%h), "h should be deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(s%barotropic%hu), "hu should be deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(s%barotropic%b), "b should be deallocated")
      if (allocated(error)) return
      call check(error,.not. allocated(s%barotropic%h0), "h0 should be deallocated")
   end subroutine test_destroy_deallocates

   subroutine test_double_init(error)
      !! Calling init twice should not leak (destroy + re-init)
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s

      call s%init(10, 10, 2, 1.0_wp, 1.0_wp, 0.45_wp)
      call s%destroy()
      call s%init(20, 15, 3, 2.0_wp, 2.0_wp, 0.3_wp)

      call check(error, s%grid%nx_phys == 20, "nx_phys should be 20 after re-init")
      if (allocated(error)) return
      call check(error, s%grid%nx_total == 26, "nx_total should be 20 + 2*3 = 26")
      if (allocated(error)) return
      call check(error, size(s%barotropic%h, 1) == 26, "h dim 1 should match new nx_total")

      call s%destroy()
   end subroutine test_double_init

   subroutine test_device_cycle_reentry(error)
      !! init -> enter_data -> step -> exit_data -> destroy, repeated.
      !! Exercises the full GPU device-mapping lifecycle three times in
      !! the same process. Catches device-memory leaks and stale OpenACC
      !! descriptors: a broken destroy/exit_data would fail on the second
      !! enter_data ("data already present") or corrupt the step output.
      type(error_type), allocatable, intent(out) :: error
      type(state_t) :: s
      integer :: cycle

      do cycle = 1, 3
         call s%init(8, 8, 2, 1.0_wp, 1.0_wp, 0.45_wp)
         s%barotropic%h = 2.0_wp
         s%barotropic%hu = 0.0_wp
         s%barotropic%hv = 0.0_wp
         s%barotropic%b = 0.0_wp

         call solver_enter_data(s)
         call solver_step(s)
         call solver_exit_data(s)

         call check(error, abs(s%barotropic%h(4, 4) - 2.0_wp) < 1.0e-10_wp, &
                    "lake-at-rest depth preserved on each cycle")
         if (allocated(error)) then
            call s%destroy()
            return
         end if

         call s%destroy()

         call check(error,.not. allocated(s%barotropic%h), &
                    "h must be deallocated after destroy on each cycle")
         if (allocated(error)) return
      end do
   end subroutine test_device_cycle_reentry

end module test_state
