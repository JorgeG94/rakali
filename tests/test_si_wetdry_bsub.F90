!! Subgrid b_sub fill: the smooth-drainage demonstration.
module test_si_wetdry_bsub
   !! Q8 Part-A gate: the analytic subgrid b_sub table (`wetdry_bsub_fill=
   !! "linear"`) gives V(eta) genuine sub-cell breakpoints so a coarse-cell
   !! tidal flat drains/floods SMOOTHLY.  Quantified against n_sub=1: the
   !! wet-fraction jump per eta step is 1 for n_sub=1 (binary) and scales as
   !! 1/n_sub for the "linear" fill, evaluated edge-wise on a 3-cell
   !! unstructured strip so the mesh-neighbour relief sweep is exercised.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp
   use rki_si_common, only: SI_BSUB_FILL_FLAT, SI_BSUB_FILL_LINEAR, &
                            parse_si_bsub_fill, si_bsub_ml_conflict
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_init_from_triangles, &
                             mesh_unstr_destroy
   use rki_si_unstr_wetdry, only: si_unstr_wetdry_fill_bsub, &
                                  si_unstr_wetdry_volume, si_unstr_wetdry_wetfrac
   use pic_logger, only: global_logger
   implicit none
   private

   public :: collect_si_wetdry_bsub_tests

contains

   subroutine collect_si_wetdry_bsub_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("bsub_parse", test_parse), &
                  new_unittest("bsub_ml_conflict_guard", test_ml_conflict_guard), &
                  new_unittest("bsub_unstructured_smoothing", test_unstructured) &
                  ]
   end subroutine collect_si_wetdry_bsub_tests

   subroutine test_ml_conflict_guard(error)
      !! The MUST-FIX conservation guard predicate (both backends' enable
      !! paths error-stop on a .true. return; the error-stop itself cannot be
      !! caught in-process, so this asserts the pure predicate directly):
      !!   linear + n_sub>=2 + multilayer(nz>1)  => CONFLICT (rejected)
      !!   flat + ML, linear + 2D, linear+ML+n_sub=1, linear+ML+nz=1 => OK.
      type(error_type), allocatable, intent(out) :: error

      ! The rejected combination (linear subgrid fill under multilayer).
      call check(error, si_bsub_ml_conflict(SI_BSUB_FILL_LINEAR, 2, .true., 4), &
                 "guard must FIRE for linear + n_sub=2 + multilayer(nz=4)")
      if (allocated(error)) return
      call check(error, si_bsub_ml_conflict(SI_BSUB_FILL_LINEAR, 8, .true., 20), &
                 "guard must FIRE for linear + n_sub=8 + multilayer(nz=20)")
      if (allocated(error)) return

      ! flat + ML: conservative (b_sub == b), allowed.
      call check(error,.not. si_bsub_ml_conflict(SI_BSUB_FILL_FLAT, 4, .true., 4), &
                 "flat + multilayer must be ALLOWED")
      if (allocated(error)) return
      ! linear + 2D barotropic: the barotropic path uses v_cell throughout.
      call check(error,.not. si_bsub_ml_conflict(SI_BSUB_FILL_LINEAR, 4, .false., 1), &
                 "linear + 2D barotropic must be ALLOWED")
      if (allocated(error)) return
      ! Degenerate subgrid (n_sub=1) == flat effectively: allowed.
      call check(error,.not. si_bsub_ml_conflict(SI_BSUB_FILL_LINEAR, 1, .true., 4), &
                 "linear + multilayer with n_sub=1 must be ALLOWED (== flat)")
      if (allocated(error)) return
      ! nz=1 multilayer degenerates to 2D: allowed.
      call check(error,.not. si_bsub_ml_conflict(SI_BSUB_FILL_LINEAR, 2, .true., 1), &
                 "linear + multilayer with nz=1 must be ALLOWED (== 2D)")
   end subroutine test_ml_conflict_guard

   subroutine test_parse(error)
      type(error_type), allocatable, intent(out) :: error
      call check(error, parse_si_bsub_fill("flat") == SI_BSUB_FILL_FLAT, "flat")
      if (allocated(error)) return
      call check(error, parse_si_bsub_fill("linear") == SI_BSUB_FILL_LINEAR, "linear")
      if (allocated(error)) return
      call check(error, parse_si_bsub_fill("") == SI_BSUB_FILL_FLAT, "empty->flat")
   end subroutine test_parse

   subroutine test_unstructured(error)
      !! A 3-quad strip (6 triangles) with a sloping bed; drive the mesh
      !! fill kernel and confirm the middle cell's wet fraction ramps in
      !! 1/n_sub steps (the same smoothing, edge-wise).
      type(error_type), allocatable, intent(out) :: error

      type(mesh_unstr_t) :: mesh
      real(wp), allocatable :: vx(:), vy(:)
      integer, allocatable :: tv(:, :)
      real(wp), allocatable :: bed(:)
      real(wp) :: mj1, mj4, mid_x
      integer :: kmid, k

      ! 3 quads in x, 1 in y => 4x2 vertex grid, 6 triangles.
      vx = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp]*1000.0_wp
      vy = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]*1000.0_wp
      allocate (tv(3, 6))
      tv(:, 1) = [1, 2, 6]
      tv(:, 2) = [1, 6, 5]
      tv(:, 3) = [2, 3, 7]
      tv(:, 4) = [2, 7, 6]
      tv(:, 5) = [3, 4, 8]
      tv(:, 6) = [3, 8, 7]
      call mesh_unstr_init_from_triangles(mesh, vx, vy, tv)

      allocate (bed(mesh%ncells_total))
      do k = 1, mesh%ncells_total
         bed(k) = -3.0_wp + 2.0_wp*(mesh%circumcentre_x(k)/3000.0_wp)
      end do

      ! Middle-x owned cell.
      kmid = 1
      mid_x = huge(1.0_wp)
      do k = 1, mesh%ncells_owned
         if (abs(mesh%circumcentre_x(k) - 1500.0_wp) < mid_x) then
            mid_x = abs(mesh%circumcentre_x(k) - 1500.0_wp)
            kmid = k
         end if
      end do

      mj1 = unstr_wetfrac_max_jump(mesh, bed, 1, SI_BSUB_FILL_FLAT, kmid)
      mj4 = unstr_wetfrac_max_jump(mesh, bed, 4, SI_BSUB_FILL_LINEAR, kmid)

      call global_logger%info("[bsub unstructured] mid-cell max wet-frac jump "// &
                              "n_sub=1/4 = "//rstr(mj1)//"/"//rstr(mj4))

      call check(error, abs(mj1 - 1.0_wp) < 1.0e-12_wp, &
                 "unstr n_sub=1 should jump 0->1, got "//rstr(mj1))
      if (allocated(error)) then
         call mesh_unstr_destroy(mesh)
         return
      end if
      call check(error, mj4 <= 0.25_wp + 1.0e-12_wp .and. mj4 < mj1, &
                 "unstr n_sub=4 linear jump "//rstr(mj4)//" not smoother (<=1/4)")
      call mesh_unstr_destroy(mesh)
   end subroutine test_unstructured

   function unstr_wetfrac_max_jump(mesh, bed, n_sub, mode, kcell) result(mj)
      type(mesh_unstr_t), intent(in) :: mesh
      real(wp), intent(in) :: bed(:)
      integer, intent(in) :: n_sub, mode, kcell
      real(wp) :: mj
      integer :: nct, nco, n, nsteps
      real(wp), allocatable :: b_sub(:, :), b_min(:), eta(:), a_cell(:)
      real(wp) :: e, elo, ehi, prev_a, cur_a

      nct = mesh%ncells_total
      nco = mesh%ncells_owned
      allocate (b_sub(n_sub, nct), b_min(nct), eta(nct), a_cell(nct))
      call si_unstr_wetdry_fill_bsub(mesh, n_sub, nct, nco, bed, mode, b_sub, b_min)
      elo = b_min(kcell) - 0.5_wp
      ehi = maxval(b_sub(:, kcell)) + 0.5_wp
      nsteps = 200
      mj = 0.0_wp
      prev_a = -1.0_wp
      do n = 0, nsteps
         e = elo + (ehi - elo)*real(n, wp)/real(nsteps, wp)
         eta = e
         call si_unstr_wetdry_wetfrac(n_sub, nct, nco, eta, b_sub, a_cell)
         cur_a = a_cell(kcell)
         if (prev_a >= 0.0_wp) mj = max(mj, abs(cur_a - prev_a))
         prev_a = cur_a
      end do
   end function unstr_wetfrac_max_jump

   pure function rstr(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write (buf, "(ES12.5)") x
      s = trim(adjustl(buf))
   end function rstr

end module test_si_wetdry_bsub
