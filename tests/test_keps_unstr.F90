!! Unit tests for the unstructured k-epsilon column kernel.
module test_keps_unstr
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, KEPS_STAB_CONSTANT, KEPS_STAB_GALPERIN, KEPS_STAB_CANUTO
   use rki_ml_keps_unstr, only: ml_keps_update_unstr_impl, &
                                ml_keps_apply_momentum_unstr_impl, &
                                ml_keps_apply_diffusion_tracer_unstr_impl
   implicit none
   private

   public :: collect_keps_unstr_tests

   ! Tiny synthetic column set: N identical wet cells, NZ layers.  The impl is
   ! cell-centric (intra-column), so a flat (nz, ncell) array suffices — no mesh.
   integer, parameter :: N = 8, NZ = 4
   real(wp), parameter :: CE2 = 1.92_wp, CMU = 0.09_wp

contains

   subroutine collect_keps_unstr_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("keps_unstr_homogeneous_decay", test_homogeneous_decay), &
                  new_unittest("keps_unstr_positivity", test_positivity), &
                  new_unittest("keps_unstr_wall_function", test_wall_function), &
                  new_unittest("keps_unstr_galperin_suppression", test_galperin_suppression), &
                  new_unittest("keps_unstr_canuto_suppression", test_canuto_suppression), &
                  new_unittest("keps_unstr_momentum_conservation", test_momentum_conservation), &
                  new_unittest("keps_unstr_momentum_shear_smoothing", test_momentum_shear_smoothing), &
                  new_unittest("keps_unstr_tracer_mixing", test_tracer_mixing) &
                  ]
   end subroutine collect_keps_unstr_tests

   subroutine test_homogeneous_decay(error)
      !! No shear, no stratification: k-eps reduces to decaying turbulence with
      !!   k(t) = k0 (1 + (ce2-1) eps0/k0 t)^(-1/(ce2-1)).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), rho(NZ, N)
      real(wp) :: tke(NZ + 1, N), eps(NZ + 1, N), nut(NZ + 1, N), kt(NZ + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      real(wp), parameter :: K0 = 1.0_wp, E0 = 0.1_wp, DT = 0.01_wp, T_END = 20.0_wp
      real(wp) :: k_exact
      integer :: step, nstep

      h = 1.0_wp
      hu = 0.0_wp                 ! no shear  -> S^2 = 0
      hv = 0.0_wp
      rho = 1025.0_wp             ! no stratification -> N^2 = 0
      tke = K0
      eps = E0
      nut = CMU*K0*K0/E0
      kt = nut
      ustar_s = 0.0_wp
      ustar_b = 0.0_wp
      height_cv = real(NZ, wp)*1.0_wp
      nstep = nint(T_END/DT)

      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      do step = 1, nstep
         call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                        ustar_s, ustar_b, height_cv, 1025.0_wp, DT, &
                                        .false., KEPS_STAB_CONSTANT)
      end do
      !$omp end target data

      k_exact = K0*(1.0_wp + (CE2 - 1.0_wp)*E0/K0*T_END)**(-1.0_wp/(CE2 - 1.0_wp))
      call check(error, abs(tke(3, 4) - k_exact)/k_exact < 0.02_wp, &
                 "k-eps homogeneous decay matches analytic power law")
   end subroutine test_homogeneous_decay

   subroutine test_positivity(error)
      !! Strong shear + large dt must keep k, eps, nu_t positive and bounded.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), rho(NZ, N)
      real(wp) :: tke(NZ + 1, N), eps(NZ + 1, N), nut(NZ + 1, N), kt(NZ + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      integer :: k, step

      h = 2.0_wp
      hv = 0.0_wp
      rho = 1025.0_wp
      do k = 1, NZ                ! strong vertical shear: u grows with layer
         hu(k, :) = h(k, :)*real(k, wp)*0.5_wp
      end do
      tke = 1.0e-6_wp
      eps = 1.0e-8_wp
      nut = 1.0e-4_wp
      kt = nut
      ustar_s = 0.0_wp
      ustar_b = 0.0_wp
      height_cv = real(NZ, wp)*2.0_wp

      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      do step = 1, 200
         call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                        ustar_s, ustar_b, height_cv, 1025.0_wp, 100.0_wp, &
                                        .false., KEPS_STAB_CONSTANT)
      end do
      !$omp end target data

      call check(error, all(tke > 0.0_wp), "k stays positive under strong shear + large dt")
      if (allocated(error)) return
      call check(error, all(eps > 0.0_wp) .and. all(nut <= 10.0_wp), &
                 "eps positive and nu_t bounded")
      if (allocated(error)) return
      call check(error, tke(3, 4) > 1.0e-4_wp, "shear production grows TKE")
   end subroutine test_positivity

   subroutine test_wall_function(error)
      !! Prescribed friction velocities pin the near-wall interfaces to the
      !! Dirichlet log-law value k = u_*^2 / sqrt(c_mu) (bed at interface 2,
      !! surface at interface NZ).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), rho(NZ, N)
      real(wp) :: tke(NZ + 1, N), eps(NZ + 1, N), nut(NZ + 1, N), kt(NZ + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      real(wp), parameter :: USB = 0.02_wp, USS = 0.01_wp
      real(wp) :: k_wb_exp, k_ws_exp
      integer :: step

      h = 2.0_wp
      hu = 0.0_wp
      hv = 0.0_wp
      rho = 1025.0_wp
      tke = 1.0e-6_wp
      eps = 1.0e-8_wp
      nut = 1.0e-4_wp
      kt = nut
      ustar_b = USB
      ustar_s = USS
      height_cv = real(NZ, wp)*2.0_wp

      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      do step = 1, 5
         call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                        ustar_s, ustar_b, height_cv, 1025.0_wp, 10.0_wp, &
                                        .true., KEPS_STAB_CONSTANT)
      end do
      !$omp end target data

      k_wb_exp = USB*USB/sqrt(CMU)
      k_ws_exp = USS*USS/sqrt(CMU)
      call check(error, abs(tke(2, 4) - k_wb_exp)/k_wb_exp < 1.0e-6_wp, &
                 "bed wall pins k = u_*b^2/sqrt(c_mu) at interface 2")
      if (allocated(error)) return
      call check(error, abs(tke(NZ, 4) - k_ws_exp)/k_ws_exp < 1.0e-6_wp, &
                 "surface wall pins k = u_*s^2/sqrt(c_mu) at interface NZ")
   end subroutine test_wall_function

   subroutine test_galperin_suppression(error)
      !! Under stable stratification Galperin suppresses nu_t below the
      !! constant-c_mu value (the scheme only enters the nu_t/kappa_t recovery).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), rho(NZ, N)
      real(wp) :: tke(NZ + 1, N), eps(NZ + 1, N), nut(NZ + 1, N), kt(NZ + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      real(wp) :: nut_c, nut_g
      integer :: k

      h = 2.0_wp
      hv = 0.0_wp
      do k = 1, NZ               ! shear + stable stratification (light at surface)
         hu(k, :) = h(k, :)*real(k, wp)*0.2_wp
         rho(k, :) = 1025.0_wp - 1.5_wp*real(k - 1, wp)
      end do
      ustar_s = 0.0_wp
      ustar_b = 0.0_wp
      height_cv = real(NZ, wp)*2.0_wp

      tke = 0.05_wp
      eps = 5.0e-3_wp
      nut = 1.0e-3_wp
      kt = 1.0e-3_wp
      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                     ustar_s, ustar_b, height_cv, 1025.0_wp, 10.0_wp, &
                                     .false., KEPS_STAB_CONSTANT)
      !$omp end target data
      nut_c = nut(3, 4)

      tke = 0.05_wp
      eps = 5.0e-3_wp
      nut = 1.0e-3_wp
      kt = 1.0e-3_wp
      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                     ustar_s, ustar_b, height_cv, 1025.0_wp, 10.0_wp, &
                                     .false., KEPS_STAB_GALPERIN)
      !$omp end target data
      nut_g = nut(3, 4)

      call check(error, nut_g < nut_c, &
                 "Galperin suppresses nu_t below constant-c_mu under stable stratification")
      if (allocated(error)) return
      call check(error, nut_g < 0.9_wp*nut_c, &
                 "Galperin nu_t suppression is meaningful (< 0.9x constant)")
   end subroutine test_galperin_suppression

   subroutine test_canuto_suppression(error)
      !! Canuto-A under stable stratification suppresses the scalar diffusivity
      !! more than the eddy viscosity, so kt < nu_t (Pr_t > 1).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), rho(NZ, N)
      real(wp) :: tke(NZ + 1, N), eps(NZ + 1, N), nut(NZ + 1, N), kt(NZ + 1, N)
      real(wp) :: ustar_s(N), ustar_b(N), height_cv(N)
      real(wp) :: nut_c, kt_c
      integer :: k

      h = 2.0_wp
      hv = 0.0_wp
      do k = 1, NZ               ! shear + strong stable stratification
         hu(k, :) = h(k, :)*real(k, wp)*0.2_wp
         rho(k, :) = 1025.0_wp - 3.0_wp*real(k - 1, wp)
      end do
      tke = 0.05_wp
      eps = 5.0e-3_wp
      nut = 1.0e-3_wp
      kt = 1.0e-3_wp
      ustar_s = 0.0_wp
      ustar_b = 0.0_wp
      height_cv = real(NZ, wp)*2.0_wp

      !$omp target data map(to: h, hu, hv, rho, ustar_s, ustar_b, height_cv) map(tofrom: tke, eps, nut, kt)
      call ml_keps_update_unstr_impl(NZ, N, h, hu, hv, rho, tke, eps, nut, kt, &
                                     ustar_s, ustar_b, height_cv, 1025.0_wp, 10.0_wp, &
                                     .false., KEPS_STAB_CANUTO)
      !$omp end target data
      nut_c = nut(3, 4)
      kt_c = kt(3, 4)

      call check(error, nut_c > 1.0e-4_wp .and. kt_c > 1.0e-4_wp, &
                 "Canuto-A nu_t, kt are positive physical values (rational form well-posed)")
      if (allocated(error)) return
      call check(error, kt_c < nut_c, &
                 "Canuto-A scalar diffusivity < eddy viscosity (Pr_t > 1) under stable stratification")
   end subroutine test_canuto_suppression

   subroutine test_momentum_conservation(error)
      !! The no-flux implicit momentum sweep redistributes momentum between
      !! interior layers but conserves the column total sum_k h_k u_k (and v).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), nut(NZ + 1, N), height_cv(N)
      real(wp), parameter :: DT = 10.0_wp
      real(wp) :: sum_u0, sum_v0, sum_u1, sum_v1
      integer :: k, step

      h = 1.0_wp
      do k = 1, NZ
         hu(k, :) = h(k, :)*real(k, wp)          ! u_k = k    (shear)
         hv(k, :) = h(k, :)*real(NZ - k, wp)     ! v_k = NZ-k (opposite shear)
      end do
      nut = 0.05_wp
      height_cv = real(NZ, wp)*1.0_wp
      sum_u0 = sum(hu(:, 4))
      sum_v0 = sum(hv(:, 4))

      !$omp target data map(to: h, nut, height_cv) map(tofrom: hu, hv)
      do step = 1, 10
         call ml_keps_apply_momentum_unstr_impl(NZ, N, h, hu, hv, nut, height_cv, DT)
      end do
      !$omp end target data

      sum_u1 = sum(hu(:, 4))
      sum_v1 = sum(hv(:, 4))
      call check(error, abs(sum_u1 - sum_u0) < 1.0e-10_wp*(abs(sum_u0) + 1.0_wp), &
                 "column x-momentum conserved under no-flux vertical diffusion")
      if (allocated(error)) return
      call check(error, abs(sum_v1 - sum_v0) < 1.0e-10_wp*(abs(sum_v0) + 1.0_wp), &
                 "column y-momentum conserved under no-flux vertical diffusion")
   end subroutine test_momentum_conservation

   subroutine test_momentum_shear_smoothing(error)
      !! Constant interior viscosity relaxes a sheared profile toward the
      !! depth-uniform mean.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hu(NZ, N), hv(NZ, N), nut(NZ + 1, N), height_cv(N)
      real(wp), parameter :: DT = 50.0_wp
      real(wp) :: u_mean, spread0, spread1
      integer :: k, step

      h = 1.0_wp
      hv = 0.0_wp
      do k = 1, NZ
         hu(k, :) = h(k, :)*real(k, wp)
      end do
      nut = 0.5_wp
      height_cv = real(NZ, wp)*1.0_wp
      u_mean = sum(hu(:, 4))/sum(h(:, 4))
      spread0 = maxval(hu(:, 4)/h(:, 4)) - minval(hu(:, 4)/h(:, 4))

      !$omp target data map(to: h, nut, height_cv) map(tofrom: hu, hv)
      do step = 1, 500
         call ml_keps_apply_momentum_unstr_impl(NZ, N, h, hu, hv, nut, height_cv, DT)
      end do
      !$omp end target data

      spread1 = maxval(hu(:, 4)/h(:, 4)) - minval(hu(:, 4)/h(:, 4))
      call check(error, spread1 < 0.01_wp*spread0, &
                 "vertical shear relaxes toward uniform under keps viscosity")
      if (allocated(error)) return
      call check(error, abs(hu(NZ, 4)/h(NZ, 4) - u_mean) < 0.05_wp, &
                 "top-layer velocity approaches the depth mean")
   end subroutine test_momentum_shear_smoothing

   subroutine test_tracer_mixing(error)
      !! Implicit keps tracer diffusion conserves the column tracer content
      !! sum_k hTr_k (no-flux) and relaxes a vertical gradient toward uniform.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NZ, N), hTr(NZ, N), kt(NZ + 1, N), height_cv(N)
      real(wp), parameter :: DT = 50.0_wp, KAPPA_BG = 1.0e-5_wp
      real(wp) :: sum0, sum1, spread0, spread1, tr_mean
      integer :: k, step

      h = 2.0_wp
      do k = 1, NZ
         hTr(k, :) = h(k, :)*(30.0_wp + real(k, wp))   ! Tr_k = 30 + k
      end do
      kt = 0.5_wp
      height_cv = real(NZ, wp)*2.0_wp
      sum0 = sum(hTr(:, 4))
      tr_mean = sum0/sum(h(:, 4))
      spread0 = maxval(hTr(:, 4)/h(:, 4)) - minval(hTr(:, 4)/h(:, 4))

      !$omp target data map(to: h, kt, height_cv) map(tofrom: hTr)
      do step = 1, 500
         call ml_keps_apply_diffusion_tracer_unstr_impl(NZ, N, h, hTr, kt, KAPPA_BG, height_cv, DT)
      end do
      !$omp end target data

      sum1 = sum(hTr(:, 4))
      spread1 = maxval(hTr(:, 4)/h(:, 4)) - minval(hTr(:, 4)/h(:, 4))
      call check(error, abs(sum1 - sum0) < 1.0e-9_wp*(abs(sum0) + 1.0_wp), &
                 "column tracer content conserved under no-flux vertical diffusion")
      if (allocated(error)) return
      call check(error, spread1 < 0.01_wp*spread0, &
                 "tracer gradient relaxes toward uniform under keps diffusion")
      if (allocated(error)) return
      call check(error, abs(hTr(NZ, 4)/h(NZ, 4) - tr_mean) < 0.05_wp, &
                 "surface-layer tracer approaches the depth mean")
   end subroutine test_tracer_mixing

end module test_keps_unstr
