!! Unit tests for the k-epsilon vertical-mixing column kernel.
module test_keps
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use rki_constants, only: wp, KEPS_STAB_CONSTANT, KEPS_STAB_GALPERIN, KEPS_STAB_CANUTO
   use rki_ml_keps, only: ml_keps_update, ml_keps_apply_momentum, &
                          ml_keps_apply_diffusion_tracer
   use rki_grid, only: hgrid_t
   implicit none
   private

   public :: collect_keps_tests

   integer, parameter :: NX = 8, NY = 8, NG = 2, NZ = 4
   real(wp), parameter :: DX = 1.0_wp, DY = 1.0_wp
   ! k-eps constants (must match the kernel) for the analytic check.
   real(wp), parameter :: CE2 = 1.92_wp, CMU = 0.09_wp

   type(hgrid_t), save :: grid
   logical, save :: grid_ready = .false.

contains

   subroutine ensure_grid()
      if (grid_ready) return
      call grid%init(NX - 2*NG, NY - 2*NG, NG, DX, DY)
      grid_ready = .true.
   end subroutine ensure_grid

   subroutine collect_keps_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      call ensure_grid()
      testsuite = [ &
                  new_unittest("keps_homogeneous_decay", test_homogeneous_decay), &
                  new_unittest("keps_positivity", test_positivity), &
                  new_unittest("keps_momentum_conservation", test_momentum_conservation), &
                  new_unittest("keps_momentum_shear_smoothing", test_momentum_shear_smoothing), &
                  new_unittest("keps_wall_function", test_wall_function), &
                  new_unittest("keps_tracer_mixing", test_tracer_mixing), &
                  new_unittest("keps_galperin_suppression", test_galperin_suppression), &
                  new_unittest("keps_canuto_suppression", test_canuto_suppression) &
                  ]
   end subroutine collect_keps_tests

   subroutine test_homogeneous_decay(error)
      !! No shear, no stratification (uniform u, rho): k-eps reduces to
      !! decaying turbulence with the analytic power law
      !!   k(t) = k0 (1 + (ce2-1) eps0/k0 t)^(-1/(ce2-1)).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ), rho(NX, NY, NZ)
      real(wp) :: tke(NX, NY, NZ + 1), eps(NX, NY, NZ + 1)
      real(wp) :: nut(NX, NY, NZ + 1), kt(NX, NY, NZ + 1)
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
      nstep = nint(T_END/DT)

      !$acc data copyin( h, hu, hv, rho) copy( tke, eps, nut, kt)
      do step = 1, nstep
         call ml_keps_update(h, hu, hv, rho, tke, eps, nut, kt, 1025.0_wp, DT, grid, NZ)
      end do
      !$acc end data

      k_exact = K0*(1.0_wp + (CE2 - 1.0_wp)*E0/K0*T_END)**(-1.0_wp/(CE2 - 1.0_wp))
      call check(error, abs(tke(4, 4, 3) - k_exact)/k_exact < 0.02_wp, &
                 "k-eps homogeneous decay matches analytic power law")
   end subroutine test_homogeneous_decay

   subroutine test_positivity(error)
      !! Strong shear with a large step must keep k, eps, nu_t positive and
      !! bounded (the semi-implicit sinks guarantee positivity).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ), rho(NX, NY, NZ)
      real(wp) :: tke(NX, NY, NZ + 1), eps(NX, NY, NZ + 1)
      real(wp) :: nut(NX, NY, NZ + 1), kt(NX, NY, NZ + 1)
      integer :: k

      h = 2.0_wp
      hv = 0.0_wp
      rho = 1025.0_wp
      ! strong vertical shear: u increases with layer index
      do k = 1, NZ
         hu(:, :, k) = h(:, :, k)*real(k, wp)*0.5_wp
      end do
      tke = 1.0e-6_wp
      eps = 1.0e-8_wp
      nut = 1.0e-4_wp
      kt = nut

      !$acc data copyin( h, hu, hv, rho) copy( tke, eps, nut, kt)
      do k = 1, 200
         call ml_keps_update(h, hu, hv, rho, tke, eps, nut, kt, 1025.0_wp, 100.0_wp, grid, NZ)
      end do
      !$acc end data

      call check(error, all(tke(NG + 1:NX - NG, NG + 1:NY - NG, :) > 0.0_wp), &
                 "k stays positive under strong shear + large dt")
      if (allocated(error)) return
      call check(error, all(eps(NG + 1:NX - NG, NG + 1:NY - NG, :) > 0.0_wp) .and. &
                 all(nut(NG + 1:NX - NG, NG + 1:NY - NG, :) <= 10.0_wp), &
                 "eps positive and nu_t bounded")
      if (allocated(error)) return
      ! shear production should have grown the TKE above the tiny seed
      call check(error, tke(4, 4, 3) > 1.0e-4_wp, &
                 "shear production grows TKE")
   end subroutine test_positivity

   subroutine test_momentum_conservation(error)
      !! The no-flux implicit momentum sweep redistributes momentum between
      !! interior layers but conserves the column total sum_k h_k u_k (and
      !! sum_k h_k v_k) to round-off, independent of nu_t.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ)
      real(wp) :: nut(NX, NY, NZ + 1)
      real(wp), parameter :: DT = 10.0_wp
      real(wp) :: sum_u0, sum_v0, sum_u1, sum_v1
      integer :: k, step

      h = 1.0_wp
      do k = 1, NZ
         hu(:, :, k) = h(:, :, k)*real(k, wp)          ! u_k = k    (shear)
         hv(:, :, k) = h(:, :, k)*real(NZ - k, wp)     ! v_k = NZ-k (opposite shear)
      end do
      nut = 0.05_wp

      sum_u0 = sum(hu(4, 4, :))
      sum_v0 = sum(hv(4, 4, :))

      !$acc data copyin( h, nut) copy( hu, hv)
      do step = 1, 10
         call ml_keps_apply_momentum(h, hu, hv, nut, grid, NZ, DT)
      end do
      !$acc end data

      sum_u1 = sum(hu(4, 4, :))
      sum_v1 = sum(hv(4, 4, :))

      call check(error, abs(sum_u1 - sum_u0) < 1.0e-10_wp*(abs(sum_u0) + 1.0_wp), &
                 "column x-momentum conserved under no-flux vertical diffusion")
      if (allocated(error)) return
      call check(error, abs(sum_v1 - sum_v0) < 1.0e-10_wp*(abs(sum_v0) + 1.0_wp), &
                 "column y-momentum conserved under no-flux vertical diffusion")
   end subroutine test_momentum_conservation

   subroutine test_momentum_shear_smoothing(error)
      !! With a constant interior viscosity the implicit sweep relaxes a
      !! sheared profile toward the depth-uniform mean.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ)
      real(wp) :: nut(NX, NY, NZ + 1)
      real(wp), parameter :: DT = 50.0_wp
      real(wp) :: u_mean, spread0, spread1
      integer :: k, step

      h = 1.0_wp
      hv = 0.0_wp
      do k = 1, NZ
         hu(:, :, k) = h(:, :, k)*real(k, wp)
      end do
      nut = 0.5_wp
      u_mean = sum(hu(4, 4, :))/sum(h(4, 4, :))
      spread0 = maxval(hu(4, 4, :)/h(4, 4, :)) - minval(hu(4, 4, :)/h(4, 4, :))

      !$acc data copyin( h, nut) copy( hu, hv)
      do step = 1, 500
         call ml_keps_apply_momentum(h, hu, hv, nut, grid, NZ, DT)
      end do
      !$acc end data

      spread1 = maxval(hu(4, 4, :)/h(4, 4, :)) - minval(hu(4, 4, :)/h(4, 4, :))
      call check(error, spread1 < 0.01_wp*spread0, &
                 "vertical shear relaxes toward uniform under k-eps viscosity")
      if (allocated(error)) return
      call check(error, abs(hu(4, 4, NZ)/h(4, 4, NZ) - u_mean) < 0.05_wp, &
                 "top-layer velocity approaches the depth mean")
   end subroutine test_momentum_shear_smoothing

   subroutine test_wall_function(error)
      !! With prescribed friction velocities the near-wall interfaces are
      !! pinned to the Dirichlet log-law value  k = u_*^2 / sqrt(c_mu)
      !! (bed at interface 2, surface at interface NZ).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ), rho(NX, NY, NZ)
      real(wp) :: tke(NX, NY, NZ + 1), eps(NX, NY, NZ + 1)
      real(wp) :: nut(NX, NY, NZ + 1), kt(NX, NY, NZ + 1)
      real(wp) :: ustar_b(NX, NY), ustar_s(NX, NY)
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

      !$acc data copyin( h, hu, hv, rho, ustar_b, ustar_s) copy( tke, eps, nut, kt)
      do step = 1, 5
         call ml_keps_update(h, hu, hv, rho, tke, eps, nut, kt, 1025.0_wp, 10.0_wp, grid, NZ, &
                             ustar_b, ustar_s)
      end do
      !$acc end data

      k_wb_exp = USB*USB/sqrt(CMU)
      k_ws_exp = USS*USS/sqrt(CMU)
      call check(error, abs(tke(4, 4, 2) - k_wb_exp)/k_wb_exp < 1.0e-6_wp, &
                 "bed wall pins k = u_*b^2/sqrt(c_mu) at interface 2")
      if (allocated(error)) return
      call check(error, abs(tke(4, 4, NZ) - k_ws_exp)/k_ws_exp < 1.0e-6_wp, &
                 "surface wall pins k = u_*s^2/sqrt(c_mu) at interface NZ")
   end subroutine test_wall_function

   subroutine test_tracer_mixing(error)
      !! Implicit k-eps tracer diffusion conserves the column tracer content
      !! sum_k hTr_k (no-flux) and relaxes a vertical gradient toward uniform.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hTr(NX, NY, NZ), kt(NX, NY, NZ + 1)
      real(wp), parameter :: DT = 50.0_wp, KAPPA_BG = 1.0e-5_wp
      real(wp) :: sum0, sum1, spread0, spread1, tr_mean
      integer :: k, step

      h = 2.0_wp
      do k = 1, NZ
         hTr(:, :, k) = h(:, :, k)*(30.0_wp + real(k, wp))   ! Tr_k = 30 + k
      end do
      kt = 0.5_wp
      sum0 = sum(hTr(4, 4, :))
      tr_mean = sum0/sum(h(4, 4, :))
      spread0 = maxval(hTr(4, 4, :)/h(4, 4, :)) - minval(hTr(4, 4, :)/h(4, 4, :))

      !$acc data copyin( h, kt) copy( hTr)
      do step = 1, 500
         call ml_keps_apply_diffusion_tracer(h, hTr, kt, KAPPA_BG, grid, NZ, DT)
      end do
      !$acc end data

      sum1 = sum(hTr(4, 4, :))
      spread1 = maxval(hTr(4, 4, :)/h(4, 4, :)) - minval(hTr(4, 4, :)/h(4, 4, :))

      call check(error, abs(sum1 - sum0) < 1.0e-9_wp*(abs(sum0) + 1.0_wp), &
                 "column tracer content conserved under no-flux vertical diffusion")
      if (allocated(error)) return
      call check(error, spread1 < 0.01_wp*spread0, &
                 "tracer gradient relaxes toward uniform under k-eps diffusion")
      if (allocated(error)) return
      call check(error, abs(hTr(4, 4, NZ)/h(4, 4, NZ) - tr_mean) < 0.05_wp, &
                 "surface-layer tracer approaches the depth mean")
   end subroutine test_tracer_mixing

   subroutine test_galperin_suppression(error)
      !! Under stable stratification the Galperin stability function suppresses
      !! the recovered eddy viscosity below the constant-c_mu value.  One step
      !! from an identical stratified, sheared state with each scheme: the k/eps
      !! update itself is scheme-independent (the scheme only enters the
      !! nu_t/kappa_t recovery), so nu_t(galperin)/nu_t(constant) =
      !! c_mu(alpha_N)/0.09 < 1 under stable stratification.
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ), rho(NX, NY, NZ)
      real(wp) :: tke(NX, NY, NZ + 1), eps(NX, NY, NZ + 1)
      real(wp) :: nut(NX, NY, NZ + 1), kt(NX, NY, NZ + 1)
      real(wp) :: nut_c, nut_g
      integer :: k

      h = 2.0_wp
      hv = 0.0_wp
      ! vertical shear + stable stratification (rho decreases upward:
      ! k=1 dense/bed, k=NZ light/surface), so N^2 > 0 in the interior.
      do k = 1, NZ
         hu(:, :, k) = h(:, :, k)*real(k, wp)*0.2_wp
         rho(:, :, k) = 1025.0_wp - 1.5_wp*real(k - 1, wp)
      end do

      ! constant scheme
      tke = 0.05_wp
      eps = 5.0e-3_wp
      nut = 1.0e-3_wp
      kt = 1.0e-3_wp
      !$acc data copyin( h, hu, hv, rho) copy( tke, eps, nut, kt)
      call ml_keps_update(h, hu, hv, rho, tke, eps, nut, kt, 1025.0_wp, 10.0_wp, &
                          grid, NZ, scheme=KEPS_STAB_CONSTANT)
      !$acc end data
      nut_c = nut(4, 4, 3)

      ! galperin scheme, identical initial state
      tke = 0.05_wp
      eps = 5.0e-3_wp
      nut = 1.0e-3_wp
      kt = 1.0e-3_wp
      !$acc data copyin( h, hu, hv, rho) copy( tke, eps, nut, kt)
      call ml_keps_update(h, hu, hv, rho, tke, eps, nut, kt, 1025.0_wp, 10.0_wp, &
                          grid, NZ, scheme=KEPS_STAB_GALPERIN)
      !$acc end data
      nut_g = nut(4, 4, 3)

      call check(error, nut_g < nut_c, &
                 "Galperin suppresses nu_t below constant-c_mu under stable stratification")
      if (allocated(error)) return
      call check(error, nut_g < 0.9_wp*nut_c, &
                 "Galperin nu_t suppression is meaningful (< 0.9x constant)")
   end subroutine test_galperin_suppression

   subroutine test_canuto_suppression(error)
      !! Canuto-A (Canuto et al. 2001 Model A) under stable stratification: it
      !! suppresses the scalar diffusivity more than the eddy viscosity, so
      !! kt < nu_t (turbulent Prandtl number > 1) -- distinct from constant
      !! (kt = nu_t) and Galperin (neutral Pr_t = 0.8, kt > nu_t).  One step from
      !! a stratified, sheared column; the rational stability functions stay
      !! positive (denominator well-posed).
      type(error_type), allocatable, intent(out) :: error
      real(wp) :: h(NX, NY, NZ), hu(NX, NY, NZ), hv(NX, NY, NZ), rho(NX, NY, NZ)
      real(wp) :: tke(NX, NY, NZ + 1), eps(NX, NY, NZ + 1)
      real(wp) :: nut(NX, NY, NZ + 1), kt(NX, NY, NZ + 1)
      real(wp) :: nut_c, kt_c
      integer :: k

      h = 2.0_wp
      hv = 0.0_wp
      ! shear + strong stable stratification (rho decreases upward), so the
      ! buoyancy number is comfortably positive after one step.
      do k = 1, NZ
         hu(:, :, k) = h(:, :, k)*real(k, wp)*0.2_wp
         rho(:, :, k) = 1025.0_wp - 3.0_wp*real(k - 1, wp)
      end do
      tke = 0.05_wp
      eps = 5.0e-3_wp
      nut = 1.0e-3_wp
      kt = 1.0e-3_wp

      !$acc data copyin( h, hu, hv, rho) copy( tke, eps, nut, kt)
      call ml_keps_update(h, hu, hv, rho, tke, eps, nut, kt, 1025.0_wp, 10.0_wp, &
                          grid, NZ, scheme=KEPS_STAB_CANUTO)
      !$acc end data
      nut_c = nut(4, 4, 3)
      kt_c = kt(4, 4, 3)

      call check(error, nut_c > 1.0e-4_wp .and. kt_c > 1.0e-4_wp, &
                 "Canuto-A nu_t, kt are positive physical values (rational form well-posed)")
      if (allocated(error)) return
      call check(error, kt_c < nut_c, &
                 "Canuto-A scalar diffusivity < eddy viscosity (Pr_t > 1) under stable stratification")
   end subroutine test_canuto_suppression

end module test_keps
