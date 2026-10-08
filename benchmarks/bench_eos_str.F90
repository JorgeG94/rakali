!! Equation-of-state kernel benchmark, structured-grid layout (single GPU)
program bench_eos_str
   !! Mirror of `bench_eos` for the structured kernel.  Same math, but
   !! the array layout is `(nx_total, ny_total, nz)` and the loop is
   !!   do concurrent (k=1:nz, j=1:ny, i=1:nx)
   !! so the fast axis is `i`, not `k` -- a different gather pattern
   !! than the unstructured `(nz, ncells)` layout the sister bench
   !! `bench_eos` measures.  Run both for an A/B on layout cost; the
   !! flop / byte counts per cell-layer are identical so the numbers
   !! are directly comparable.
   !!
   !! `nghost` is set to 0 here so the kernel does exactly `nx*ny*nz`
   !! work -- matches the `ncells*nz` work of the unstr bench.
   !!
   !! Usage:
   !!   bench_eos_str                       Default sweep (matches bench_eos)
   !!   bench_eos_str <ncells>              Single ncells, default nz
   !!   bench_eos_str <ncells> <nz>         Single (ncells, nz)
   !!   bench_eos_str <ncells> <nz> <reps>  Custom rep count
   !!
   !! `ncells` is split into `nx * ny` with `nx ~ sqrt(ncells)` and
   !! `ny = ceil(ncells/nx)`, so the actual cell count may exceed the
   !! requested value by a few rows -- the CSV reports what was used.
   !!
   !! Human-readable lines on stderr (logger).  CSV on stdout.
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_ml_eos, only: ml_equation_of_state
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_timer, only: timer_type
   implicit none

   character(len=*), parameter :: BACKEND = &
#if defined(RAKALI_GPU_OFFLOAD)
                                  "openacc"
#else
   "cpu_serial"
#endif

   real(wp), parameter :: RHO_0 = 1025.0_wp
   real(wp), parameter :: BETA_S = 0.78_wp
   real(wp), parameter :: ALPHA_T = 0.17_wp
   real(wp), parameter :: S_REF = 0.0_wp
   real(wp), parameter :: T_REF = 15.0_wp

   integer, allocatable :: ncells_sweep(:), nz_sweep(:)
   integer :: timed_reps, warmup_reps, val
   integer :: i_n, i_nz
   character(len=32) :: arg

   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      read (arg, *) val
      if (val > 0) then
         allocate (ncells_sweep(1))
         ncells_sweep(1) = val
      end if
   end if
   if (.not. allocated(ncells_sweep)) then
#if defined(RAKALI_GPU_OFFLOAD)
      ncells_sweep = [10000, 100000, 1000000, 4000000]
#else
      ncells_sweep = [1000, 10000, 100000]
#endif
   end if

   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      read (arg, *) val
      if (val > 0) then
         allocate (nz_sweep(1))
         nz_sweep(1) = val
      end if
   end if
   if (.not. allocated(nz_sweep)) then
      nz_sweep = [5, 10, 20, 40]
   end if

   timed_reps = 0
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      read (arg, *) timed_reps
   end if
   if (timed_reps <= 0) then
#if defined(RAKALI_GPU_OFFLOAD)
      timed_reps = 500
#else
      timed_reps = 100
#endif
   end if
   warmup_reps = max(5, timed_reps/10)

   write (*, '(a)') "backend,nx,ny,ncells,nz,warmup,reps,total_s,"// &
      "cell_layers_per_s,bytes_per_s_GB,flops_per_s_Gflops,max_rho_diff"

   do i_n = 1, size(ncells_sweep)
      do i_nz = 1, size(nz_sweep)
         call run_eos_str(ncells_sweep(i_n), nz_sweep(i_nz), &
                          warmup_reps, timed_reps)
      end do
   end do

contains

   subroutine run_eos_str(ncells_target, nz, n_warmup, n_timed)
      integer, intent(in) :: ncells_target, nz, n_warmup, n_timed

      type(hgrid_t) :: grid
      real(wp), allocatable :: h_layer(:, :, :), hS_layer(:, :, :), hT_layer(:, :, :)
      real(wp), allocatable :: rho_layer(:, :, :), rho_ref(:, :, :)
      type(timer_type) :: tmr
      integer :: r, i, j, k, nx_phys, ny_phys, nx, ny, ncells
      real(wp) :: elapsed, cl_rate, gbps, gflops, drift
      real(wp) :: bytes_per_step, flops_per_step
      character(len=512) :: csv

      ! Split ncells_target into a near-square (nx, ny) grid.  The
      ! kernel touches every cell in [1:nx, 1:ny, 1:nz] so the actual
      ! work is nx*ny*nz; report it accurately.
      nx_phys = nint(sqrt(real(ncells_target, wp)))
      if (nx_phys < 1) nx_phys = 1
      ny_phys = (ncells_target + nx_phys - 1)/nx_phys

      ! `nghost = 0` so nx_total = nx_phys, ny_total = ny_phys -- the
      ! kernel work matches the unstr bench's `ncells * nz` exactly.
      call grid%init(nx_phys, ny_phys, 0, 1.0_wp, 1.0_wp)
      nx = grid%nx_total
      ny = grid%ny_total
      ncells = nx*ny

      allocate (h_layer(nx, ny, nz))
      allocate (hS_layer(nx, ny, nz))
      allocate (hT_layer(nx, ny, nz))
      allocate (rho_layer(nx, ny, nz))
      allocate (rho_ref(nx, ny, nz))

      call build_columns_str(h_layer, hS_layer, hT_layer, nx, ny, nz)

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               rho_ref(i, j, k) = RHO_0 &
                                  + BETA_S*(hS_layer(i, j, k)/h_layer(i, j, k) - S_REF) &
                                  - ALPHA_T*(hT_layer(i, j, k)/h_layer(i, j, k) - T_REF)
            end do
         end do
      end do

      !$acc enter data copyin(h_layer, hS_layer, hT_layer) create(rho_layer)

      do r = 1, n_warmup
         call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                   RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, &
                                   grid, nz)
      end do
      call gpu_sync()

      call tmr%start()
      do r = 1, n_timed
         call ml_equation_of_state(h_layer, hS_layer, hT_layer, rho_layer, &
                                   RHO_0, BETA_S, S_REF, ALPHA_T, T_REF, &
                                   grid, nz)
      end do
      call gpu_sync()
      call tmr%stop()
      elapsed = tmr%get_elapsed_time()

      !$acc update self(rho_layer)

      drift = 0.0_wp
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               drift = max(drift, abs(rho_layer(i, j, k) - rho_ref(i, j, k)))
            end do
         end do
      end do

      !$acc exit data delete(h_layer, hS_layer, hT_layer, rho_layer)

      cl_rate = real(ncells, wp)*real(nz, wp)*real(n_timed, wp)/elapsed
      bytes_per_step = real(ncells, wp)*real(nz, wp)*4.0_wp*8.0_wp
      gbps = bytes_per_step*real(n_timed, wp)/elapsed/1.0e9_wp
      flops_per_step = real(ncells, wp)*real(nz, wp)*6.0_wp
      gflops = flops_per_step*real(n_timed, wp)/elapsed/1.0e9_wp

      call logger%info("--- EOS-str nx="//to_string(nx)//" ny="//to_string(ny) &
                       //" ("//to_string(ncells)//" cells) nz="//to_string(nz) &
                       //" ["//BACKEND//"] " &
                       //to_string(n_warmup)//"w+"//to_string(n_timed)//"t ---")
      call logger%info("  Total time:         "//to_string(elapsed)//" s")
      call logger%info("  Cell-layers/s:      "//to_string(cl_rate))
      call logger%info("  Effective BW:       "//to_string(gbps)//" GB/s")
      call logger%info("  Achieved Gflops/s:  "//to_string(gflops))
      call logger%info("  Max abs(rho - ref): "//to_string(drift))

      write (csv, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",'// &
             'es12.5,",",es12.5,",",es12.5,",",es12.5,",",es10.3)') &
         BACKEND, nx, ny, ncells, nz, n_warmup, n_timed, &
         elapsed, cl_rate, gbps, gflops, drift
      write (*, '(a)') trim(csv)

      deallocate (h_layer, hS_layer, hT_layer, rho_layer, rho_ref)
   end subroutine run_eos_str

   subroutine build_columns_str(h, hS, hT, nx, ny, nz)
      !! Stratified IC mirroring `bench_eos`: S 35 -> 30 PSU, T 15 -> 20
      !! degC, h gently varying with depth.  The horizontal axes are
      !! filled identically across all (i, j) -- the bench is about the
      !! kernel's per-cell cost, not horizontal patterns.
      integer, intent(in) :: nx, ny, nz
      real(wp), intent(out) :: h(nx, ny, nz), hS(nx, ny, nz), hT(nx, ny, nz)

      real(wp) :: frac, hk, S, T
      integer :: i, j, k

      do k = 1, nz
         frac = real(k - 1, wp)/real(max(nz - 1, 1), wp)
         hk = 2.0_wp + 0.5_wp*frac
         S = 35.0_wp - 5.0_wp*frac
         T = 15.0_wp + 5.0_wp*frac
         do j = 1, ny
            do i = 1, nx
               h(i, j, k) = hk
               hS(i, j, k) = hk*S
               hT(i, j, k) = hk*T
            end do
         end do
      end do
   end subroutine build_columns_str

   subroutine gpu_sync()
      !! No-op (acc directives synchronous under stdpar).
   end subroutine gpu_sync

end program bench_eos_str
