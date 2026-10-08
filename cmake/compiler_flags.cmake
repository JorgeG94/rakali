# Compiler-specific flags for Rakali.
#
# Mutates CMAKE_Fortran_FLAGS / CMAKE_Fortran_FLAGS_{DEBUG,RELEASE} and may
# append to RAKALI_COMPILE_DEFS.  Kept as a global-flag include (rather than
# target-scoped INTERFACE libraries) because GPU offload (-stdpar=gpu / -acc)
# must be applied uniformly at compile *and* link time across every Fortran TU
# and dependent target — easiest to guarantee via the global flag set.

# Compilers that expose the `openacc` Fortran module + acc_set_device_num API
# regardless of whether OpenACC is the active codegen backend. Used by
# rki_mem_report to query device memory even on the OpenMP backend (do
# concurrent on these compilers still dispatches through the OpenACC runtime).
# Extend the list as new toolchains gain support (AMD flang OpenACC, etc).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID
                                                 STREQUAL "Cray")
  list(APPEND RAKALI_COMPILE_DEFS RAKALI_HAS_OPENACC_RUNTIME)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  # Common NVHPC flags
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfree -Mbackslash")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -Mbounds -Mchkptr -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")
  if(RAKALI_VERBOSE_COMPILE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Minfo=all")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minfo=all")
    set(CMAKE_Fortran_FLAGS_RELEASE
        "${CMAKE_Fortran_FLAGS_RELEASE} -Minfo=accel,opt")
  endif()

  if(RAKALI_ENABLE_GPU)
    # NVHPC consumes RAKALI_GPU_ARCH as `-gpu=ccNN`; the default (gfx90a, for
    # AMD/LLVM-flang) is invalid here, so fail fast with guidance rather than
    # letting NVHPC emit a cryptic codegen error.
    if(NOT RAKALI_GPU_ARCH MATCHES "^cc[0-9]+$")
      message(
        FATAL_ERROR
          "NVHPC GPU build needs RAKALI_GPU_ARCH=ccNN (e.g. cc70/cc80/cc90); got "
          "'${RAKALI_GPU_ARCH}'. The default is gfx90a for AMD/LLVM-flang — pass "
          "-DRAKALI_GPU_ARCH=cc70 for an NVIDIA NVHPC build.")
    endif()
    if(RAKALI_PARALLEL_BACKEND STREQUAL "openacc")
      # -acc=gpu is needed so OpenACC directives like `!$acc routine seq` on the
      # inlined `*_impl` pure helpers (e.g. in `rki_kernel_flux_unstr` and
      # `rki_kernel_extrapolate_unstr`) are explicitly enabled. Without it NVHPC
      # may parse the directives but not generate device-callable versions of
      # the helpers, leading to runtime CUDA_ERROR_ILLEGAL_ADDRESS inside `do
      # concurrent` GPU kernels that span module boundaries.
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=gpu -acc=gpu -gpu=${RAKALI_GPU_ARCH},mem:separate"
      )
      message(
        STATUS
          "GPU enabled (do concurrent + OpenACC, NVHPC): targeting ${RAKALI_GPU_ARCH}"
      )
    else() # openmp
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=gpu -mp=gpu -gpu=${RAKALI_GPU_ARCH},mem:separate"
      )
      message(
        STATUS
          "GPU enabled (do concurrent + OpenMP target, NVHPC): targeting ${RAKALI_GPU_ARCH}"
      )
    endif()
    # NVHPC 26.5 reads loop trip-counts from the device copy, inserting an extra
    # per-kernel data refresh in multi-loop `!$acc kernels` regions (the
    # barotropic fast loop) — ~2x on the barotropic solver. Revert to host
    # trip-counts until the compiler regression (NVIDIA TPR #38714) is fixed.
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 26.5)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -gpu=tripcount:host")
    endif()
    list(APPEND RAKALI_COMPILE_DEFS RAKALI_GPU_OFFLOAD)
  elseif(RAKALI_ENABLE_THREADS)
    if(RAKALI_PARALLEL_BACKEND STREQUAL "openacc")
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=multicore -acc=multicore")
      message(STATUS "CPU multicore (do concurrent + OpenACC threaded, NVHPC)")
    else() # openmp
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -stdpar=multicore -mp=multicore")
      message(STATUS "CPU multicore (do concurrent + OpenMP host, NVHPC)")
    endif()
  else()
    message(STATUS "CPU serial (no threading)")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
      "${CMAKE_Fortran_FLAGS} -ffree-form -std=gnu -O3 -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fcheck=all -fbacktrace -Wall -Wextra")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops")
  if(RAKALI_ENABLE_COVERAGE)
    # gcov instrumentation: disable optimisation so line counts map cleanly,
    # apply --coverage to both compile and link (gcc driver expands it to
    # -fprofile-arcs -ftest-coverage and links libgcov).
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    message(STATUS "Coverage enabled: -O0 -g --coverage (GNU/gfortran)")
  endif()
  if(RAKALI_VERBOSE_COMPILE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopt-info-loop")
    message(
      STATUS
        "GNU verbose compile: loop optimisation diagnostics enabled (-fopt-info-loop)"
    )
  endif()
  if(RAKALI_ENABLE_GPU)
    # Fail loud: RAKALI_ENABLE_GPU=ON here is always an explicit user request
    # (the auto-detect default is OFF for GNU) — a silent CPU-only downgrade
    # would misreport what was built.
    message(
      FATAL_ERROR
        "RAKALI_ENABLE_GPU=ON but gfortran has no `do concurrent` GPU offload "
        "path. Use the NVHPC toolchain (nvfortran) for GPU builds, or "
        "configure with -DRAKALI_ENABLE_GPU=OFF for a gfortran CPU build.")
  endif()
  if(RAKALI_ENABLE_THREADS)
    if(RAKALI_PARALLEL_BACKEND STREQUAL "openmp")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
      message(STATUS "CPU multicore (do concurrent + OpenMP host, gfortran)")
    else() # openacc / default
      set(RAKALI_NPROC
          "0"
          CACHE STRING
                "Thread count for -ftree-parallelize-loops (0 = auto-detect)")
      if(RAKALI_NPROC GREATER 0)
        set(_gnu_nproc ${RAKALI_NPROC})
      else()
        include(ProcessorCount)
        ProcessorCount(_gnu_nproc)
        if(_gnu_nproc EQUAL 0)
          set(_gnu_nproc 4)
        endif()
      endif()
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -ftree-parallelize-loops=${_gnu_nproc}")
      message(
        STATUS
          "CPU multicore (do concurrent, gfortran): -ftree-parallelize-loops=${_gnu_nproc}"
      )
    endif()
  else()
    message(STATUS "CPU serial (no threading)")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Classic ifort (CMAKE_Fortran_COMPILER_ID == "Intel"). Kept distinct from ifx
  # ("IntelLLVM") with STREQUAL — MATCHES would conflate the two. No OpenMP here
  # on purpose: ifort's `do concurrent` codegen is mature and the analytical
  # suite passes serially at -O3.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -stand f18")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model=precise")
  if(RAKALI_ENABLE_GPU)
    message(
      FATAL_ERROR
        "RAKALI_ENABLE_GPU=ON but classic ifort has no GPU offload path. Use "
        "NVHPC (nvfortran) for GPU builds, or -DRAKALI_ENABLE_GPU=OFF.")
  endif()
  message(STATUS "ifort host CPU (serial do concurrent)")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  # -qopenmp is mandatory for ifx, not optional. At -O2+ ifx miscompiles
  # multi-index `do concurrent` bodies that write via `intent(out)` scalars
  # through pure-subroutine calls — the writes silently become zero. The warning
  # is "Locality information is ignored without one of these command line
  # qualifiers '-qopenmp or -parallel'", and the analytical tests fail wholesale
  # (compute_flux_hll returns zeros for non-trivial gradients). See
  # tmp/ifx_repro/flux_repro.f90 for the standalone repro. Since `do concurrent`
  # is the project's primary parallel construct, we add -qopenmp
  # unconditionally; RAKALI_PARALLEL_BACKEND only gates GPU target offload here.
  #
  # -fp-model=precise pins away from the default fast-math at -O2+ (FMA
  # contraction + reassociation + reciprocal approx) so symmetry and
  # mass-conservation tests match gfortran / nvhpc IEEE behavior.
  #
  # -xHost is intentionally omitted: Gadi login and compute nodes can have
  # different ISA, so a login-node build with -xHost may emit instructions
  # unsupported on the compute node.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -stand f18 -qopenmp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model=precise")
  if(RAKALI_PARALLEL_BACKEND STREQUAL "openmp" AND RAKALI_ENABLE_GPU)
    # -fopenmp-target-do-concurrent lowers `do concurrent` onto the OpenMP
    # target region — ifx's do-concurrent GPU-offload path. Without it the DC
    # loops run on the host even with the target backend selected.
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -fopenmp-targets=spir64 -fopenmp-target-do-concurrent"
    )
    message(STATUS "ifx OpenMP target offload enabled "
                   "(-fopenmp-targets=spir64 -fopenmp-target-do-concurrent)")
  elseif(RAKALI_ENABLE_GPU)
    message(
      FATAL_ERROR
        "RAKALI_ENABLE_GPU=ON with ifx requires the OpenMP backend "
        "(-DRAKALI_PARALLEL_BACKEND=openmp -> -fopenmp-targets=spir64); the "
        "'${RAKALI_PARALLEL_BACKEND}' backend has no ifx GPU path. Either add "
        "-DRAKALI_PARALLEL_BACKEND=openmp or drop -DRAKALI_ENABLE_GPU=ON.")
  else()
    message(STATUS "ifx host CPU (-qopenmp, do concurrent)")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ")
  if(RAKALI_PARALLEL_BACKEND STREQUAL "openmp")
    if(RAKALI_ENABLE_GPU)
      # OpenMP target offload to an AMD/NVIDIA GPU.  `--offload-arch` selects
      # the device (RAKALI_GPU_ARCH, default gfx90a = MI200/MI250); applied at
      # compile AND link (global flags) so the offload images are produced +
      # bundled.
      #
      # -fPIC is REQUIRED on the AMD offload path: LLVM Flang emits
      # position-dependent relocations (R_X86_64_32/32S) by default, which
      # ld.lld rejects ("cannot be used against local symbol; recompile with
      # -fPIC") when the offload objects are bundled/device-linked. Scoped to
      # the GPU-offload build only (this branch is already gated on
      # RAKALI_ENABLE_GPU) so host-only builds are unaffected.
      #
      # -fdo-concurrent-to-openmp=device lowers `do concurrent` onto an OpenMP
      # target region (the GPU); =host would map to host multicore, =none leaves
      # it serial. This is the Flang DC-offload switch, mirroring ifx's
      # -fopenmp-target-do-concurrent above.
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -fopenmp --offload-arch=${RAKALI_GPU_ARCH} -fPIC -fdo-concurrent-to-openmp=device"
      )
      message(
        STATUS
          "LLVM Flang OpenMP target offload "
          "(--offload-arch=${RAKALI_GPU_ARCH} -fdo-concurrent-to-openmp=device)"
      )
      list(APPEND RAKALI_COMPILE_DEFS RAKALI_GPU_OFFLOAD)
    else()
      # Host multicore: map `do concurrent` onto host OpenMP threads.
      set(CMAKE_Fortran_FLAGS
          "${CMAKE_Fortran_FLAGS} -fopenmp -fdo-concurrent-to-openmp=host")
      message(STATUS "LLVM Flang OpenMP host "
                     "(-fopenmp -fdo-concurrent-to-openmp=host)")
    endif()
  elseif(RAKALI_ENABLE_GPU)
    message(
      FATAL_ERROR
        "RAKALI_ENABLE_GPU=ON with LLVM Flang requires the OpenMP backend "
        "(-DRAKALI_PARALLEL_BACKEND=openmp -> --offload-arch); the "
        "'${RAKALI_PARALLEL_BACKEND}' backend has no Flang GPU path. Either "
        "add -DRAKALI_PARALLEL_BACKEND=openmp or drop -DRAKALI_ENABLE_GPU=ON.")
  else()
    # No threading: leave `do concurrent` serial (no OpenMP mapping).
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -fdo-concurrent-to-openmp=none")
    message(STATUS "LLVM Flang host CPU "
                   "(serial do concurrent, -fdo-concurrent-to-openmp=none)")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 --cpp --std=f23")
  set(CMAKE_Fortran_FLAGS_DEBUG "")
  set(CMAKE_Fortran_FLAGS_RELEASE "")
endif()

# Applied LAST, via add_compile_options rather than by appending to
# CMAKE_Fortran_FLAGS: CMake emits CMAKE_Fortran_FLAGS, then
# CMAKE_Fortran_FLAGS_<CONFIG>, then COMPILE_OPTIONS, so only this form can
# override something set on the Release line (e.g. `-march=native`). Directory
# scope, set before any target or fetched dependency exists, so the subprojects
# compile with it too.
if(RAKALI_EXTRA_FORTRAN_FLAGS)
  separate_arguments(_rki_extra_flags NATIVE_COMMAND
                     "${RAKALI_EXTRA_FORTRAN_FLAGS}")
  foreach(_rki_flag IN LISTS _rki_extra_flags)
    add_compile_options("$<$<COMPILE_LANGUAGE:Fortran>:${_rki_flag}>")
  endforeach()
  unset(_rki_flag)
  unset(_rki_extra_flags)
endif()
