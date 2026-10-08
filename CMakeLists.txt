cmake_minimum_required(VERSION 3.25)
project(
  rakali
  VERSION 0.1.0
  LANGUAGES Fortran)

set(main_lib core_rakali)
set(exe_name rakali)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Build options

# GPU offload is OPT-IN (default OFF). Enable with -DRAKALI_ENABLE_GPU=ON on the
# NVHPC toolchain (the production `do concurrent` + OpenACC GPU path). An
# explicit ON with a compiler that has no GPU offload path is a configure-time
# FATAL_ERROR (see cmake/compiler_flags.cmake) — never a silent CPU-only
# downgrade. The value is reported on the "GPU offload:" line in the build
# summary below.
option(RAKALI_ENABLE_GPU "Build with GPU offloading (do concurrent + OpenACC)"
       OFF)

option(RAKALI_ENABLE_THREADS
       "Enable multicore threading (do concurrent + OpenACC multicore)" OFF)

# Parallel backend choice. Most builds on `main` use OpenACC (NVHPC path) — the
# OpenMP variant is consumed via the auto-regenerated `auto/dc-openmp` branch,
# which sets this to `openmp` and applies overlay patches under patches/openmp/.
# Setting `openmp` on `main` directly only changes compiler flags — the source
# still carries `!$acc` directives until acc_to_omp.py runs.
set(RAKALI_PARALLEL_BACKEND
    "openacc"
    CACHE STRING "Parallel backend: openacc | openmp (advanced; leave default)")
set_property(CACHE RAKALI_PARALLEL_BACKEND PROPERTY STRINGS openacc openmp)
mark_as_advanced(RAKALI_PARALLEL_BACKEND)
if(NOT RAKALI_PARALLEL_BACKEND MATCHES "^(openacc|openmp)$")
  message(FATAL_ERROR "RAKALI_PARALLEL_BACKEND must be one of: openacc, openmp "
                      "(got '${RAKALI_PARALLEL_BACKEND}')")
endif()

option(RAKALI_ENABLE_DOUBLE "Use double precision (real64)" ON)
option(RAKALI_ENABLE_TESTING "Build test suite" ON)
# Escape hatch for flags compiler_flags.cmake has no opinion about. Motivating
# case: `-tp` on NVHPC / `-march` on GNU+ifx. Both target the BUILD HOST's CPU
# by default, so a binary built on one machine can die with SIGILL on another.
# CI pins a baseline; a cluster user building on a login node wants the same.
set(RAKALI_EXTRA_FORTRAN_FLAGS
    ""
    CACHE STRING
          "Extra Fortran flags, appended after every flag this project sets")
option(
  RAKALI_BUILD_BENCHMARKS
  "Build benchmark executables (dev-only; some drive core GPU device routines directly, so they cannot link against a shared core — keep OFF with RAKALI_BUILD_SHARED)"
  OFF)
option(RAKALI_VERBOSE_COMPILE "Show detailed compiler diagnostics (-Minfo=all)"
       OFF)
option(RAKALI_ENABLE_NVTX
       "Enable NVTX profiling annotations for Nsight Systems" OFF)
option(RAKALI_BUILD_SHARED
       "Build librakali as a shared (.so) library instead of a static archive"
       OFF)
option(
  RAKALI_TESTS_LINK_SHARED
  "Link tests against a shared core (.so) for faster test builds. ON by default: only the tests in tests/CMakeLists.txt's RAKALI_STATIC_CORE_TESTS list (device code in their own source) keep the static core; everything else links the .so. Set OFF to force every test back to the static core."
  ON)
option(
  RAKALI_ENABLE_PYTHON
  "Build the C-binding shared library (librakali.so) for the Python/Julia API"
  OFF)
option(
  RAKALI_ENABLE_NETCDF
  "Build the NetCDF-backed I/O subsystem (driver, output, gauges, restart, forcing, nesting, public API). Disable for portability testing on AMD/Intel where bringing along NetCDF is inconvenient — kernels, solver and benchmarks still build."
  ON)
option(
  RAKALI_ENABLE_COVERAGE
  "Build with gcov instrumentation (-O0 -g --coverage). GNU/gfortran only; requires lcov + genhtml on PATH for the `coverage` target."
  OFF)

# LFortran 0.64 (alpha) needs a handful of source-level workarounds
# (module-local parameter hoists for the imported-parameter purity bug; a real
# 2.0**k in place of scale() for the negative-exponent bug).  They are all
# guarded by the LFORTRAN_PASSING macro so the default build compiles the REAL
# code — pass -DRAKALI_LFORTRAN_PASSING=ON to build lfortran 0.64, and leave it
# OFF to test a newer lfortran against the unmodified code (the guards then
# vanish).  MREs + details in tmp_local_artifacts/lf_nz_mre/.
option(
  RAKALI_LFORTRAN_PASSING
  "Activate the LFortran 0.64 source workarounds (defines LFORTRAN_PASSING)"
  OFF)
if(RAKALI_LFORTRAN_PASSING)
  add_compile_definitions(LFORTRAN_PASSING)
endif()

# Max vertical layers for per-column fixed-size stack workspace inside
# do-concurrent kernels (NZ_STACK_MAX in rki_constants).  Lower it on the
# AMD/flang offload path when the device link reports "stack frame size exceeds
# limit" — the per-work-item frame of the heavy column solvers (e.g.
# kappa-shear) scales with this.  Must stay >= nz + 1 (the model's vertical
# layer count plus one) — `validate_config` refuses the run fail-loud otherwise.
# NOT 2*nz+2: the Redi neutral-surface locals are declared 2*NZ_STACK_MAX+2 and
# scale with this value.  128 covers the largest shipped case (nz=90) and test
# (nz=100).
set(RAKALI_NZ_STACK_MAX
    "128"
    CACHE
      STRING
      "Max vertical layers for per-column stack workspace (lower for AMD/flang if the device link hits stack-frame limits)"
)

# Unstructured-multilayer per-column stack bound (BPG_NZ_MAX in
# rki_ml_common_unstr) — DISTINCT from RAKALI_NZ_STACK_MAX and deliberately
# small: the unstructured path is thread-per-cell with many such locals per
# kernel frame, so raising this multiplies per-thread local memory across every
# unstructured column kernel.  The largest supported layer count is BPG_NZ_MAX -
# 1 (the k-eps interface tridiagonal indexes to nz + 1); validate_config
# enforces that fail-loud.  Raise only if you need deeper unstructured columns.
set(RAKALI_BPG_NZ_MAX
    "20"
    CACHE
      STRING
      "Max unstructured-multilayer stack bound (largest legal nz_layers is this minus 1; raise only for deeper unstructured columns)"
)

# Python/Julia API depends on the NetCDF I/O modules transitively (rki_api ->
# rki_output, rki_handle -> rki_output_unstr).  Reject the combination up front
# instead of failing later with a confusing missing- module error.
if(RAKALI_ENABLE_PYTHON AND NOT RAKALI_ENABLE_NETCDF)
  message(
    FATAL_ERROR "RAKALI_ENABLE_PYTHON requires RAKALI_ENABLE_NETCDF=ON "
                "(the C API uses the NetCDF-backed output/restart modules).")
endif()

# LFortran (alpha) can't link a macOS shared library: CMake's Apple platform
# rules inject the SONAME flag `-install_name @rpath/lib….dylib`, which
# lfortran's CLI rejects ("argument not expected: -install_name").  The shared
# core is only a build-speed / FFI convenience, and the static-core path is
# fully supported (tests fall back to it per-target), so force it off here — no
# .dylib is emitted and every test/app links the static archive.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LFortran")
  set(RAKALI_BUILD_SHARED OFF)
  set(RAKALI_TESTS_LINK_SHARED OFF)
  message(
    STATUS
      "LFortran: forcing the static core (no shared .dylib — lfortran lacks the macOS -install_name flag)"
  )
endif()

# Shared library requires all code (including dependencies) built with -fPIC.
# The shared-core test/benchmark path (RAKALI_TESTS_LINK_SHARED) needs it too.
if(RAKALI_BUILD_SHARED
   OR RAKALI_ENABLE_PYTHON
   OR (RAKALI_TESTS_LINK_SHARED AND (RAKALI_ENABLE_TESTING
                                     OR RAKALI_BUILD_BENCHMARKS)))
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type (Debug/Release)" FORCE)
endif()

# GPU architecture — compiler-dependent default.  AMD/LLVM flang targets gfx90a
# (MI200/MI250 OpenMP-target offload, `--offload-arch`); every other compiler
# (notably NVHPC, which consumes it as `-gpu=ccNN`) defaults to cc70. Override
# explicitly with -DRAKALI_GPU_ARCH=<arch>.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  set(_rakali_default_gpu_arch "gfx90a")
else()
  set(_rakali_default_gpu_arch "cc70")
endif()
set(RAKALI_GPU_ARCH
    "${_rakali_default_gpu_arch}"
    CACHE STRING "GPU arch: AMD gfxNNN (LLVM flang) or NVIDIA ccNN (NVHPC)")

# Precision — applied to main_lib target after it's created (see below) Collect
# compile definitions in a list, apply them after add_library()
set(RAKALI_COMPILE_DEFS)
if(RAKALI_ENABLE_DOUBLE)
  list(APPEND RAKALI_COMPILE_DEFS RAKALI_DOUBLE_PRECISION)
endif()

# NVTX profiling annotations
if(RAKALI_ENABLE_NVTX)
  list(APPEND RAKALI_COMPILE_DEFS USE_NVTX)
  message(STATUS "NVTX profiling enabled")
endif()

# Compiler-specific flags (mutates CMAKE_Fortran_FLAGS* and RAKALI_COMPILE_DEFS)
include(compiler_flags)

if(RAKALI_ENABLE_COVERAGE AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(
    WARNING "RAKALI_ENABLE_COVERAGE=ON requires GNU/gfortran (current: "
            "${CMAKE_Fortran_COMPILER_ID}). Coverage flags NOT applied; the "
            "`coverage` target will not be created.")
endif()

# Print compile flags summary
message(STATUS "")
message(STATUS "=== Rakali Build Configuration ===")
message(
  STATUS
    "Compiler:       ${CMAKE_Fortran_COMPILER_ID} (${CMAKE_Fortran_COMPILER})")
message(STATUS "Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "Fortran flags:  ${CMAKE_Fortran_FLAGS}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug flags:    ${CMAKE_Fortran_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  message(STATUS "Release flags:  ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()
if(RAKALI_EXTRA_FORTRAN_FLAGS)
  message(STATUS "Extra flags:    ${RAKALI_EXTRA_FORTRAN_FLAGS}")
endif()
message(STATUS "GPU offload:    ${RAKALI_ENABLE_GPU}")
message(STATUS "Double prec:    ${RAKALI_ENABLE_DOUBLE}")
message(STATUS "NetCDF I/O:     ${RAKALI_ENABLE_NETCDF}")
message(STATUS "Tests:          ${RAKALI_ENABLE_TESTING}")
message(STATUS "Benchmarks:     ${RAKALI_BUILD_BENCHMARKS}")
message(STATUS "=================================")
message(STATUS "")

# Dependencies
if(RAKALI_ENABLE_TESTING)
  enable_testing()
  # Expose unit-test-only module symbols (e.g. the unsplit reference routines +
  # directionally-split building blocks in rki_continuity) to the test suite via
  # `#ifdef RAKALI_ENABLE_TESTING`, without widening the production API.
  list(APPEND RAKALI_COMPILE_DEFS RAKALI_ENABLE_TESTING)
  if(NOT TARGET test-drive::test-drive)
    find_package("test-drive" REQUIRED)
  endif()
endif()

if(NOT TARGET pic::pic)
  find_package("pic" REQUIRED)
endif()

# NetCDF Fortran for I/O Try CMake config first (modern installs); fall back to
# pkg-config (autotools installs)
if(RAKALI_ENABLE_NETCDF)
  find_package(netCDF-Fortran QUIET HINTS "$ENV{NETCDFF_DIR}"
               "$ENV{NETCDF_FORTRAN_DIR}")

  if(netCDF-Fortran_FOUND)
    message(
      STATUS "NetCDF-Fortran found via CMake config: ${netCDF-Fortran_DIR}")
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(NETCDFF REQUIRED IMPORTED_TARGET netcdf-fortran)
    add_library(netCDF::netcdff ALIAS PkgConfig::NETCDFF)
    message(
      STATUS "NetCDF-Fortran found via pkg-config: ${NETCDFF_LINK_LIBRARIES}")
  endif()
else()
  list(APPEND RAKALI_COMPILE_DEFS RAKALI_NO_NETCDF)
  message(
    STATUS "NetCDF disabled: skipping I/O subsystem (driver, output, gauges, "
           "restart, forcing, nesting, public API, main executable). "
           "Kernels, solver and benchmarks still build.")
endif()

# Benchmark-only API surface.  Kernels guarded by `#ifdef RAKALI_BENCHMARK_API`
# expose internal helpers (e.g. `ml_equation_of_state_impl`) only when the
# benchmark suite is being built, so production deployments don't get a wider
# public ABI than they need.
if(RAKALI_BUILD_BENCHMARKS)
  list(APPEND RAKALI_COMPILE_DEFS RAKALI_BENCHMARK_API)
endif()

# Per-column stack-workspace bound (NZ_STACK_MAX in rki_constants). Plumbed as a
# preprocessor define so it is build-configurable without editing source.
list(APPEND RAKALI_COMPILE_DEFS RKI_NZ_STACK_MAX=${RAKALI_NZ_STACK_MAX})

# Unstructured-multilayer stack bound (BPG_NZ_MAX in rki_ml_common_unstr). Same
# plumbing pattern: build-configurable without editing source.
list(APPEND RAKALI_COMPILE_DEFS RKI_BPG_NZ_MAX=${RAKALI_BPG_NZ_MAX})

# Flang/AMD offload advisory: the per-work-item device stack frame of the heavy
# per-column solvers (e.g. kappa-shear) scales with NZ_STACK_MAX, and on some
# targets (MI250X / gfx90a) the default can exceed the AMDGPU static stack-frame
# limit, making the device link fail with "stack frame size (...) exceeds
# limit".
if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang") # matches "Flang" and "LLVMFlang"
  message(
    STATUS
      "Flang detected: NZ_STACK_MAX=${RAKALI_NZ_STACK_MAX}. If a GPU-offload device "
      "link reports 'stack frame size exceeds limit', reconfigure with a smaller "
      "-DRAKALI_NZ_STACK_MAX=<n> (must stay >= nz + 1 for your vertical layer count)."
  )
endif()

# Core library — built via an OBJECT library so the same compiled objects back
# BOTH an always-static core AND, on request, a shared core, without a second
# compile.
#
# Why the split: tests and benchmarks have their own `do concurrent` kernels
# that call core `!$acc routine seq` device functions.  Relocatable device code
# (RDC) can only be device-linked WITHIN a single executable/archive — it cannot
# cross a shared-library boundary.  So tests/benches must link a STATIC core
# (device code device-links into each exe).  The app and the Python/Julia FFI
# have no device kernels of their own into the core, so they can link the SHARED
# core (.so) when RAKALI_BUILD_SHARED=ON.
set(core_objs ${main_lib}_objs)
add_library(${core_objs} OBJECT)
set_target_properties(${core_objs} PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions(${core_objs} PUBLIC ${RAKALI_COMPILE_DEFS})
target_include_directories(${core_objs} PUBLIC ${CMAKE_BINARY_DIR}/modules)

set(libraries_to_link)

# Source directories (populate the OBJECT library)
add_subdirectory(src)

# Module output directory — compilation happens in the OBJECT library
set_target_properties(${core_objs} PROPERTIES Fortran_MODULE_DIRECTORY
                                              ${CMAKE_BINARY_DIR}/modules)

# Link dependencies (PUBLIC on the objects ⇒ propagate to every core consumer)
list(APPEND libraries_to_link pic::pic)
if(RAKALI_ENABLE_NETCDF)
  list(APPEND libraries_to_link netCDF::netcdff)
endif()
target_link_libraries(${core_objs} PUBLIC ${libraries_to_link})
if(RAKALI_ENABLE_NVTX AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  target_compile_options(${core_objs} PUBLIC "-cudalib=nvtx")
  target_link_options(${core_objs} PUBLIC "-cudalib=nvtx")
endif()

# STATIC core — ALWAYS built.  This is what tests + benchmarks link, so their
# device kernels device-link within the executable.
add_library(${main_lib} STATIC)
target_link_libraries(${main_lib} PUBLIC ${core_objs})
set_target_properties(${main_lib} PROPERTIES POSITION_INDEPENDENT_CODE ON
                                             LINKER_LANGUAGE Fortran)

# SHARED core (.so) — built from the SAME object library as the static core, so
# no second compile.  Two independent consumers want it: * the app +
# Python/Julia FFI, when RAKALI_BUILD_SHARED is on (they cross no device
# boundary into the core); * the test suite + benchmarks, when
# RAKALI_TESTS_LINK_SHARED is on — every test/bench except the few with an
# on-device call into a core `!$acc routine seq` links this instead of the
# static archive, skipping the per-exe static pull + RDC device link for a much
# faster build.
if(RAKALI_BUILD_SHARED
   OR (RAKALI_TESTS_LINK_SHARED AND (RAKALI_ENABLE_TESTING
                                     OR RAKALI_BUILD_BENCHMARKS)))
  add_library(${main_lib}_shared SHARED)
  target_link_libraries(${main_lib}_shared PUBLIC ${core_objs})
  set_target_properties(
    ${main_lib}_shared
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE Fortran
               OUTPUT_NAME rakali_core)

  # NVHPC's `pgacclnk` two-stage device link picks up a PRE-EXISTING
  # librakali_core.so from the build tree and hands it to `ld` as an input,
  # which then refuses with "input file 'librakali_core.so' is the same as
  # output file" and exits 2 -- while leaving the stale .so in place, so a test
  # run immediately afterwards silently executes old code.  That has produced
  # both a phantom all-green and a phantom mass-failure in one session.
  #
  # The .so is nowhere on the link line (verified: inputs are objects, libpic
  # and the NetCDF shared libs), so this is internal to the NVHPC link wrapper
  # and not a CMake dependency cycle.  It is NOT the Makefiles generator either
  # -- it reproduces identically under Ninja.  Deleting the output before
  # linking is the reliable fix, and it makes the manual `rm librakali_core.so`
  # workaround unnecessary.
  add_custom_command(
    TARGET ${main_lib}_shared
    PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E rm -f $<TARGET_FILE:${main_lib}_shared>
    COMMENT "Removing stale librakali_core.so before relink (NVHPC pgacclnk)")
endif()

# Distributable core linked by the app + FFI: a SHARED .so when requested, else
# the static core.  (Tests pick their core per-target in tests/CMakeLists.txt.)
if(RAKALI_BUILD_SHARED)
  set(core_dist ${main_lib}_shared)
  message(
    STATUS
      "librakali_core: SHARED (.so) for app/FFI; STATIC core for tests/benchmarks"
  )
else()
  set(core_dist ${main_lib})
endif()

# Executable — depends on rki_driver, which is part of the NetCDF I/O subsystem.
# Without NetCDF the namelist-driven simulation has no way to write output, so
# we skip the executable entirely; the library + benchmarks are still useful for
# portability testing.
if(RAKALI_ENABLE_NETCDF)
  add_executable(${exe_name} ${PROJECT_SOURCE_DIR}/app/main.F90)
  target_include_directories(${exe_name} PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(${exe_name} PRIVATE ${core_dist})
endif()

# Namelist knob-reference generator (host-only; no NetCDF / GPU needed). Builds
# a default config, registers the full strict schema, and renders every knob to
# Markdown.  Regenerate docs/generated_nml_knobs.md with: ./rakali_nml_doc
# docs/generated_nml_knobs.md
add_executable(rakali_nml_doc ${PROJECT_SOURCE_DIR}/app/nml_doc.F90)
target_include_directories(rakali_nml_doc PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(rakali_nml_doc PRIVATE ${core_dist})

# Python/Julia C API shared library.  The API modules are compiled into the core
# objects as well (via src/CMakeLists.txt's API_SOURCES) so the static tests can
# link them; the shared lib just re-exports those symbols.
if(RAKALI_ENABLE_PYTHON)
  add_library(
    rakali_shared SHARED ${PROJECT_SOURCE_DIR}/src/api/rki_api.F90
                         ${PROJECT_SOURCE_DIR}/src/api/rki_api_unstr.F90)
  set_target_properties(
    rakali_shared
    PROPERTIES OUTPUT_NAME "rakali"
               Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
               POSITION_INDEPENDENT_CODE ON)
  target_include_directories(rakali_shared
                             PRIVATE "${CMAKE_BINARY_DIR}/modules")
  target_link_libraries(rakali_shared PRIVATE ${core_dist})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" AND RAKALI_ENABLE_GPU)
    # `nvfortran -shared` does not pull in libnvomp, but libaccdevice needs its
    # __hxSetThreadAffinityToNumaNode. An executable gets it through the
    # driver's -acc runtime set; a library dlopen'd by a plain Python process
    # would otherwise fail with "undefined symbol" at load time.
    target_link_libraries(rakali_shared PRIVATE -Wl,--no-as-needed nvomp)
  endif()
  message(STATUS "Building Python/Julia C API shared library: librakali.so")
endif()

if(RAKALI_ENABLE_TESTING)
  add_subdirectory(tests)
endif()
if(RAKALI_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# Coverage report target: build with -DRAKALI_ENABLE_COVERAGE=ON, then `cmake
# --build build --target coverage` runs the test suite, captures .gcda counters
# with lcov, filters out system + build-tree paths, and renders an HTML report
# under build/coverage_report/index.html.
if(RAKALI_ENABLE_COVERAGE AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  find_program(LCOV_EXECUTABLE lcov)
  find_program(GENHTML_EXECUTABLE genhtml)
  if(LCOV_EXECUTABLE AND GENHTML_EXECUTABLE)
    # `--ignore-errors mismatch,unused` were added in lcov 2.x; on lcov 1.x
    # those keys are unknown arguments and the whole capture step fails. Detect
    # the major version once and gate the flags accordingly.
    execute_process(
      COMMAND ${LCOV_EXECUTABLE} --version
      OUTPUT_VARIABLE _lcov_version_raw
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCH "([0-9]+)\\.([0-9]+)" _lcov_version_match
                 "${_lcov_version_raw}")
    set(_lcov_major "${CMAKE_MATCH_1}")
    if(_lcov_major GREATER_EQUAL 2)
      set(_lcov_capture_ignore --ignore-errors mismatch,gcov,source,unused)
      set(_lcov_remove_ignore --ignore-errors unused)
      message(STATUS "Coverage: lcov ${CMAKE_MATCH_0} (2.x flags)")
    else()
      set(_lcov_capture_ignore --ignore-errors gcov,source)
      set(_lcov_remove_ignore "")
      message(STATUS "Coverage: lcov ${CMAKE_MATCH_0} (1.x flags)")
    endif()

    add_custom_target(
      coverage
      COMMAND ${CMAKE_COMMAND} -E remove -f coverage.info coverage_filtered.info
      COMMAND ${LCOV_EXECUTABLE} --directory . --zerocounters
      COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
      COMMAND ${LCOV_EXECUTABLE} --directory . --capture --output-file
              coverage.info ${_lcov_capture_ignore}
      COMMAND
        ${LCOV_EXECUTABLE} --remove coverage.info '/usr/*' '/opt/*' '/Library/*'
        '*/build/*' '*/tests/*' --output-file coverage_filtered.info
        ${_lcov_remove_ignore}
      COMMAND ${GENHTML_EXECUTABLE} coverage_filtered.info --output-directory
              coverage_report --ignore-errors source
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
      COMMENT "Running tests and generating gcov/lcov coverage report")
    message(STATUS "Coverage target enabled: build target `coverage`")
  else()
    message(
      WARNING
        "RAKALI_ENABLE_COVERAGE=ON but lcov/genhtml not found on PATH "
        "(lcov=${LCOV_EXECUTABLE}, genhtml=${GENHTML_EXECUTABLE}). "
        "Install with `brew install lcov` (macOS) or `apt install lcov` (Linux)."
    )
  endif()
endif()
