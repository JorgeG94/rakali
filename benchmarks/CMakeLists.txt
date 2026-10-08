# Benchmark suite for tidal solver

# Benchmarks link the SHARED core by default (fast link), mirroring the test
# suite (RAKALI_TESTS_LINK_SHARED, on by default). A self-contained `do
# concurrent` — e.g. bench_single_gpu's array copy — links against the .so fine;
# benches that CALL a core `!$acc routine seq` from their own `do concurrent` (a
# device→device reference that can't cross a .so boundary) were dropped rather
# than kept on the static core. When the shared core wasn't built, `bench_core`
# falls back to static for everything.
if(TARGET ${main_lib}_shared)
  set(bench_core ${main_lib}_shared)
else()
  set(bench_core ${main_lib})
endif()

# Single-GPU kernel-level benchmark
add_executable(bench_single_gpu bench_single_gpu.F90)
target_include_directories(bench_single_gpu
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(bench_single_gpu PRIVATE ${bench_core})

# Single-GPU time-integrator benchmark (performance + convergence)
add_executable(benchmark_timestepping_single_gpu
               benchmark_timestepping_single_gpu.F90)
target_include_directories(benchmark_timestepping_single_gpu
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(benchmark_timestepping_single_gpu PRIVATE ${bench_core})

# Unstructured KNP solver benchmark (regular triangulated mesh)
add_executable(bench_unstr bench_unstr.F90)
target_include_directories(bench_unstr PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(bench_unstr PRIVATE ${bench_core})

# Structured-layout EOS kernel bench — drives the structured
# `ml_equation_of_state` against a (nx, ny, nz) grid.
add_executable(bench_eos_str bench_eos_str.F90)
target_include_directories(bench_eos_str PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(bench_eos_str PRIVATE ${bench_core})

# Full-physics multilayer bench — drives `solver_unstr_step` in multilayer mode
# on a regular triangulated mesh with a stratified IC. Exercises EOS + BPG +
# vmix + tracer advection + remap chain together.
add_executable(bench_multilayer bench_multilayer.F90)
target_include_directories(bench_multilayer
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(bench_multilayer PRIVATE ${bench_core})

# Structured sister of bench_multilayer — drives `solver_step` in multilayer
# mode on a Cartesian grid.  Bonus: prints a per-kernel `profiler_report` after
# the run.
add_executable(bench_multilayer_str bench_multilayer_str.F90)
target_include_directories(bench_multilayer_str
                           PRIVATE "${CMAKE_BINARY_DIR}/modules")
target_link_libraries(bench_multilayer_str PRIVATE ${bench_core})
