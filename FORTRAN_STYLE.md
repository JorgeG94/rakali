# Tidal Fortran Style Guide

This document defines the coding conventions for all Fortran code within the tidal project. External dependencies (pic, etc.) may follow their own conventions.

## Naming Conventions

### Files
- Module files: `rki_<name>.f90` (lowercase with underscores)
- Preprocessed files: `rki_<name>.F90` (uppercase extension)
- Test files: `test_rki_<name>.f90`

### Modules
- Module names match filenames: `rki_<name>`
- One module per file (except for submodules)

```fortran
! Good
module rki_physical_fragment

! Bad
module PhysicalFragment
module physical_fragment  ! missing rki_ prefix
```

### Derived Types
- All derived types use `_t` suffix
- Use lowercase with underscores

```fortran
! Good
type :: calculation_result_t
type :: physical_fragment_t
type :: energy_t

! Bad
type :: CalculationResult
type :: calculation_result   ! missing _t suffix
type :: TFragment
```

### Variables and Procedures
- Use lowercase with underscores (snake_case)
- Use descriptive names - avoid single letters except for loop indices or if it is a very local, easy to deduce variable

```fortran
! Good
integer :: n_atoms
real(dp) :: total_energy
subroutine calculate_gradient(fragment, result)

! Bad
integer :: nAtoms      ! camelCase
integer :: na          ! too short
real(dp) :: E          ! too short, unclear
```

### Constants
- Use UPPERCASE with underscores for parameters

```fortran
! Good
integer, parameter :: MAX_ITERATIONS = 100
real(dp), parameter :: BOHR_TO_ANGSTROM = 0.529177210903_dp

! Bad
integer, parameter :: maxIterations = 100
```

## Required Practices

### Use Statements
- **Always use `only` clause** - enforced by `fortitude` linter

```fortran
! Good
use pic_types, only: dp, int32
use rki_result_types, only: calculation_result_t, energy_t

! Bad - will fail linter
use pic_types
use rki_result_types
```

### Implicit None
- **Always include `implicit none`** in modules and programs

```fortran
module rki_example
   use pic_types, only: dp
   implicit none
   private
   ! ...
end module rki_example
```

### Intent Declarations
- **Always declare intent** for all procedure arguments

```fortran
! Good
subroutine compute_energy(fragment, result)
   type(physical_fragment_t), intent(in) :: fragment
   type(calculation_result_t), intent(out) :: result

! Bad
subroutine compute_energy(fragment, result)
   type(physical_fragment_t) :: fragment  ! missing intent
```

### Private by Default (linter enforced)
- Modules should be `private` by default
- Explicitly declare `public` entities

```fortran
module rki_example
   implicit none
   private

   public :: my_public_type_t
   public :: my_public_subroutine

   type :: my_public_type_t
      ! ...
   end type

   type :: internal_helper_t  ! stays private
      ! ...
   end type
end module
```

### Limit Procedure Arguments
- Public subroutines/functions should have **6 or fewer arguments**
- If more are needed, group related arguments into a derived type
- This improves readability, maintainability, and makes API changes easier

```fortran
! Bad - too many arguments
subroutine run_calculation(coords, elements, charge, multiplicity, &
                           method, basis, max_iter, tolerance, &
                           energy, gradient, error)
   real(dp), intent(in) :: coords(:,:)
   integer, intent(in) :: elements(:)
   integer, intent(in) :: charge, multiplicity
   character(*), intent(in) :: method, basis
   integer, intent(in) :: max_iter
   real(dp), intent(in) :: tolerance
   real(dp), intent(out) :: energy
   real(dp), intent(out) :: gradient(:,:)
   type(error_t), intent(out) :: error
   ! 11 arguments - hard to use and maintain

! Good - grouped into logical types
subroutine run_calculation(system, config, result, error)
   type(system_t), intent(in) :: system           ! coords, elements, charge, mult
   type(calc_config_t), intent(in) :: config      ! method, basis, max_iter, tol
   type(calculation_result_t), intent(out) :: result  ! energy, gradient
   type(error_t), intent(out), optional :: error
   ! 4 arguments - clean and extensible
```

- **Exception**: Simple utility functions (e.g., `to_bohr(value)`) can have minimal arguments
- Internal/private procedures have more flexibility but should still aim for clarity

## Forbidden Practices

### No GOTO Statements
- `goto` is forbidden in this project under any circumstance.  Use
  structured control flow: `do`, `if`, `select case`, `exit`,
  `cycle`, named `block ... end block`.
- Linter (fortitude) catches `goto` and `continue` labels.
- **Do not** reintroduce `goto` even for early-return / cleanup
  cascades.  Canonical replacements are below.

```fortran
! Bad — bail with goto:
   if (error) goto 999
   ! ...
999 continue
   print *, "Error occurred"
   return

! Good — early return:
   if (error) then
      print *, "Error occurred"
      return
   end if
```

**Replacement for `testdrive`-style tests** (the most common goto
pattern this project saw): wrap the body in a named `block` and
`exit` to skip remaining checks while still hitting the cleanup
after the block.

```fortran
! Bad — labeled cleanup target:
   call check(error, c1, "m1")
   if (allocated(error)) goto 99
   call check(error, c2, "m2")
   if (allocated(error)) goto 99
   ! ...
99 call cor%destroy()
   call ms%destroy()
   end subroutine

! Good — named block, structured exit, single cleanup:
   checks: block
      call check(error, c1, "m1")
      if (allocated(error)) exit checks
      call check(error, c2, "m2")
      if (allocated(error)) exit checks
      ! ...
   end block checks
   call cor%destroy()
   call ms%destroy()
   end subroutine
```

For **phased cleanup cascades** (where each label rolls back a
different set of resources — the F77 idiom of `goto 99` / `goto 88`
/ `goto 77` for nested acquire/release pairs), nest blocks around
each acquire/release pair:

```fortran
   call nc_open_read(FN, ncid)
   ncid_scope: block
      call check(error, ...)
      if (allocated(error)) exit ncid_scope
      allocate (read_buf(...))
      buf_scope: block
         call check(error, ...)
         if (allocated(error)) exit buf_scope
         ! ... more checks ...
      end block buf_scope
      deallocate (read_buf)
   end block ncid_scope
   call nc_close(ncid)
```

Each block's cleanup (`deallocate`, `nc_close`, etc.) sits
immediately after the matching `end block`; the structure mirrors
RAII / stack-unwinding.

### No Arithmetic IF
- Use `if-then-else` or `select case`

```fortran
! Bad (arithmetic IF - Fortran 77)
   if (x) 10, 20, 30

! Good
   if (x < 0) then
      ! negative case
   else if (x == 0) then
      ! zero case
   else
      ! positive case
   end if
```

### No Implicit SAVE
- Avoid module-level variables that retain state
- If state is needed, use explicit `save` attribute and document why

```fortran
! Bad - implicit save behavior
module rki_bad_example
   integer :: counter = 0  ! implicitly saved, retains value between calls
end module

! Good - use derived types to manage state
module rki_good_example
   type :: counter_t
      integer :: value = 0
   end type
end module
```

### No COMMON Blocks
- Use module variables or derived types instead

```fortran
! Bad
common /shared_data/ x, y, z

! Good
module rki_shared_data
   use pic_types, only: dp
   implicit none
   type :: shared_data_t
      real(dp) :: x, y, z
   end type
end module
```

### No Naked Print Statements
- **NEVER** use `print *` or `write(*,*)` for output
- **ALWAYS** use `pic_logger` (`global_logger`)
- This is enforced - naked prints will be rejected in code review

```fortran
! Bad - forbidden
print *, "Starting calculation"
write(*,*) "Energy:", energy
write(6,*) "Done"

! Good - use logger
use pic_logger, only: logger => global_logger
call logger%info("Starting calculation")
call logger%info("Energy: " // to_char(energy))
call logger%info("Done")
```

### No Emojis in Fortran Code
- No emojis in Fortran source files (`.f90`, `.F90`)
- This includes comments, strings, and documentation
- Keep output professional and portable

```fortran
! Bad
call logger%info("Calculation complete! 🎉")
!! 🚀 Fast implementation of MBE

! Good
call logger%info("Calculation complete")
!! Fast implementation of MBE
```

- Python scripts (validation, tooling) may use emojis if desired

### No EQUIVALENCE
- Use proper type conversions or `transfer()` if absolutely necessary

### No Fixed-Form Source
- All code must be free-form (`.f90` / `.F90`)
- No `.f` or `.F` files

### No Assumed-Size Arrays
- Use assumed-shape arrays with explicit interface

```fortran
! Bad
subroutine process_array(arr, n)
   real(dp) :: arr(*)  ! assumed-size
   integer :: n

! Good
subroutine process_array(arr)
   real(dp), intent(in) :: arr(:)  ! assumed-shape
```

### No External Statements for Internal Procedures
- Use `contains` for internal procedures
- Use explicit interfaces via modules

```fortran
! Bad
external :: my_function

! Good - use modules
use rki_my_module, only: my_function
```

## Recommended Practices

### Use pic Library Utilities
- **Always prefer pic functionality** over implementing your own
- Ensures consistency, reduces code duplication, and leverages tested implementations

| pic Module | Functionality | Use Instead Of |
|------------|---------------|----------------|
| `pic_types` | Kind parameters (`dp`, `int32`, etc.) | Literal kinds (`real(8)`) |
| `pic_logger` | Logging (`global_logger`) | `print *` / `write(*,*)` (**forbidden**) |
| `pic_timer` | Performance timing | Manual `cpu_time` calls |
| `pic_strings` | String utilities | Custom string manipulation |
| `pic_sorting` | Sorting algorithms | Hand-rolled sorts |
| `pic_math` | Math utilities | Reimplementing common math |
| `pic_test_helpers` | Test utilities (`is_equal`) | Custom comparison functions |

```fortran
! Good - use pic utilities
use pic_logger, only: logger => global_logger
use pic_timer, only: timer_t
use pic_sorting, only: argsort

call logger%info("Starting calculation")
call logger%debug("Processing fragment " // to_char(i))
call logger%warning("Large system detected, may be slow")
call logger%error("Invalid input: negative charge")

type(timer_t) :: timer
call timer%start()
! ... work ...
call timer%stop()
call logger%info("Elapsed: " // timer%elapsed_string())

indices = argsort(energies)  ! sorted indices

! Bad - NEVER use naked print statements
print *, "Starting calculation"  ! no log levels
write(*,*) "Debug info"          ! can't be silenced, clutters output
```

**Why logger over print?**
- Log levels (debug/info/warning/error) - control verbosity
- Consistent formatting across the codebase
- Can be redirected to files or silenced entirely

- Check pic documentation before implementing common functionality
- If pic is missing something useful, consider contributing it upstream

### No Naked BLAS Calls

#### pic-blas
- **Always use pic-blas** for BLAS/LAPACK - never call BLAS/LAPACK directly
- Provides cleaner interface and handles different BLAS implementations

```fortran
! Good - use pic-blas wrappers
use pic_blas, only: pic_gemm, pic_dot

call pic_gemm(A, B, C)
result = pic_dot(x, y)

! Bad - naked BLAS calls
call dgemm('N', 'N', m, n, k, 1.0d0, A, lda, B, ldb, 0.0d0, C, ldc)
result = ddot(n, x, 1, y, 1)
```

- Exceptions: Only if pic-blas lacks needed functionality (then consider contributing)

### Kind Parameters
- Use `pic_types` for portable kind definitions
- Never use literal kind numbers

```fortran
! Good
use pic_types, only: dp, int32
real(dp) :: energy
integer(int32) :: count

! Bad
real(8) :: energy      ! non-portable
real*8 :: energy       ! obsolete syntax
double precision :: e  ! obsolete
```

### Array Operations
- Prefer whole-array operations over explicit loops when clear

```fortran
! Good
gradient = 0.0_dp
total = sum(energies)

! Also fine for complex operations
do i = 1, n_atoms
   gradient(:, i) = gradient(:, i) + contribution(:, i)
end do
```

### Block Constructs for Limiting Scope
- Use `block` to limit variable scope and improve readability
- Useful for temporary variables needed only in a small section
- Helps prevent accidental reuse of variables

```fortran
! Good - temporary variables scoped to where they're needed
subroutine process_fragments(fragments, total_energy)
   type(fragment_t), intent(in) :: fragments(:)
   real(dp), intent(out) :: total_energy

   integer :: i

   total_energy = 0.0_dp
   do i = 1, size(fragments)
      block
         real(dp) :: frag_energy
         real(dp) :: correction

         call compute_energy(fragments(i), frag_energy)
         call compute_correction(fragments(i), correction)
         total_energy = total_energy + frag_energy + correction
      end block
   end do
end subroutine

! Bad - temporaries pollute the entire procedure scope
subroutine process_fragments(fragments, total_energy)
   type(fragment_t), intent(in) :: fragments(:)
   real(dp), intent(out) :: total_energy

   integer :: i
   real(dp) :: frag_energy    ! visible everywhere, might be misused later
   real(dp) :: correction     ! visible everywhere

   total_energy = 0.0_dp
   do i = 1, size(fragments)
      call compute_energy(fragments(i), frag_energy)
      call compute_correction(fragments(i), correction)
      total_energy = total_energy + frag_energy + correction
   end do
end subroutine
```

### Associate Construct for Readability
- Use `associate` to create short aliases for long expressions
- Improves readability without runtime cost

```fortran
! Good
associate(coords => fragment%coordinates, &
          n => fragment%n_atoms)
   do i = 1, n
      distance = norm2(coords(:, i) - origin)
   end do
end associate

! Bad - repeated long expressions
do i = 1, fragment%n_atoms
   distance = norm2(fragment%coordinates(:, i) - origin)
end do
```

### Error Handling
- Use the `error_t` type from `rki_error` module
- Check and propagate errors

```fortran
use rki_error, only: error_t, create_error, has_error

subroutine my_subroutine(input, output, error)
   type(input_t), intent(in) :: input
   type(output_t), intent(out) :: output
   type(error_t), intent(out), optional :: error

   if (invalid_input) then
      if (present(error)) then
         call create_error(error, "Invalid input: reason")
      end if
      return
   end if
end subroutine
```

### Documentation
- Use `!!` for FORD documentation comments
- Document public interfaces

```fortran
type :: calculation_result_t
   !! Container for quantum chemistry calculation results
   !!
   !! Holds energy, gradient, Hessian, and associated metadata
   !! from a single-point or property calculation.

   type(energy_t) :: energy
      !! Total and component energies (SCF, correlation, etc.)
   real(dp), allocatable :: gradient(:,:)
      !! Nuclear gradient (3, n_atoms) in Hartree/Bohr
end type
```

### Citing algorithms (cite the paper, not another codebase)

When a kernel implements a published method, cite the **paper** in the
docstring — not another model's source file. Write `!! Wright (1997) nonlinear
EOS`, not `!! ported from <other-model>/MOM_EOS_Wright.F90`. Source-file
references rot (the other project renames/moves things), read as "we copied
this," and bury the actual provenance — the paper — that a reader needs.

It is fine to **read a non-GPL reference implementation** (SCHISM, ROMS,
…) to understand an algorithm — then implement it ourselves from the maths and
cite the paper *they* cite. **GPL code is off-limits**: don't read, port, or
line-by-line reference it (e.g. GOTM, GPLv2). Rakali is not GPL and must stay
that way, so any closure with only a GPL reference implementation is written
clean-room from the published equations. The closure *interface* (what inputs a
scheme consumes / returns) is a standard functional shape and fine to mirror;
only the GPL *source code* is protected expression.

### Memory Management
- Provide `destroy` procedures for types with allocatable components
- Clean up allocatable arrays when no longer needed

```fortran
type :: my_type_t
   real(dp), allocatable :: data(:)
contains
   procedure :: destroy => my_type_destroy
end type

subroutine my_type_destroy(this)
   class(my_type_t), intent(inout) :: this
   if (allocated(this%data)) deallocate(this%data)
end subroutine
```

### Prefer Allocatable Over Pointer
- Use `allocatable` instead of `pointer` when possible
- Allocatable arrays are automatically deallocated when out of scope
- Less risk of memory leaks and dangling pointers
- Compiler can optimize allocatable better

```fortran
! Good - automatic cleanup, no leak risk
type :: data_container_t
   real(dp), allocatable :: values(:)
   real(dp), allocatable :: matrix(:,:)
end type

! Bad - manual cleanup required, leak risk
type :: data_container_t
   real(dp), pointer :: values(:) => null()
   real(dp), pointer :: matrix(:,:) => null()
end type
```

- **When to use pointer**: Only when you need pointer semantics (aliasing, linked structures, polymorphic returns)

### Pure and Elemental Procedures
- **Default to `pure`.** When you add a new function or subroutine, make it
  `pure` unless it genuinely needs a side effect (I/O, logging, mutating
  module state, an `allocate` + `!$acc enter data`, or calling a non-pure
  callee). The bar is "can this be pure?", not "should this be pure?" — reach
  for the side effect only when the procedure's job requires it. Pure is the
  norm in this codebase (the kernels, EOS, pressure, tracer, and ALE leaves
  are pure); a new non-pure leaf helper should be the exception you can justify.
- Use `elemental` for scalar operations that should work on arrays.
- Enables compiler optimizations, documents intent, and keeps the
  `do concurrent` GPU path clean (loop bodies must call pure code).
- If a procedure can't be pure only because it lazily allocates a scratch
  workspace, that's a signal to lift the workspace onto a state slot
  (allocate once at init, free at destroy) — see the persistent-workspace
  pattern under *Do Concurrent + OpenACC* — which then unblocks `pure`.

```fortran
! Good - pure function, no side effects
pure function kinetic_energy(mass, velocity) result(energy)
   real(dp), intent(in) :: mass, velocity
   real(dp) :: energy
   energy = 0.5_dp * mass * velocity**2
end function

! Good - elemental, works on scalars and arrays
elemental function to_bohr(angstrom) result(bohr)
   real(dp), intent(in) :: angstrom
   real(dp) :: bohr
   bohr = angstrom * ANGSTROM_TO_BOHR
end function

! Usage: works on scalar or array
r_bohr = to_bohr(1.5_dp)           ! scalar
coords_bohr = to_bohr(coords_ang)   ! array
```

### Avoid Deep Nesting
- Maximum 3-4 levels of indentation
- Use early returns, `cycle`, `exit` to reduce nesting
- Extract deeply nested code to separate subroutines

```fortran
! Bad - deeply nested, hard to follow
do i = 1, n
   if (condition1) then
      if (condition2) then
         do j = 1, m
            if (condition3) then
               ! actual work buried here
            end if
         end do
      end if
   end if
end do

! Good - early cycle, flat structure
do i = 1, n
   if (.not. condition1) cycle
   if (.not. condition2) cycle

   do j = 1, m
      if (.not. condition3) cycle
      ! actual work at reasonable depth
   end do
end do

! Also good - extract to subroutine
do i = 1, n
   call process_item(items(i), result)
end do
```

### No Magic Numbers
- Use named constants for any non-obvious literal values
- Makes code self-documenting and easier to maintain

```fortran
! Bad - what do these numbers mean?
if (n_atoms > 50) then
   cutoff = 4.5
end if
tolerance = 1.0e-8

! Good - self-documenting
integer, parameter :: LARGE_SYSTEM_THRESHOLD = 50
real(dp), parameter :: DEFAULT_DISTANCE_CUTOFF = 4.5_dp  ! Angstrom
real(dp), parameter :: SCF_CONVERGENCE_TOLERANCE = 1.0e-8_dp

if (n_atoms > LARGE_SYSTEM_THRESHOLD) then
   cutoff = DEFAULT_DISTANCE_CUTOFF
end if
tolerance = SCF_CONVERGENCE_TOLERANCE
```

### Allocatable Character Strings
- Use `character(len=:), allocatable` for dynamic strings
- Avoid fixed-length strings that waste memory or truncate

```fortran
! Good - allocatable string
character(len=:), allocatable :: method_name
character(len=:), allocatable :: error_message

method_name = "GFN2-xTB"  ! automatically sized
error_message = "Error in " // trim(filename) // ": " // trim(reason)

! Bad - fixed length, may truncate or waste space
character(len=32) :: method_name   ! what if name is longer?
character(len=256) :: error_message  ! wastes memory for short messages
```

### Do Concurrent + OpenACC (Rakali's GPU Path)

Rakali's GPU offloading story is **`do concurrent` + OpenACC**, compiled with NVHPC `-stdpar=gpu -acc=gpu`. This gives us:
- Data-parallel loops written in portable Fortran (`do concurrent`) — same source works on CPU, multicore, and GPU.
- OpenACC directives (`!$acc enter data`, `!$acc parallel loop reduction(...)`) for the parts `do concurrent` can't express: explicit device data management and reductions on GPU.

NVHPC has matured significantly — many older folklore cautions no longer apply. This section captures what's real in the current (26.x) compiler.

#### What works well
- **`!$acc enter/exit data` on allocatable arrays.** Rock-solid since ~25.x. Put them at the start/end of a long-lived scope (e.g., `solver_enter_data`), map state arrays once, keep them resident, only transfer back for diagnostic output or user callbacks. This is the foundation of Rakali's performance.
- **`do concurrent` on large k,j,i loops with `local(...)` clauses.** Compiles cleanly and coalesces well. Put the contiguous (leading) array index LAST in the index list — `do concurrent(k,j,i)` for `arr(i,j,k)` (see the index-order quirk below for why the spread is 1× vs 16–23× depending on backend).
- **`!$acc update self/device` for in-place host↔device sync.** No need for full exit/enter round trips on mid-simulation callbacks.
- **Reductions** via `!$acc parallel loop collapse(2) reduction(max:...)`, written directly with no `#ifdef` wrapping. On non-OpenACC compilers the directive is an inert comment and the loop runs sequentially.

```fortran
! Good - independent iterations, local scratch via local() clause
do concurrent (j=ng+1:ny-ng, i=ng+1:nx-ng) &
   local(tmp, flux_x, flux_y)
   tmp = 0.5_wp * (u(i,j) + u(i+1,j))
   flux_x = ...
   flux_y = ...
   out(i,j) = flux_x + flux_y
end do
```

#### Required pattern for every `*_enter_data` / `*_exit_data` routine

Every derived type whose components are referenced from `do concurrent` (or
`!$acc parallel loop`) bodies MUST attach the parent struct with
`!$acc enter data copyin(this)` **before** any component attaches.
Reverse the order on `exit_data` (components first, parent last).

```fortran
subroutine foo_enter_data(this)
   class(foo_t), intent(inout) :: this
   !$acc enter data copyin(this)              ! PARENT FIRST
   !$acc enter data copyin(this%arr1, this%arr2)
   call this%sub_struct%enter_data()          ! recursive: same rule
end subroutine

subroutine foo_exit_data(this)
   class(foo_t), intent(inout) :: this
   call this%sub_struct%exit_data()
   !$acc exit data delete(this%arr1, this%arr2)
   !$acc exit data delete(this)               ! PARENT LAST
end subroutine
```

**Why:** When a DC body references `foo%arr1(i, j)`, NVHPC needs the parent
struct's descriptor on the device to resolve the component pointer.
Without `copyin(this)` the compiler falls back to UVM / managed memory and
page-faults the descriptor on every kernel launch.  `nsys` profiles show
this as "memcpy" eating 90%+ of kernel time.  Discovered 2026-05-20 in a
solver fast loop — adding the missing `copyin(this)` lines to ~13
`enter_data` routines (and a missing sub-state `enter_data()` call) turned
a 60 s smoke run into a 5 s smoke run (12× speedup).  The lesson: audit
*every* `enter_data` routine for the parent `copyin(this)`; a path added
later can miss it everywhere.

The coastal `solver_enter_data` in `src/core/coastal/solver/rki_solver.F90`
is the canonical reference implementation.

#### What remains genuinely painful
- **Nested derived types with allocatable-array components mapped to the device.** Parent-before-components ordering in `!$acc enter data` is mandatory (map the outer struct first, then each array component separately). Reverse the order on `exit data`. Omitting a component — even one no kernel reads from — produces cryptic runtime failures or silent corruption (an unattached allocatable leaves an invalid pointer in the parent's device descriptor; the first kernel to dereference *anything* on it trips). See the parent-before-components pattern above.
- **Indirected access through a derived-type-array component from inside a kernel.** Even with the parent + per-element attach in place, writing `state%multilayer%tracers(it)%hTr(i,j,k)` directly inside an `!$acc parallel loop` or `do concurrent` body hits `CUDA_ERROR_ILLEGAL_ADDRESS` on NVHPC. Two levels of indirection — the array-of-derived-types descriptor *plus* the allocatable component descriptor — aren't reachable through the device's dereferencing chain. **Fix:** outer-shim + inner-impl. The outer subroutine takes the parent struct, dereferences `state%multilayer%tracers(idx)%hTr` on the *host*, and passes flat 3D arrays as actual arguments to an inner subroutine whose dummies are `real(wp), intent(...) :: hTr_layer(grid%nx_total, grid%ny_total, nz)`. Fortran argument-passing collapses the indirection at the call site; the inner kernel sees a plain array.

  Every `ml_*_tracer` kernel in `kernels/structured/multilayer/rki_ml_tracers.F90` follows this pattern. `compute_layer_sums` / `compute_layer_sums_impl` in `rki_init.F90` is the canonical reduction example. On the unstructured side: `ml_equation_of_state_unstr` → `_impl`, `ml_distribute_with_remap_unstr` → `_impl`, `halo_exchange_layers_unstr` → `_impl`.

  ```fortran
  ! BAD - kernel chases tracer-registry indirection on device, fails at runtime
  !$acc parallel loop collapse(3) reduction(+:acc)
  do k = 1, nz
     do j = jlo, jhi
        do i = ilo, ihi
           acc = acc + state%multilayer%tracers(state%multilayer%idx_salinity)%hTr(i, j, k)
        end do
     end do
  end do

  ! GOOD - outer shim hands flat arrays to the device kernel
  call sum_tracer_impl(state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                       state%grid, state%multilayer%nz_ml, acc)

  subroutine sum_tracer_impl(hTr_layer, grid, nz, acc)
     type(hgrid_t), intent(in) :: grid
     integer, intent(in) :: nz
     real(wp), intent(in) :: hTr_layer(grid%nx_total, grid%ny_total, nz)
     real(wp), intent(out) :: acc
     ! ...
     !$acc parallel loop collapse(3) reduction(+:acc)
     do k = 1, nz
        do j = jlo, jhi
           do i = ilo, ihi
              acc = acc + hTr_layer(i, j, k)
           end do
        end do
     end do
  end subroutine
  ```
- **Cross-module pure helpers called from `do concurrent` on GPU.** Same-module pure calls work fine — `compute_flux` calls `flux_cell` from inside `do concurrent` and lowers correctly on NVHPC. The genuine restriction is *cross-module*: even with `!$acc routine seq` on the callee, NVHPC's stdpar GPU path does not currently lower a cross-compilation-unit pure subroutine call from inside `do concurrent` (`NVFORTRAN-S-1061` at compile time, or wrong-result silent failure depending on the NVHPC release). Two workarounds, in order of preference:
  1. Move the callee into the same module as the caller (often the cleanest fix when the helper is small and only used in one place).
  2. Drop the loop to `!$acc parallel loop` with explicit `present(...)` and `private(...)` clauses — that lowering path *does* handle cross-module routines marked `!$acc routine seq`. Canonical example: `ml_distribute_with_remap` in `kernels/structured/multilayer/rki_ml_dynamics.F90` calling `vcoord_target_dz_column` from `rki_vcoord`.

- **Stack-size limits on per-thread automatic arrays inside `do concurrent`.** NVHPC's default GPU thread stack is small (tens of KB). Fixed-size workspace arrays like `p_int(BPG_NZ_MAX)` in the BPG kernels are fine at modest sizes; unbounded local arrays can fail silently. See `BPG_NZ_MAX = 20` in `rki_constants` and the validation check in `validate_config`.
- **OpenMP target offload for `do concurrent`.** Not our path. We compile with `-stdpar=gpu -acc=gpu`, not with `-mp=gpu`. The OpenMP-target + `do concurrent` combination has historically been rougher on NVHPC; OpenACC + stdpar is the well-trodden path and what all our performance work is built around. **Don't mix OpenMP target directives with `do concurrent` in new code.**

- **Cross-iteration read/write on the same array (Jacobi vs Gauss-Seidel).** A `do concurrent` that reads a *neighbour* index and writes its *own* index on the **same** array gives different answers per backend: GPU lowers to one launch (threads read the pre-update image → **Jacobi**); CPU / `-stdpar=multicore` runs sequentially (later iterations see earlier writes → **Gauss-Seidel**) — a small systematic CPU↔GPU divergence (~1e-6/step, accumulating) the bit-compare CI catches but GPU-only testing misses. **Fix:** snapshot the read-source into a scratch buffer on entry, read from scratch, write to the live array. Examples: `ml_advect_tracer_unstr` (`tracer_advect_scratch`) and the Smagorinsky momentum snapshot (`horvisc_hu_pre`).

#### Other NVHPC / gfortran codegen quirks

- **`do concurrent` index order — put the contiguous (leading) array index
  LAST.** For column-major `arr(i,j,k)` write `do concurrent(k,j,i)` (i
  rightmost). The penalty for getting it wrong is backend-dependent and large:
  - **NVHPC `-stdpar=gpu`** (production): order-immune. stdpar auto-collapses
    the nest and remaps the contiguous index to the warp lane itself (`-Minfo`
    prints `auto-collapsed-innermost`). All permutations measured identical
    (~780 GB/s AXPY; ties on a 4-array C-grid stencil).
  - **OpenMP-target (`-mp=gpu`, the `auto/openmp` backend) + CPU cache**:
    load-bearing. `!$omp … collapse(3)` is literal — the innermost loop maps
    to the fastest thread, so a non-contiguous innermost index strides global
    memory. Measured on the same stencil: i-innermost 0.11 s vs k-innermost
    **1.8–2.6 s (16–23× slower)**, bit-identical results.

  So `(k,j,i)` is free on stdpar and essential everywhere else — never place
  the contiguous index anywhere but last. MRE: `local_archive/omp_stencil_mre.F90`.
- **`-gpu=managed` is banned.** Host-allocate a scratch array, write it from a
  `do concurrent`, then read it back with the `sum()` intrinsic → it silently
  returns 0 under managed memory. Use `!$acc parallel loop reduction(...)` for
  reductions over fresh scratch, and keep state explicitly mapped.
- **Don't split a `present(optional)` test into two `do concurrent` loops.**
  NVHPC compiled the no-optional branch of a split remap ~25× slower than the
  with-optional one. Gate the optional *inside* one loop (as the vdiff kernels
  do), not `if (present(x)) then <DC> else <DC>`.
- **Keep a loop-invariant `select case` inside one `do concurrent`.** Splitting
  into one DC per case doubles the GPU launches (~1.8× slower); a uniform branch
  every iteration takes is essentially free. MRE-measured.
- **Don't loop-swap or hoist invariants on bandwidth-bound kernels.** NVHPC
  already does it; the extra scratch arrays just add register spills (a MUSCL
  rewrite came out ~20% *slower*). Micro-opt only when a profile says the kernel
  is compute-bound, not bandwidth-bound.
- **CUDA-graph capture.** A bare `do concurrent` can't be CUDA-graph-captured
  (stdpar inserts a sync per loop); wrapping the loop body in
  `!$acc kernels async(q)` (body unchanged) makes it capturable (~5× on
  launch-bound sequences). MRE: `tools/mre_dc_cuda_graph.F90`.
- **gfortran `do concurrent local(x)` reassignment (gfortran 15.1).** gfortran
  corrupts a `local(x)` scalar when it's reassigned across multiple if/else
  branches in the body (NVHPC is unaffected). Use one local per write site
  (`h_face_E/W/N/S`, not a single reused `h_face`).
- **Assumed-shape dummy arguments inside `do concurrent` kernels.** NVHPC
  walks the Fortran array descriptor on every kernel launch when a dummy is
  declared `arr(:,:,:)`. The vmix incident produced 1.4 M per-launch
  descriptor-walk memcpys over 52 s of runtime until the hot kernels were
  converted to explicit-shape. **Rule: any array dummy referenced inside a
  `do concurrent` body must use explicit-shape** — pass integer dimensions
  and declare `arr(nx, ny, nz)`. Proven pattern: `kappa_shear_merge_into_kv_kt`
  and `epbl_merge_into_kv_kt`.

  **Waiver** for cadence-bounded paths (diag fills that fire once per output
  frame, outer-shim impls that receive registry-dereferenced allocatables of
  varying size, init-only routines): annotate the declaration line or the
  line immediately above it with `! assumed-shape-ok: <reason>`. The reason
  should identify why conversion isn't applicable (e.g. face-sized array with
  mismatched first dim, TBP polymorphic dispatch, tracer-registry outer-shim).

  Enforced diff-aware by `tools/dc_assumed_shape_lint.py` (pre-commit hook
  `dc-assumed-shape`): only newly-added offenders in `src/`/`app/` fail.

- **`map(alloc:/delete:)` workspace arrays can't be unit-tested by host
  pre-fill** — the host buffer isn't the device buffer. Assert via
  run-twice-and-compare instead.

- **A GPU kernel behind a type-bound procedure must not be reached through
  *overridable* (polymorphic) dispatch.** If a `do concurrent` kernel — or a
  routine that launches one / maps state to the device — is invoked as a TBP
  on a `class(...)` actual, the call goes through a runtime vtable. The
  compiler then can't devirtualize or inline it, and NVHPC's stdpar /
  OpenMP-target ABI does not inline polymorphic dispatch into device code: at
  best a missed inline, at worst the **`class(...)` box gets mapped to the
  device and the run crashes** (the AMD libomptarget overlap crash; same root
  as the `!$omp declare target`-doesn't-inline footgun). **Rule: any TBP that
  contains or launches a GPU kernel is `non_overridable`, OR the polymorphic
  TBP is a thin wrapper that immediately delegates (via `select type`) to a
  non-polymorphic `type(...)` `*_impl` that does the device work.** Keep the
  *device-touching* code on the concrete type; let only host orchestration be
  polymorphic.

  This is forward-looking, not a retrofit: a repo-wide sweep adding
  `non_overridable` to all 174 host-side TBPs benchmarked as a measured null
  (−0.13%, inside run noise; byte-identical binary) on the GPU seamount bench
  precisely *because* no kernel currently sits behind a `class(...)` TBP — the
  `_impl` delegation already keeps them off it. The rule exists so a future
  kernel doesn't quietly regress (or crash on AMD) by being bound directly to a
  polymorphic type.

#### Restrictions that apply to any `do concurrent`
- No `exit`, `cycle`, `return`, or `goto` inside the loop body (Fortran-language rule, independent of backend).
- Side effects (I/O, logging, pure-function violations) are undefined behavior.
- Reductions need a separate OpenACC directive; `do concurrent` alone can't express them portably.

#### Rule of thumb
- **Default to `do concurrent` + `local(...)` for data-parallel kernels.** Portable across CPU (serial or multicore via `-stdpar=multicore`) and GPU.
- **Drop to `!$acc parallel loop` when you need reductions or explicit collapse.** Write the directive directly, no `#ifdef` needed — non-OpenACC compilers treat it as a comment.
- **Never rely on `do concurrent` for reductions or ordered operations.** It's a parallelism hint, not a primitive.
  - *Yes, F2023 added a `reduce(+:s)` locality spec, so `do concurrent (…) reduce(+:s)` is legal — we still don't use it.* The `!$acc parallel loop reduction(...)` form degrades gracefully: it's an inert comment to any compiler not built with `-acc`, so the identical source runs as a correct **serial** reduction on gfortran/ifx and a GPU reduction on NVHPC — one source, no `#ifdef`. `do concurrent reduce` instead needs *uniform* F2023 `reduce` support across NVHPC + gfortran + ifx (new and uneven; gfortran's `do concurrent` doesn't GPU-offload at all). And the failure mode is silent: a plain `do concurrent` with `s = s + …` and **no** `reduce` clause is a data race on the shared scalar (the `sum()`-returns-0 footgun above). So reductions stay on `!$acc parallel loop reduction(...)`.

### Submodules for Large Modules
- Use submodules to separate interface from implementation
- Reduces recompilation when only implementation changes
- Useful for large modules with many procedures

```fortran
! In rki_heavy_module.f90 - interface only
module rki_heavy_module
   implicit none
   private
   public :: compute_expensive_thing

   interface
      module subroutine compute_expensive_thing(input, output)
         real(dp), intent(in) :: input(:)
         real(dp), intent(out) :: output(:)
      end subroutine
   end interface
end module

! In rki_heavy_module_impl.f90 - implementation
submodule (rki_heavy_module) rki_heavy_module_impl
contains
   module subroutine compute_expensive_thing(input, output)
      real(dp), intent(in) :: input(:)
      real(dp), intent(out) :: output(:)
      ! Long implementation here
      ! Changes here don't trigger recompilation of modules that use rki_heavy_module
   end subroutine
end submodule
```

## Units

- **Internal units**: Bohr (length), Hartree (energy)
- **Conversion**: Use `to_bohr()` / `to_angstrom()` from `rki_physical_fragment`
- **Document units** in comments when not obvious

```fortran
real(dp) :: bond_length      ! in Bohr
real(dp) :: energy           ! in Hartree
real(dp) :: frequency_cm1    ! in cm^-1 (document non-atomic units)
```

## File Structure Template

```fortran
!! Brief module description
module rki_example
   !! Extended module documentation
   !!
   !! More details about the module purpose and usage.
   use pic_types, only: dp, int32
   use rki_other_module, only: needed_type_t
   implicit none
   private

   public :: example_type_t
   public :: example_subroutine

   type :: example_type_t
      !! Type documentation
      integer :: n_items
         !! Number of items
      real(dp), allocatable :: values(:)
         !! Array of values in Hartree
   contains
      procedure :: compute => example_compute
      procedure :: destroy => example_destroy
   end type example_type_t

contains

   subroutine example_compute(this, input, output)
      !! Subroutine documentation
      class(example_type_t), intent(inout) :: this
      real(dp), intent(in) :: input
      real(dp), intent(out) :: output

      ! Implementation
   end subroutine example_compute

   subroutine example_destroy(this)
      !! Clean up allocated memory
      class(example_type_t), intent(inout) :: this

      if (allocated(this%values)) deallocate(this%values)
   end subroutine example_destroy

end module rki_example
```

## Enforcement

- **Linter**: Run `fortitude check` before committing
- **CI**: GitHub Actions will fail if linter errors exist
- **Code Review**: Verify conventions are followed in PRs
