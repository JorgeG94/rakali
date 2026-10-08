!! Solution state for the unstructured central-upwind solver
module rki_state_unstr
   !! Composition wrapper around the unstructured solver state.
   !! Mirrors the structured `state_t` decomposition: this type is a
   !! thin shell that owns mode flags + timing + BC config + the four
   !! slope-limiter knobs, and composes the sub-state objects:
   !!
   !!   `state%barotropic`  -- depth-integrated CV + edge reconstruction
   !!                          + RK2 saves + Chapman state + per-edge
   !!                          mass flux + barotropic physics scalars.
   !!   `state%multilayer`  -- per-layer dynamics + tracer registry +
   !!                          vcoord + workspaces + EOS / mixing
   !!                          / BPG scalar config.
   !!
   !! The mesh metadata stays as a separate `mesh_unstr_t` argument
   !! threaded through the kernels (not embedded on `state`) -- same
   !! convention as before the split.
   !!
   !! Variable conventions:
   !!   stage = w  = h + b           (free-surface elevation)
   !!   xmom  = uh                   (x-momentum)
   !!   ymom  = vh                   (y-momentum)
   !!   bed   = b                    (bathymetry, static after init)
   !!   height = h                   (= stage - bed, recomputed each step)
   !!
   !! All centroid arrays are sized to `mesh%ncells_total` so halo
   !! cells are available to the kernels without special-casing.
   use rki_constants, only: wp
   use rki_mesh_unstr, only: mesh_unstr_t, mesh_unstr_si_quality_gate
   use rki_boundary_unstr, only: BC_WALL, MAX_BC_GROUPS
   use rki_barotropic_state_unstr, only: unstr_barotropic_state_t
   use rki_multilayer_state_unstr, only: unstr_multilayer_state_t
   use rki_surface_forcing_unstr, only: unstr_surface_forcing_t
   use rki_si_unstr_state, only: semi_implicit_state_unstr_t, SI_PRECOND_MG
   use rki_si_common, only: si_bsub_ml_conflict
   use rki_timestepper, only: TIME_RK2
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none
   private

   public :: unstr_state_t
   public :: LOW_FROUDE_OFF, LOW_FROUDE_1, LOW_FROUDE_2

   integer, parameter :: EDGES_PER_TRI = 3
      !! Number of edges per cell. Triangles only for now.

   ! Low-Froude correction modes for the central-upwind flux. Mirrored
   ! by the inlined `compute_local_froude_impl` copies in
   ! `rki_kernel_flux_unstr`; the kernel-side copies are intentional
   ! duplicates so NVHPC can inline the helper into `flux_one_cell`.
   integer, parameter :: LOW_FROUDE_OFF = 0
   integer, parameter :: LOW_FROUDE_1 = 1
   integer, parameter :: LOW_FROUDE_2 = 2

   type :: unstr_state_t
      !! Container for the solution state, factored into sub-objects.

      ! ----- sizes (cached from the source mesh) -----
      integer :: ncells_owned = 0
         !! Cells owned by this rank.
      integer :: ncells_total = 0
         !! ncells_owned + ncells_halo. All centroid/edge arrays are this long.

      ! ----- sub-state objects -----
      type(unstr_barotropic_state_t) :: barotropic
         !! Centroid CV + edge reconstruction + per-step accumulators
         !! + RK2 saves + Chapman BC state + per-edge mass flux +
         !! barotropic scalar physics.  See `rki_barotropic_state_unstr`.

      ! Coupled hydrostatic multilayer
      logical :: use_multilayer = .false.
         !! When true, `multilayer` is initialised and mapped.
      type(unstr_multilayer_state_t) :: multilayer
         !! Per-layer arrays, vcoord, tracer registry, workspaces, and
         !! all multilayer scalar config.  See `rki_multilayer_state_unstr`.

      ! Surface forcing fields (heat flux + virtual freshwater flux).
      ! Always allocated and zero-filled so the surface-flux kernel can
      ! run unconditionally; non-zero values are pushed in via the C API.
      type(unstr_surface_forcing_t) :: forcing
         !! Surface tracer forcing.  See `rki_surface_forcing_unstr`.

      ! Semi-implicit (Casulli) free-surface mode.
      logical :: use_semi_implicit = .false.
         !! When true, `si` is initialised + mapped and the solver dispatches
         !! to the unstructured SI step (`time_scheme="semi_implicit"`).
      type(semi_implicit_state_unstr_t) :: si
         !! CG workspace + edge face-depths for the TPFA Helmholtz eta-solve.
         !! Allocated only when `use_semi_implicit`.  See `rki_si_unstr_state`.
      real(wp) :: dt_fixed = 0.0_wp
         !! Fixed timestep (s) for the semi-implicit mode (bypasses the
         !! gravity-wave CFL).  Set from cfg%dt_fixed at init; 0 on the
         !! explicit path.

      ! ----- physics knobs (config-loaded later) -----
      real(wp) :: minimum_allowed_height = 1.0e-3_wp
         !! Minimum allowed depth. Cells below this are treated as dry.
      real(wp) :: cfl = 0.5_wp
         !! Default CFL for the unstructured solver.

      ! ----- slope-limiter beta values -----
      ! Wet/dry blending uses `hfactor` to interpolate between the dry
      ! value (small, conservative) and the wet value (close to 1).
      real(wp) :: beta_w = 1.0_wp
      real(wp) :: beta_w_dry = 0.2_wp
      real(wp) :: beta_uh = 0.8_wp
      real(wp) :: beta_uh_dry = 0.2_wp
      real(wp) :: beta_vh = 0.8_wp
      real(wp) :: beta_vh_dry = 0.2_wp

      logical :: extrapolate_velocity_second_order = .true.
         !! When true, the extrapolation kernel converts centroid momenta
         !! to velocities before reconstruction and back at the edge
         !! midpoints.  Much better behaved at wet/dry fronts.

      integer :: low_froude_mode = LOW_FROUDE_OFF
         !! Low-Froude correction mode for the central-upwind flux
         !! (LOW_FROUDE_OFF / _1 / _2).

      ! ----- boundary condition lookup table -----
      ! Fixed-size integer array: lives inside the struct itself, so
      ! `map(to: state)` carries it across to the device automatically.
      integer  :: bc_type(MAX_BC_GROUPS) = BC_WALL
         !! BC type for each boundary group `g = 1..MAX_BC_GROUPS`.
      real(wp) :: bc_eta_ext(MAX_BC_GROUPS) = 0.0_wp
         !! Per-group external stage (m) for `BC_TIDAL`.
      real(wp) :: bc_un_ext(MAX_BC_GROUPS) = 0.0_wp
         !! Per-group external normal velocity (m/s) for `BC_INFLOW`.
      real(wp) :: bc_h_ext(MAX_BC_GROUPS) = 0.0_wp
         !! Per-group external depth (m) for `BC_OUTFLOW`.
      real(wp) :: bc_Qn_ext(MAX_BC_GROUPS) = 0.0_wp
         !! Per-group specific discharge (m^2/s) for `BC_DISCHARGE`.
      real(wp) :: bc_sponge_strength(MAX_BC_GROUPS) = 0.0_wp
         !! Sponge relaxation rate [1/s] per BC group (0 = inactive).
         !! Consumed by `apply_sponge_unstr`. Rides the struct copyin
         !! like `bc_eta_ext` — no separate device map needed.
      real(wp) :: bc_sponge_h_ref(MAX_BC_GROUPS) = 0.0_wp
         !! Sponge reference depth h_ref [m] per BC group.
         !! `apply_sponge_unstr` relaxes `h = stage - bed` toward this value.

      ! ----- timing + diagnostics -----
      real(wp) :: dt = 0.0_wp
         !! Current timestep (s).
      real(wp) :: time = 0.0_wp
         !! Current simulation time (s).
      integer  :: step_count = 0
         !! Total steps taken since init.
      real(wp) :: mass_error = 0.0_wp
         !! Cumulative negative-mass clipped by `protect`.
      integer  :: time_scheme = TIME_RK2
         !! Explicit time integrator (`TIME_*` from rki_timestepper).
         !! Default RK2 reproduces the historical SSP-RK2 path bit-for-bit.
         !! 2D barotropic path only; multilayer keeps RK2.
   end type unstr_state_t

   public :: state_unstr_init
   public :: state_unstr_enable_si
   public :: state_unstr_destroy
   public :: state_unstr_enter_data
   public :: state_unstr_exit_data

contains

   subroutine state_unstr_init(state, mesh, nz_ml)
      !! Allocate every array in `state` to match the cell counts in
      !! `mesh`. All arrays are zero-initialised.  Calling on an
      !! already-initialised state is safe -- it destroys first.
      !!
      !! When `nz_ml` is present and >= 1, the multilayer sub-object
      !! is initialised: `use_multilayer` is set and
      !! `state%multilayer` allocates per-layer arrays + tracer registry.
      !! If `nz_ml` is absent the state stays 2D depth-averaged.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in)    :: mesh
      integer, intent(in), optional :: nz_ml

      integer :: n

      call state_unstr_destroy(state)

      n = mesh%ncells_total
      state%ncells_owned = mesh%ncells_owned
      state%ncells_total = n

      call state%barotropic%init(mesh)
      call state%forcing%init(mesh)

      state%dt = 0.0_wp
      state%time = 0.0_wp
      state%step_count = 0
      state%mass_error = 0.0_wp

      if (present(nz_ml)) then
         if (nz_ml < 1) then
            call logger%error("state_unstr: nz_ml must be >= 1 (got "// &
                              to_char(nz_ml)//")")
            return
         end if
         state%use_multilayer = .true.
         call state%multilayer%init(mesh, nz_ml)
         call logger%info("state_unstr: multilayer enabled with nz_ml = "// &
                          to_char(nz_ml))
      end if

      call logger%info("state_unstr: allocated for "//to_char(state%ncells_owned)// &
                       " owned cells ("//to_char(n)//" total)")
   end subroutine state_unstr_init

   subroutine state_unstr_enter_data(state)
      !! Map the unstructured state arrays into device memory.
      !!
      !! Called once after `state_unstr_init` and bathymetry/IC setup,
      !! before any kernel runs.  All subsequent kernel work happens
      !! on the device until `state_unstr_exit_data` is called.
      !! Mirrors `rki_solver::solver_enter_data`.
      type(unstr_state_t), intent(inout) :: state

      ! Map the parent struct descriptor first so the device knows the
      ! shape of `unstr_state_t`.  Components (and the value-type
      ! sub-objects' allocatables) attach to it below.  Without the
      ! parent map first, NVHPC stdpar can't dereference component
      ! arrays from inside a do-concurrent body.
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data

      call state%barotropic%enter_data()
      call state%forcing%enter_data()

      if (state%use_semi_implicit) then
         call state%si%enter_data()
      end if

      if (state%use_multilayer) then
         call state%multilayer%enter_data()

         ! Generic tracer registry: two-step attach pattern (parent
         ! array descriptor first, then per-element hTr/hTr0).  Mirrors
         ! the structured `solver_tracers_enter_data`.  Without this,
         ! NVHPC stdpar can't dereference `tracers(it)%hTr` from a
         ! do-concurrent body and the kernel hits CUDA_ERROR_ILLEGAL_ADDRESS.
         if (allocated(state%multilayer%tracers)) then
            call tracers_enter_data(state)
         end if
      end if
   end subroutine state_unstr_enter_data

   subroutine tracers_enter_data(state)
      !! Per-element copyin for the tracer registry.  Split out so the
      !! enclosing routine stays linear.
      type(unstr_state_t), intent(inout) :: state

      integer :: it

      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
      do it = 1, size(state%multilayer%tracers)
         ! [acc->omp-cpu] dropped device directive (host no-op): enter data
      end do

      ! The device map above is an element-by-element snapshot: growing
      ! the registry after this point would strand every mapping.
      ! register_passive_tracer refuses while locked.
      state%multilayer%registry_locked = .true.
   end subroutine tracers_enter_data

   subroutine state_unstr_exit_data(state)
      !! Copy state arrays back from the device and release device
      !! memory.  Counterpart to `state_unstr_enter_data`.
      type(unstr_state_t), intent(inout) :: state

      ! Components first, parent struct last (mirror of enter_data order).
      if (state%use_multilayer) then
         ! Generic tracer registry: per-element components first, then
         ! the array descriptor (mirror of enter_data order).
         if (allocated(state%multilayer%tracers)) then
            call tracers_exit_data(state)
         end if
         call state%multilayer%exit_data()
      end if

      if (state%use_semi_implicit) then
         call state%si%exit_data()
      end if

      call state%forcing%exit_data()
      call state%barotropic%exit_data()

      ! Parent struct descriptor last.
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
   end subroutine state_unstr_exit_data

   subroutine tracers_exit_data(state)
      !! Per-element copyout/delete for the tracer registry.
      type(unstr_state_t), intent(inout) :: state

      integer :: it

      do it = 1, size(state%multilayer%tracers)
         ! [acc->omp-cpu] dropped device directive (host no-op): exit data
         ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      end do
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data

      state%multilayer%registry_locked = .false.
   end subroutine tracers_exit_data

   subroutine state_unstr_destroy(state)
      !! Release every allocatable array and reset counters.
      type(unstr_state_t), intent(inout) :: state

      call state%barotropic%destroy()
      call state%multilayer%destroy()
      call state%forcing%destroy()
      call state%si%destroy()

      state%ncells_owned = 0
      state%ncells_total = 0
      state%dt = 0.0_wp
      state%time = 0.0_wp
      state%step_count = 0
      state%mass_error = 0.0_wp
      state%use_multilayer = .false.
      state%use_semi_implicit = .false.
      state%dt_fixed = 0.0_wp
   end subroutine state_unstr_destroy

   subroutine state_unstr_enable_si(state, mesh, precond, rtol, dt_fixed, ok, &
                                    wetdry, n_sub, bsub_fill)
      !! Enable the semi-implicit (Casulli) free-surface mode on this state:
      !! run the fail-loud TPFA mesh-quality gate, allocate the `si` slot
      !! (with the P3.2 multilayer column slots when `use_multilayer`, AND/OR
      !! the P3.3 Casulli wet/dry slots when `wetdry` — the two slot groups
      !! gate independently, so 3D + wet/dry coexist), and stash the CG
      !! controls.  Call AFTER `state_unstr_init` (so the multilayer
      !! sub-state exists and `nz_ml` is known) and BEFORE
      !! `state_unstr_enter_data`.  `ok=.false.` (gate failure or a bad
      !! precond having already error-stopped) leaves `use_semi_implicit`
      !! off so the caller can abort cleanly.
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh
      integer, intent(in) :: precond
         !! `SI_PRECOND_*` (already parsed; MG = aggregation multigrid here).
      real(wp), intent(in) :: rtol
         !! CG relative residual tolerance.
      real(wp), intent(in) :: dt_fixed
         !! Fixed timestep (s); the SI mode bypasses the gravity-wave CFL.
      logical, intent(out) :: ok
      logical, intent(in), optional :: wetdry
         !! Casulli wetting/drying (2D barotropic + 3D sigma multilayer;
         !! the linear b_sub fill is 2D-only — see si_bsub_ml_conflict).
         !! Default .false.
      integer, intent(in), optional :: n_sub
         !! Subgrid sub-bottoms per cell (>= 1).  Default 1 (flat bottom).
      integer, intent(in), optional :: bsub_fill
         !! Sub-bottom fill mode (`SI_BSUB_FILL_*`).  Default 0 (FLAT).

      logical :: wd
      integer :: ns, nzml, bf

      wd = .false.
      if (present(wetdry)) wd = wetdry
      ns = 1
      if (present(n_sub)) ns = max(1, n_sub)
      bf = 0
      if (present(bsub_fill)) bf = bsub_fill

      ! CONSERVATION GUARD (fail-loud): the `linear` subgrid b_sub fill
      ! (n_sub >= 2) is INCOMPATIBLE with the multilayer (nz > 1) wet/dry
      ! step.  The free-surface eta-solve conserves the SUBGRID wetted volume
      ! v_cell(eta) = mean_s max(eta - b_sub_s, 0) (si_unstr_wetdry_volume),
      ! but `si_unstr_update_thickness_ml` and the recovery Dirichlet ghost
      ! consume max(eta - b, 0) on the PLAIN bed.  Under `flat` (b_sub == b)
      ! these are identical so mass conserves; under `linear` they diverge
      ! across the partial-flooding band => mass leaks at the wet/dry front —
      ! the exact defect the port exists to prevent.  A v_cell-consistent ML
      ! thickness update is a documented deferred limitation; until then
      ! `linear + ML` is rejected rather than run silently.  `flat + ML`
      ! (conservative, tested) and `linear + 2D` (barotropic uses v_cell
      ! throughout) both stay allowed.
      if (si_bsub_ml_conflict(bf, ns, state%use_multilayer, state%multilayer%nz_ml)) then
         call logger%error("state_unstr_enable_si: wetdry_bsub_fill='linear' "// &
                           "(n_sub>=2) is not supported with the multilayer "// &
                           "(nz>1) wet/dry step — the ML thickness update "// &
                           "consumes the plain bed, not the subgrid v_cell, so "// &
                           "mass leaks at the wet/dry front.  Use "// &
                           "wetdry_bsub_fill='flat' with multilayer, or 'linear' "// &
                           "only with the 2D barotropic path.")
         error stop "state_unstr_enable_si: linear b_sub incompatible with multilayer wet/dry"
      end if

      call mesh_unstr_si_quality_gate(mesh, ok)
      if (.not. ok) return

      ! The multilayer column slots (nz>1) and the wet/dry slots gate
      ! INDEPENDENTLY in `si%init`, so 3D + wet/dry coexist: the ML wet/dry
      ! step (`si_unstr_step_ml`, wet/dry branch) reuses the mesh-agnostic
      ! V(eta)/a/Newton primitives + the dsig scale-together column drain.
      if (state%use_multilayer .and. state%multilayer%nz_ml > 1) then
         nzml = state%multilayer%nz_ml
         call state%si%init(mesh, nz=nzml, wetdry=wd, n_sub=ns)
      else
         call state%si%init(mesh, wetdry=wd, n_sub=ns)
      end if
      state%si%bsub_fill = bf
      state%si%precond = precond
      state%si%rtol = rtol
      state%dt_fixed = dt_fixed
      state%use_semi_implicit = .true.

      ! Aggregation-MG preconditioner: v1 is the 2D barotropic path only
      ! (wet/dry Newton + multilayer column solves keep plain/Jacobi CG).
      ! Fail loud rather than silently ignore an MG request there.
      if (precond == SI_PRECOND_MG) then
         if (wd .or. (state%use_multilayer .and. state%multilayer%nz_ml > 1)) then
            call logger%error("state_unstr_enable_si: precond='mg' (aggregation "// &
                              "multigrid) is 2D-barotropic only in v1 — the wet/dry "// &
                              "Newton and multilayer column solves are not yet wired "// &
                              "to the MG hierarchy.  Use precond='jacobi' with "// &
                              "wet/dry or multilayer.")
            error stop "state_unstr_enable_si: mg precond not supported with wetdry/multilayer (v1)"
         end if
         ! Build the (static) level hierarchy now (geometric aggregation;
         ! operator coefficients are filled per step by the refresh in the SI
         ! step).  Host-only; runs before enter_data.
         call state%si%build_mg(mesh)
      end if
   end subroutine state_unstr_enable_si

end module rki_state_unstr
