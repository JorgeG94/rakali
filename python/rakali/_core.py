"""Core Rakali wrapper: ctypes bindings to the librakali shared library."""
import ctypes
import os
import re
import sys
import numpy as np
from pathlib import Path

from ._lifecycle import (_SimBase, _owned_view, _require_file,
                         _require_writable_dir)


_TIME_EPS = 1e-10  # Must match 1.0e-10_wp in rki_advance (rki_api.F90)

# Platform-specific shared-library suffixes.  macOS produces .dylib,
# Linux/BSD produces .so, Windows produces .dll.  We probe all of them
# so the same package works across platforms without users having to
# set RAKALI_LIB by hand.
if sys.platform == "darwin":
    _LIB_SUFFIXES = (".dylib", ".so")
elif sys.platform == "win32":
    _LIB_SUFFIXES = (".dll",)
else:
    _LIB_SUFFIXES = (".so",)


def _find_library():
    """Find the librakali shared library via env var or common locations.

    Search order (first hit wins):
      1. ``RAKALI_LIB`` environment variable (explicit override)
      2. ``rakali/lib/librakali.{dylib,so,dll}`` next to the installed
         package (forward-compat for binary wheels that bundle the lib)
      3. any ``build*/`` directory at the repo root whose name starts
         with ``build`` and that contains a ``librakali.{dylib,so,dll}``
         -- covers ``build``, ``build_nvhpc``, ``build_gcc``,
         ``build_validate/<backend>``, etc.

    On macOS both ``.dylib`` and ``.so`` suffixes are accepted (CMake
    defaults to ``.dylib`` for ``SHARED`` libraries, but some build
    scripts explicitly set ``.so``).
    """
    # 1. Explicit path
    env = os.environ.get("RAKALI_LIB")
    if env:
        if not os.path.isfile(env):
            raise FileNotFoundError(
                f"RAKALI_LIB is set to {env!r}, which is not a file")
        return env

    pkg_dir = Path(__file__).resolve().parent

    def _with_suffixes(base: Path):
        """Yield ``base`` with each platform shared-library suffix."""
        for suffix in _LIB_SUFFIXES:
            yield base.with_suffix(suffix)

    # 2. Bundled inside the installed package (future binary wheel layout)
    for bundled in _with_suffixes(pkg_dir / "lib" / "librakali.so"):
        if bundled.is_file():
            return str(bundled)

    # 3. Any build*/ directory at the repo root. Recursively walk at
    # depth 1-2 so both flat (``build/librakali.so``) and nested
    # (``build_validate/<backend>/librakali.so``) layouts work.
    repo_dir = pkg_dir.parent.parent  # python/rakali -> repo root
    for suffix in _LIB_SUFFIXES:
        name = f"librakali{suffix}"
        # Direct children: build*/librakali.so
        for cand in sorted(repo_dir.glob(f"build*/{name}")):
            if cand.is_file():
                return str(cand)
        # One level deeper: build*/something/librakali.so
        for cand in sorted(repo_dir.glob(f"build*/*/{name}")):
            if cand.is_file():
                return str(cand)

    tried_suffixes = "/".join(s.lstrip(".") for s in _LIB_SUFFIXES)
    raise FileNotFoundError(
        f"Cannot find librakali.{{{tried_suffixes}}}. "
        f"Set the RAKALI_LIB environment variable to the absolute path, "
        f"or build the shared library with "
        f"`cmake -B build -S . -DRAKALI_ENABLE_PYTHON=ON && "
        f"cmake --build build`."
    )


def _load_lib():
    """Load the shared library and set up function signatures.

    Every entry point takes the opaque ``void* handle`` as its first
    argument (by value, ``c_void_p``).  Lifecycle entries that mutate the
    handle itself (``rki_create``, ``rki_create_from_namelist``,
    ``rki_destroy``) take it via ``POINTER(c_void_p)`` so the callee can
    write back into the caller's slot.
    """
    lib = ctypes.CDLL(_find_library())

    HANDLE = ctypes.c_void_p           # void* (by value) for ordinary calls
    HANDLE_OUT = ctypes.POINTER(HANDLE) # void** for create / destroy

    p_int = ctypes.POINTER(ctypes.c_int)

    # Solver working precision (bytes): pick the numpy dtype + ctypes
    # scalar used for the zero-copy field views so a single-precision
    # build (RAKALI_ENABLE_DOUBLE=OFF, wp=real32) is read as float32
    # instead of being mis-interpreted as float64.  Process-global, so
    # cache it on the lib once.  Setters stay c_double (their Fortran ABI
    # converts to wp internally).
    lib.rki_working_precision.argtypes = []
    lib.rki_working_precision.restype = ctypes.c_int
    _wp_bytes = lib.rki_working_precision()
    lib._rki_wp_ctype = ctypes.c_float if _wp_bytes == 4 else ctypes.c_double
    lib._rki_wp_dtype = np.float32 if _wp_bytes == 4 else np.float64

    # --- Lifecycle ---
    lib.rki_create.argtypes = [ctypes.c_char_p, ctypes.c_int, HANDLE_OUT]
    lib.rki_create.restype = ctypes.c_int

    lib.rki_finalize.argtypes = [HANDLE]
    lib.rki_finalize.restype = ctypes.c_int

    lib.rki_create_from_namelist.argtypes = [ctypes.c_char_p, ctypes.c_int, HANDLE_OUT]
    lib.rki_create_from_namelist.restype = ctypes.c_int

    # In-memory create (no temp file): config passed as a namelist buffer.
    lib.rki_create_from_string.argtypes = [ctypes.c_char_p, ctypes.c_int, HANDLE_OUT]
    lib.rki_create_from_string.restype = ctypes.c_int
    lib.rki_create_from_namelist_string.argtypes = [
        ctypes.c_char_p, ctypes.c_int, HANDLE_OUT]
    lib.rki_create_from_namelist_string.restype = ctypes.c_int

    lib.rki_destroy.argtypes = [HANDLE_OUT]
    lib.rki_destroy.restype = None

    # --- Stepping ---
    lib.rki_step.argtypes = [HANDLE, ctypes.c_int]
    lib.rki_step.restype = ctypes.c_int

    lib.rki_advance.argtypes = [HANDLE, ctypes.c_double]
    lib.rki_advance.restype = ctypes.c_int

    # --- Queries ---
    lib.rki_get_time.argtypes = [HANDLE]
    lib.rki_get_time.restype = ctypes.c_double

    lib.rki_get_dt.argtypes = [HANDLE]
    lib.rki_get_dt.restype = ctypes.c_double

    lib.rki_get_step_count.argtypes = [HANDLE]
    lib.rki_get_step_count.restype = ctypes.c_int

    lib.rki_get_nz_ml.argtypes = [HANDLE]
    lib.rki_get_nz_ml.restype = ctypes.c_int

    _ptr_int_int = [
        HANDLE,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
    ]
    _ptr_int_int_int = _ptr_int_int + [ctypes.POINTER(ctypes.c_int)]

    lib.rki_get_grid_info.argtypes = [
        HANDLE,
        ctypes.POINTER(ctypes.c_int),     # nx
        ctypes.POINTER(ctypes.c_int),     # ny
        ctypes.POINTER(ctypes.c_int),     # nghost
        ctypes.POINTER(ctypes.c_double),  # dx
        ctypes.POINTER(ctypes.c_double),  # dy
    ]
    lib.rki_get_grid_info.restype = None

    # --- 2D field pointers ---
    for name in ['rki_get_h_ptr', 'rki_get_hu_ptr', 'rki_get_hv_ptr', 'rki_get_b_ptr',
                 'rki_get_h0_ptr', 'rki_get_hu0_ptr', 'rki_get_hv0_ptr']:
        getattr(lib, name).argtypes = _ptr_int_int
        getattr(lib, name).restype = ctypes.c_int

    for name in ['rki_get_clamped_elevation', 'rki_get_manning_n',
                 'rki_get_coriolis_f']:
        getattr(lib, name).argtypes = [HANDLE]
        getattr(lib, name).restype = ctypes.c_double

    # --- 3D field pointers ---
    for name in ['rki_get_h_layer_ptr', 'rki_get_hu_layer_ptr',
                 'rki_get_hv_layer_ptr', 'rki_get_salinity_ptr',
                 'rki_get_temperature_ptr']:
        getattr(lib, name).argtypes = _ptr_int_int_int
        getattr(lib, name).restype = ctypes.c_int

    # --- Setters ---
    lib.rki_set_bathymetry.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_bathymetry.restype = ctypes.c_int

    lib.rki_set_initial_h.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_initial_h.restype = ctypes.c_int

    # --- Forcing setters ---
    lib.rki_set_wind.argtypes = [
        HANDLE,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_wind.restype = ctypes.c_int

    lib.rki_set_atm_pressure.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_atm_pressure.restype = ctypes.c_int

    lib.rki_set_precipitation.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_precipitation.restype = ctypes.c_int

    lib.rki_set_heat_flux.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_heat_flux.restype = ctypes.c_int

    lib.rki_set_evap_minus_precip.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
    ]
    lib.rki_set_evap_minus_precip.restype = ctypes.c_int

    # --- DA / state setters ---
    for name in ['rki_set_h', 'rki_set_hu', 'rki_set_hv']:
        getattr(lib, name).argtypes = [
            HANDLE, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int
        ]
        getattr(lib, name).restype = ctypes.c_int

    lib.rki_set_salinity.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ]
    lib.rki_set_salinity.restype = ctypes.c_int

    lib.rki_set_temperature.argtypes = [
        HANDLE, ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ]
    lib.rki_set_temperature.restype = ctypes.c_int

    # --- GPU sync ---
    lib.rki_sync_host.argtypes = [HANDLE]
    lib.rki_sync_host.restype = None

    lib.rki_sync_from_device.argtypes = [HANDLE]
    lib.rki_sync_from_device.restype = None

    # --- Restart save / load ---
    lib.rki_save_restart.argtypes = [HANDLE, ctypes.c_char_p, ctypes.c_int]
    lib.rki_save_restart.restype = ctypes.c_int

    lib.rki_load_restart.argtypes = [HANDLE, ctypes.c_char_p, ctypes.c_int]
    lib.rki_load_restart.restype = ctypes.c_int

    # --- Boundary condition setters (mid-sim forcing updates) ---
    for name in ['rki_set_inflow_velocity', 'rki_set_inflow_salinity',
                 'rki_set_inflow_temperature',
                 'rki_set_discharge_flux', 'rki_set_clamped_elevation']:
        getattr(lib, name).argtypes = [HANDLE, ctypes.c_double]
        getattr(lib, name).restype = ctypes.c_int

    lib.rki_set_use_ppm_tracer.argtypes = [HANDLE, ctypes.c_int]
    lib.rki_set_use_ppm_tracer.restype = ctypes.c_int
    lib.rki_set_use_hllc.argtypes = [HANDLE, ctypes.c_int]
    lib.rki_set_use_hllc.restype = ctypes.c_int

    lib.rki_set_tidal_constituents.argtypes = [
        HANDLE,
        ctypes.POINTER(ctypes.c_double),  # amp
        ctypes.POINTER(ctypes.c_double),  # phase
        ctypes.POINTER(ctypes.c_double),  # omega
        ctypes.c_int,                     # n
    ]
    lib.rki_set_tidal_constituents.restype = ctypes.c_int

    return lib


# Module-level library handle (loaded once)
_lib = None


def _get_lib():
    global _lib
    if _lib is None:
        _lib = _load_lib()
    return _lib


class Rakali(_SimBase):
    """Python interface to the Rakali shallow water solver.

    Usage (from namelist file):
        sim = Rakali("river_plume.nml")

    Usage (programmatic):
        sim = Rakali(
            nx=300, ny=150, dx=10.0, dy=10.0,
            t_end=3600.0,
            manning_n=0.025,
            bathymetry_file="bathymetry.nc",
            bc_west="inflow", bc_east="open",
            inflow_velocity=0.5,
            use_multilayer=True, nz_layers=5,
            initial_salinity=35.0,
        )

    The most common namelist parameters can be passed as keyword arguments
    (see ``Rakali._PARAM_MAP``); anything else needs a namelist file. See
    docs/REFERENCE.md for the full parameter list.

    Notes
    -----
    Arrays returned by properties (h, hu, hv, b, etc.) are zero-copy
    views into Fortran memory. They use Fortran (column-major) ordering.
    The arrays include ghost cells; use h_phys, etc. for the physical domain.
    """

    # Map Python kwargs to (namelist_group, fortran_name)
    _PARAM_MAP = {
        # grid_nml
        'nx': ('grid_nml', 'nx'), 'ny': ('grid_nml', 'ny'),
        'dx': ('grid_nml', 'dx'), 'dy': ('grid_nml', 'dy'),
        'nghost': ('grid_nml', 'nghost'),
        # time_nml
        't_end': ('time_nml', 't_end'), 'cfl': ('time_nml', 'cfl'),
        'dt_max': ('time_nml', 'dt_max'), 'dt_fixed': ('time_nml', 'dt_fixed'),
        'cfl_interval': ('time_nml', 'cfl_interval'),
        'time_scheme': ('time_nml', 'time_scheme'),
        # physics_nml
        'manning_n': ('physics_nml', 'manning_n'),
        'dry_threshold': ('physics_nml', 'dry_threshold'),
        'wind_stress_x': ('physics_nml', 'wind_stress_x'),
        'wind_stress_y': ('physics_nml', 'wind_stress_y'),
        'coriolis_f': ('physics_nml', 'coriolis_f'),
        'use_hllc': ('physics_nml', 'use_hllc'),
        # output_nml
        'output_to_file': ('output_nml', 'output_to_file'),
        'output_dir': ('output_nml', 'output_dir'),
        'output_interval': ('output_nml', 'output_interval'),
        'output_format': ('output_nml', 'output_format'),
        'restart_interval': ('output_nml', 'restart_interval'),
        'async_output': ('output_nml', 'async_output'),
        'async_buffer_slots': ('output_nml', 'async_buffer_slots'),
        'bathymetry_file': ('output_nml', 'bathymetry_file'),
        'restart_file': ('output_nml', 'restart_file'),
        'forcing_file': ('output_nml', 'forcing_file'),
        'forcing_x_origin': ('output_nml', 'forcing_x_origin'),
        'forcing_y_origin': ('output_nml', 'forcing_y_origin'),
        'gauge_file': ('output_nml', 'gauge_file'),
        'gauge_output_interval': ('output_nml', 'gauge_output_interval'),
        'compress_output': ('output_nml', 'compress_output'),
        'compress_level': ('output_nml', 'compress_level'),
        # boundary_nml
        'bc_west': ('boundary_nml', 'bc_west'),
        'bc_east': ('boundary_nml', 'bc_east'),
        'bc_south': ('boundary_nml', 'bc_south'),
        'bc_north': ('boundary_nml', 'bc_north'),
        'tidal_amplitude': ('boundary_nml', 'tidal_amplitude'),
        'tidal_period': ('boundary_nml', 'tidal_period'),
        'n_tidal_constituents': ('boundary_nml', 'n_tidal_constituents'),
        'tidal_amp': ('boundary_nml', 'tidal_amp'),
        'tidal_phase': ('boundary_nml', 'tidal_phase'),
        'tidal_omega': ('boundary_nml', 'tidal_omega'),
        'inflow_velocity': ('boundary_nml', 'inflow_velocity'),
        'inflow_salinity': ('boundary_nml', 'inflow_salinity'),
        'inflow_temperature': ('boundary_nml', 'inflow_temperature'),
        'discharge_flux': ('boundary_nml', 'discharge_flux'),
        'clamped_elevation': ('boundary_nml', 'clamped_elevation'),
        'sponge_width': ('boundary_nml', 'sponge_width'),
        'sponge_strength': ('boundary_nml', 'sponge_strength'),
        'sponge_h_ref': ('boundary_nml', 'sponge_h_ref'),
        'nesting_file': ('boundary_nml', 'nesting_file'),
        'nesting_x_origin': ('boundary_nml', 'nesting_x_origin'),
        'nesting_y_origin': ('boundary_nml', 'nesting_y_origin'),
        # nonhydrostatic_nml
        'use_nonhydrostatic': ('nonhydrostatic_nml', 'use_nonhydrostatic'),
        'use_multilayer': ('nonhydrostatic_nml', 'use_multilayer'),
        'nz_layers': ('nonhydrostatic_nml', 'nz_layers'),
        'nh_tolerance': ('nonhydrostatic_nml', 'nh_tolerance'),
        'nh_max_iterations': ('nonhydrostatic_nml', 'nh_max_iterations'),
        'nh_preconditioner': ('nonhydrostatic_nml', 'nh_preconditioner'),
        'nh_solve_interval': ('nonhydrostatic_nml', 'nh_solve_interval'),
        'bpg_method': ('nonhydrostatic_nml', 'bpg_method'),
        'ml_min_depth': ('nonhydrostatic_nml', 'ml_min_depth'),
        'rho_0': ('nonhydrostatic_nml', 'rho_0'),
        'vmix_nu_max': ('nonhydrostatic_nml', 'vmix_nu_max'),
        'vmix_nu_bg': ('nonhydrostatic_nml', 'vmix_nu_bg'),
        'vmix_ri_alpha': ('nonhydrostatic_nml', 'vmix_ri_alpha'),
        'hdiff_kappa': ('nonhydrostatic_nml', 'hdiff_kappa'),
        'smagorinsky_cs': ('nonhydrostatic_nml', 'smagorinsky_cs'),
        'split_ratio': ('nonhydrostatic_nml', 'split_ratio'),
        'split_averaging': ('nonhydrostatic_nml', 'split_averaging'),
        'out_h_layer': ('output_nml', 'out_h_layer'),
        'out_velocity': ('output_nml', 'out_velocity'),
        'out_salinity': ('output_nml', 'out_salinity'),
        'out_temperature': ('output_nml', 'out_temperature'),
        'out_density': ('output_nml', 'out_density'),
        'out_w': ('output_nml', 'out_w'),
        # tracer_nml (new — salinity + temperature tracer configuration)
        'initial_salinity': ('tracer_nml', 'initial_salinity'),
        'S_ref': ('tracer_nml', 'S_ref'),
        'beta_S': ('tracer_nml', 'beta_S'),
        'S_min': ('tracer_nml', 'S_min'),
        'S_max': ('tracer_nml', 'S_max'),
        'kappa_S_bg': ('tracer_nml', 'kappa_S_bg'),
        'initial_temperature': ('tracer_nml', 'initial_temperature'),
        'T_ref': ('tracer_nml', 'T_ref'),
        'alpha_T': ('tracer_nml', 'alpha_T'),
        'T_min': ('tracer_nml', 'T_min'),
        'T_max': ('tracer_nml', 'T_max'),
        'kappa_T_bg': ('tracer_nml', 'kappa_T_bg'),
        'use_ppm_tracer': ('tracer_nml', 'use_ppm_tracer'),
        # initial_condition_nml
        'init_condition': ('initial_condition_nml', 'init_condition'),
        'eta_initial': ('initial_condition_nml', 'eta_initial'),
        'h0': ('initial_condition_nml', 'h0'),
        'hump_amplitude': ('initial_condition_nml', 'hump_amplitude'),
        'hump_sigma': ('initial_condition_nml', 'hump_sigma'),
        'hump_x_offset': ('initial_condition_nml', 'hump_x_offset'),
        'hump_y_offset': ('initial_condition_nml', 'hump_y_offset'),
        'init_u_shear': ('initial_condition_nml', 'init_u_shear'),
        # vcoord_nml (vertical coordinate configuration)
        'vcoord_type': ('vcoord_nml', 'vcoord_type'),
        'remap_method': ('vcoord_nml', 'remap_method'),
        'vcoord_depth_transition': ('vcoord_nml', 'vcoord_depth_transition'),
        'vcoord_blend_width': ('vcoord_nml', 'vcoord_blend_width'),
        'zstar_h_surf_target': ('vcoord_nml', 'zstar_h_surf_target'),
        'zstar_stretching': ('vcoord_nml', 'zstar_stretching'),
        'zstar_h_min': ('vcoord_nml', 'zstar_h_min'),
        'zstar_n_surf': ('vcoord_nml', 'zstar_n_surf'),
        # logging_nml
        'log_level': ('logging_nml', 'log_level'),
        'status_interval': ('logging_nml', 'status_interval'),
        # metadata_nml
        'reference_time': ('metadata_nml', 'reference_time'),
        'institution': ('metadata_nml', 'institution'),
    }

    def __init__(self, namelist=None, defer_init=False, **kwargs):
        """Create a Rakali simulation.

        Parameters
        ----------
        namelist : str, optional
            Path to a Fortran namelist file.
        defer_init : bool, optional
            If True, only allocate state (phase 1). You must call
            set_bathymetry() / set_initial_h() then finalize() before
            stepping. Default False (one-shot init).
        **kwargs
            Namelist parameters as keyword arguments; see
            ``Rakali._PARAM_MAP`` for the supported names.

        Examples
        --------
        One-shot (existing behaviour):

        >>> sim = Rakali("config.nml")
        >>> sim = Rakali(nx=200, ny=100, dx=10.0, dy=10.0, t_end=3600.0)

        Two-phase (programmatic bathymetry):

        >>> sim = Rakali(nx=200, ny=100, dx=10.0, dy=10.0, t_end=3600.0,
        ...              defer_init=True)
        >>> sim.set_bathymetry(my_bathy_array)
        >>> sim.set_initial_h(my_depth_array)
        >>> sim.finalize()
        """
        self._lib = _get_lib()
        self._handle = ctypes.c_void_p()
        self._alive = False  # set True after create; gates teardown

        if namelist is not None and kwargs:
            raise ValueError("Pass either a namelist file or keyword arguments, not both")

        if namelist is not None:
            # A missing input file would end in a Fortran ERROR STOP that
            # kills the interpreter, so fail early with a Python exception.
            self._precheck_input_files(
                self._scan_namelist_inputs(namelist),
                Path(namelist).resolve().parent)
            # Existing namelist file: create from the path, chdir so any
            # relative paths inside it (bathymetry_file, …) resolve.
            nml_path = str(Path(namelist).resolve())
            nml_dir = str(Path(namelist).resolve().parent)
            old_cwd = os.getcwd()
            os.chdir(nml_dir)
            try:
                path = nml_path.encode()
                if defer_init:
                    rc = self._lib.rki_create(path, len(path),
                                              ctypes.byref(self._handle))
                    self._finalised = False
                else:
                    rc = self._lib.rki_create_from_namelist(
                        path, len(path), ctypes.byref(self._handle))
                    self._finalised = True
                if rc != 0:
                    raise RuntimeError(f"Failed to create simulation (rc={rc})")
            finally:
                os.chdir(old_cwd)
        elif kwargs:
            # Programmatic config: build the namelist text and hand it to
            # the C side as an in-memory buffer — NO temp file, NO chdir.
            # Relative paths (bathymetry_file, …) resolve against the
            # caller's cwd, which is left untouched.
            self._precheck_input_files(
                {k: v for k, v in kwargs.items() if k in self._INPUT_FILE_KEYS},
                Path.cwd())
            buf = self._kwargs_to_namelist(kwargs).encode()
            if defer_init:
                rc = self._lib.rki_create_from_string(
                    buf, len(buf), ctypes.byref(self._handle))
                self._finalised = False
            else:
                rc = self._lib.rki_create_from_namelist_string(
                    buf, len(buf), ctypes.byref(self._handle))
                self._finalised = True
            if rc != 0:
                raise RuntimeError(f"Failed to create simulation (rc={rc})")
        else:
            raise ValueError("Pass a namelist file or keyword arguments")

        self._alive = True
        self._query_grid()

    # Namelist keys naming files the library API READS at create time. (Only
    # bathymetry: forcing / nesting / gauge / restart / output files are
    # consumed by the Fortran driver executable, never by this API.)
    _INPUT_FILE_KEYS = ('bathymetry_file',)

    # Accepted by the namelist schema but silently ignored through the
    # library API -- rejecting them beats a run that quietly writes nothing.
    _UNSUPPORTED_KWARGS = frozenset((
        'output_to_file', 'output_dir', 'output_interval', 'output_format',
        'restart_interval', 'async_output', 'async_buffer_slots',
        'restart_file', 'forcing_file', 'forcing_x_origin',
        'forcing_y_origin', 'gauge_file', 'gauge_output_interval',
        'compress_output', 'compress_level', 'out_h_layer', 'out_velocity',
        'out_salinity', 'out_temperature', 'out_density', 'out_w',
        'nesting_file', 'nesting_x_origin', 'nesting_y_origin',
    ))

    # Enumerated values. The Fortran parsers silently fall back to a default
    # on a typo (or stop the process), so reject unknown names up front.
    _ENUM_VALUES = {
        'bc_west': ('wall', 'open', 'tidal', 'nested', 'inflow', 'discharge',
                    'clamped', 'sponge', 'chapman'),
        'time_scheme': ('euler', 'rk2', 'rk3', 'rk43'),
        'remap_method': ('pcm', 'plm', 'ppm', 'ppm_h4', 'pqm'),
        'vcoord_type': ('sigma', 'zsigma', 'zstar', 'zstar_lite',
                        'zstar_full', 'zstar_sigma'),
    }
    for _k in ('bc_east', 'bc_south', 'bc_north'):
        _ENUM_VALUES[_k] = _ENUM_VALUES['bc_west']
    del _k

    # Keys that must be strictly positive (nonpositive values reach a Fortran
    # ERROR STOP in validate_config).
    _POSITIVE = ('nx', 'ny', 'dx', 'dy', 't_end', 'cfl', 'nz_layers')

    @classmethod
    def _scan_namelist_inputs(cls, namelist):
        """Return {key: path} for input-file entries in a namelist file."""
        nml = Path(namelist)
        if not nml.is_file():
            raise FileNotFoundError(f"Namelist file not found: {nml}")
        found = {}
        pat = re.compile(r'^\s*(%s)\s*=\s*["\']([^"\']*)["\']'
                         % '|'.join(cls._INPUT_FILE_KEYS),
                         re.IGNORECASE | re.MULTILINE)
        for m in pat.finditer(nml.read_text()):
            found[m.group(1).lower()] = m.group(2)
        return found

    @staticmethod
    def _precheck_input_files(entries, base_dir):
        """Raise FileNotFoundError for any named input file that is absent."""
        for key, value in entries.items():
            if not value:
                continue
            path = Path(str(value))
            if not path.is_absolute():
                path = Path(base_dir) / path
            if not path.is_file():
                raise FileNotFoundError(
                    f"{key} = {value!r} not found (looked for {path})")

    @staticmethod
    def _format_value(v):
        """Format a Python value as a Fortran namelist value."""
        if isinstance(v, np.generic):
            v = v.item()  # np.bool_ / np.float64 / np.int64 -> builtin
        if v is None:
            raise ValueError("None is not a valid namelist value")
        if isinstance(v, float) and not np.isfinite(v):
            raise ValueError(f"non-finite value {v!r} is not allowed")
        if isinstance(v, bool):
            return '.true.' if v else '.false.'
        elif isinstance(v, str):
            return f'"{v}"'
        elif isinstance(v, (list, tuple)):
            return ', '.join(Rakali._format_value(x) for x in v)
        elif isinstance(v, float):
            return repr(v)
        else:
            return str(v)

    @classmethod
    def _kwargs_to_namelist(cls, kwargs):
        """Convert keyword arguments to a Fortran namelist string."""
        # Group parameters by namelist group
        groups = {}
        for key, value in kwargs.items():
            if key in cls._UNSUPPORTED_KWARGS:
                raise ValueError(
                    f"'{key}' is read by the Fortran driver executable, not "
                    f"the library API, so it would be silently ignored here. "
                    f"Use save_restart() for restarts, or run the namelist "
                    f"with the `rakali` executable for file output/forcing.")
            if key not in cls._PARAM_MAP:
                raise ValueError(
                    f"Unknown parameter '{key}'. "
                    f"Valid parameters: {sorted(cls._PARAM_MAP.keys())}"
                )
            if isinstance(value, np.generic):
                value = value.item()
            if key in cls._ENUM_VALUES:
                if (not isinstance(value, str)
                        or value.lower() not in cls._ENUM_VALUES[key]):
                    raise ValueError(
                        f"{key}={value!r} is not valid; choose one of "
                        f"{cls._ENUM_VALUES[key]}")
            if key in cls._POSITIVE and not (
                    isinstance(value, (int, float))
                    and not isinstance(value, bool) and value > 0):
                raise ValueError(f"{key} must be a positive number, "
                                 f"got {value!r}")
            group, fname = cls._PARAM_MAP[key]
            if group not in groups:
                groups[group] = []
            groups[group].append((fname, value))

        # Build namelist string
        lines = ["! Auto-generated by Rakali Python API\n"]
        for group, params in groups.items():
            lines.append(f"&{group}")
            for fname, value in params:
                lines.append(f"  {fname} = {cls._format_value(value)}")
            lines.append("/\n")
        return '\n'.join(lines)

    def _query_grid(self):
        nx = ctypes.c_int()
        ny = ctypes.c_int()
        ng = ctypes.c_int()
        dx = ctypes.c_double()
        dy = ctypes.c_double()
        self._lib.rki_get_grid_info(
            self._handle,
            ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(ng),
            ctypes.byref(dx), ctypes.byref(dy),
        )
        self._nx_phys = nx.value
        self._ny_phys = ny.value
        self._nghost = ng.value
        self._nx_total = nx.value + 2 * ng.value
        self._ny_total = ny.value + 2 * ng.value
        self._dx = dx.value
        self._dy = dy.value

    def _check_alive(self):
        if not self._alive:
            raise RuntimeError("Simulation has been destroyed")

    def _check_ready(self):
        self._check_alive()
        if not self._finalised:
            raise RuntimeError(
                "Simulation not finalised. Call finalize() after "
                "set_bathymetry() / set_initial_h()."
            )

    # --- Field accessors ---

    def _get_2d_field(self, getter_name):
        """Get a zero-copy numpy view of a 2D Fortran array."""
        self._check_alive()
        ptr = ctypes.c_void_p()
        nx = ctypes.c_int()
        ny = ctypes.c_int()
        rc = getattr(self._lib, getter_name)(
            self._handle,
            ctypes.byref(ptr), ctypes.byref(nx), ctypes.byref(ny)
        )
        if rc != 0:
            raise RuntimeError(f"{getter_name} failed (rc={rc})")
        n = nx.value * ny.value
        arr = _owned_view(ptr.value, self._lib._rki_wp_ctype, n, self)
        # Fortran column-major: first index is fast
        return arr.reshape((nx.value, ny.value), order='F')

    def _get_3d_field(self, getter_name):
        """Get a zero-copy numpy view of a 3D Fortran array."""
        self._check_alive()
        ptr = ctypes.c_void_p()
        nx = ctypes.c_int()
        ny = ctypes.c_int()
        nz = ctypes.c_int()
        rc = getattr(self._lib, getter_name)(
            self._handle,
            ctypes.byref(ptr), ctypes.byref(nx), ctypes.byref(ny), ctypes.byref(nz)
        )
        if rc != 0:
            raise RuntimeError(f"{getter_name} failed (rc={rc})")
        n = nx.value * ny.value * nz.value
        arr = _owned_view(ptr.value, self._lib._rki_wp_ctype, n, self)
        return arr.reshape((nx.value, ny.value, nz.value), order='F')

    def _phys_slice(self, arr):
        """Slice a 2D or 3D array to the physical domain (no ghost cells)."""
        ng = self._nghost
        if arr.ndim == 2:
            return arr[ng:-ng, ng:-ng]
        else:
            return arr[ng:-ng, ng:-ng, :]

    @property
    def h(self):
        """Water depth (full domain including ghost cells)."""
        return self._get_2d_field('rki_get_h_ptr')

    @property
    def hu(self):
        """x-momentum (full domain)."""
        return self._get_2d_field('rki_get_hu_ptr')

    @property
    def hv(self):
        """y-momentum (full domain)."""
        return self._get_2d_field('rki_get_hv_ptr')

    @property
    def b(self):
        """Bathymetry (full domain)."""
        return self._get_2d_field('rki_get_b_ptr')

    @property
    def h0(self):
        """RK2 save of h at start of stage 1 (full domain). Valid between
        stage 1 and stage 2 of the current step; overwritten next step."""
        return self._get_2d_field('rki_get_h0_ptr')

    @property
    def hu0(self):
        """RK2 save of hu at start of stage 1 (full domain)."""
        return self._get_2d_field('rki_get_hu0_ptr')

    @property
    def hv0(self):
        """RK2 save of hv at start of stage 1 (full domain)."""
        return self._get_2d_field('rki_get_hv0_ptr')

    @property
    def clamped_elevation(self):
        """Current scalar eta value on the BC_CLAMPED edge (m)."""
        self._check_alive()
        return float(self._lib.rki_get_clamped_elevation(self._handle))

    @property
    def manning_n(self):
        """Manning roughness coefficient (s/m^{1/3})."""
        self._check_alive()
        return float(self._lib.rki_get_manning_n(self._handle))

    @property
    def coriolis_f(self):
        """Coriolis parameter f (1/s)."""
        self._check_alive()
        return float(self._lib.rki_get_coriolis_f(self._handle))

    @property
    def h_phys(self):
        """Water depth (physical domain only, no ghost cells)."""
        return self._phys_slice(self.h)

    @property
    def hu_phys(self):
        """x-momentum (physical domain)."""
        return self._phys_slice(self.hu)

    @property
    def hv_phys(self):
        """y-momentum (physical domain)."""
        return self._phys_slice(self.hv)

    @property
    def b_phys(self):
        """Bathymetry (physical domain)."""
        return self._phys_slice(self.b)

    @property
    def h_layer(self):
        """Layer thickness (3D, physical domain). Requires use_multilayer."""
        return self._phys_slice(self._get_3d_field('rki_get_h_layer_ptr'))

    @property
    def salinity(self):
        """Layer salinity in PSU (3D, physical domain).

        Computed from ``hS_layer / h_layer``.  Cells with effectively
        zero water depth (``h_layer < 1e-6``) are returned as
        ``np.nan`` rather than ``0.0`` so downstream operations like
        ``.min()``, ``.max()``, and matplotlib's ``imshow`` don't
        report or display the dry land cells as if they were
        freshwater.  Use ``np.nanmin``/``np.nanmax`` to skip dry
        cells, or ``np.ma.masked_invalid(sim.salinity)`` if you need
        a masked array.
        """
        hS = self._get_3d_field('rki_get_salinity_ptr')
        h_lay = self._get_3d_field('rki_get_h_layer_ptr')
        ng = self._nghost
        hS_phys = hS[ng:-ng, ng:-ng, :]
        h_phys = h_lay[ng:-ng, ng:-ng, :]
        with np.errstate(invalid='ignore', divide='ignore'):
            S = np.where(h_phys > 1e-6, hS_phys / h_phys, np.nan)
        return S

    @property
    def temperature(self):
        """Layer temperature in degC (3D, physical domain).

        Computed from ``hT_layer / h_layer``.  Dry cells are returned
        as ``np.nan``; use ``np.nanmin``/``np.nanmax`` to skip them.
        """
        hT = self._get_3d_field('rki_get_temperature_ptr')
        h_lay = self._get_3d_field('rki_get_h_layer_ptr')
        ng = self._nghost
        hT_phys = hT[ng:-ng, ng:-ng, :]
        h_phys = h_lay[ng:-ng, ng:-ng, :]
        with np.errstate(invalid='ignore', divide='ignore'):
            T = np.where(h_phys > 1e-6, hT_phys / h_phys, np.nan)
        return T

    # --- Properties ---

    @property
    def time(self):
        """Current simulation time (seconds)."""
        self._check_alive()
        return self._lib.rki_get_time(self._handle)

    @property
    def dt(self):
        """Last timestep used (seconds)."""
        self._check_alive()
        return self._lib.rki_get_dt(self._handle)

    @property
    def step_count(self):
        """Total number of solver steps taken."""
        self._check_alive()
        return self._lib.rki_get_step_count(self._handle)

    @property
    def nz(self):
        """Number of sigma layers (0 if multilayer is inactive)."""
        self._check_alive()
        return self._lib.rki_get_nz_ml(self._handle)

    @property
    def dx(self):
        """Grid spacing in x (metres)."""
        return self._dx

    @property
    def dy(self):
        """Grid spacing in y (metres)."""
        return self._dy

    @property
    def nx(self):
        """Number of physical cells in x."""
        return self._nx_phys

    @property
    def ny(self):
        """Number of physical cells in y."""
        return self._ny_phys

    # --- Methods ---

    def finalize(self):
        """Complete two-phase initialisation.

        Fills bathymetry ghost cells, sets initial conditions (if
        set_initial_h was not called), initialises boundary conditions,
        and uploads state to GPU. Must be called before step/advance.
        """
        self._check_alive()
        if self._finalised:
            raise RuntimeError("Simulation already finalised")
        rc = self._lib.rki_finalize(self._handle)
        if rc != 0:
            raise RuntimeError(f"Finalize failed (rc={rc})")
        self._finalised = True

    def prepare(self):
        """Complete setup and lock in the run (canonical name).

        Identical to :meth:`finalize`; the ``prepare`` spelling matches
        :class:`~rakali.unstr.UnstructuredSimulation` so the two front
        ends read the same (create -> set_* -> prepare -> step/advance).
        """
        self.finalize()

    def step(self, n=1):
        """Take n solver steps.

        The library clamps each step to the configured ``t_end``; once the
        run has reached ``t_end`` further steps would make no progress, so
        this raises instead of silently counting empty steps.
        """
        self._check_ready()
        t0 = self.time
        rc = self._lib.rki_step(self._handle, n)
        if rc != 0:
            raise RuntimeError(f"Solver step failed (rc={rc})")
        if n > 0 and not self.time > t0:
            raise RuntimeError(
                f"step() made no progress at t={self.time:g} s: the run has "
                f"reached its t_end. Increase t_end when creating the "
                f"simulation.")

    def advance(self, target_time, callback=None, callback_interval=None,
                progress=False):
        """Advance simulation to target_time (seconds).

        Parameters
        ----------
        target_time : float
            Simulation time to advance to.
        callback : callable, optional
            Function called periodically with ``callback(sim)``.
            The simulation object is synced to host before each call,
            so fields (h_phys, salinity, etc.) are readable.
        callback_interval : float, optional
            Time between callback (and progress-bar) updates, in
            seconds of sim time. If omitted, defaults to
            ``(target_time - current_time) / 10`` when a ``callback``
            is supplied (10 invocations across the run), or
            ``/100`` when only ``progress=True`` is used (smooth bar).
            Pass an explicit value if you want a smoother bar *and*
            a less frequent callback (e.g. ``callback_interval=1.0``).
        progress : bool, optional
            If True, show a tqdm progress bar covering the time range
            from ``self.time`` to ``target_time``. Falls back silently
            if tqdm isn't installed. Composes with ``callback`` — both
            fire on the same interval.

        Examples
        --------
        >>> sim.advance(3600.0, callback=lambda s: print(f"t={s.time:.0f}s"),
        ...             callback_interval=60.0)

        >>> sim.advance(86400.0, progress=True)   # tqdm bar over 24 h
        """
        self._check_ready()

        # Set up tqdm if requested
        pbar = None
        if progress:
            try:
                from tqdm.auto import tqdm
                pbar = tqdm(
                    total=max(target_time - self.time, 0.0),
                    desc=f"advance to t={target_time:.1f}s",
                    unit='s', unit_scale=False, leave=True,
                )
            except ImportError:
                import warnings
                warnings.warn(
                    "progress=True requires tqdm; running without progress bar. "
                    "Install with `pip install tqdm`.",
                    RuntimeWarning, stacklevel=2,
                )
                pbar = None

        # If neither callback nor progress, use the fast path
        if callback is None and pbar is None:
            rc = self._lib.rki_advance(self._handle, target_time)
            if rc != 0:
                raise RuntimeError(f"Advance failed (rc={rc})")
            return

        # Callback loop: advance in intervals, sync, call back, update bar
        if callback_interval is None:
            if callback is None:
                # Progress bar only — fine intervals for a smooth bar
                callback_interval = max(
                    (target_time - self.time) / 100.0, _TIME_EPS)
            else:
                # Has a user callback — coarser default (backward compat
                # with test 5 and the long-standing 10-callback default).
                callback_interval = max(
                    (target_time - self.time) / 10.0, _TIME_EPS)
        if callback_interval <= 0:
            raise ValueError("callback_interval must be positive")

        try:
            t_now = self.time
            last_t = t_now
            while t_now < target_time - _TIME_EPS:
                t_next = min(t_now + callback_interval, target_time)
                rc = self._lib.rki_advance(self._handle, t_next)
                if rc != 0:
                    raise RuntimeError(f"Advance failed (rc={rc})")
                # Lightweight GPU->host sync (data stays on GPU)
                self._lib.rki_sync_host(self._handle)
                if callback is not None:
                    callback(self)
                t_now = self.time
                if pbar is not None:
                    pbar.update(t_now - last_t)
                    last_t = t_now
        finally:
            if pbar is not None:
                pbar.close()

    def _check_shape(self, arr, label):
        """Validate that ``arr`` matches the ``(nx, ny)`` grid shape."""
        if arr.shape != (self._nx_phys, self._ny_phys):
            raise ValueError(
                f"{label}: expected shape "
                f"({self._nx_phys}, {self._ny_phys}), got {arr.shape}"
            )
        return arr

    def set_bathymetry(self, b):
        """Set bathymetry from a numpy array.

        ``b`` must have shape ``(nx, ny)``.
        """
        self._check_alive()
        b_arr = np.asarray(b, dtype=np.float64)
        b_local = self._check_shape(b_arr, "set_bathymetry")
        b_f = np.asfortranarray(b_local)
        rc = self._lib.rki_set_bathymetry(
            self._handle,
            b_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set bathymetry (rc={rc})")

    def set_initial_h(self, h):
        """Set initial water depth from a numpy array.

        ``h`` must have shape ``(nx, ny)``.
        """
        self._check_alive()
        h_arr = np.asarray(h, dtype=np.float64)
        h_local = self._check_shape(h_arr, "set_initial_h")
        h_f = np.asfortranarray(h_local)
        rc = self._lib.rki_set_initial_h(
            self._handle,
            h_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set initial h (rc={rc})")

    def set_wind(self, wx, wy):
        """Set wind stress from two numpy arrays (nx_phys, ny_phys).

        Both components (Pa) are set together. On GPU builds the device
        copy is updated in-place — no full exit/enter cycle needed.
        Can be called mid-simulation to change forcing.
        """
        self._check_alive()
        wx_f = np.asfortranarray(wx, dtype=np.float64)
        wy_f = np.asfortranarray(wy, dtype=np.float64)
        expected = (self._nx_phys, self._ny_phys)
        if wx_f.shape != expected or wy_f.shape != expected:
            raise ValueError(
                f"Expected shape {expected}, got wx={wx_f.shape}, wy={wy_f.shape}"
            )
        rc = self._lib.rki_set_wind(
            self._handle,
            wx_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            wy_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set wind (rc={rc})")

    def set_atm_pressure(self, p):
        """Set atmospheric pressure from a numpy array (nx_phys, ny_phys).

        Units: Pa. On GPU builds the device copy is updated in-place.
        Can be called mid-simulation.
        """
        self._check_alive()
        p_f = np.asfortranarray(p, dtype=np.float64)
        if p_f.shape != (self._nx_phys, self._ny_phys):
            raise ValueError(
                f"Expected shape ({self._nx_phys}, {self._ny_phys}), "
                f"got {p_f.shape}"
            )
        rc = self._lib.rki_set_atm_pressure(
            self._handle,
            p_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set atm_pressure (rc={rc})")

    def set_precipitation(self, p):
        """Set precipitation rate from a numpy array (nx_phys, ny_phys).

        Units: m/s. On GPU builds the device copy is updated in-place.
        Can be called mid-simulation.
        """
        self._check_alive()
        p_f = np.asfortranarray(p, dtype=np.float64)
        if p_f.shape != (self._nx_phys, self._ny_phys):
            raise ValueError(
                f"Expected shape ({self._nx_phys}, {self._ny_phys}), "
                f"got {p_f.shape}"
            )
        rc = self._lib.rki_set_precipitation(
            self._handle,
            p_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set precipitation (rc={rc})")

    def set_heat_flux(self, q):
        """Set net surface heat flux from a numpy array (nx_phys, ny_phys).

        Units: W/m^2, positive into the ocean.  Applied to the top
        multilayer layer as hT(nz) += dt * Q / (rho_0 * cp). Does not
        change h_layer — companion `set_precipitation` handles volume.
        Multilayer only. On GPU builds the device copy is updated
        in-place. Can be called mid-simulation.
        """
        self._check_alive()
        q_f = np.asfortranarray(q, dtype=np.float64)
        if q_f.shape != (self._nx_phys, self._ny_phys):
            raise ValueError(
                f"Expected shape ({self._nx_phys}, {self._ny_phys}), "
                f"got {q_f.shape}"
            )
        rc = self._lib.rki_set_heat_flux(
            self._handle,
            q_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set heat_flux (rc={rc})")

    def set_evap_minus_precip(self, ep):
        """Set virtual surface freshwater flux from (nx_phys, ny_phys).

        Units: m/s, positive for net evaporation. Applied as a virtual
        salt flux on the top layer: hS(nz) += dt * (E-P) * S_top.
        Does not change h_layer. Multilayer only. Can be called
        mid-simulation.
        """
        self._check_alive()
        ep_f = np.asfortranarray(ep, dtype=np.float64)
        if ep_f.shape != (self._nx_phys, self._ny_phys):
            raise ValueError(
                f"Expected shape ({self._nx_phys}, {self._ny_phys}), "
                f"got {ep_f.shape}"
            )
        rc = self._lib.rki_set_evap_minus_precip(
            self._handle,
            ep_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set evap_minus_precip (rc={rc})")

    # --- Data assimilation setters (mid-simulation state injection) ---

    def _set_2d_field(self, getter_name, label, arr):
        """Copy a (nx_phys, ny_phys) numpy array into a Fortran state slice."""
        self._check_alive()
        a_f = np.asfortranarray(arr, dtype=np.float64)
        if a_f.shape != (self._nx_phys, self._ny_phys):
            raise ValueError(
                f"Expected shape ({self._nx_phys}, {self._ny_phys}), "
                f"got {a_f.shape}"
            )
        rc = getattr(self._lib, getter_name)(
            self._handle,
            a_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set {label} (rc={rc})")

    def set_h(self, h):
        """Set water depth from a numpy array (nx_phys, ny_phys).

        Works mid-simulation: GPU data is refreshed in-place via target
        update. For pre-finalize use, prefer set_initial_h() which carries
        the same semantics but is named to make the IC intent explicit.
        """
        self._set_2d_field('rki_set_h', 'h', h)

    def set_hu(self, hu):
        """Set x-momentum (h*u) from a numpy array (nx_phys, ny_phys)."""
        self._set_2d_field('rki_set_hu', 'hu', hu)

    def set_hv(self, hv):
        """Set y-momentum (h*v) from a numpy array (nx_phys, ny_phys)."""
        self._set_2d_field('rki_set_hv', 'hv', hv)

    def set_salinity(self, S):
        """Set per-layer salinity from a numpy array (nx_phys, ny_phys, nz).

        Units: PSU. Internally we store hS = h_layer * S; the setter
        multiplies by the current layer thickness on the host. Requires
        a multilayer simulation (raises RuntimeError otherwise).
        """
        self._check_alive()
        nz = self.nz
        if nz <= 0:
            raise RuntimeError(
                "set_salinity requires multilayer mode (nz_layers > 0)"
            )
        s_f = np.asfortranarray(S, dtype=np.float64)
        expected = (self._nx_phys, self._ny_phys, nz)
        if s_f.shape != expected:
            raise ValueError(
                f"Expected shape {expected}, got {s_f.shape}"
            )
        rc = self._lib.rki_set_salinity(
            self._handle,
            s_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys, nz,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set salinity (rc={rc})")

    def set_temperature(self, T):
        """Set per-layer temperature from a numpy array (nx_phys, ny_phys, nz).

        Units: degC. Internally we store hT = h_layer * T; the setter
        multiplies by the current layer thickness on the host. Requires
        a multilayer simulation (raises RuntimeError otherwise).
        """
        self._check_alive()
        nz = self.nz
        if nz <= 0:
            raise RuntimeError(
                "set_temperature requires multilayer mode (nz_layers > 0)"
            )
        t_f = np.asfortranarray(T, dtype=np.float64)
        expected = (self._nx_phys, self._ny_phys, nz)
        if t_f.shape != expected:
            raise ValueError(
                f"Expected shape {expected}, got {t_f.shape}"
            )
        rc = self._lib.rki_set_temperature(
            self._handle,
            t_f.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            self._nx_phys, self._ny_phys, nz,
        )
        if rc != 0:
            raise RuntimeError(f"Failed to set temperature (rc={rc})")

    # --- Boundary condition setters (mid-sim forcing updates) ---

    def _set_bc_scalar(self, fn_name, label, value):
        """Call a C scalar BC setter and translate the rc into a Python error."""
        self._check_alive()
        rc = getattr(self._lib, fn_name)(self._handle, float(value))
        if rc != 0:
            raise RuntimeError(f"Failed to set {label} (rc={rc})")

    def set_inflow_velocity(self, v):
        """Set the prescribed normal velocity for ``BC_INFLOW`` edges (m/s).

        Can be called mid-simulation (e.g., from a callback) to drive
        time-varying inflow from observed data. The next solver step
        picks up the new value automatically — no GPU sync needed
        because the boundary kernel re-reads ``state%bc`` into locals
        at each invocation.
        """
        self._set_bc_scalar('rki_set_inflow_velocity', 'inflow_velocity', v)

    def set_inflow_salinity(self, s):
        """Set the prescribed inflow salinity (PSU) for ``BC_INFLOW`` edges.

        Negative values revert to zero-gradient extrapolation. Only
        meaningful in multilayer mode.
        """
        self._set_bc_scalar('rki_set_inflow_salinity', 'inflow_salinity', s)

    def set_inflow_temperature(self, t):
        """Set the prescribed inflow temperature (degC) for ``BC_INFLOW`` edges.

        Negative values revert to zero-gradient extrapolation. Only
        meaningful in multilayer mode.
        """
        self._set_bc_scalar('rki_set_inflow_temperature',
                            'inflow_temperature', t)

    def set_discharge_flux(self, q):
        """Set the prescribed volume flux for ``BC_DISCHARGE`` edges (m^2/s).

        Drives a time-varying river discharge schedule when called from
        a callback loop.
        """
        self._set_bc_scalar('rki_set_discharge_flux', 'discharge_flux', q)

    def set_clamped_elevation(self, eta):
        """Set the prescribed surface elevation for ``BC_CLAMPED`` edges (m)."""
        self._set_bc_scalar('rki_set_clamped_elevation',
                            'clamped_elevation', eta)

    def set_use_ppm_tracer(self, enable: bool = True):
        """Enable/disable high-order tracer advection.

        When True (default), uses PPM (Piecewise Parabolic Method) for
        third-order tracer advection.  When False, falls back to
        first-order upwind.

        Can be called mid-simulation.  Also settable at construction via
        the ``use_ppm_tracer`` keyword argument.
        """
        self._check_alive()
        rc = self._lib.rki_set_use_ppm_tracer(self._handle, int(bool(enable)))
        if rc != 0:
            raise RuntimeError(f"Failed to set use_ppm_tracer (rc={rc})")

    def set_use_hllc(self, enable: bool = True):
        """Switch the structured barotropic Riemann solver.

        When True, the flux kernel uses HLLC (contact-wave upwinding for
        transverse momentum); when False, plain HLL.  Default off.

        Structured-only — the unstructured solver
        (``rakali.unstr.UnstructuredSimulation``) does not expose this
        knob because KNP central-upwind already resolves contacts.

        Can be called mid-simulation.  Also settable at construction via
        the ``use_hllc`` keyword argument.
        """
        self._check_alive()
        rc = self._lib.rki_set_use_hllc(self._handle, int(bool(enable)))
        if rc != 0:
            raise RuntimeError(f"Failed to set use_hllc (rc={rc})")

    def set_tidal_constituents(self, amp, phase, omega):
        """Replace the tidal constituent set for ``BC_TIDAL`` edges.

        Parameters
        ----------
        amp : array-like
            Amplitudes (m), one per constituent.
        phase : array-like
            Phases (radians), one per constituent.
        omega : array-like
            Angular frequencies (rad/s), one per constituent.

        All three arrays must have the same length, between 0 and
        ``MAX_TIDAL_CONSTITUENTS`` (10). Passing zero-length arrays
        disables multi-constituent forcing entirely.

        The kernel evaluates ``eta(t) = sum_k amp_k * cos(omega_k * t - phase_k)``
        on the host before each boundary kernel call.
        """
        self._check_alive()
        amp_a = np.ascontiguousarray(amp, dtype=np.float64)
        phase_a = np.ascontiguousarray(phase, dtype=np.float64)
        omega_a = np.ascontiguousarray(omega, dtype=np.float64)
        if amp_a.ndim != 1 or phase_a.ndim != 1 or omega_a.ndim != 1:
            raise ValueError("amp, phase, omega must be 1-D arrays")
        if not (amp_a.size == phase_a.size == omega_a.size):
            raise ValueError(
                f"amp/phase/omega must have the same length "
                f"(got {amp_a.size}, {phase_a.size}, {omega_a.size})"
            )
        n = amp_a.size
        rc = self._lib.rki_set_tidal_constituents(
            self._handle,
            amp_a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            phase_a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            omega_a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            n,
        )
        if rc != 0:
            if rc == -2:
                raise ValueError(
                    f"Number of constituents out of range "
                    f"(got {n}, max is MAX_TIDAL_CONSTITUENTS=10)"
                )
            raise RuntimeError(f"Failed to set tidal constituents (rc={rc})")

    # --- Restart save / load ---

    def save_restart(self, filename):
        """Write the current state to a NetCDF restart file.

        Captures the barotropic state (h, hu, hv, b, time) and, if the
        simulation is multilayer, also h_layer, hu_layer, hv_layer, and
        hS_layer. The current step count and configuration are NOT
        saved — restoring requires the same namelist/config used to
        create the simulation.
        """
        self._check_ready()
        _require_writable_dir(filename)
        path = str(filename).encode()
        rc = self._lib.rki_save_restart(self._handle, path, len(path))
        if rc != 0:
            raise RuntimeError(f"save_restart failed (rc={rc})")

    def save_checkpoint(self, filename):
        """Write a single-step mid-RK checkpoint as a NetCDF file.

        Snapshots everything needed to replay the next RK stage from
        scratch in the 2D barotropic pathway:

        * Full-grid (ghost-included) arrays: ``h, hu, hv, b`` (current
          state) and ``h0, hu0, hv0`` (RK2 saves written by stage 1).
        * Scalars: ``time, step_count, dt, manning_n, coriolis_f``,
          grid metadata (``nx_total, ny_total, nghost, dx, dy``), and
          the south-edge ``clamped_elevation``.

        Multilayer / NH state is **not** saved by this method — it is
        intended for the 2D path used by `python/examples/real_bathymetry/lakes_entrance`.
        Use ``save_restart`` for the multilayer-aware NetCDF restart.
        """
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "save_checkpoint() requires the xarray package."
            ) from e

        self._check_ready()
        self.sync_host()

        nx_total = self._nx_phys + 2 * self._nghost
        ny_total = self._ny_phys + 2 * self._nghost

        # Copies, not zero-copy views — the caller may keep stepping.
        h   = self.h.copy()
        hu  = self.hu.copy()
        hv  = self.hv.copy()
        b   = self.b.copy()
        h0  = self.h0.copy()
        hu0 = self.hu0.copy()
        hv0 = self.hv0.copy()

        dims = ('x_total', 'y_total')
        data_vars = {
            'h':   (dims, h,   {'units': 'm',     'long_name': 'water depth (ghost-included)'}),
            'hu':  (dims, hu,  {'units': 'm2 s-1','long_name': 'x-momentum (ghost-included)'}),
            'hv':  (dims, hv,  {'units': 'm2 s-1','long_name': 'y-momentum (ghost-included)'}),
            'b':   (dims, b,   {'units': 'm',     'long_name': 'bottom elevation',
                                'positive': 'up'}),
            'h0':  (dims, h0,  {'units': 'm',     'long_name': 'RK2 save of h at start of stage 1'}),
            'hu0': (dims, hu0, {'units': 'm2 s-1','long_name': 'RK2 save of hu at start of stage 1'}),
            'hv0': (dims, hv0, {'units': 'm2 s-1','long_name': 'RK2 save of hv at start of stage 1'}),
        }

        attrs = {
            'title': 'Rakali per-step checkpoint (2D barotropic)',
            'source': 'Rakali Python API (rakali.Rakali.save_checkpoint)',
            'time': float(self.time),
            'step_count': int(self.step_count),
            'dt': float(self.dt),
            'nx_total': int(nx_total),
            'ny_total': int(ny_total),
            'nx_phys': int(self._nx_phys),
            'ny_phys': int(self._ny_phys),
            'nghost': int(self._nghost),
            'dx': float(self._dx),
            'dy': float(self._dy),
            'manning_n': float(self.manning_n),
            'coriolis_f': float(self.coriolis_f),
            'clamped_elevation': float(self.clamped_elevation),
        }

        ds = xr.Dataset(data_vars=data_vars, attrs=attrs)
        enc = {v: {'zlib': True, 'complevel': 1} for v in ds.data_vars}
        ds.to_netcdf(str(filename), encoding=enc)

    @classmethod
    def from_restart(cls, restart_file, namelist=None, **kwargs):
        """Create a Rakali simulation by loading state from a restart file.

        The namelist (or kwargs) provides the configuration — grid,
        boundary conditions, physics parameters — exactly as for a fresh
        Rakali() call. The restart file provides the state arrays
        (h, hu, hv, b, time, plus multilayer fields if applicable).
        The grid in the restart file must match the configured grid.

        Equivalent to creating the simulation with ``defer_init=True``,
        loading the restart state, then calling ``finalize()``.

        Parameters
        ----------
        restart_file : str
            Path to a NetCDF restart file written by ``save_restart``
            (or by the Fortran ``restart_write`` routine).
        namelist : str, optional
            Path to a Fortran namelist file (mutually exclusive with kwargs).
        **kwargs
            Namelist parameters as keyword arguments — same set as the
            ``Rakali()`` constructor.

        Examples
        --------
        >>> sim = Rakali("config.nml")
        >>> sim.advance(3600.0)
        >>> sim.save_restart("checkpoint.nc")
        >>> sim.destroy()
        ...
        >>> sim = Rakali.from_restart("checkpoint.nc", "config.nml")
        >>> sim.time
        3600.0
        """
        # Resolve the restart path before any cwd changes inside __init__
        _require_file(restart_file)
        restart_path = os.path.abspath(str(restart_file))
        # from_restart always uses two-phase init internally
        kwargs.pop('defer_init', None)
        sim = cls(namelist=namelist, defer_init=True, **kwargs)
        try:
            path_b = restart_path.encode()
            rc = sim._lib.rki_load_restart(sim._handle, path_b, len(path_b))
            if rc != 0:
                raise RuntimeError(
                    f"rki_load_restart('{restart_path}') failed (rc={rc})"
                )
            sim.finalize()
        except Exception:
            sim.destroy()
            raise
        return sim

    # --- xarray export ---

    def to_xarray(self, time_dim=True):
        """Return the current simulation state as an ``xarray.Dataset``.

        Builds a snapshot at the current ``sim.time`` with cell-centred
        spatial coordinates and a length-1 time dimension (so multiple
        snapshots can be stacked with ``xr.concat(..., dim='time')``).

        Always-on variables (CF-style ``units`` and ``long_name`` attrs):

        =========  ========  ===============================================
        Variable   Units     Description
        =========  ========  ===============================================
        ``h``      m         water depth (0 on dry cells)
        ``b``      m         bathymetry (positive up)
        ``eta``    m         surface elevation = h + b
        ``hu``     m^2 s^-1  x-momentum (0 on dry cells)
        ``hv``     m^2 s^-1  y-momentum (0 on dry cells)
        ``u``      m s^-1    x-velocity (hu/h, NaN on dry cells)
        ``v``      m s^-1    y-velocity (hv/h, NaN on dry cells)
        =========  ========  ===============================================

        Multilayer mode adds 3-D variables on a ``z`` layer-index dimension
        (bottom-up: ``z=1`` is the bed layer, ``z=nz`` the surface layer): ``h_layer`` (m), ``u_layer`` (m/s, NaN on dry),
        ``v_layer`` (m/s, NaN on dry), ``salinity`` (PSU, NaN on dry).

        Derived quantities (``u``, ``v``, ``u_layer``, ``v_layer``,
        ``salinity``) use ``np.nan`` as the dry-cell sentinel rather
        than ``0.0``.  NaN propagates naturally through ``.min()``,
        ``.max()``, ``.std()``, matplotlib, and xarray, so downstream
        analysis won't accidentally treat the land strip as
        zero-velocity / fresh water.  Use ``np.nanmin``/``np.nanmax``
        if you need explicit NaN-skipping.

        Parameters
        ----------
        time_dim : bool, optional
            If True (default), include a length-1 ``time`` dimension and
            coordinate.  Set to False to drop it (useful when you only
            ever care about a single snapshot and don't plan to concat).

        Returns
        -------
        xarray.Dataset
            Snapshot of the simulation at the current time.  Numpy
            arrays are *copied* from the Fortran-side buffers — the
            returned dataset is independent of the live simulation
            state and is safe to keep around after further stepping.

        Raises
        ------
        ImportError
            If xarray is not installed in the Python environment.

        Examples
        --------
        Single snapshot::

            sim = Rakali("config.nml")
            sim.advance(3600.0)
            ds = sim.to_xarray()
            ds.h.plot()                # via xarray's matplotlib hook

        Time series via callback::

            snapshots = []
            sim.advance(7200.0,
                        callback=lambda s: snapshots.append(s.to_xarray()),
                        callback_interval=600.0)
            ds = xr.concat(snapshots, dim='time')
            ds.salinity.isel(z=0).plot()
        """
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "to_xarray() requires the xarray package. "
                "Install with `pip install xarray`."
            ) from e

        self._check_ready()
        # Lightweight GPU->host sync so the property accessors below
        # don't trigger the heavy ensure_host (exit_data) path.
        self.sync_host()

        nx, ny = self._nx_phys, self._ny_phys
        nz = self.nz
        dx, dy = self._dx, self._dy

        x_coord = (np.arange(nx) + 0.5) * dx
        y_coord = (np.arange(ny) + 0.5) * dy

        coords = {
            'x': ('x', x_coord, {
                'units': 'm',
                'long_name': 'distance east of west boundary',
                'axis': 'X',
            }),
            'y': ('y', y_coord, {
                'units': 'm',
                'long_name': 'distance north of south boundary',
                'axis': 'Y',
            }),
        }
        if nz > 0:
            coords['z'] = ('z', np.arange(1, nz + 1), {
                'long_name': 'layer index (1 = bed, nz = surface)',
                'positive': 'up',
            })
        if time_dim:
            coords['time'] = ('time', np.array([self.time]), {
                'units': 'seconds since simulation start',
                'long_name': 'simulation time',
                'axis': 'T',
            })

        # Snapshot the barotropic fields (copies, not views — safe to
        # keep across further stepping)
        h = self.h_phys.copy()
        hu = self.hu_phys.copy()
        hv = self.hv_phys.copy()
        b = self.b_phys.copy()

        # Derived velocities: NaN on dry cells so downstream analysis
        # treats land as missing data, not "still water moving at 0 m/s".
        with np.errstate(invalid='ignore', divide='ignore'):
            u = np.where(h > 1e-6, hu / h, np.nan)
            v = np.where(h > 1e-6, hv / h, np.nan)
        eta = h + b

        def _wrap2d(arr):
            if time_dim:
                return (('time', 'x', 'y'), arr[np.newaxis, :, :])
            return (('x', 'y'), arr)

        def _wrap3d(arr):
            if time_dim:
                return (('time', 'x', 'y', 'z'), arr[np.newaxis, :, :, :])
            return (('x', 'y', 'z'), arr)

        data_vars = {
            'h':   (*_wrap2d(h),   {'units': 'm',
                                     'long_name': 'water depth'}),
            'b':   (*_wrap2d(b),   {'units': 'm',
                                     'long_name': 'bottom elevation',
                                     'positive': 'up'}),
            'eta': (*_wrap2d(eta), {'units': 'm',
                                     'long_name': 'surface elevation (h + b)'}),
            'hu':  (*_wrap2d(hu),  {'units': 'm2 s-1',
                                     'long_name': 'x-momentum (h*u)'}),
            'hv':  (*_wrap2d(hv),  {'units': 'm2 s-1',
                                     'long_name': 'y-momentum (h*v)'}),
            'u':   (*_wrap2d(u),   {'units': 'm s-1',
                                     'long_name': 'x-velocity'}),
            'v':   (*_wrap2d(v),   {'units': 'm s-1',
                                     'long_name': 'y-velocity'}),
        }

        if nz > 0:
            h_lay = self.h_layer.copy()
            sal = self.salinity.copy()      # already NaN-on-dry from the property
            tmp = self.temperature.copy()   # already NaN-on-dry from the property

            # Per-layer momentum: read raw, slice physical, divide by h_layer.
            # NaN on dry cells, same convention as the barotropic u, v above.
            hu_lay_full = self._get_3d_field('rki_get_hu_layer_ptr')
            hv_lay_full = self._get_3d_field('rki_get_hv_layer_ptr')
            ng = self._nghost
            hu_lay = hu_lay_full[ng:-ng, ng:-ng, :].copy()
            hv_lay = hv_lay_full[ng:-ng, ng:-ng, :].copy()
            with np.errstate(invalid='ignore', divide='ignore'):
                u_lay = np.where(h_lay > 1e-6, hu_lay / h_lay, np.nan)
                v_lay = np.where(h_lay > 1e-6, hv_lay / h_lay, np.nan)

            # Density (two-tracer EOS) is intentionally NOT exposed
            # here — rho_0, beta_S, alpha_T aren't in the C API
            # surface, and recomputing from salinity + temperature
            # downstream is trivial.

            data_vars['h_layer'] = (*_wrap3d(h_lay), {
                'units': 'm', 'long_name': 'layer thickness'})
            data_vars['u_layer'] = (*_wrap3d(u_lay), {
                'units': 'm s-1', 'long_name': 'layer x-velocity'})
            data_vars['v_layer'] = (*_wrap3d(v_lay), {
                'units': 'm s-1', 'long_name': 'layer y-velocity'})
            data_vars['salinity'] = (*_wrap3d(sal), {
                'units': 'PSU', 'long_name': 'salinity'})
            data_vars['temperature'] = (*_wrap3d(tmp), {
                'units': 'degC', 'long_name': 'temperature'})

        attrs = {
            'title': 'Rakali simulation snapshot',
            'source': 'Rakali Python API (rakali.Rakali.to_xarray)',
            'Conventions': 'CF-1.11',
            'dx': float(dx),
            'dy': float(dy),
            'nx': int(nx),
            'ny': int(ny),
            'time': float(self.time),
            'step_count': int(self.step_count),
        }
        if nz > 0:
            attrs['nz'] = int(nz)

        return xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)

    def sync_host(self):
        """Lightweight GPU-to-host sync that keeps data on GPU.

        Copies solution fields to host so they can be read via
        properties (h_phys, salinity, etc.) without removing them
        from GPU memory. The next solver step proceeds without
        re-uploading. Preferred over sync_from_device() for
        mid-simulation diagnostics.
        """
        self._check_alive()
        self._lib.rki_sync_host(self._handle)

    def sync_from_device(self):
        """Full GPU-to-host sync that removes data from GPU.

        Use at end of simulation or when GPU memory must be freed.
        For mid-simulation callbacks, prefer sync_host() instead.
        """
        self._check_alive()
        self._lib.rki_sync_from_device(self._handle)

    def _destroy_handle(self):
        """Free the C handle (idempotent via `_alive`).

        Teardown entry points — `close()`, `destroy()`, `with` exit and
        `__del__` — come from the shared :class:`_SimBase` mixin so the
        structured and unstructured wrappers behave identically.
        """
        if getattr(self, '_alive', False):
            self._lib.rki_destroy(ctypes.byref(self._handle))
            self._alive = False

    def __repr__(self):
        if not self._alive:
            return "Rakali(destroyed)"
        nz = self.nz
        ml = f", nz={nz}" if nz > 0 else ""
        return (
            f"Rakali({self._nx_phys}x{self._ny_phys}, "
            f"dx={self._dx}m{ml}, t={self.time:.1f}s)"
        )

    def _repr_html_(self):
        """Rich HTML representation for Jupyter notebooks."""
        if not self._alive:
            return "<pre>Rakali(destroyed)</pre>"

        # Try to read live field ranges, but don't crash the repr if
        # the simulation isn't finalised yet.
        try:
            self._check_ready()
            self.sync_host()
            h = self.h_phys
            h_range = f"{h.min():.4f} to {h.max():.4f}"
            nz = self.nz
            sal_range = ""
            temp_range = ""
            if nz > 0:
                S = self.salinity
                sal_range = f"{S.min():.3f} to {S.max():.3f}"
                T = self.temperature
                temp_range = f"{T.min():.3f} to {T.max():.3f}"
            field_rows = [
                ("h (m)", h_range),
            ]
            if nz > 0:
                field_rows.append(("salinity (PSU)", sal_range))
                field_rows.append(("temperature (degC)", temp_range))
        except Exception:
            field_rows = []
            nz = self.nz

        rows = [
            ("Grid", f"{self._nx_phys} x {self._ny_phys} cells"),
            ("Spacing", f"dx = {self._dx} m, dy = {self._dy} m"),
            ("Domain", f"{self._nx_phys * self._dx:.0f} m "
                       f"x {self._ny_phys * self._dy:.0f} m"),
            ("Layers", str(nz) if nz > 0 else "barotropic"),
            ("Time", f"{self.time:.2f} s"),
            ("Steps", str(self.step_count)),
            ("dt (last)", f"{self.dt:.4g} s"),
        ]
        rows.extend(field_rows)

        body = "".join(
            f"<tr><td style='text-align:right;padding-right:1em;color:#666'>"
            f"<b>{k}</b></td><td>{v}</td></tr>"
            for k, v in rows
        )
        return (
            "<div style='font-family:monospace'>"
            "<b>Rakali simulation</b>"
            f"<table style='border-collapse:collapse;margin-top:0.3em'>{body}</table>"
            "</div>"
        )
