"""Tests for ``rakali.yaml_loader``.

Pure-Python tests over the example YAML files in ``python/examples/``
plus a handful of synthetic snippets (written via ``tmp_path``) for the
error paths. The loader doesn't touch the solver singleton, so these
tests can run alongside the rest of the suite without contention.
"""
from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import pytest

from rakali.yaml_loader import (
    Bathymetry,
    BoundaryConditions,
    Bowl,
    Coastline,
    Config,
    DamBreak,
    Domain,
    Forcing,
    Funnel,
    GaussianBump,
    InitialCondition,
    Island,
    Mode,
    Physics,
    Salinity,
    Solitary,
    Tide,
    Tilt,
    load_config,
    to_rakali_kwargs,
)


EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


# ---------------------------------------------------------------------------
# Round-trip: every shipped example YAML must load
# ---------------------------------------------------------------------------

EXAMPLE_FILES = sorted(EXAMPLES_DIR.glob("*.yaml"))


@pytest.mark.parametrize("path", EXAMPLE_FILES, ids=lambda p: p.name)
def test_examples_load(path):
    """Every shipped example YAML parses without error and produces a Config."""
    cfg = load_config(path)
    assert isinstance(cfg, Config)
    assert cfg.domain.type in ("structured", "unstructured")


def test_examples_dir_is_not_empty():
    """Guards against the parametrise above silently producing zero tests."""
    assert EXAMPLE_FILES, f"no YAML files found in {EXAMPLES_DIR}"


# ---------------------------------------------------------------------------
# Spot-checks against specific examples
# ---------------------------------------------------------------------------

def test_structured_dam_break_values():
    cfg = load_config(EXAMPLES_DIR / "live_dam_break.yaml")

    assert cfg.domain == Domain(
        type="structured", nx=200, ny=100, dx=5.0, dy=5.0, t_end=6000.0,
    )
    assert cfg.bathymetry == Bathymetry(shape="flat", depth=5.0)

    assert cfg.initial_condition.type == "dam_break"
    assert cfg.initial_condition.dam_break == DamBreak(
        geometry="axis", axis="x", split_fraction=0.5, h_high=5.0, h_low=1.0,
    )

    assert cfg.physics.cfl == 0.45
    assert cfg.physics.manning_n == 0.025
    assert cfg.physics.smagorinsky_cs is None  # not set in this example

    assert cfg.boundary_conditions.west == "wall"
    assert cfg.boundary_conditions.east == "wall"
    assert cfg.boundary_conditions.south == "wall"
    assert cfg.boundary_conditions.north == "wall"
    assert cfg.boundary_conditions.groups == {}

    # Only the mode fields actually present in the YAML propagate;
    # use_nonhydrostatic isn't set here, so it remains None.
    assert cfg.mode == Mode(
        use_multilayer=False, use_nonhydrostatic=None,
        log_level="warning", defer_init=True,
    )
    # Sections that don't apply to this example are absent.
    assert cfg.forcing is None
    assert cfg.salinity is None


def test_kelvin_wave_gaussian_bump():
    cfg = load_config(EXAMPLES_DIR / "live_kelvin_wave.yaml")

    assert cfg.initial_condition.type == "gaussian_bump"
    assert cfg.initial_condition.bump == GaussianBump(
        amplitude=1.0, x0_fraction=0.15, sigma_x=30000.0, decay_from="south",
    )
    assert cfg.forcing.coriolis_f == 1.0e-4
    assert cfg.physics.smagorinsky_cs == 0.05


def test_lock_exchange_salinity_dam_break():
    """Lock exchange has type=still for mass; the dam-break is in salinity."""
    cfg = load_config(EXAMPLES_DIR / "live_lock_exchange.yaml")

    assert cfg.domain.nz == 12
    assert cfg.initial_condition.type == "still"
    assert cfg.initial_condition.dam_break is None
    assert cfg.salinity == Salinity(initial=33.0, s_left=36.0, s_right=30.0)
    assert cfg.mode.use_multilayer is True
    assert cfg.physics.bpg_method == "jacobian"


def test_river_plume_coastline_subblock():
    cfg = load_config(EXAMPLES_DIR / "live_river_plume.yaml")

    assert cfg.bathymetry.shape == "flat"
    assert cfg.bathymetry.land_height == 1.5
    assert cfg.bathymetry.coastline == Coastline(
        west_strip_cells=6,
        river_mouth_j_fraction=0.65,
        river_mouth_halfwidth_cells=5,
    )
    assert cfg.forcing == Forcing(
        inflow_velocity=0.6, inflow_salinity=2.0, coriolis_f=4.0e-4, tide=None,
    )


def test_tidal_estuary_funnel_and_tide():
    cfg = load_config(EXAMPLES_DIR / "live_tidal_estuary.yaml")

    assert cfg.bathymetry.shape == "funnel"
    assert cfg.bathymetry.depth is None  # funnel uses depth_head/_mouth instead
    assert cfg.bathymetry.funnel == Funnel(
        depth_head=5.0, depth_mouth=22.0,
        width_head_cells=8, width_mouth_cells=38,
    )
    assert cfg.forcing.tide == Tide(amplitude=1.5, period=44712.0)
    assert cfg.boundary_conditions.east == "tidal"


def test_solitary_wave_nonhydrostatic():
    cfg = load_config(EXAMPLES_DIR / "live_solitary_wave.yaml")

    assert cfg.initial_condition.type == "solitary"
    assert cfg.initial_condition.solitary == Solitary(amplitude=1.0, x0=80.0)
    assert cfg.mode.use_nonhydrostatic is True
    assert cfg.mode.use_multilayer is False
    assert cfg.domain.nz == 6


def test_vorticity_island_subblock():
    cfg = load_config(EXAMPLES_DIR / "live_vorticity_2d.yaml")

    assert cfg.bathymetry.island == Island(
        cx_fraction=0.25, cy_fraction=0.5, radius_cells=12,
    )
    assert cfg.forcing.coriolis_f == 0.0


# ---------------------------------------------------------------------------
# Unstructured examples
# ---------------------------------------------------------------------------

def test_unstructured_disk_radial_dam_break():
    cfg = load_config(EXAMPLES_DIR / "circular_dam_break_unstr.yaml")

    assert cfg.domain == Domain(
        type="unstructured", mesh_shape="disk",
        radius=20.0, n_across=80, n_boundary=128, t_end=1.0,
    )
    assert cfg.initial_condition.dam_break == DamBreak(
        geometry="radial", centre=[0.0, 0.0], radius=2.5,
        h_high=2.0, h_low=0.5,
    )
    assert cfg.boundary_conditions.groups == {1: "wall"}
    # Structured-only fields are unset.
    assert cfg.boundary_conditions.west is None


def test_unstructured_rectangle_two_groups():
    cfg = load_config(EXAMPLES_DIR / "dam_break_unstr.yaml")

    assert cfg.domain.mesh_shape == "rectangle"
    assert cfg.domain.length == 100.0
    assert cfg.domain.width == 2.0
    assert cfg.boundary_conditions.groups == {1: "wall", 2: "open"}


def test_thacker_bowl_tilted_ic():
    cfg = load_config(EXAMPLES_DIR / "thacker_bowl_unstr.yaml")

    assert cfg.bathymetry.shape == "bowl"
    assert cfg.bathymetry.bowl == Bowl(L=5.0, h0=1.0)
    assert cfg.initial_condition.type == "tilted"
    assert cfg.initial_condition.tilt == Tilt(z_lake=-0.3, amplitude=0.1)
    assert cfg.domain.t_end is None  # script computes it from theory


# ---------------------------------------------------------------------------
# Defaults: missing optional sections fall through to None / defaults
# ---------------------------------------------------------------------------

def _write(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(dedent(body))
    return path


def test_minimal_config(tmp_path):
    """Only ``domain`` is required; everything else is optional."""
    path = _write(tmp_path, """
        domain:
          type: structured
          nx: 10
          ny: 10
          dx: 1.0
          dy: 1.0
          t_end: 1.0
    """)
    cfg = load_config(path)
    assert cfg.domain.nx == 10
    assert cfg.bathymetry is None
    assert cfg.initial_condition is None
    assert cfg.forcing is None
    assert cfg.physics is None
    assert cfg.boundary_conditions is None
    assert cfg.salinity is None
    assert cfg.mode is None


def test_mode_unset_fields_stay_none(tmp_path):
    """Fields omitted from a present block stay None — they must not be
    forced into ``to_rakali_kwargs`` and must not override Rakali's
    own defaults."""
    path = _write(tmp_path, """
        domain:
          type: structured
        mode:
          use_multilayer: true
    """)
    cfg = load_config(path)
    assert cfg.mode == Mode(
        use_multilayer=True, use_nonhydrostatic=None,
        log_level=None, defer_init=None,
    )

    # Crucially, only the explicitly-set field survives the kwargs adapter.
    kw = to_rakali_kwargs(cfg)
    assert kw["use_multilayer"] is True
    assert "use_nonhydrostatic" not in kw
    assert "log_level" not in kw
    assert "defer_init" not in kw


# ---------------------------------------------------------------------------
# Boundary conditions: structured edges, group_N tags, and mixing
# ---------------------------------------------------------------------------

def test_bc_partial_structured(tmp_path):
    path = _write(tmp_path, """
        domain: {type: structured}
        boundary_conditions:
          west: inflow
          east: open
    """)
    cfg = load_config(path)
    assert cfg.boundary_conditions == BoundaryConditions(
        west="inflow", east="open",
    )


def test_bc_groups_only(tmp_path):
    path = _write(tmp_path, """
        domain: {type: unstructured}
        boundary_conditions:
          groups:
            1: wall
            3: open
    """)
    cfg = load_config(path)
    assert cfg.boundary_conditions.groups == {1: "wall", 3: "open"}
    assert cfg.boundary_conditions.west is None


def test_bc_groups_rejected_for_structured(tmp_path):
    """Structured domains must use named edges, not tag groups."""
    path = _write(tmp_path, """
        domain: {type: structured}
        boundary_conditions:
          west: wall
          groups: {2: open}
    """)
    with pytest.raises(ValueError,
                       match="groups is invalid for structured domains"):
        load_config(path)


def test_bc_edges_rejected_for_unstructured(tmp_path):
    """Unstructured domains must use tag groups, not named edges."""
    path = _write(tmp_path, """
        domain: {type: unstructured}
        boundary_conditions:
          west: wall
    """)
    with pytest.raises(ValueError,
                       match="edges .* invalid for unstructured domains"):
        load_config(path)


def test_use_hllc_rejected_for_unstructured(tmp_path):
    """`use_hllc` is structured-only — the unstructured solver always
    uses KNP central-upwind, so accepting the flag silently would be
    misleading."""
    path = _write(tmp_path, """
        domain: {type: unstructured, mesh_shape: rectangle, length: 1.0, width: 1.0, n_across: 4}
        physics: {use_hllc: true}
    """)
    with pytest.raises(ValueError,
                       match="use_hllc is a structured-grid-only knob"):
        load_config(path)


def test_use_hllc_accepted_for_structured(tmp_path):
    """Sanity check: physics.use_hllc round-trips through the loader for
    structured domains."""
    path = _write(tmp_path, """
        domain: {type: structured, nx: 10, ny: 10, dx: 1.0, dy: 1.0, t_end: 1.0}
        physics: {use_hllc: true}
    """)
    cfg = load_config(path)
    assert cfg.physics.use_hllc is True


def test_bc_groups_keys_must_be_int(tmp_path):
    path = _write(tmp_path, """
        domain: {type: unstructured}
        boundary_conditions:
          groups:
            'one': wall
    """)
    with pytest.raises(ValueError, match="groups keys must be integers"):
        load_config(path)


def test_bc_unknown_key(tmp_path):
    """Unknown BC keys are caught by the standard ``_build`` unknown-key
    check (no more bespoke parser)."""
    path = _write(tmp_path, """
        domain: {type: structured}
        boundary_conditions:
          northwest: wall
    """)
    with pytest.raises(
        ValueError, match="unknown key.*'northwest'.*BoundaryConditions"
    ):
        load_config(path)


# ---------------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------------

def test_missing_domain_block(tmp_path):
    path = _write(tmp_path, """
        physics: {cfl: 0.4}
    """)
    with pytest.raises(ValueError, match="missing required 'domain' block"):
        load_config(path)


def test_top_level_not_a_mapping(tmp_path):
    path = _write(tmp_path, """
        - just
        - a
        - list
    """)
    with pytest.raises(ValueError, match="top-level YAML must be a mapping"):
        load_config(path)


def test_unknown_top_level_key(tmp_path):
    path = _write(tmp_path, """
        domain: {type: structured}
        nonsense: 42
    """)
    with pytest.raises(ValueError, match="unknown key.*'nonsense'"):
        load_config(path)


def test_unknown_nested_key(tmp_path):
    path = _write(tmp_path, """
        domain:
          type: structured
          spurious_field: 1
    """)
    with pytest.raises(ValueError, match="unknown key.*'spurious_field'.*Domain"):
        load_config(path)


def test_unknown_subblock_key(tmp_path):
    path = _write(tmp_path, """
        domain: {type: structured}
        bathymetry:
          shape: flat
          funnel: {depth_head: 1.0, depth_mouth: 2.0, width_head_cells: 1, width_mouth_cells: 2, extra: 3}
    """)
    with pytest.raises(ValueError, match="'extra'.*Funnel"):
        load_config(path)


def test_block_must_be_mapping(tmp_path):
    path = _write(tmp_path, """
        domain: [not, a, mapping]
    """)
    with pytest.raises(ValueError, match="expected dict for Domain"):
        load_config(path)


# ---------------------------------------------------------------------------
# load_config accepts both str and Path
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# to_rakali_kwargs adapter
# ---------------------------------------------------------------------------

def test_rakali_kwargs_dam_break():
    cfg = load_config(EXAMPLES_DIR / "live_dam_break.yaml")
    kw = to_rakali_kwargs(cfg)

    assert kw["nx"] == 200 and kw["ny"] == 100
    assert kw["dx"] == 5.0 and kw["dy"] == 5.0
    assert kw["t_end"] == 6000.0
    assert kw["cfl"] == 0.45
    assert kw["manning_n"] == 0.025
    assert kw["bc_west"] == "wall"
    assert kw["use_multilayer"] is False
    assert kw["log_level"] == "warning"
    assert kw["defer_init"] is True
    # Fields not set in this YAML must NOT appear in kwargs.
    assert "nz_layers" not in kw
    assert "smagorinsky_cs" not in kw
    assert "inflow_velocity" not in kw


def test_rakali_kwargs_river_plume_full_set():
    cfg = load_config(EXAMPLES_DIR / "live_river_plume.yaml")
    kw = to_rakali_kwargs(cfg)

    assert kw["nz_layers"] == 3
    assert kw["inflow_velocity"] == 0.6
    assert kw["inflow_salinity"] == 2.0
    assert kw["coriolis_f"] == 4.0e-4
    assert kw["initial_salinity"] == 35.0
    assert kw["bpg_method"] == "jacobian"
    assert kw["bc_west"] == "inflow"
    assert kw["bc_east"] == "open"
    assert kw["use_multilayer"] is True
    # No manning_n in this example's YAML, so it must not be forced into kwargs.
    assert "manning_n" not in kw


def test_rakali_kwargs_rejects_unstructured():
    cfg = load_config(EXAMPLES_DIR / "circular_dam_break_unstr.yaml")
    with pytest.raises(ValueError, match="domain.type='structured'"):
        to_rakali_kwargs(cfg)


def test_load_config_accepts_string_path():
    path = EXAMPLES_DIR / "live_dam_break.yaml"
    cfg_from_str = load_config(str(path))
    cfg_from_path = load_config(path)
    assert cfg_from_str.domain == cfg_from_path.domain
