"""YAML config loader for the example simulations.

Reads the unified-schema YAML files in ``python/examples/*.yaml`` and
returns a typed ``Config`` object. The schema is documented in those
example files; in short:

    domain:               type (structured|unstructured) + grid/mesh fields
    bathymetry:           shape (flat|funnel|bowl) + sub-blocks
    initial_condition:    type (still|dam_break|gaussian_bump|solitary|tilted)
    forcing:              inflow_velocity, inflow_salinity, coriolis_f, tide
    physics:              cfl, manning_n, smagorinsky_cs, hdiff_kappa, bpg_method, use_hllc, time_scheme
    boundary_conditions:  west/east/south/north OR group_N
    salinity:             initial, s_left, s_right (multilayer only)
    mode:                 use_multilayer, use_nonhydrostatic, log_level, defer_init

Usage::

    from rakali.yaml_loader import load_config
    cfg = load_config("python/examples/live_dam_break.yaml")
    print(cfg.domain.nx, cfg.physics.cfl)
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any, Optional, Union

import yaml


# ---------------------------------------------------------------------------
# domain
# ---------------------------------------------------------------------------

@dataclass
class Domain:
    type: str  # "structured" | "unstructured"
    # structured grid:
    nx: Optional[int] = None
    ny: Optional[int] = None
    dx: Optional[float] = None
    dy: Optional[float] = None
    # unstructured mesh:
    mesh_shape: Optional[str] = None  # "disk" | "rectangle"
    radius: Optional[float] = None
    length: Optional[float] = None
    width: Optional[float] = None
    n_across: Optional[int] = None
    n_boundary: Optional[int] = None
    # both:
    nz: Optional[int] = None
    t_end: Optional[float] = None


# ---------------------------------------------------------------------------
# bathymetry
# ---------------------------------------------------------------------------

@dataclass
class Funnel:
    depth_head: float
    depth_mouth: float
    width_head_cells: int
    width_mouth_cells: int


@dataclass
class Bowl:
    L: float
    h0: float


@dataclass
class Island:
    cx_fraction: float
    cy_fraction: float
    radius_cells: int


@dataclass
class Coastline:
    west_strip_cells: int
    river_mouth_j_fraction: float
    river_mouth_halfwidth_cells: int


@dataclass
class Bathymetry:
    shape: str  # "flat" | "funnel" | "bowl"
    depth: Optional[float] = None
    land_height: Optional[float] = None
    funnel: Optional[Funnel] = None
    bowl: Optional[Bowl] = None
    island: Optional[Island] = None
    coastline: Optional[Coastline] = None


# ---------------------------------------------------------------------------
# initial condition
# ---------------------------------------------------------------------------

@dataclass
class DamBreak:
    geometry: str  # "axis" | "radial"
    axis: Optional[str] = None
    split_fraction: Optional[float] = None
    centre: Optional[list] = None
    radius: Optional[float] = None
    h_high: Optional[float] = None
    h_low: Optional[float] = None


@dataclass
class GaussianBump:
    amplitude: float
    x0_fraction: float
    sigma_x: float
    decay_from: str


@dataclass
class Solitary:
    amplitude: float
    x0: float


@dataclass
class Tilt:
    z_lake: float
    amplitude: float


@dataclass
class InitialCondition:
    type: str  # "still" | "dam_break" | "gaussian_bump" | "solitary" | "tilted"
    dam_break: Optional[DamBreak] = None
    bump: Optional[GaussianBump] = None
    solitary: Optional[Solitary] = None
    tilt: Optional[Tilt] = None


# ---------------------------------------------------------------------------
# forcing
# ---------------------------------------------------------------------------

@dataclass
class Tide:
    amplitude: float
    period: float


@dataclass
class Forcing:
    inflow_velocity: Optional[float] = None
    inflow_salinity: Optional[float] = None
    coriolis_f: Optional[float] = None
    tide: Optional[Tide] = None


# ---------------------------------------------------------------------------
# physics, boundaries, salinity, mode
# ---------------------------------------------------------------------------

@dataclass
class Physics:
    cfl: Optional[float] = None
    manning_n: Optional[float] = None
    smagorinsky_cs: Optional[float] = None
    hdiff_kappa: Optional[float] = None
    bpg_method: Optional[str] = None
    use_hllc: Optional[bool] = None
    time_scheme: Optional[str] = None


@dataclass
class BoundaryConditions:
    # structured (named edges):
    west: Optional[str] = None
    east: Optional[str] = None
    south: Optional[str] = None
    north: Optional[str] = None
    # unstructured (tag-id → bc-type mapping):
    groups: dict = field(default_factory=dict)


@dataclass
class Salinity:
    initial: Optional[float] = None
    s_left: Optional[float] = None
    s_right: Optional[float] = None


@dataclass
class Mode:
    # All optional — only fields explicitly set in the YAML propagate
    # through ``to_rakali_kwargs``. Fields left at None fall through to
    # the Rakali constructor's own defaults.
    use_multilayer: Optional[bool] = None
    use_nonhydrostatic: Optional[bool] = None
    log_level: Optional[str] = None
    defer_init: Optional[bool] = None


# ---------------------------------------------------------------------------
# top-level
# ---------------------------------------------------------------------------

@dataclass
class Config:
    domain: Domain
    bathymetry: Optional[Bathymetry] = None
    initial_condition: Optional[InitialCondition] = None
    forcing: Optional[Forcing] = None
    physics: Optional[Physics] = None
    boundary_conditions: Optional[BoundaryConditions] = None
    salinity: Optional[Salinity] = None
    mode: Optional[Mode] = None


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def _build(cls, raw: Optional[dict]):
    """Construct a dataclass from a raw YAML dict, recursing into nested
    dataclass-typed fields. Unknown keys raise ValueError."""
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError(
            f"expected dict for {cls.__name__}, got {type(raw).__name__}"
        )

    known = {f.name: f for f in fields(cls)}
    unknown = set(raw) - set(known)
    if unknown:
        raise ValueError(
            f"unknown key(s) {sorted(unknown)} in '{cls.__name__}' block"
        )

    kwargs: dict[str, Any] = {}
    for name, f in known.items():
        if name not in raw:
            continue
        value = raw[name]
        nested = _nested_dataclass(f.type)
        if nested is not None and value is not None:
            kwargs[name] = _build(nested, value)
        else:
            kwargs[name] = value
    return cls(**kwargs)


def _nested_dataclass(annotation):
    """Return the nested dataclass type if ``annotation`` is a dataclass
    or Optional[dataclass]; else None. Annotations may be strings under
    ``from __future__ import annotations`` — resolve those via ``_RESOLVED``."""
    if isinstance(annotation, str):
        # Strip the Optional[...] wrapper textually if present.
        text = annotation.strip()
        if text.startswith("Optional[") and text.endswith("]"):
            text = text[len("Optional["):-1].strip()
        return _RESOLVED.get(text)

    origin = getattr(annotation, "__origin__", None)
    if origin is Union:
        args = [a for a in annotation.__args__ if a is not type(None)]
        if len(args) == 1:
            annotation = args[0]
    if is_dataclass(annotation):
        return annotation
    return None


# Forward-ref → dataclass map (populated after all classes are defined).
_RESOLVED: dict[str, type] = {
    cls.__name__: cls
    for cls in (
        Domain, Funnel, Bowl, Island, Coastline, Bathymetry,
        DamBreak, GaussianBump, Solitary, Tilt, InitialCondition,
        Tide, Forcing, Physics, BoundaryConditions, Salinity, Mode,
        Config,
    )
}


_EDGE_FIELDS = ("west", "east", "south", "north")


def _validate_physics(cfg: Config) -> None:
    """Reject physics knobs that only apply to one domain type."""
    if cfg.physics is None:
        return
    if cfg.domain.type == "unstructured" and cfg.physics.use_hllc is not None:
        raise ValueError(
            "physics.use_hllc is a structured-grid-only knob; the "
            "unstructured solver always uses Kurganov-Noelle-Petrova "
            "central-upwind"
        )


def _validate_boundary_conditions(cfg: Config) -> None:
    """Reject BC keys that don't match the domain type, and check that
    ``groups`` keys are integers."""
    bc = cfg.boundary_conditions
    if bc is None:
        return

    has_edges = any(getattr(bc, e) is not None for e in _EDGE_FIELDS)
    has_groups = bool(bc.groups)

    if cfg.domain.type == "structured" and has_groups:
        raise ValueError(
            "boundary_conditions.groups is invalid for structured domains; "
            "use west/east/south/north instead"
        )
    if cfg.domain.type == "unstructured" and has_edges:
        raise ValueError(
            "boundary_conditions edges (west/east/south/north) are invalid "
            "for unstructured domains; use 'groups: {tag: bc_type}' instead"
        )

    for key in bc.groups:
        if not isinstance(key, int):
            raise ValueError(
                f"boundary_conditions.groups keys must be integers, "
                f"got {key!r}"
            )


def load_config(path: Union[str, Path]) -> Config:
    """Load a YAML config file and return a populated ``Config``."""
    path = Path(path)
    with path.open("r") as f:
        raw = yaml.safe_load(f) or {}

    if not isinstance(raw, dict):
        raise ValueError(f"{path}: top-level YAML must be a mapping")

    if "domain" not in raw:
        raise ValueError(f"{path}: missing required 'domain' block")

    cfg = _build(Config, raw)
    _validate_physics(cfg)
    _validate_boundary_conditions(cfg)
    return cfg


# ---------------------------------------------------------------------------
# Adapter to the structured ``Rakali(...)`` constructor
# ---------------------------------------------------------------------------

def to_rakali_kwargs(cfg: Config) -> dict:
    """Return kwargs dict for the structured ``Rakali(...)`` constructor.

    Only fields that are set in the config are included, so callers that
    omit (e.g.) ``physics.manning_n`` get the Rakali default rather than
    a forced None.
    """
    if cfg.domain.type != "structured":
        raise ValueError(
            f"to_rakali_kwargs requires domain.type='structured', "
            f"got '{cfg.domain.type}'"
        )

    # Omit unset fields: a literal None would be written into the namelist
    # (a fatal parse error) instead of letting Rakali's defaults apply.
    kw: dict = {
        k: v for k, v in (
            ("nx", cfg.domain.nx),
            ("ny", cfg.domain.ny),
            ("dx", cfg.domain.dx),
            ("dy", cfg.domain.dy),
            ("t_end", cfg.domain.t_end),
        ) if v is not None
    }
    if cfg.domain.nz is not None:
        kw["nz_layers"] = cfg.domain.nz

    bc = cfg.boundary_conditions
    if bc is not None:
        for side in ("west", "east", "south", "north"):
            v = getattr(bc, side)
            if v is not None:
                kw[f"bc_{side}"] = v

    if cfg.physics is not None:
        for src, dst in (
            ("cfl", "cfl"),
            ("manning_n", "manning_n"),
            ("smagorinsky_cs", "smagorinsky_cs"),
            ("hdiff_kappa", "hdiff_kappa"),
            ("bpg_method", "bpg_method"),
            ("use_hllc", "use_hllc"),
            ("time_scheme", "time_scheme"),
        ):
            v = getattr(cfg.physics, src)
            if v is not None:
                kw[dst] = v

    if cfg.forcing is not None:
        for src, dst in (
            ("inflow_velocity", "inflow_velocity"),
            ("inflow_salinity", "inflow_salinity"),
            ("coriolis_f", "coriolis_f"),
        ):
            v = getattr(cfg.forcing, src)
            if v is not None:
                kw[dst] = v

    if cfg.salinity is not None and cfg.salinity.initial is not None:
        kw["initial_salinity"] = cfg.salinity.initial

    if cfg.mode is not None:
        for name in ("use_multilayer", "use_nonhydrostatic",
                     "log_level", "defer_init"):
            v = getattr(cfg.mode, name)
            if v is not None:
                kw[name] = v

    return kw


__all__ = [
    "load_config",
    "to_rakali_kwargs",
    "Config",
    "Domain",
    "Bathymetry", "Funnel", "Bowl", "Island", "Coastline",
    "InitialCondition", "DamBreak", "GaussianBump", "Solitary", "Tilt",
    "Forcing", "Tide",
    "Physics",
    "BoundaryConditions",
    "Salinity",
    "Mode",
]
