"""Phase 1b smoke tests for the multilayer C API + Python wrapper.

These exercise the unstructured multilayer plumbing from Python:

* :meth:`enable_multilayer` allocates the layer arrays on the Fortran
  side and ``nz_ml`` reports the expected depth.
* Right after :meth:`finalize`, ``h_layer`` / ``hu_layer`` carry
  ``dsig_k`` times the barotropic state (uniform vertical profile).
* After stepping, the load-bearing Phase 1 invariants hold on the
  Python side:

    sum_k h_layer(cell, k)  == height(cell)
    sum_k hu_layer(cell, k) == xmom(cell)
    sum_k hv_layer(cell, k) == ymom(cell)

  to machine precision because the solver's end-of-step resync
  explicitly overwrites the barotropic fields with the layer sums.

* With a barotropic-only physics (Phase 1) and a shear-free initial
  condition, every layer stays at the depth-averaged velocity, so the
  velocity spread across layers per cell is zero.

These are smoke / invariant checks, not physics validation -- the
analytical side will move in with Phase 2+.
"""
from __future__ import annotations

import numpy as np
import pytest

from rakali.unstr import UnstructuredMesh, UnstructuredSimulation


NZ = 5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_square(max_area: float = 0.25):
    """A 10x2 m channel with walls on all sides."""
    return UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=10.0, y1=2.0, max_area=max_area,
        bc_tags={"left": 1, "right": 1, "bottom": 1, "top": 1},
    )


def _flat_lake(sim, depth: float = 1.0):
    bed = np.zeros(sim.ncells, dtype=np.float64)
    stage = np.full(sim.ncells, depth, dtype=np.float64)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


def _tilted_stage(sim, slope: float = 0.05, mean: float = 1.0):
    """Perturbed stage that drives motion without dry-out."""
    cx = sim.centroid_x
    bed = np.zeros_like(cx)
    stage = mean + slope * cx
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_nz_ml_reports_expected_value():
    """enable_multilayer sets the layer count; reading back matches."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        assert sim.nz_ml == 0, "no-multilayer default should report 0"
        sim.enable_multilayer(nz=NZ)
        assert sim.nz_ml == NZ


def test_enable_multilayer_rejected_after_prepare():
    """The wrapper raises if multilayer is switched on mid-lifecycle --
    the underlying C entrypoint does not tolerate a rekeyed state after
    the device mapping has been taken."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        _flat_lake(sim)
        sim.prepare()
        with pytest.raises(RuntimeError, match="enable_multilayer"):
            sim.enable_multilayer(nz=NZ)


def test_layer_shapes_after_finalize():
    """Layer views come back as (ncells_owned, nz_ml) numpy arrays and
    carry dsig_k * barotropic after finalize populates them from the
    IC."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        sim.enable_multilayer(nz=NZ)
        _flat_lake(sim, depth=2.0)
        sim.prepare()

        hl = sim.h_layer
        assert hl.shape == (sim.ncells, NZ)
        # Uniform dsig = 1/nz, so every slot carries H/nz = 2/5
        assert np.allclose(hl, 2.0 / NZ, atol=1e-14)
        # Zero momentum IC -> zero everywhere in layers
        assert np.allclose(sim.hu_layer, 0.0, atol=1e-14)
        assert np.allclose(sim.hv_layer, 0.0, atol=1e-14)


def test_layer_sum_invariants_after_step():
    """The Phase 1 load-bearing invariants: after stepping, layer sums
    must equal the barotropic fields to machine precision (the resync
    at the end of solver_unstr_step writes these by construction)."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        sim.enable_multilayer(nz=NZ)
        _tilted_stage(sim, slope=0.05)
        sim.advance(0.1)

        hl = sim.h_layer
        hul = sim.hu_layer
        hvl = sim.hv_layer

        # The invariants the resync guarantees
        h_sum = hl.sum(axis=1)
        hu_sum = hul.sum(axis=1)
        hv_sum = hvl.sum(axis=1)

        assert np.max(np.abs(h_sum - sim.height)) < 1e-12
        assert np.max(np.abs(hu_sum - sim.xmom)) < 1e-12
        assert np.max(np.abs(hv_sum - sim.ymom)) < 1e-12


def test_uniform_ic_stays_shear_free():
    """With barotropic-only physics and an IC that has no vertical
    shear, layer velocities u_k = hu_k / h_k must stay uniform across
    layers at every cell. Disable the baroclinic path (beta_s=0,
    alpha_t=0) so sigma BPG + vertical exchange cannot inject
    spurious shear from the tilted-stage IC -- matches the fix
    applied on the Fortran side in test_multilayer_unstr.F90."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        sim.enable_multilayer(nz=NZ)
        sim.set_eos(beta_s=0.0, alpha_t=0.0)
        _tilted_stage(sim, slope=0.05)
        sim.advance(0.05)

        hl = sim.h_layer
        hul = sim.hu_layer
        hvl = sim.hv_layer

        # Skip cells that dried out in the step
        wet = hl.min(axis=1) > 1e-10
        assert wet.any(), "no wet cells -- test setup is broken"

        uk = hul[wet] / hl[wet]
        vk = hvl[wet] / hl[wet]
        # spread = max - min across layers, per cell
        u_spread = uk.max(axis=1) - uk.min(axis=1)
        v_spread = vk.max(axis=1) - vk.min(axis=1)
        assert np.max(u_spread) < 1e-10
        assert np.max(v_spread) < 1e-10


def test_nonuniform_dsig_default_is_uniform():
    """The Phase 1 default is uniform sigma layers: dsig_k = 1/nz.
    Pressed into the h_layer array via dsig*H, so checking one cell's
    layer thicknesses recovers that."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        sim.enable_multilayer(nz=NZ)
        _flat_lake(sim, depth=3.0)
        sim.prepare()

        hl = sim.h_layer
        # every cell should carry five layers of 3/5 = 0.6
        assert np.allclose(hl, 0.6, atol=1e-14)
        # layer thickness fractions sum to 1 (=> sum h_k = H)
        assert np.allclose(hl.sum(axis=1), 3.0, atol=1e-14)


# ---------------------------------------------------------------------------
# VCOORD_ZSTAR_FULL — F4 Python API exposure
# ---------------------------------------------------------------------------


def _shelf_lake(sim, h_shallow: float = 10.0, h_deep: float = 200.0):
    """Hard step at x=5: shallow on left, deep on right.  Wraps the
    Python wrapper's bathymetry / IC setters."""
    cx = sim.centroid_x
    bed = np.where(cx < 5.0, -h_shallow, -h_deep).astype(np.float64)
    stage = np.zeros_like(bed)
    sim.set_bathymetry(bed)
    sim.set_initial_stage(stage)


def test_zstar_full_set_vcoord_accepts():
    """`set_vcoord('zstar_full', ...)` is accepted and reaches the C API
    without raising.  Smoke check for the Python plumbing."""
    mesh = _build_square()
    with UnstructuredSimulation(mesh) as sim:
        sim.enable_multilayer(nz=NZ)
        sim.set_vcoord(
            coord_type="zstar_full",
            zstar_h_surf_target=0.5,
            zstar_n_surf=2,
            zstar_stretching="uniform",
            zstar_h_min=1e-4,
        )
        _flat_lake(sim, depth=3.0)
        sim.prepare()
        # Solver init must run cleanly with ZSTAR_FULL.
        hl = sim.h_layer
        # Mass per column: sum(h_layer) == H == 3.0
        assert np.allclose(hl.sum(axis=1), 3.0, atol=1e-12)


def test_zstar_full_shelf_surface_target():
    """Headline ZSTAR_FULL property: surface ROMS layer thickness ≈
    H_SURF_TARGET on both sides of a hard shelf-break, regardless of
    total column depth.  z*-lite cannot deliver this (surface would
    scale with H)."""
    mesh = UnstructuredMesh.rectangle(
        x0=0.0, y0=0.0, x1=10.0, y1=2.0, max_area=0.5,
        bc_tags={"left": 1, "right": 1, "bottom": 1, "top": 1},
    )
    NZ_FULL = 10
    H_SURF_TARGET = 2.0
    with UnstructuredSimulation(mesh) as sim:
        sim.enable_multilayer(nz=NZ_FULL)
        sim.set_vcoord(
            coord_type="zstar_full",
            zstar_h_surf_target=H_SURF_TARGET,
            zstar_n_surf=3,
            zstar_stretching="uniform",
        )
        _shelf_lake(sim, h_shallow=10.0, h_deep=200.0)
        sim.prepare()

        hl = sim.h_layer  # (ncells, nz), axis 1 in ROMS order: k=nz-1 is surface
        cx = sim.centroid_x[:hl.shape[0]]

        # ROMS surface = last index along axis 1
        dz_surf = hl[:, -1]
        # On the shallow side (cx < 5) and the deep side, surface
        # layer should be near H_SURF_TARGET (between H/2 and 2H)
        shallow_mask = cx < 5.0
        deep_mask = cx >= 5.0
        assert np.all(dz_surf[shallow_mask] > 0.5 * H_SURF_TARGET)
        assert np.all(dz_surf[shallow_mask] < 2.0 * H_SURF_TARGET)
        assert np.all(dz_surf[deep_mask] > 0.5 * H_SURF_TARGET)
        assert np.all(dz_surf[deep_mask] < 2.0 * H_SURF_TARGET)
