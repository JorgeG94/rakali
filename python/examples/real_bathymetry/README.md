# Real-bathymetry case studies

These scripts drive the structured solver (`Rakali`) on real coastal bathymetry.
**The bathymetry files are not distributed with Rakali**, so the scripts will not run
until you supply your own data in the expected format. They are kept as worked examples
of how to set up a real-domain run from Python: loading a raster, tidal/surge boundary
forcing from a callback, wind, and snapshot output.

| Script | Expects | Notes |
|---|---|---|
| `lakes_entrance/run_le.py` | `lakes_entrance/lakes_entrance_30m.nc` (variable `b`, shape `(nx, ny)`) | 72 h surge on a clamped south boundary |
| `sydney_basin/parramatta/run_river.py` | `sydney_basin/bathymetry_25m_merged.nc` | river flood, snapshots to NetCDF |
| `sydney_basin/parramatta/run_river_minimal.py` | same | minimal version of the above |

Run from the repository root after building `librakali.so`
(`-DRAKALI_ENABLE_PYTHON=ON`), for example:

```bash
PYTHONPATH=python python python/examples/real_bathymetry/lakes_entrance/run_le.py
```

For self-contained examples that need no external data, see the parent directory.
