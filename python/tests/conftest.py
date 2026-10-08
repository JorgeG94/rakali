"""Pytest configuration for the Rakali Python wrapper tests.

These tests load librakali.so via ctypes and drive the unstructured
solver from Python. They require:

* a built shared library (cmake -DRAKALI_ENABLE_PYTHON=ON, or RAKALI_LIB set)
* meshpy in the conda env (provided by environments/rakali_env_3.13.yml)

Run from the repo root::

    PYTHONPATH=python pytest python/tests/

Tests must not run in parallel (they share one process and non-thread-safe
HDF5) and each test must call ``destroy()`` (or use the
context-manager form) before the next one creates a new simulation.
"""
import sys
from pathlib import Path

# Make `import rakali` work without an editable install when running
# pytest from the repo root.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_PY_DIR = _REPO_ROOT / "python"
if str(_PY_DIR) not in sys.path:
    sys.path.insert(0, str(_PY_DIR))

# Quiet the pic_logger to WARNING for the whole test session: each test
# rebuilds a mesh + state and the INFO lines drown out pytest output.
# Set RAKALI_TEST_LOG_LEVEL=8 (info) to restore the default if needed.
import os as _os

try:
    from rakali.unstr import set_log_level, LOG_WARNING  # noqa: E402

    set_log_level(int(_os.environ.get("RAKALI_TEST_LOG_LEVEL", LOG_WARNING)))
except FileNotFoundError:
    # No shared library: pure-Python tests (mesh builder, YAML loader) still
    # run; solver tests fail with the library-not-found message.
    pass
