#!/usr/bin/env python3
"""Drift gate for docs/CLOSURE_MATRIX.md.

Fails (exit 1) when the closure documentation has drifted from the code:

  1. Every closure-selector knob in src/core/rki_config.F90 must appear in
     docs/CLOSURE_MATRIX.md.  The `vmix_use_*` family is *discovered* from the
     config, so a new vertical-mixing scheme cannot be added without a matrix
     row.  The other closure/forcing toggles are checked by name.
  2. Every `test_*` named in the matrix must be registered in
     tests/CMakeLists.txt.
  3. The two knob-driven tracer-advection / remap defaults in the config
     (`use_ppm_tracer = .true.`, `remap_method = "ppm"`) must match the default
     the matrix states.  These guard the Tracer Advection & Reconstruction
     section the same way the vmix checks guard the vertical-mixing section.

Wire into pre-commit + CI.  The matrix exists because prose drifts silently:
CAPABILITIES_AND_LIMITATIONS.md still claimed PP81 was the only coastal mixing
scheme three closures after KPP landed.  This makes that fail the build.

Run from anywhere; paths are resolved relative to this file.
"""
import re
import sys
from pathlib import Path
from typing import List

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "src" / "core" / "rki_config.F90"
MATRIX = ROOT / "docs" / "CLOSURE_MATRIX.md"
TESTS_CMAKE = ROOT / "tests" / "CMakeLists.txt"

# Named coastal closure / forcing selector knobs (besides the auto-discovered
# vmix_use_* family).  Each must exist in the config AND have a matrix cell.
CURATED_KNOBS = [
    "smagorinsky_cs",
    "biharmonic_nu4",
    "use_linear_bottom_drag",
    "wind_surface_layer",
    "bpg_method",
]


def check_tracer_defaults(config: str, matrix: str) -> List[str]:
    """Guard the two knob-driven tracer-advection / ALE-remap defaults.

    Same mechanism as the vmix default checks: grep the default straight out
    of the config_t declaration, then assert the doc states the same default.
    If the code and the doc disagree, the doc is the bug (per CLAUDE.md: code
    is authority) — fail loudly so it's fixed in the same PR.
    """
    errors: List[str] = []

    # use_ppm_tracer — logical default in config_t (the `= .value.` form).
    m = re.search(r"logical\s*::\s*use_ppm_tracer\s*=\s*(\.\w+\.)", config)
    if not m:
        errors.append(
            "use_ppm_tracer declaration not found in rki_config.F90 — "
            "renamed/removed in code, update this script."
        )
    else:
        default = m.group(1)
        if default != ".true.":
            errors.append(
                f"use_ppm_tracer default in rki_config.F90 is {default}, "
                "expected .true. — update the matrix and this script together."
            )
        # The doc must state the same default for this knob.
        if not re.search(r"use_ppm_tracer.*default\s*`?\.true\.`?", matrix):
            errors.append(
                "docs/CLOSURE_MATRIX.md does not state use_ppm_tracer's "
                f"default ({default}) in its tracer-advection section."
            )

    # remap_method — character default in config_t (the `= "value"` form).
    m = re.search(r'character\(len=\d+\)\s*::\s*remap_method\s*=\s*"(\w+)"', config)
    if not m:
        errors.append(
            "remap_method declaration not found in rki_config.F90 — "
            "renamed/removed in code, update this script."
        )
    else:
        default = m.group(1)
        if default != "ppm":
            errors.append(
                f'remap_method default in rki_config.F90 is "{default}", '
                'expected "ppm" — update the matrix and this script together.'
            )
        # The doc must state the same default for this knob.
        if f'`"{default}"`' not in matrix:
            errors.append(
                'docs/CLOSURE_MATRIX.md does not state remap_method\'s default '
                f'("{default}") as `"{default}"`.'
            )

    return errors


def main() -> int:
    for path in (CONFIG, MATRIX, TESTS_CMAKE):
        if not path.exists():
            print(f"check_closure_matrix: FAIL — missing {path}", file=sys.stderr)
            return 1

    config = CONFIG.read_text()
    matrix = MATRIX.read_text()
    cmake = TESTS_CMAKE.read_text()

    # 1a. discover vmix_use_* toggles *declared* in config_t (the `= default`
    #     form, so comments and local `logical :: a, b` decls don't match).
    vmix_knobs = sorted(set(re.findall(r"logical\s*::\s*(vmix_use_\w+)\s*=", config)))
    knobs = vmix_knobs + CURATED_KNOBS

    errors: List[str] = []

    # the curated names must actually exist in the config (catch a rename here)
    phantom = [k for k in CURATED_KNOBS if k not in config]
    if phantom:
        errors.append(
            "curated knob(s) not found in rki_config.F90 — renamed/removed in "
            f"code, update this script: {phantom}"
        )

    # 1b. every selector knob must have a row in the matrix
    missing_in_matrix = [k for k in knobs if k not in matrix]
    if missing_in_matrix:
        errors.append(
            "closure knob(s) in the config with no row in docs/CLOSURE_MATRIX.md: "
            f"{missing_in_matrix}"
        )

    # 2. every test named in the matrix must be registered
    matrix_tests = sorted(set(re.findall(r"\btest_[a-z0-9_]+", matrix)))
    missing_tests = [t for t in matrix_tests if t not in cmake]
    if missing_tests:
        errors.append(
            "test(s) named in docs/CLOSURE_MATRIX.md but not in "
            f"tests/CMakeLists.txt: {missing_tests}"
        )

    # 3. the two knob-driven tracer-advection / ALE-remap defaults must match
    errors.extend(check_tracer_defaults(config, matrix))

    if errors:
        for e in errors:
            print(f"check_closure_matrix: FAIL — {e}", file=sys.stderr)
        return 1

    print(
        f"check_closure_matrix: OK — {len(knobs)} closure knobs "
        f"({len(vmix_knobs)} vmix discovered) + {len(matrix_tests)} tests "
        "reconciled + tracer-advection defaults (use_ppm_tracer, remap_method) verified"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
