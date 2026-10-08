#!/usr/bin/env python3
"""Pre-commit guard: every rakali ctest entry must carry a regime LABEL.

The coastal path is being extracted into a separate open-source repo (ocean +
sea-ice stay closed), so every test registered in tests/CMakeLists.txt must be
attributable to a regime suite (`ctest -L coastal` / `-L ocean` / `-L core`).
Two things are enforced, all in tests/CMakeLists.txt:

  1. Every RAKALI_TESTS row has exactly 5 `|`-separated fields and the 5th
     (regime) is one of: coastal / ocean / core / ocean+core.
  2. Every literally-named `add_test(NAME rakali_...)` appears in some
     `set_tests_properties(... LABELS ...)`.

The pic / test-drive dependency self-tests are registered by their own
(vendored) build systems and are deliberately NOT covered here — they are not
rakali tests.  CMake enforces 1-2 fail-loud at configure time as well; this
hook catches it before a configure is ever run.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS_CMAKE = ROOT / "tests" / "CMakeLists.txt"

VALID_REGIMES = {"coastal", "ocean", "core", "ocean+core"}


def main() -> int:
    text = TESTS_CMAKE.read_text(encoding="utf-8", errors="replace")
    errors: list[str] = []

    # ------------------------------------------------------------------ 1.
    # RAKALI_TESTS rows: every quoted |-row must have 5 fields, valid regime.
    rows = [m for m in re.findall(r'"(test_[^"]*)"', text) if "|" in m]
    if not rows:
        errors.append("no RAKALI_TESTS rows found — parser broken or list renamed")
    for row in rows:
        fields = row.split("|")
        if len(fields) != 5:
            errors.append(
                f"RAKALI_TESTS row '{row}' has {len(fields)} fields; expected 5 "
                "(<test-name>|<module-name>|<collect-fn>|<suite-display>|<regime>)"
            )
        elif fields[4] not in VALID_REGIMES:
            errors.append(
                f"RAKALI_TESTS row '{fields[0]}': invalid regime '{fields[4]}' "
                f"(must be one of: {', '.join(sorted(VALID_REGIMES))})"
            )

    # ------------------------------------------------------------------ 2.
    # Literally-named add_test entries must be labelled via an explicit set_tests_properties(... LABELS ...).
    labelled: set[str] = set()
    for props in re.findall(r"set_tests_properties\(([^)]*)\)", text):
        if "LABELS" not in props:
            continue
        labelled.update(re.findall(r"\brakali_\w+", props))
    for name in re.findall(r"add_test\(\s*NAME\s+(rakali_\w+)", text):
        # names containing ${...} are handled by the loops above; the regex
        # only matches fully literal names.
        if name not in labelled:
            errors.append(
                f"ctest entry '{name}' is registered with add_test but never "
                "given a regime LABEL via set_tests_properties(... LABELS ...)"
            )

    if errors:
        print("check_test_regime_labels: FAIL")
        for e in errors:
            print(f"  - {e}")
        print(
            "\nEvery rakali ctest entry must carry a regime label so the "
            "coastal / ocean / core suites stay separable (ctest -L <regime>)."
        )
        return 1
    print(
        f"check_test_regime_labels: OK "
        f"({len(rows)} RAKALI_TESTS rows)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
