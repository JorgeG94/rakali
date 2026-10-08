#!/usr/bin/env python3
"""Count lines of code in src/ — total, code, comment, blank — per subdirectory and overall.

Counts Fortran (.F90, .f90) sources. A line is a "comment" if its first non-whitespace
character is '!' (this also covers OpenACC/OMP directives, which is fine for a rough LoC
count; tweak the classifier if you want them as code).

Pass --tests to additionally tally tests/ — same LoC breakdown plus a count of
`new_unittest(` invocations per file as a per-suite test count.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

FORTRAN_SUFFIXES = {".F90", ".f90"}
UNITTEST_RE = re.compile(r"\bnew_unittest\(")


def classify(line: str) -> str:
    stripped = line.strip()
    if not stripped:
        return "blank"
    if stripped.startswith("!"):
        return "comment"
    return "code"


def count_file(path: Path) -> dict[str, int]:
    counts = {"code": 0, "comment": 0, "blank": 0, "total": 0, "tests": 0}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            counts[classify(line)] += 1
            counts["total"] += 1
            counts["tests"] += len(UNITTEST_RE.findall(line))
    return counts


def render_table(root: Path, files: list[Path], *, include_tests: bool) -> None:
    per_dir: dict[Path, dict[str, int]] = defaultdict(
        lambda: {"code": 0, "comment": 0, "blank": 0, "total": 0, "tests": 0, "files": 0}
    )
    grand = {"code": 0, "comment": 0, "blank": 0, "total": 0, "tests": 0, "files": 0}

    for path in files:
        c = count_file(path)
        rel_dir = path.parent.relative_to(root)
        bucket = per_dir[rel_dir]
        for k in ("code", "comment", "blank", "total", "tests"):
            bucket[k] += c[k]
            grand[k] += c[k]
        bucket["files"] += 1
        grand["files"] += 1

    cols = ["files", "code", "comment", "blank", "total"]
    if include_tests:
        cols.append("tests")
    header = f"{'directory':<40} " + " ".join(f"{c:>8}" for c in cols)
    print(header)
    print("-" * len(header))
    for d in sorted(per_dir):
        b = per_dir[d]
        label = str(d) if str(d) != "." else "(root)"
        row = " ".join(f"{b[c]:>8}" for c in cols)
        print(f"{label:<40} {row}")
    print("-" * len(header))
    total_row = " ".join(f"{grand[c]:>8}" for c in cols)
    print(f"{'TOTAL':<40} {total_row}")


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("src", nargs="?", default=str(repo_root / "src"),
                        help="source directory (default: <repo>/src)")
    parser.add_argument("--tests", nargs="?", const=str(repo_root / "tests"), default=None,
                        metavar="TESTS_DIR",
                        help="also tally tests/ with a `tests` column (default dir: <repo>/tests)")
    args = parser.parse_args()

    src = Path(args.src)
    if not src.is_dir():
        print(f"error: {src} is not a directory", file=sys.stderr)
        return 1

    src_files = sorted(p for p in src.rglob("*") if p.is_file() and p.suffix in FORTRAN_SUFFIXES)
    if not src_files:
        print(f"no Fortran sources under {src}")
    else:
        render_table(src, src_files, include_tests=False)

    if args.tests:
        tests = Path(args.tests)
        if not tests.is_dir():
            print(f"\nerror: {tests} is not a directory", file=sys.stderr)
            return 1
        test_files = sorted(p for p in tests.rglob("*") if p.is_file() and p.suffix in FORTRAN_SUFFIXES)
        if not test_files:
            print(f"\nno Fortran sources under {tests}")
        else:
            print()
            print(f"=== tests ({tests.relative_to(repo_root) if tests.is_relative_to(repo_root) else tests}) ===")
            render_table(tests, test_files, include_tests=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
