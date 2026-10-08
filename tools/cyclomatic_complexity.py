#!/usr/bin/env python3
"""cyclomatic_complexity.py — McCabe complexity per Fortran procedure.

Diagnostic, not a hard gate by default: it ranks subroutines / functions /
programs by cyclomatic complexity so you can spot refactoring candidates (deeply
branched dispatchers, god-routines). Cyclomatic complexity here is

    CC = 1 + (# decision points)

where decision points are counted from the structured-control keywords:

    if (...) / else if (...)         each condition           +1
    do / do while / do concurrent    each loop                +1
    case (...)                       each case label          +1   (case default: not counted)
    where (...) / elsewhere          each masked branch       +1
    forall (...)                     each                      +1
    go to / goto                     each                      +1
    .and. / .or.                     only with --logical-ops  +1

`else`, `end if`, `select case`, `contains`, and interface-block bodies are not
counted. Internal procedures (after `contains`) are scored as their own units —
a parent routine's score excludes its nested procedures' bodies.

This is a line-based approximation (string- and continuation-aware), not a real
Fortran parser; it's deliberately conservative and consistent rather than
exact. Good enough to rank routines.

Usage:
    python tools/cyclomatic_complexity.py                 # scans src/
    python tools/cyclomatic_complexity.py src app
    python tools/cyclomatic_complexity.py --top 40 --threshold 20
    python tools/cyclomatic_complexity.py --all
    python tools/cyclomatic_complexity.py --max 60        # exit 1 if any CC > 60
    python tools/cyclomatic_complexity.py --csv

Options:
    --threshold N   Mark routines with CC >= N (default 15).
    --top N         Show the N highest-CC routines (default 30).
    --all           Show every routine (overrides --top).
    --max N         Exit non-zero if any routine has CC > N (gate). Default: off.
    --logical-ops   Also count .and./.or. (extended McCabe).
    --csv           Emit machine-readable CSV (cc,file,line,name) instead.
    --exclude RE    Comma-separated name regexes to exclude. Default excludes
                    lifecycle boilerplate: *_enter_data / *_exit_data / *_destroy.
    --no-exclude    Include lifecycle/boilerplate routines in the ranking.

Lifecycle routines (`_init`/`_destroy`/`_enter_data`/`_exit_data`) are mostly
mechanical allocate/deallocate/`!$acc` chains — a destructor with 40
`if (allocated(x)) deallocate(x)` lines scores CC 41 without being a refactor
target. The `_destroy`/`_enter_data`/`_exit_data` trio is excluded by default so
the ranking surfaces genuine cognitive complexity; `_init` is kept (it can carry
real setup logic).
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

SOURCE_SUFFIXES = (".F90", ".f90")

# ---- decision-point patterns (case-insensitive, applied to label-stripped code) ----
RE_IF = re.compile(r"^\s*(?:\w+\s*:\s*)?if\s*\(", re.IGNORECASE)
RE_ELSEIF = re.compile(r"^\s*else\s*if\s*\(", re.IGNORECASE)
RE_DO = re.compile(r"^\s*(?:\w+\s*:\s*)?do\b", re.IGNORECASE)
RE_CASE = re.compile(r"^\s*case\s*\(", re.IGNORECASE)
RE_WHERE = re.compile(r"^\s*(?:\w+\s*:\s*)?where\s*\(", re.IGNORECASE)
RE_ELSEWHERE = re.compile(r"^\s*else\s*where\b", re.IGNORECASE)
RE_FORALL = re.compile(r"^\s*(?:\w+\s*:\s*)?forall\s*\(", re.IGNORECASE)
RE_GOTO = re.compile(r"\bgo\s*to\b", re.IGNORECASE)
RE_LOGICAL_OP = re.compile(r"\.(and|or)\.", re.IGNORECASE)

# ---- program-unit boundaries ----
RE_OPEN_SUB = re.compile(r"\bsubroutine\s+(\w+)", re.IGNORECASE)
RE_OPEN_FUNC = re.compile(r"\bfunction\s+(\w+)", re.IGNORECASE)
RE_OPEN_PROG = re.compile(r"^\s*program\s+(\w+)", re.IGNORECASE)
RE_END_ANY = re.compile(r"^\s*end\s*$", re.IGNORECASE)
RE_END_UNIT = re.compile(r"^\s*end\s+(subroutine|function|program)\b", re.IGNORECASE)
RE_INTERFACE = re.compile(r"^\s*(?:abstract\s+)?interface\b", re.IGNORECASE)
RE_END_INTERFACE = re.compile(r"^\s*end\s+interface\b", re.IGNORECASE)
RE_LABEL = re.compile(r"^\s*\d+\s+")

# Lifecycle / boilerplate routines: allocate / deallocate / acc-enter-exit
# chains whose complexity is mechanical (e.g. `if (allocated(x)) deallocate(x)`
# repeated per field), not cognitive. Excluded from the ranking by default —
# pass --no-exclude to include them, or --exclude to override the patterns.
# `_init` is intentionally NOT here: it can carry real setup logic.
DEFAULT_EXCLUDE_PATTERNS = (r"_enter_data$", r"_exit_data$", r"_destroy$")


@dataclass
class Proc:
    name: str
    path: Path
    line: int
    cc: int = 1


def strip_string_literals(line: str) -> str:
    out = []
    quote = None
    for ch in line:
        if quote:
            out.append(ch if ch == quote else " ")
            if ch == quote:
                quote = None
        else:
            out.append(ch)
            if ch in ("'", '"'):
                quote = ch
    return "".join(out)


def code_part(line: str) -> str:
    s = strip_string_literals(line)
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def logical_lines(lines: list[str]) -> Iterable[tuple[int, str]]:
    """Yield (start_lineno, joined_code) folding `&` continuations, comments
    stripped, leading numeric labels removed."""
    i = 0
    n = len(lines)
    while i < n:
        start = i
        parts: list[str] = []
        while i < n:
            c = code_part(lines[i]).rstrip()
            if c.endswith("&"):
                parts.append(c[:-1])
                i += 1
                if i < n:
                    nxt = code_part(lines[i]).lstrip()
                    if nxt.startswith("&"):
                        lines[i] = lines[i].replace("&", " ", 1)
            else:
                parts.append(c)
                i += 1
                break
        joined = " ".join(parts)
        yield start + 1, RE_LABEL.sub("", joined)


def decision_points(code: str, count_logical_ops: bool) -> int:
    d = 0
    if RE_ELSEIF.match(code):
        d += 1
    elif RE_IF.match(code):
        d += 1
    if RE_DO.match(code):
        d += 1
    if RE_CASE.match(code):
        d += 1
    if RE_ELSEWHERE.match(code):
        d += 1
    elif RE_WHERE.match(code):
        d += 1
    if RE_FORALL.match(code):
        d += 1
    d += len(RE_GOTO.findall(code))
    if count_logical_ops:
        d += len(RE_LOGICAL_OP.findall(code))
    return d


def scan_file(path: Path, count_logical_ops: bool) -> list[Proc]:
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return []
    lines = text.splitlines()
    procs: list[Proc] = []
    stack: list[Proc] = []
    iface_depth = 0

    for lineno, code in logical_lines(lines):
        if not code.strip():
            continue

        # Interface blocks: their subroutine/function lines are declarations.
        if RE_INTERFACE.match(code):
            iface_depth += 1
            continue
        if RE_END_INTERFACE.match(code):
            iface_depth = max(0, iface_depth - 1)
            continue
        if iface_depth > 0:
            continue

        # Unit close (bare `end`, or `end subroutine/function/program`).
        if RE_END_UNIT.match(code) or RE_END_ANY.match(code):
            if stack:
                procs.append(stack.pop())
            continue

        # Unit open. `end ...` already handled above, so these are real opens.
        m = RE_OPEN_SUB.search(code)
        if m and not code.lstrip().lower().startswith("end"):
            stack.append(Proc(m.group(1), path, lineno))
            continue
        m = RE_OPEN_PROG.match(code)
        if m:
            stack.append(Proc(m.group(1), path, lineno))
            continue
        m = RE_OPEN_FUNC.search(code)
        if m and not code.lstrip().lower().startswith("end"):
            stack.append(Proc(m.group(1), path, lineno))
            continue

        if "contains" == code.strip().lower():
            continue

        if stack:
            stack[-1].cc += decision_points(code, count_logical_ops)

    # Any unclosed units (malformed / truncated) still get reported.
    procs.extend(reversed(stack))
    return procs


def walk_files(paths: list[Path]) -> Iterable[Path]:
    for p in paths:
        if p.is_file() and p.suffix in SOURCE_SUFFIXES:
            yield p
        elif p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.suffix in SOURCE_SUFFIXES:
                    yield f


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", default=["src"], type=Path)
    ap.add_argument("--threshold", type=int, default=15)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--max", type=int, default=None,
                    help="exit non-zero if any routine CC exceeds this")
    ap.add_argument("--logical-ops", action="store_true")
    ap.add_argument("--csv", action="store_true")
    ap.add_argument("--exclude", default=None,
                    help="comma-separated name regexes to exclude (default: "
                         "lifecycle routines *_enter_data / *_exit_data / *_destroy)")
    ap.add_argument("--no-exclude", action="store_true",
                    help="include lifecycle/boilerplate routines in the ranking")
    args = ap.parse_args()

    procs: list[Proc] = []
    for f in walk_files(args.paths):
        procs.extend(scan_file(f, args.logical_ops))

    if args.no_exclude:
        patterns: list[str] = []
    elif args.exclude:
        patterns = [p.strip() for p in args.exclude.split(",") if p.strip()]
    else:
        patterns = list(DEFAULT_EXCLUDE_PATTERNS)
    excl = [re.compile(p, re.IGNORECASE) for p in patterns]
    n_excluded = sum(1 for p in procs if any(r.search(p.name) for r in excl))
    procs = [p for p in procs if not any(r.search(p.name) for r in excl)]

    if not procs:
        msg = "cyclomatic_complexity: no procedures found"
        if n_excluded:
            msg += f" ({n_excluded} excluded by name pattern)"
        print(msg)
        return 0

    procs.sort(key=lambda p: p.cc, reverse=True)

    if args.csv:
        print("cc,file,line,name")
        for p in procs:
            print(f"{p.cc},{p.path},{p.line},{p.name}")
        return 0

    shown = procs if args.all else procs[:args.top]
    print(f"  {'CC':>4}  {'routine':<44}  location")
    print(f"  {'-'*4}  {'-'*44}  {'-'*30}")
    for p in shown:
        flag = f"   *** >= threshold ({args.threshold})" if p.cc >= args.threshold else ""
        loc = f"{p.path}:{p.line}"
        print(f"  {p.cc:>4}  {p.name:<44.44}  {loc}{flag}")

    ccs = [p.cc for p in procs]
    n = len(ccs)
    mean = sum(ccs) / n
    srt = sorted(ccs)
    median = srt[n // 2] if n % 2 else (srt[n // 2 - 1] + srt[n // 2]) / 2
    over = sum(1 for c in ccs if c >= args.threshold)
    print()
    summary = (f"cyclomatic_complexity: {n} procedures, max {max(ccs)}, "
               f"mean {mean:.1f}, median {median:g}, "
               f"{over} over threshold {args.threshold}")
    if n_excluded:
        summary += f"  ({n_excluded} lifecycle routines excluded; --no-exclude to include)"
    print(summary)

    if args.max is not None:
        worst = max(ccs)
        if worst > args.max:
            print(f"[cyclomatic_complexity] --max {args.max}: worst is {worst}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
