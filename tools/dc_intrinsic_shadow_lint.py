#!/usr/bin/env python3
"""dc_intrinsic_shadow_lint.py — catch Fortran-intrinsic-shadowing locals in
device code.

NVHPC 26.3 (and, defensively, ifx) has a code-generation bug: naming a local
variable after a Fortran intrinsic and listing it in the `local(...)` clause of
a `do concurrent` loop silently produces *approximate* math at the use site
(~1e-7 relative drift, even at `-O0 -Kieee -Mnofma`). Renaming the variable
(e.g. `scale -> loc_scale`) fixes it. The same hazard applies to a local named
after an intrinsic inside an `!$acc routine seq` / `!$omp declare target`
helper called from a `do concurrent`.

This linter flags two unambiguous shapes:

  ERROR  dc-local-intrinsic-shadow
         An intrinsic name appears in a `local(...)` / `local_init(...)` clause
         of a `do concurrent`. Names there are always variables, so this is the
         documented bug. Always exits non-zero.

  WARN   device-decl-intrinsic-shadow
         A non-parameter variable is *declared* with an intrinsic name in a file
         that contains device code (`do concurrent`, `!$acc routine`, or
         `!$omp declare target`). Catches the helper-routine variant
         (e.g. `ppm_limit_pos`). `parameter` constants are skipped — they are
         compile-time-substituted and not affected by the codegen bug.

Usage:
    python tools/dc_intrinsic_shadow_lint.py PATH [PATH ...]
    python tools/dc_intrinsic_shadow_lint.py src app --fail-on-warn
    python tools/dc_intrinsic_shadow_lint.py --diff origin/main src app
    python tools/dc_intrinsic_shadow_lint.py --list-names

Options:
    --fail-on-warn   Exit non-zero on WARN hits too (errors always exit 1).
    --diff REF       Only consider lines added vs REF (for PR gating).
    --names a,b,c    Replace the default trap-intrinsic set.
    --extra-names    Add names to the default set (comma-separated).
    --list-names     Print the active trap set and exit.

Related: claude memory `feedback_nvhpc_local_intrinsic_shadow`; reproducer
`tools/mre_nvhpc_dc_local_intrinsic_shadow.F90` (which keeps `scale` ON PURPOSE
— pass it explicitly only if you want to confirm the linter trips).
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable

SOURCE_SUFFIXES = (".F90", ".f90")

# Default trap set: intrinsic names empirically grouped with the `scale` bug
# (see claude memory). Deliberately NOT the full intrinsic list — that would be
# noisy and unverified. Extend with --extra-names if a new trap is found.
DEFAULT_TRAP_INTRINSICS = (
    "scale", "dim", "count", "sum", "product", "dot_product",
    "nint", "ceiling", "floor", "aint", "anint", "mod", "modulo", "sign",
)

DC_RE = re.compile(r"\bdo\s+concurrent\b", re.IGNORECASE)
ACC_ROUTINE_RE = re.compile(r"^\s*!\$acc\s+routine\b", re.IGNORECASE)
OMP_DECLARE_TARGET_RE = re.compile(r"^\s*!\$omp\s+declare\s+target\b", re.IGNORECASE)
# A line that opens a variable declaration: a type-spec followed (eventually) by `::`.
DECL_RE = re.compile(
    r"^\s*(real|integer|logical|character|complex|double\s+precision|"
    r"double\s+complex|type\s*\(|class\s*\()",
    re.IGNORECASE,
)
LOCAL_CLAUSE_RE = re.compile(r"\blocal(?:_init)?\s*\(", re.IGNORECASE)


# ---- Fortran line helpers (mirrors tools/dc_audit.py) ----

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
    """Executable portion of a line: trailing comment stripped (string-aware)."""
    s = strip_string_literals(line)
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def join_continuations(lines: list[str], start: int) -> tuple[str, int]:
    """Join `&`-continued code starting at `start`; return (joined, last_idx)."""
    parts: list[str] = []
    i = start
    while i < len(lines):
        c = code_part(lines[i]).rstrip()
        if c.endswith("&"):
            parts.append(c[:-1])
            i += 1
            if i < len(lines):
                nxt = code_part(lines[i]).lstrip()
                if nxt.startswith("&"):
                    lines[i] = lines[i].replace("&", " ", 1)
        else:
            parts.append(c)
            return " ".join(parts), i
    return " ".join(parts), i


def split_top_level(s: str, sep: str = ",") -> list[str]:
    """Split on `sep` ignoring nested parentheses."""
    out, buf, depth = [], [], 0
    for ch in s:
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            buf.append(ch)
        elif ch == sep and depth == 0:
            out.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf).strip())
    return out


def paren_group(s: str, open_idx: int) -> tuple[str, int]:
    """Given `s` and the index of a `(`, return (inside, index_of_matching_close)."""
    depth = 0
    for k in range(open_idx, len(s)):
        if s[k] == "(":
            depth += 1
        elif s[k] == ")":
            depth -= 1
            if depth == 0:
                return s[open_idx + 1:k], k
    return s[open_idx + 1:], len(s)


def local_clause_names(joined: str) -> list[str]:
    """Extract every name in `local(...)` / `local_init(...)` clauses."""
    names: list[str] = []
    for m in LOCAL_CLAUSE_RE.finditer(joined):
        inside, _ = paren_group(joined, joined.index("(", m.start()))
        for item in split_top_level(inside):
            name = item.strip().split("(")[0].strip().lower()
            if name:
                names.append(name)
    return names


def declared_names(joined_code: str) -> list[str]:
    """Names declared on a `... :: a, b(=init), c(dims)` line. Empty if not a
    declaration or if it carries the `parameter` attribute (constants are safe)."""
    if "::" not in joined_code or not DECL_RE.match(joined_code):
        return []
    head, _, rhs = joined_code.partition("::")
    if re.search(r"\bparameter\b", head, re.IGNORECASE):
        return []
    names: list[str] = []
    for item in split_top_level(rhs):
        # strip array spec `(...)`, char-len `*n`, init `= ...` / pointer `=> ...`
        name = re.split(r"[(*=]", item, maxsplit=1)[0].strip().lower()
        if name:
            names.append(name)
    return names


# ---- file walking + diff filter (mirrors openmp_portability_lint.py) ----

def walk_files(paths: list[Path]) -> Iterable[Path]:
    for p in paths:
        if p.is_file() and p.suffix in SOURCE_SUFFIXES:
            yield p
        elif p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.suffix in SOURCE_SUFFIXES:
                    yield f


def changed_lines(ref: str) -> dict[Path, set[int]]:
    out = subprocess.run(
        ["git", "diff", "--unified=0", "--no-color", ref, "--"],
        capture_output=True, text=True, check=True,
    ).stdout
    result: dict[Path, set[int]] = {}
    current: Path | None = None
    for line in out.splitlines():
        if line.startswith("+++ b/"):
            current = Path(line[6:])
            result.setdefault(current, set())
        elif line.startswith("@@") and current is not None:
            m = re.match(r"@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@", line)
            if m:
                start, n = int(m.group(1)), int(m.group(2) or "1")
                for i in range(start, start + n):
                    result[current].add(i)
    return result


def file_has_device_code(lines: list[str]) -> bool:
    for ln in lines:
        if DC_RE.search(code_part(ln)) or ACC_ROUTINE_RE.match(ln) \
           or OMP_DECLARE_TARGET_RE.match(ln):
            return True
    return False


def lint(paths: list[Path], trap: set[str],
         diff_lines: dict[Path, set[int]] | None
         ) -> list[tuple[str, Path, int, str, str, str]]:
    """Returns (severity, path, line, rule_id, snippet, name)."""
    hits: list[tuple[str, Path, int, str, str, str]] = []
    for f in walk_files(paths):
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        lines = text.splitlines()
        line_filter = None
        if diff_lines is not None:
            rel = f.resolve().relative_to(Path.cwd().resolve()) if f.is_absolute() else f
            line_filter = diff_lines.get(rel)
            if not line_filter:
                continue
        device = file_has_device_code(lines)
        i = 0
        while i < len(lines):
            lineno = i + 1
            code = code_part(lines[i])
            joined, last = join_continuations(lines, i)
            in_diff = line_filter is None or any(
                n in line_filter for n in range(lineno, last + 2))

            if DC_RE.search(code):
                if in_diff:
                    for name in local_clause_names(joined):
                        if name in trap:
                            hits.append(("error", f, lineno,
                                         "dc-local-intrinsic-shadow",
                                         joined.strip(), name))
                i = last + 1
                continue

            if device and in_diff:
                for name in declared_names(joined):
                    if name in trap:
                        hits.append(("warn", f, lineno,
                                     "device-decl-intrinsic-shadow",
                                     joined.strip(), name))
            i = last + 1
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", default=["src", "app"], type=Path)
    ap.add_argument("--fail-on-warn", action="store_true")
    ap.add_argument("--diff", metavar="REF", help="only scan lines added vs REF")
    ap.add_argument("--names", help="replace default trap set (comma-separated)")
    ap.add_argument("--extra-names", help="add to default trap set (comma-separated)")
    ap.add_argument("--list-names", action="store_true")
    args = ap.parse_args()

    if args.names:
        trap = {n.strip().lower() for n in args.names.split(",") if n.strip()}
    else:
        trap = set(DEFAULT_TRAP_INTRINSICS)
    if args.extra_names:
        trap |= {n.strip().lower() for n in args.extra_names.split(",") if n.strip()}

    if args.list_names:
        print("trap intrinsics:", ", ".join(sorted(trap)))
        return 0

    diff_lines = changed_lines(args.diff) if args.diff else None
    hits = lint(args.paths, trap, diff_lines)

    n_err = sum(1 for s, *_ in hits if s == "error")
    n_warn = sum(1 for s, *_ in hits if s == "warn")

    for sev, path, ln, rule, snippet, name in hits:
        prefix = "ERROR" if sev == "error" else "WARN "
        print(f"{prefix}  {path}:{ln}  [{rule}]  shadows intrinsic `{name}`")
        print(f"        {snippet}")
        if rule == "dc-local-intrinsic-shadow":
            print(f"        FIX: rename `{name}` (e.g. `loc_{name}`) in the local() clause + body.")
        else:
            print(f"        FIX: rename `{name}` (e.g. `loc_{name}`) — shadows an intrinsic in device-reachable code.")
        print()

    if not hits:
        print("dc_intrinsic_shadow_lint: clean")
        return 0
    print(f"dc_intrinsic_shadow_lint: {n_err} error(s), {n_warn} warning(s)")
    if n_err > 0:
        return 1
    if n_warn > 0 and args.fail_on_warn:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
