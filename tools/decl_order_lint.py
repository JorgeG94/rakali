#!/usr/bin/env python3
"""decl_order_lint.py — flag forward references to not-yet-declared identifiers
in array dimension specifications.

ifx warns #8586 ("Implicit type is given to allow out-of-order declaration.
Non-standard extension") when an explicit-shape array declaration's bounds
reference a dummy or local integer that is declared on a *later* line of the
same specification part, e.g.::

    real(wp), intent(in) :: h_layer(nx, ny, nz)   ! uses nx, ny, nz
    integer,  intent(in) :: nx, ny, nz             ! declared after — ifx #8586

This is non-standard Fortran but accepted silently by gfortran and nvfortran.
ifx enforces the standard.

Rule ID: ``decl-order-forward-ref``

Fix: move integer dimension declarations above the first array that uses them.

Usage::

    python tools/decl_order_lint.py PATH [PATH ...]
    python tools/decl_order_lint.py src app tests benchmarks
    python tools/decl_order_lint.py --diff origin/main src app

Options::

    --diff REF   Only consider files/lines changed vs REF (for PR gating).

Exit codes: 0 = clean, 1 = one or more ERROR hits.

Known limitations
-----------------
* Spec-part detection is line-based: the spec part ends at the first executable
  statement (assignment, call, do, if, …) or ``contains``.  It does *not* parse
  full Fortran grammar, so extremely unusual layouts (``entry``, ``block``
  constructs, statement functions) may mis-detect boundaries — but these are
  vanishingly rare in this codebase.
* Identifier extraction from dimension specs uses a simple regex over the
  parenthesised shape expression.  Complex parameter expressions (``max(a,b)``)
  are handled correctly for the names inside.
* Module-level ``use``-associated names and module ``parameter`` constants that
  are visible from outside the subroutine/function are intentionally *not*
  flagged — the lint only flags identifiers that appear as local/dummy
  declarations in the same spec part.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable

SOURCE_SUFFIXES = (".F90", ".f90")

# ---------------------------------------------------------------------------
# Fortran line utilities (shared style with dc_intrinsic_shadow_lint.py)
# ---------------------------------------------------------------------------

def strip_string_literals(line: str) -> str:
    out: list[str] = []
    quote: str | None = None
    for ch in line:
        if quote:
            out.append(" ")
            if ch == quote:
                quote = None
        else:
            if ch in ("'", '"'):
                quote = ch
                out.append(" ")
            else:
                out.append(ch)
    return "".join(out)


def code_part(line: str) -> str:
    """Return code text with trailing Fortran comment stripped."""
    s = strip_string_literals(line)
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def join_continuations(lines: list[str], start: int) -> tuple[str, int]:
    """Join &-continued lines; return (joined_code, last_physical_line_index)."""
    parts: list[str] = []
    i = start
    while i < len(lines):
        c = code_part(lines[i]).rstrip()
        if c.endswith("&"):
            parts.append(c[:-1].rstrip())
            i += 1
            if i < len(lines):
                nxt = code_part(lines[i]).lstrip()
                if nxt.startswith("&"):
                    # Replace the leading & so the merge doesn't leave a stray &
                    lines[i] = lines[i].replace("&", " ", 1)
        else:
            parts.append(c)
            return " ".join(parts), i
    return " ".join(parts), i


# ---------------------------------------------------------------------------
# Spec-part / procedure parsing
# ---------------------------------------------------------------------------

# Matches the opening of a subroutine or function (not interface bodies,
# module procedure declarations, or end statements).
PROC_START_RE = re.compile(
    r"^\s*(?:(?:pure|elemental|recursive|impure)\s+)*"
    r"(?:subroutine|(?:(?:real|integer|logical|character|complex|double\s+precision"
    r"|type\s*\([^)]+\))\s+)?function)\s+\w",
    re.IGNORECASE,
)
END_PROC_RE = re.compile(
    r"^\s*end\s*(?:subroutine|function)\b", re.IGNORECASE,
)
# End of specification part: first executable or CONTAINS.
# We conservatively detect common executable openers.
EXEC_START_RE = re.compile(
    r"^\s*(?:"
    r"contains\b"
    r"|call\s+\w"
    r"|(?:do|if|else|end\s*if|end\s*do|select\s+case|end\s+select|where|"
    r"forall|end\s+forall|associate|end\s+associate)\b"
    r"|return\b"
    r"|\w+\s*(?:\([^)]*\)\s*)?="  # assignment (incl. indexed)
    r")",
    re.IGNORECASE,
)

# Declarations: lines that start with a type-spec and contain ::
DECL_RE = re.compile(
    r"^\s*(?:real|integer|logical|character|complex|double\s+precision"
    r"|double\s+complex|type\s*\(|class\s*\()",
    re.IGNORECASE,
)

# Dimension attribute: dimension(...)
DIM_ATTR_RE = re.compile(r"\bdimension\s*\(", re.IGNORECASE)

# Identifier used in expressions (inside dimension specs)
IDENT_RE = re.compile(r"\b([a-zA-Z]\w*)\b")

# Fortran built-in functions / keywords that appear in dimension specs but are
# not user-declared identifiers.  We skip these so we don't falsely require
# them to be declared.
BUILTIN_SKIP = frozenset([
    "max", "min", "int", "real", "kind", "size", "lbound", "ubound",
    "merge", "mod", "abs", "sign", "trim", "len", "len_trim",
    "selected_real_kind", "selected_int_kind",
    "wp", "dp", "sp", "i4", "i8",  # common kind aliases treated as constants
    "true", "false",
])


def paren_groups(text: str) -> list[str]:
    """Return list of top-level (…) contents in text."""
    results: list[str] = []
    depth = 0
    buf: list[str] = []
    for ch in text:
        if ch == "(":
            depth += 1
            if depth == 1:
                buf = []
            else:
                buf.append(ch)
        elif ch == ")":
            depth -= 1
            if depth == 0:
                results.append("".join(buf))
                buf = []
            elif depth > 0:
                buf.append(ch)
        elif depth > 0:
            buf.append(ch)
    return results


def dim_spec_names(joined: str) -> list[str]:
    """Extract identifiers used in array dimension specs of a declaration line.

    Handles both:
      - ``real(wp) :: arr(nx, ny, nz)``       (shape after entity name)
      - ``real(wp), dimension(nx,ny) :: arr``  (dimension attribute)

    Returns the raw identifier strings (lowercase) that appear inside shape
    specifications, excluding obvious builtins and numeric literals.
    """
    if "::" not in joined or not DECL_RE.match(joined):
        return []

    head, _, rhs = joined.partition("::")
    names: list[str] = []

    # 1. dimension(...) attribute in the head
    for m in DIM_ATTR_RE.finditer(head):
        # find the matching paren group
        start = joined.index("(", m.start())
        groups = paren_groups(joined[start:])
        if groups:
            for ident in IDENT_RE.findall(groups[0]):
                il = ident.lower()
                if il not in BUILTIN_SKIP:
                    names.append(il)

    # 2. Per-entity shape specs in the rhs: entity(dims)
    # Each entity in rhs is separated by top-level commas.
    for entity in _split_top_level(rhs):
        entity = entity.strip()
        # entity name is up to first ( or = or *
        m = re.match(r"([a-zA-Z]\w*)\s*\(", entity)
        if m:
            # find the paren group for this entity
            paren_start = entity.index("(")
            groups = paren_groups(entity[paren_start:])
            if groups:
                for ident in IDENT_RE.findall(groups[0]):
                    il = ident.lower()
                    if il not in BUILTIN_SKIP:
                        names.append(il)
    return names


def _split_top_level(s: str, sep: str = ",") -> list[str]:
    """Split string on sep ignoring nested parens."""
    out: list[str] = []
    buf: list[str] = []
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            buf.append(ch)
        elif ch == sep and depth == 0:
            out.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf))
    return out


def declared_names_in_line(joined: str) -> list[str]:
    """Return all variable names declared on a ``:: name1, name2(dims), ...`` line.

    Skips ``parameter`` constants (those are visible globally and any order is fine).
    """
    if "::" not in joined or not DECL_RE.match(joined):
        return []
    head, _, rhs = joined.partition("::")
    if re.search(r"\bparameter\b", head, re.IGNORECASE):
        return []
    names: list[str] = []
    for entity in _split_top_level(rhs):
        entity = entity.strip()
        # name is up to first ( = * whitespace
        m = re.match(r"([a-zA-Z]\w*)", entity)
        if m:
            names.append(m.group(1).lower())
    return names


# ---------------------------------------------------------------------------
# Core analysis: scan a single procedure's spec part
# ---------------------------------------------------------------------------

def check_spec_part(
    spec_lines: list[tuple[int, str]],  # (1-based lineno, joined_code)
) -> list[tuple[int, int, str, str]]:
    """Analyse one procedure's spec part.

    Returns list of (array_decl_lineno, integer_decl_lineno, array_name_used, dim_ident).
    """
    hits: list[tuple[int, int, str, str]] = []

    # Build: for each line, what names does it DECLARE?
    # Map: name -> lineno of its declaration
    declared_at: dict[str, int] = {}
    for lineno, joined in spec_lines:
        for name in declared_names_in_line(joined):
            if name not in declared_at:
                declared_at[name] = lineno

    # Now scan array declarations and check if their dim specs reference
    # names declared LATER.
    for lineno, joined in spec_lines:
        dim_names = dim_spec_names(joined)
        if not dim_names:
            continue
        # Which array entity names are on this line?
        arr_names = declared_names_in_line(joined)
        for dim_name in dim_names:
            if dim_name not in declared_at:
                continue  # not a local/dummy — skip (module param, use-assoc, etc.)
            if declared_at[dim_name] > lineno:
                # Forward reference
                for arr in arr_names:
                    hits.append((lineno, declared_at[dim_name], arr, dim_name))
                break  # one hit per declaration line is enough
    return hits


# ---------------------------------------------------------------------------
# File scanner
# ---------------------------------------------------------------------------

def scan_file(path: Path) -> list[tuple[int, int, str, str, str]]:
    """Scan one file; return (array_lineno, int_lineno, array_name, dim_id, proc_name)."""
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return []

    raw_lines = text.splitlines()
    results: list[tuple[int, int, str, str, str]] = []

    # Build a list of (lineno, joined_code, physical_last_lineno) for the file.
    logical_lines: list[tuple[int, str, int]] = []
    i = 0
    while i < len(raw_lines):
        lineno = i + 1
        joined, last = join_continuations(raw_lines, i)
        logical_lines.append((lineno, joined, last + 1))
        i = last + 1

    # Walk procedures: subroutine / function ... end subroutine/function
    in_proc = False
    proc_name = ""
    spec_lines: list[tuple[int, str]] = []
    spec_done = False
    depth = 0  # nesting of blocks inside spec part (interface etc.)
    contains_depth = 0  # for nested contains

    # Track interface blocks so we don't parse their internals as real procs
    in_interface = False
    interface_depth = 0

    for lineno, joined, _ in logical_lines:
        code = joined.strip()
        code_lower = code.lower()

        # Interface block tracking (skip interior procedures)
        if re.match(r"^\s*interface\b", code, re.IGNORECASE):
            in_interface = True
            interface_depth += 1
            continue
        if in_interface:
            if re.match(r"^\s*end\s*interface\b", code, re.IGNORECASE):
                interface_depth -= 1
                if interface_depth == 0:
                    in_interface = False
            continue

        if not in_proc:
            m = PROC_START_RE.match(code)
            if m:
                in_proc = True
                spec_done = False
                spec_lines = []
                depth = 0
                contains_depth = 0
                # Extract procedure name
                nm = re.search(r"(?:subroutine|function)\s+(\w+)", code, re.IGNORECASE)
                proc_name = nm.group(1) if nm else "?"
        else:
            if END_PROC_RE.match(code):
                if not spec_done and spec_lines:
                    for ah in check_spec_part(spec_lines):
                        results.append((*ah, proc_name))
                in_proc = False
                spec_lines = []
                spec_done = False
                continue

            if spec_done:
                continue

            # Check for end of spec part
            if re.match(r"^\s*contains\b", code, re.IGNORECASE):
                if not spec_done and spec_lines:
                    for ah in check_spec_part(spec_lines):
                        results.append((*ah, proc_name))
                spec_done = True
                continue

            # Detect executable statements that end the spec part.
            # But be careful: some lines look like assignments but are
            # :: declarations.  Only trigger if no '::' in the line.
            if "::" not in code and EXEC_START_RE.match(code):
                if spec_lines:
                    for ah in check_spec_part(spec_lines):
                        results.append((*ah, proc_name))
                spec_done = True
                continue

            # Accumulate spec-part lines (only those with ::)
            if "::" in code and DECL_RE.match(code):
                spec_lines.append((lineno, joined))

    return results


# ---------------------------------------------------------------------------
# File walking + diff filter
# ---------------------------------------------------------------------------

def walk_files(paths: list[Path]) -> Iterable[Path]:
    for p in paths:
        if p.is_file() and p.suffix in SOURCE_SUFFIXES:
            yield p
        elif p.is_dir():
            for f in sorted(p.rglob("*")):
                if f.suffix in SOURCE_SUFFIXES:
                    yield f


def changed_files(ref: str) -> set[Path]:
    out = subprocess.run(
        ["git", "diff", "--name-only", ref, "--"],
        capture_output=True, text=True, check=True,
    ).stdout
    return {Path(p.strip()) for p in out.splitlines() if p.strip()}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "paths", nargs="*", default=["src", "app"], type=Path,
        help="Files or directories to scan (default: src app)",
    )
    ap.add_argument(
        "--diff", metavar="REF",
        help="Only scan files changed vs REF (for PR gating).",
    )
    args = ap.parse_args()

    diff_files: set[Path] | None = None
    if args.diff:
        diff_files = changed_files(args.diff)

    all_hits: list[tuple[Path, int, int, str, str, str]] = []
    file_counts: dict[Path, int] = {}

    for f in walk_files(args.paths):
        if diff_files is not None:
            # compare relative path
            try:
                rel = f.resolve().relative_to(Path.cwd().resolve())
            except ValueError:
                rel = f
            if rel not in diff_files and f not in diff_files:
                continue

        hits = scan_file(f)
        if hits:
            file_counts[f] = len(hits)
            for arr_ln, int_ln, arr_name, dim_id, proc in hits:
                all_hits.append((f, arr_ln, int_ln, arr_name, dim_id, proc))

    if not all_hits:
        print("decl_order_lint: clean")
        return 0

    for path, arr_ln, int_ln, arr_name, dim_id, proc in all_hits:
        print(
            f"ERROR  {path}:{arr_ln}  [decl-order-forward-ref]  "
            f"in `{proc}`: array `{arr_name}` uses dim `{dim_id}` "
            f"which is declared at line {int_ln} (after this line)"
        )
        print(
            f"       FIX: declare the integer dims before the arrays that use them "
            f"(move line {int_ln} above line {arr_ln})."
        )
        print()

    total = len(all_hits)
    print(f"decl_order_lint: {total} error(s) in {len(file_counts)} file(s)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
