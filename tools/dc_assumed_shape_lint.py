#!/usr/bin/env python3
"""dc_assumed_shape_lint.py — catch assumed-shape dummy arguments referenced
inside `do concurrent` kernels.

NVHPC walks the Fortran array descriptor on every kernel launch when an
assumed-shape dummy (`arr(:,:,:)`) is referenced inside a `do concurrent`
loop — a past vertical-mixing kernel incident triggered
1.4 M descriptor-walk memcpys over 52 s of runtime until all hot kernels
were converted to explicit-shape (`arr(nx,ny,nz)`).  The fix is enforced
here for all new code; the conversion pattern is well-proven in
`kappa_shear_merge_into_kv_kt` and `epbl_merge_into_kv_kt`.

Rule enforced:
  ERROR  dc-assumed-shape-dummy
         An assumed-shape dummy argument (declared with a shape spec that
         contains ONLY colons, commas, and spaces, e.g. `(:)`, `(:,:)`,
         `(:,:,:)`) is referenced inside a `do concurrent` region of the
         same procedure.  `pointer` and `allocatable` dummies are skipped
         (they can't be converted to explicit-shape).  Assumed-SIZE `(*)`
         dummies are also skipped.

Waiver:
  Annotate the DECLARATION LINE or the line IMMEDIATELY ABOVE it with the
  comment token `assumed-shape-ok` to waive a specific dummy.  Recommended
  form:
      ! assumed-shape-ok: <reason>
      real(wp), intent(in) :: arr(:, :, :)
  or on the same line (trailing comment):
      real(wp), intent(in) :: arr(:, :, :)  ! assumed-shape-ok: diag cadence-bounded

  Cadence-bounded paths (diag fills that fire once per output frame, EPBL
  column kernel's tracer args that run at thermo cadence) qualify for the
  waiver.  Hot per-stage kernels (every RK stage, every time-step) do NOT.

Usage:
    python3 tools/dc_assumed_shape_lint.py src app
    python3 tools/dc_assumed_shape_lint.py --diff origin/main src app
    python3 tools/dc_assumed_shape_lint.py src app            # full scan

Options:
    --diff REF   Only report findings whose declaration line or in-DC
                 reference line is among lines added vs REF.  When REF is
                 not resolvable (shallow clone, fresh worktree) the script
                 exits 0 with an informational message — the dedicated CI
                 workflow is the authoritative gate.

References:
  - FORTRAN_STYLE.md §GPU / do concurrent — explicit-shape rule
  - CLAUDE.md Gotchas — "Assumed-shape dummies in do concurrent kernels"
  - Memory: feedback_nvhpc_explicit_shape_kernels
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
# Line helpers (mirror dc_intrinsic_shadow_lint.py)
# ---------------------------------------------------------------------------

def strip_string_literals(line: str) -> str:
    out: list[str] = []
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
    """Return the executable portion of a line (trailing comment stripped)."""
    s = strip_string_literals(line)
    bang = s.find("!")
    return line[:bang] if bang >= 0 else line


def join_continuations(lines: list[str], start: int) -> tuple[str, int]:
    """Join &-continued code starting at start; return (joined, last_idx)."""
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


# ---------------------------------------------------------------------------
# Fortran analysis patterns
# ---------------------------------------------------------------------------

PROC_OPEN_RE = re.compile(
    r"^\s*(?:pure\s+|elemental\s+|recursive\s+)*"
    r"(?:subroutine|function)\s+(\w+)\s*(\([^)]*\))?",
    re.IGNORECASE,
)
PROC_END_RE = re.compile(
    r"^\s*end\s+(?:subroutine|function)\b",
    re.IGNORECASE,
)

# Declaration that might carry a dummy: type-spec followed by :: name(shape)
DECL_RE = re.compile(
    r"^\s*(?:real|integer|logical|character|complex|double\s+precision|"
    r"double\s+complex)\b",
    re.IGNORECASE,
)

# Assumed-shape dimension spec: contains ONLY colons, commas, spaces, digits
# (digits for assumed-rank ".."), but NOT "*" (assumed-size) and not "=>" / "allocatable"
# We detect it on the declared name's parenthesised spec.
# Strategy: after "::", split names, check each one's paren group.

POINTER_ATTR_RE = re.compile(r"\bpointer\b", re.IGNORECASE)
ALLOCATABLE_ATTR_RE = re.compile(r"\ballocatable\b", re.IGNORECASE)

DC_RE = re.compile(r"\bdo\s+concurrent\b", re.IGNORECASE)
DO_RE = re.compile(r"^\s*do\b", re.IGNORECASE)
END_DO_RE = re.compile(r"^\s*end\s*do\b", re.IGNORECASE)

WAIVER_TOKEN = "assumed-shape-ok"


def is_assumed_shape_spec(spec: str) -> bool:
    """Return True if spec (without outer parens) is an assumed-shape spec.

    Fortran assumed-shape: every dimension has NO explicit upper bound, i.e.
    each dimension looks like  ':'  or  'lower:'  (the upper bound is absent,
    signalled by a trailing colon with nothing after it).  Examples:
      (:)        → assumed-shape
      (:, :)     → assumed-shape
      (:, :, :)  → assumed-shape
      (1:)       → assumed-shape (explicit lower, no upper)

    Explicit-shape: every dimension has explicit bounds on both sides, e.g.:
      (n)        → explicit (single bound → size n, lower=1)
      (0:nz_ref) → explicit (both bounds given)
      (nx, ny)   → explicit

    Assumed-size `(*)` is skipped separately.
    """
    if not spec.strip():
        return False
    if "*" in spec:
        return False
    if "=" in spec:
        return False
    # Must contain at least one ':'
    if ":" not in spec:
        return False
    # Each comma-separated dimension must have an absent upper bound
    # (trailing colon with nothing after it, possibly just whitespace).
    dims = [d.strip() for d in spec.split(",")]
    for d in dims:
        if ":" not in d:
            # No colon at all → explicit-size dimension (e.g. just 'n')
            return False
        parts = d.split(":")
        # parts[0] = lower bound (may be empty → default 1)
        # parts[1] = upper bound (must be empty for assumed-shape)
        if len(parts) != 2:
            return False
        upper = parts[1].strip()
        if upper:
            # Upper bound present → explicit-shape dimension
            return False
    return True


def paren_content(s: str, open_idx: int) -> tuple[str, int]:
    """Return (content, close_idx) for the ( at open_idx."""
    depth = 0
    for k in range(open_idx, len(s)):
        if s[k] == "(":
            depth += 1
        elif s[k] == ")":
            depth -= 1
            if depth == 0:
                return s[open_idx + 1:k], k
    return s[open_idx + 1:], len(s) - 1


def split_top_level_comma(s: str) -> list[str]:
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
        elif ch == "," and depth == 0:
            out.append("".join(buf).strip())
            buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf).strip())
    return out


def assumed_shape_dummies_from_decl(joined: str) -> list[str]:
    """Return list of dummy names declared with assumed-shape in this line.

    Skips pointer/allocatable dummies (different semantics; can't be
    converted to explicit-shape).
    """
    if "::" not in joined or not DECL_RE.match(joined):
        return []
    head, _, rhs = joined.partition("::")
    # Skip pointer or allocatable dummies
    if POINTER_ATTR_RE.search(head) or ALLOCATABLE_ATTR_RE.search(head):
        return []
    names: list[str] = []
    for item in split_top_level_comma(rhs.strip()):
        # item looks like:  name  OR  name(spec)  OR  name = init
        # strip init
        item = item.split("=")[0].strip()
        paren_pos = item.find("(")
        if paren_pos < 0:
            continue  # scalar dummy — not an array
        name = item[:paren_pos].strip().lower()
        spec_raw, _ = paren_content(item, paren_pos)
        if is_assumed_shape_spec(spec_raw):
            names.append(name)
    return names


def proc_dummy_args(sig: str) -> set[str]:
    """Extract dummy argument names from a procedure signature line."""
    paren_pos = sig.find("(")
    if paren_pos < 0:
        return set()
    inside, _ = paren_content(sig, paren_pos)
    # result(...) clause at the end — strip
    inside = re.sub(r"\bresult\s*\([^)]*\)", "", inside, flags=re.IGNORECASE)
    names = {n.strip().lower() for n in inside.split(",") if n.strip()}
    return names


def has_waiver(lines: list[str], decl_lineno: int) -> bool:
    """Return True if the declaration or any of the 3 preceding lines carries the
    waiver token.  This handles multi-line comment blocks placed immediately above
    the declaration."""
    idx = decl_lineno - 1  # 0-based
    for look in range(4):  # declaration line + 3 lines above
        check_idx = idx - look
        if check_idx < 0:
            break
        if WAIVER_TOKEN in lines[check_idx]:
            return True
        # Stop looking back once we hit a non-comment, non-blank, non-docstring line
        # that isn't the declaration itself.
        if look > 0:
            stripped = lines[check_idx].strip()
            if stripped and not stripped.startswith("!"):
                break
    return False


# ---------------------------------------------------------------------------
# Word-boundary reference check
# ---------------------------------------------------------------------------

def _word_re(name: str) -> re.Pattern[str]:
    return re.compile(r"\b" + re.escape(name) + r"\b", re.IGNORECASE)


# ---------------------------------------------------------------------------
# File walker + diff filter (mirror dc_intrinsic_shadow_lint.py)
# ---------------------------------------------------------------------------

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
                start = int(m.group(1))
                n = int(m.group(2) or "1")
                for ln in range(start, start + n):
                    result[current].add(ln)
    return result


# ---------------------------------------------------------------------------
# Core lint logic
# ---------------------------------------------------------------------------

def lint_file(
    path: Path,
    diff_filter: set[int] | None,
) -> list[tuple[int, int, str, str]]:
    """Scan one file; return list of (decl_lineno, ref_lineno, dummy_name, context).

    ref_lineno is the line number of the first reference found in the DC region.
    """
    try:
        raw_text = path.read_text(errors="replace")
    except OSError:
        return []

    lines = raw_text.splitlines()
    hits: list[tuple[int, int, str, str]] = []

    i = 0
    while i < len(lines):
        # Detect procedure open
        joined_sig, last_sig = join_continuations(list(lines), i)
        m = PROC_OPEN_RE.match(joined_sig)
        if not m:
            i += 1
            continue

        # Collect procedure body: scan for end subroutine/function
        proc_start = i
        proc_name = m.group(1).lower()
        dummy_args = proc_dummy_args(joined_sig)

        # Find end of procedure
        proc_end = last_sig + 1
        depth = 0
        j = last_sig + 1
        while j < len(lines):
            code = code_part(lines[j])
            if PROC_OPEN_RE.match(code):
                depth += 1
            elif PROC_END_RE.match(code):
                if depth == 0:
                    proc_end = j
                    break
                depth -= 1
            j += 1

        proc_lines = lines[proc_start:proc_end + 1]

        # --- Pass 1: collect assumed-shape dummies in this procedure ---
        # Maps: dummy_name → declaration_lineno (1-based in file)
        assumed_shape: dict[str, int] = {}

        k = 0
        while k < len(proc_lines):
            abs_lineno = proc_start + k + 1  # 1-based
            joined, last_local = join_continuations(list(proc_lines), k)
            for name in assumed_shape_dummies_from_decl(joined):
                if name in dummy_args:
                    if not has_waiver(lines, abs_lineno):
                        assumed_shape[name] = abs_lineno
            k = last_local + 1

        if not assumed_shape:
            i = proc_end + 1
            continue

        # --- Pass 2: find do concurrent regions, check references ---
        k = 0
        while k < len(proc_lines):
            code = code_part(proc_lines[k])
            joined_dc, last_dc = join_continuations(list(proc_lines), k)

            if DC_RE.search(joined_dc):
                # Found a do concurrent — find its matching end do
                dc_start_local = k
                dc_end_local = k
                nest = 0
                m2 = k
                while m2 < len(proc_lines):
                    c2 = code_part(proc_lines[m2])
                    if DO_RE.match(c2):
                        nest += 1
                    if END_DO_RE.match(c2):
                        nest -= 1
                        if nest == 0:
                            dc_end_local = m2
                            break
                    m2 += 1

                # Scan the DC body for references to assumed-shape dummies
                for body_k in range(dc_start_local, dc_end_local + 1):
                    abs_body_lineno = proc_start + body_k + 1
                    body_code = code_part(proc_lines[body_k])
                    for name, decl_ln in assumed_shape.items():
                        if _word_re(name).search(body_code):
                            # Apply diff filter: either decl or ref must be in diff
                            if diff_filter is not None:
                                if (decl_ln not in diff_filter
                                        and abs_body_lineno not in diff_filter):
                                    break
                            hits.append((decl_ln, abs_body_lineno, name,
                                         proc_lines[body_k].rstrip()))
                            # Only report first reference per dummy per DC
                            break

                k = dc_end_local + 1
                continue

            k = last_dc + 1

        i = proc_end + 1

    return hits


def lint(
    paths: list[Path],
    diff_lines: dict[Path, set[int]] | None,
) -> list[tuple[Path, int, int, str, str]]:
    """Return (path, decl_lineno, ref_lineno, dummy_name, context)."""
    results: list[tuple[Path, int, int, str, str]] = []
    for f in walk_files(paths):
        rel = f.resolve().relative_to(Path.cwd().resolve()) if f.is_absolute() else f
        filt = None
        if diff_lines is not None:
            filt = diff_lines.get(rel)
            if not filt:
                # also try absolute path key
                filt = diff_lines.get(f)
            if filt is None:
                continue
        for decl_ln, ref_ln, name, ctx in lint_file(f, filt):
            results.append((f, decl_ln, ref_ln, name, ctx))
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", default=["src", "app"], type=Path,
                    help="files or directories to scan (default: src app)")
    ap.add_argument("--diff", metavar="REF",
                    help="only report findings on lines added vs REF")
    args = ap.parse_args()

    diff_lines: dict[Path, set[int]] | None = None
    if args.diff:
        try:
            diff_lines = changed_lines(args.diff)
        except subprocess.CalledProcessError as exc:
            print(f"dc_assumed_shape_lint: cannot resolve diff base '{args.diff}' "
                  f"(exit {exc.returncode}) — skipping cleanly.")
            print("  (the dc-assumed-shape CI workflow is the authoritative gate)")
            return 0

    hits = lint([Path(p) for p in args.paths], diff_lines)

    if not hits:
        print("dc_assumed_shape_lint: clean")
        return 0

    for path, decl_ln, ref_ln, name, ctx in hits:
        print(f"ERROR  {path}:{decl_ln}  [dc-assumed-shape-dummy]  "
              f"assumed-shape dummy `{name}` referenced in `do concurrent` "
              f"(first ref at line {ref_ln})")
        print(f"       ref: {ctx.strip()}")
        print(f"  FIX: add explicit-shape integer dims to the procedure signature "
              f"and change `{name}(:,...)` → `{name}(nx,ny,...)` (or equivalent).")
        print(f"  WAIVER: add `! assumed-shape-ok: <reason>` on or above the "
              f"declaration at line {decl_ln} for cadence-bounded paths.")
        print()

    print(f"dc_assumed_shape_lint: {len(hits)} error(s)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
