#!/usr/bin/env python3
"""Diagnostic scan: classify every `do concurrent` loop in the source tree.

Reports per-loop:
  - file:line
  - number of induction indices
  - clauses present (local / local_init / shared / default / reduce / mask)
  - what sits directly above: !$acc, !$omp, plain comment, or nothing

Aggregates stats so we can decide the transformation strategy before writing it.

Usage:
    python tools/dc_audit.py            # scans src/
    python tools/dc_audit.py src tests  # scan multiple roots
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from _parallel import pmap, resolve_jobs

DC_RE = re.compile(r"\bdo\s+concurrent\b", re.IGNORECASE)
ACC_RE = re.compile(r"^\s*!\$acc\b", re.IGNORECASE)
OMP_RE = re.compile(r"^\s*!\$omp\b", re.IGNORECASE)
COMMENT_RE = re.compile(r"^\s*!")
KNOWN_CLAUSES = ("local_init", "local", "shared", "default", "reduce")


@dataclass
class Loop:
    path: Path
    line: int
    header: str
    n_indices: int
    clauses: list[str] = field(default_factory=list)
    has_mask: bool = False
    above: str = "nothing"  # nothing | acc | omp | comment | code


def strip_string_literals(line: str) -> str:
    """Replace string literal contents with spaces so `!` inside strings doesn't
    look like a comment marker. Keeps column positions stable."""
    out = []
    i = 0
    quote = None
    while i < len(line):
        ch = line[i]
        if quote:
            if ch == quote:
                quote = None
                out.append(ch)
            else:
                out.append(" ")
        else:
            if ch in ("'", '"'):
                quote = ch
                out.append(ch)
            else:
                out.append(ch)
        i += 1
    return "".join(out)


def code_part(line: str) -> str:
    """Return the executable portion (strip trailing comment, ignoring strings)."""
    s = strip_string_literals(line)
    bang = s.find("!")
    if bang >= 0:
        return line[:bang]
    return line


def join_continuations(lines: list[str], start: int) -> tuple[str, int]:
    """Starting at `start`, follow `&` continuations and return joined code +
    the index of the last consumed line."""
    parts = []
    i = start
    while i < len(lines):
        c = code_part(lines[i]).rstrip()
        if c.endswith("&"):
            parts.append(c[:-1])
            i += 1
            # skip leading `&` on next line if present
            if i < len(lines):
                nxt = code_part(lines[i]).lstrip()
                if nxt.startswith("&"):
                    lines[i] = lines[i].replace("&", " ", 1)
        else:
            parts.append(c)
            return (" ".join(parts), i)
    return (" ".join(parts), i)


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


def parse_dc_header(joined: str) -> tuple[int, list[str], bool] | None:
    """Parse a joined `do concurrent (...)` line.
    Returns (n_indices, clauses, has_mask) or None if it can't be parsed."""
    m = DC_RE.search(joined)
    if not m:
        return None
    rest = joined[m.end():].lstrip()
    if not rest.startswith("("):
        return None
    # find matching close paren
    depth = 0
    end = -1
    for idx, ch in enumerate(rest):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                end = idx
                break
    if end < 0:
        return None
    inside = rest[1:end]
    trailing = rest[end + 1:].strip()
    # Heuristic split: header items are comma-separated at depth 0.
    # Items look like:
    #   <var> = <lo>:<hi>[:<step>]   (induction)
    #   <logical-expr>               (mask, must contain no '=' at top level except .eq./.ne./etc.)
    #   local(...) / shared(...) / default(none) / reduce(op:var) / local_init(...)
    items = split_top_level(inside, ",")
    n_indices = 0
    clauses: list[str] = []
    has_mask = False
    for it in items:
        low = it.lower().lstrip()
        matched_clause = None
        for c in KNOWN_CLAUSES:
            if low.startswith(c) and (len(low) == len(c) or low[len(c)] in "( "):
                matched_clause = c
                break
        if matched_clause:
            clauses.append(matched_clause)
            continue
        # induction triplet?  pattern: name = expr : expr [ : expr ]
        # Use ':' presence at depth 0 outside the lhs '=' as the signal.
        if "=" in it and ":" in it:
            # crude but adequate
            n_indices += 1
            continue
        # otherwise treat as mask
        if it.strip():
            has_mask = True
    # Locality specs sit after the triplet parens: e.g. ` local(a,b) shared(c)`
    if trailing:
        t = trailing
        while t:
            tlow = t.lower()
            matched = None
            for c in KNOWN_CLAUSES:
                if tlow.startswith(c) and (len(t) == len(c) or t[len(c)] in "( "):
                    matched = c
                    break
            if not matched:
                break
            clauses.append(matched)
            rest_after = t[len(matched):].lstrip()
            if rest_after.startswith("("):
                depth = 0
                j = -1
                for k, ch in enumerate(rest_after):
                    if ch == "(":
                        depth += 1
                    elif ch == ")":
                        depth -= 1
                        if depth == 0:
                            j = k
                            break
                if j < 0:
                    break
                t = rest_after[j + 1:].lstrip()
            else:
                t = rest_after
    return (n_indices, clauses, has_mask)


def classify_above(lines: list[str], dc_line_idx: int) -> str:
    """Look at the nearest non-blank line above; classify."""
    j = dc_line_idx - 1
    while j >= 0 and lines[j].strip() == "":
        j -= 1
    if j < 0:
        return "nothing"
    prev = lines[j]
    if ACC_RE.match(prev):
        return "acc"
    if OMP_RE.match(prev):
        return "omp"
    if COMMENT_RE.match(prev):
        return "comment"
    return "code"


def scan_file(path: Path) -> list[Loop]:
    try:
        text = path.read_text()
    except UnicodeDecodeError:
        return []
    lines = text.splitlines()
    loops: list[Loop] = []
    i = 0
    while i < len(lines):
        raw = lines[i]
        code = code_part(raw)
        if DC_RE.search(code):
            joined, last = join_continuations(lines, i)
            parsed = parse_dc_header(joined)
            if parsed is not None:
                n_idx, clauses, mask = parsed
                loops.append(Loop(
                    path=path,
                    line=i + 1,
                    header=joined.strip(),
                    n_indices=n_idx,
                    clauses=clauses,
                    has_mask=mask,
                    above=classify_above(lines, i),
                ))
            i = last + 1
        else:
            i += 1
    return loops


def _files(roots: list[Path]) -> list[Path]:
    files: list[Path] = []
    for root in roots:
        if root.is_file():
            files.append(root)
            continue
        files.extend(sorted(root.rglob("*.F90")))
        files.extend(sorted(root.rglob("*.f90")))
    return files


def walk(roots: list[Path], jobs: int = 1) -> list[Loop]:
    loops: list[Loop] = []
    for got in pmap(scan_file, _files(roots), jobs, "dc_audit"):
        loops.extend(got)
    return loops


def report(loops: list[Loop]) -> None:
    if not loops:
        print("No `do concurrent` loops found.")
        return
    n = len(loops)
    print(f"# do concurrent audit\n")
    print(f"Total loops: {n}\n")

    print("## What sits directly above\n")
    above = Counter(l.above for l in loops)
    for k in ("acc", "omp", "comment", "code", "nothing"):
        v = above.get(k, 0)
        pct = 100 * v / n
        print(f"  {k:<8} {v:5d}  ({pct:5.1f}%)")

    print("\n## Index count\n")
    idx = Counter(l.n_indices for l in loops)
    for k in sorted(idx):
        print(f"  {k} indices: {idx[k]}")

    print("\n## Clauses\n")
    clause_count = Counter()
    for l in loops:
        for c in l.clauses:
            clause_count[c] += 1
    if clause_count:
        for c, v in clause_count.most_common():
            print(f"  {c:<11} {v}")
    else:
        print("  (none)")
    masks = sum(1 for l in loops if l.has_mask)
    print(f"  mask        {masks}")

    print("\n## Top files by loop count\n")
    by_file = Counter(str(l.path) for l in loops)
    for f, v in by_file.most_common(15):
        print(f"  {v:4d}  {f}")

    print("\n## Loops needing attention (mask, or unusual clauses)\n")
    flagged = [l for l in loops if l.has_mask or "default" in l.clauses or "shared" in l.clauses]
    if not flagged:
        print("  (none)")
    for l in flagged[:50]:
        print(f"  {l.path}:{l.line}  mask={l.has_mask} clauses={l.clauses}")
    if len(flagged) > 50:
        print(f"  ... and {len(flagged) - 50} more")

    print("\n## Loops with no acc/omp directive above (would lose parallelism)\n")
    bare = [l for l in loops if l.above in ("nothing", "code", "comment")]
    print(f"  {len(bare)} of {n} ({100*len(bare)/n:.1f}%)")


def main() -> int:
    raw = sys.argv[1:]
    strict = False
    if "--strict" in raw:
        strict = True
        raw = [a for a in raw if a != "--strict"]
    jobs = 1
    for flag in ("-j", "--jobs", "-t"):
        if flag in raw:
            i = raw.index(flag)
            try:
                jobs = int(raw[i + 1])
            except (IndexError, ValueError):
                print(f"error: {flag} needs an integer argument", file=sys.stderr)
                return 1
            raw = raw[:i] + raw[i + 2:]
    if "--help" in raw or "-h" in raw:
        print("usage: dc_audit.py [--strict] [-j N] [PATH ...]")
        print("  --strict   exit non-zero if any loop is flagged (mask or odd clauses)")
        print("  -j N       scan across N worker processes (0 = one per core)")
        return 0
    args = raw or ["src"]
    roots = [Path(a) for a in args]
    for r in roots:
        if not r.exists():
            print(f"error: {r} does not exist", file=sys.stderr)
            return 1
    loops = walk(roots, resolve_jobs(jobs))
    report(loops)
    if strict:
        flagged = [l for l in loops if l.has_mask or "default" in l.clauses or "shared" in l.clauses]
        if flagged:
            print(f"\n[dc_audit] --strict: {len(flagged)} loop(s) flagged", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
