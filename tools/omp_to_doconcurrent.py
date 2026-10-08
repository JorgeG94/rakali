#!/usr/bin/env python3
"""Convert `!$omp target teams distribute parallel do` (without reduction) to
`do concurrent`. Loops with `reduction(...)` clauses are left as-is.

Usage: python tools/omp_to_doconcurrent.py [--check] file [file ...]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


OMP_OPEN_RE = re.compile(
    r"^(\s*)(!\$omp\s+target\s+teams\s+distribute\s+parallel\s+do)\b(.*)$",
    re.IGNORECASE,
)
OMP_CONT_RE = re.compile(r"^\s*!\$omp&\s*(.*)$", re.IGNORECASE)
OMP_END_RE = re.compile(
    r"^\s*!\$omp\s+end\s+target\s+teams\s+distribute\s+parallel\s+do\b.*$",
    re.IGNORECASE,
)
DO_RE = re.compile(
    r"^(\s*)do\s+(\w+)\s*=\s*(.+?)\s*,\s*(.+?)\s*$",
    re.IGNORECASE,
)
DO_CONCURRENT_RE = re.compile(r"^\s*do\s+concurrent\b", re.IGNORECASE)
END_DO_RE = re.compile(r"^\s*end\s+do\s*$", re.IGNORECASE)
BLANK_OR_COMMENT_RE = re.compile(r"^\s*(!.*)?$")


def split_top_level_commas(s: str) -> list[str]:
    """Split on commas that are not inside parentheses."""
    parts = []
    cur = []
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        parts.append(tail)
    return parts


def find_clause(text: str, name: str) -> tuple[int, int, str] | None:
    """Find a `name(...)` clause in text. Return (start, end, body) or None."""
    pat = re.compile(r"\b" + re.escape(name) + r"\s*\(", re.IGNORECASE)
    m = pat.search(text)
    if not m:
        return None
    start = m.start()
    # Find matching close paren
    i = m.end() - 1  # at the '('
    depth = 0
    while i < len(text):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                body = text[m.end():i]
                return (start, i + 1, body.strip())
        i += 1
    return None


def convert(text: str) -> tuple[str, int, list[str]]:
    """Return (new_text, n_converted, warnings)."""
    lines = text.split("\n")
    out: list[str] = []
    warnings: list[str] = []
    i = 0
    n_converted = 0

    while i < len(lines):
        line = lines[i]
        m = OMP_OPEN_RE.match(line)
        if not m:
            out.append(line)
            i += 1
            continue

        indent = m.group(1)
        rest = m.group(3)

        # Collect continuation lines
        directive_idxs = [i]
        # Build flattened directive text
        flat = rest.rstrip()
        # Strip trailing & for line continuation joining
        if flat.endswith("&"):
            flat = flat[:-1].rstrip()
        j = i + 1
        while j < len(lines):
            cm = OMP_CONT_RE.match(lines[j])
            if not cm:
                break
            cont = cm.group(1).rstrip()
            had_amp = cont.endswith("&")
            if had_amp:
                cont = cont[:-1].rstrip()
            flat = flat + " " + cont
            directive_idxs.append(j)
            j += 1

        # Reduction? leave as-is.
        if re.search(r"\breduction\s*\(", flat, re.IGNORECASE):
            for idx in directive_idxs:
                out.append(lines[idx])
            i = j
            continue

        # Parse collapse(N)
        collapse_n = 1
        cm2 = re.search(r"\bcollapse\s*\(\s*(\d+)\s*\)", flat, re.IGNORECASE)
        if cm2:
            collapse_n = int(cm2.group(1))

        # Parse private(...)
        priv_clause = find_clause(flat, "private")
        private_vars: list[str] = []
        if priv_clause is not None:
            private_vars = split_top_level_commas(priv_clause[2])

        # NOTE: any map(...) clauses are intentionally dropped — data is
        # managed via separate !$omp target enter/exit data directives.

        # Find next collapse_n `do <var> = <lo>, <hi>` lines
        loop_idxs: list[int] = []
        loop_specs: list[tuple[str, str, str]] = []  # (var, lo, hi)
        k = j
        while len(loop_specs) < collapse_n and k < len(lines):
            do_m = DO_RE.match(lines[k])
            if do_m:
                loop_idxs.append(k)
                loop_specs.append(
                    (do_m.group(2), do_m.group(3).strip(), do_m.group(4).strip())
                )
                k += 1
            elif BLANK_OR_COMMENT_RE.match(lines[k]):
                # tolerate blank/comment between directive and first do
                loop_idxs.append(k)  # keep position but don't consume as a loop
                # Actually we shouldn't auto-skip these. For safety, abort.
                break
            else:
                break

        if len(loop_specs) != collapse_n:
            warnings.append(
                f"line {i+1}: could not parse {collapse_n} `do` headers"
            )
            for idx in directive_idxs:
                out.append(lines[idx])
            i = j
            continue

        # Walk the body tracking nesting depth to find where the outer N
        # loops close. The `!$omp end target teams distribute parallel do`
        # directive is OPTIONAL in OpenMP, so we don't rely on finding it
        # — we infer the loop end from `do/end do` balancing.
        end_do_outer: list[int] = []  # indices of the N outer `end do` lines
        depth = collapse_n
        scan = k
        while scan < len(lines):
            ls = lines[scan]
            if DO_CONCURRENT_RE.match(ls):
                depth += 1
            elif DO_RE.match(ls):
                depth += 1
            elif END_DO_RE.match(ls):
                if depth <= collapse_n:
                    end_do_outer.append(scan)
                depth -= 1
                if depth == 0:
                    break
            scan += 1

        if depth != 0 or len(end_do_outer) != collapse_n:
            warnings.append(
                f"line {i+1}: could not balance {collapse_n} loops (depth={depth})"
            )
            for idx in directive_idxs:
                out.append(lines[idx])
            i = j
            continue

        end_do_idxs = sorted(end_do_outer)
        last_end_do = end_do_idxs[-1]

        # Optional `!$omp end target teams distribute parallel do` after
        # the outermost end-do — if present, drop it.
        end_idx: int | None = None
        nxt = last_end_do + 1
        while nxt < len(lines) and BLANK_OR_COMMENT_RE.match(lines[nxt]):
            # only skip blank lines, not comments — but a directive line is
            # also a "comment" by our regex; check explicitly first
            if OMP_END_RE.match(lines[nxt]):
                break
            if lines[nxt].strip() == "":
                nxt += 1
                continue
            break
        if nxt < len(lines) and OMP_END_RE.match(lines[nxt]):
            end_idx = nxt

        # Build the new `do concurrent(...)` line(s)
        spec_str = ", ".join(f"{v}={lo}:{hi}" for v, lo, hi in loop_specs)
        first_do_indent = re.match(r"^(\s*)", lines[loop_idxs[0]]).group(1)
        if private_vars:
            new_do_lines = [
                f"{first_do_indent}do concurrent({spec_str}) &",
                f"{first_do_indent}   local({', '.join(private_vars)})",
            ]
        else:
            new_do_lines = [f"{first_do_indent}do concurrent({spec_str})"]

        skip = set(directive_idxs)
        skip.update(loop_idxs[1:])
        skip.update(end_do_idxs[:-1])
        if end_idx is not None:
            skip.add(end_idx)

        last_to_emit = end_idx if end_idx is not None else last_end_do

        # Emit from i through last_to_emit
        for emit_idx in range(i, last_to_emit + 1):
            if emit_idx in skip:
                continue
            if emit_idx == loop_idxs[0]:
                out.extend(new_do_lines)
                continue
            out.append(lines[emit_idx])

        n_converted += 1
        i = last_to_emit + 1

    return "\n".join(out), n_converted, warnings


def main() -> int:
    args = sys.argv[1:]
    check_only = False
    if args and args[0] == "--check":
        check_only = True
        args = args[1:]
    if not args:
        print("usage: omp_to_doconcurrent.py [--check] file [file ...]", file=sys.stderr)
        return 2

    total_converted = 0
    any_warnings = False
    for path_str in args:
        path = Path(path_str)
        text = path.read_text()
        new_text, n_converted, warnings = convert(text)
        if warnings:
            any_warnings = True
            for w in warnings:
                print(f"{path}:{w}", file=sys.stderr)
        if n_converted == 0 and new_text == text:
            continue
        if check_only:
            print(f"{path}: would convert {n_converted} loops")
        else:
            path.write_text(new_text)
            print(f"{path}: converted {n_converted} loops")
        total_converted += n_converted

    print(f"Total: {total_converted} loops {'would be ' if check_only else ''}converted",
          file=sys.stderr)
    return 1 if any_warnings else 0


if __name__ == "__main__":
    sys.exit(main())
