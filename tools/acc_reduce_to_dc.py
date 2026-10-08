#!/usr/bin/env python3
"""Convert `!$acc parallel loop ... reduction(op:vars)` to F2023
`do concurrent(...) reduce(op:vars)`.

Why: a reduction expressed as an OpenACC directive forces the enclosing
procedure to carry a compute directive, which (a) makes the OpenMP-target
translation put an `!$omp target` region there -- illegal in a `pure`
procedure and fatal for any `do concurrent` that calls it -- and (b) needs a
translator at all.  Expressed as `do concurrent ... reduce`, the same loop is
plain Fortran: portable to nvfortran / gfortran / ifx / flang, `pure`-safe,
and needs no translation for the OpenMP variant.

Deliberately NOT converted:
  * directives carrying `async(...)` -- those are the CUDA-graph batching
    sites; async is an OpenACC-only accelerator we keep on NVIDIA and strip
    elsewhere, so they must stay `!$acc`.
  * anything the loop-structure parser cannot balance (reported, left alone).

`present(...)` / `default(...)` clauses are dropped: they are OpenACC data
assertions with no `do concurrent` equivalent (residency is the caller's job
via enter/exit data).  `private(...)` becomes `local(...)`.

Usage: python3 tools/acc_reduce_to_dc.py [--write] [paths...]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from omp_to_doconcurrent import (  # noqa: E402
    split_top_level_commas, find_clause, DO_RE, END_DO_RE,
    DO_CONCURRENT_RE, BLANK_OR_COMMENT_RE,
)

ACC_OPEN_RE = re.compile(r"^(\s*)(!\$acc\s+parallel\s+loop)\b(.*)$", re.I)
# OpenACC free-form continuation: the next line may start `!$acc&` OR bare
# `!$acc` (both legal).  Only consume it when the PREVIOUS line actually
# ended in `&`, else a following directive would be swallowed.
ACC_CONT_RE = re.compile(r"^\s*!\$acc&?\s*(.*)$", re.I)
REDUCTION_RE = re.compile(r"\breduction\s*\(\s*([^:)]+?)\s*:\s*([^)]*)\)", re.I)


def convert(lines: list[str]) -> tuple[list[str], int, list[str]]:
    out: list[str] = []
    warnings: list[str] = []
    n_conv = 0
    i = 0
    while i < len(lines):
        m = ACC_OPEN_RE.match(lines[i])
        if not m:
            out.append(lines[i]); i += 1; continue

        # ---- flatten directive + continuations ----
        directive_idxs = [i]
        raw = m.group(3).rstrip()
        cont_expected = raw.endswith("&")
        flat = raw[:-1].rstrip() if cont_expected else raw
        j = i + 1
        while cont_expected and j < len(lines):
            cm = ACC_CONT_RE.match(lines[j])
            if not cm:
                break
            cont = cm.group(1).rstrip()
            cont_expected = cont.endswith("&")
            if cont_expected:
                cont = cont[:-1].rstrip()
            flat += " " + cont
            directive_idxs.append(j)
            j += 1

        def keep(reason: str = ""):
            if reason:
                warnings.append(f"line {i+1}: {reason}")
            for idx in directive_idxs:
                out.append(lines[idx])

        reds = REDUCTION_RE.findall(flat)
        if not reds:
            keep(); i = j; continue
        if re.search(r"\basync\s*\(", flat, re.I):
            keep(); i = j; continue          # async stays OpenACC by design

        collapse_n = 1
        cm2 = re.search(r"\bcollapse\s*\(\s*(\d+)\s*\)", flat, re.I)
        if cm2:
            collapse_n = int(cm2.group(1))

        priv = find_clause(flat, "private")
        private_vars = split_top_level_commas(priv[2]) if priv else []

        # ---- locate the collapse_n `do` headers ----
        loop_specs: list[tuple[str, str, str]] = []
        loop_idxs: list[int] = []
        k = j
        while len(loop_specs) < collapse_n and k < len(lines):
            dm = DO_RE.match(lines[k])
            if not dm:
                break
            loop_idxs.append(k)
            loop_specs.append((dm.group(2), dm.group(3).strip(), dm.group(4).strip()))
            k += 1
        if len(loop_specs) != collapse_n:
            keep(f"could not parse {collapse_n} `do` headers"); i = j; continue

        # ---- balance the loop nest ----
        end_outer: list[int] = []
        depth = collapse_n
        scan = k
        while scan < len(lines):
            ls = lines[scan]
            if DO_CONCURRENT_RE.match(ls) or DO_RE.match(ls):
                depth += 1
            elif END_DO_RE.match(ls):
                if depth <= collapse_n:
                    end_outer.append(scan)
                depth -= 1
                if depth == 0:
                    break
            scan += 1
        if depth != 0 or len(end_outer) != collapse_n:
            keep(f"could not balance {collapse_n} loops (depth={depth})")
            i = j; continue

        # ---- emit ----
        spec = ", ".join(f"{v}={lo}:{hi}" for v, lo, hi in loop_specs)
        red_str = " ".join(f"reduce({op.strip()}:{vs.strip()})" for op, vs in reds)
        indent = re.match(r"^(\s*)", lines[loop_idxs[0]]).group(1)
        head = f"{indent}do concurrent({spec})"
        clauses = ([f"local({', '.join(private_vars)})"] if private_vars else []) + \
                  ([red_str] if red_str else [])
        oneline = head + ("".join(" " + c for c in clauses))
        LIMIT = 110          # keep well under fortitude S001 (178)
        if len(oneline) <= LIMIT:
            new_do = [oneline]
        else:
            # Wrap clause-by-clause, and if a single clause (a long reduce
            # list) still overflows, break it at commas.
            new_do, cont = [head + " &"], f"{indent}   "
            parts: list[str] = []
            for c in clauses:
                if len(cont) + len(c) <= LIMIT:
                    parts.append(c)
                    continue
                bits, cur = c.split(", "), ""
                for b in bits:
                    cand = (cur + ", " + b) if cur else b
                    if len(cont) + len(cand) > LIMIT and cur:
                        parts.append(cur + ", &"); cur = b
                    else:
                        cur = cand
                if cur:
                    parts.append(cur)
            for n_, pc in enumerate(parts):
                sep = "" if (pc.endswith("&") or n_ == len(parts) - 1) else " &"
                new_do.append(f"{cont}{pc}{sep}")

        body_start, body_end = k, end_outer[0]
        dedent = "   " * (collapse_n - 1)
        body = []
        for idx in range(body_start, body_end):
            ln = lines[idx]
            body.append(ln[len(dedent):] if dedent and ln.startswith(dedent + " ") else ln)

        out.extend(new_do)
        out.extend(body)
        out.append(f"{indent}end do")
        n_conv += 1
        i = end_outer[-1] + 1
    return out, n_conv, warnings


def main() -> int:
    args = sys.argv[1:]
    write = "--write" in args
    roots = [a for a in args if not a.startswith("--")] or ["src"]
    files: list[Path] = []
    for r in roots:
        p = Path(r)
        files += sorted(p.rglob("*.F90")) if p.is_dir() else [p]
    total, nfiles, allw = 0, 0, []
    for f in files:
        text = f.read_text()
        if "!$acc parallel loop" not in text:
            continue
        new, n, warns = convert(text.splitlines())
        allw += [f"{f}: {w}" for w in warns]
        if n:
            total += n; nfiles += 1
            print(f"  {n:3d}  {f}")
            if write:
                f.write_text("\n".join(new) + ("\n" if text.endswith("\n") else ""))
    for w in allw:
        print(f"  SKIP {w}", file=sys.stderr)
    print(f"\n{total} reductions converted across {nfiles} files"
          f"{'' if write else '  (dry run -- pass --write)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
