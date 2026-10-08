#!/usr/bin/env python3
"""Forbid `transfer()` inside a `do concurrent` body.

nvfortran 26.5 MISCOMPILES it on `-stdpar=gpu`: a plain elementwise

    do concurrent(i=1:N)
       got(i) = transfer(x(i), 0_int64)
    end do

gets ~95% of elements wrong (94645/100000), reading a SHIFTED element --
`got 0x4008000000000000` (3.0) where `0x4000000000000000` (2.0) was expected.
No reduction is involved; it is an index/scoping miscompile in the
do-concurrent lowering.  The identical expression under `!$acc parallel loop`
is CORRECT, which is why chksum-style kernels must stay on OpenACC.

Silent wrong answers, no crash, no warning -- exactly the class a lint has to
catch, because nothing downstream will.

Fix: keep the loop as `!$acc parallel loop`, or move the bit manipulation
host-side.  `popcnt` alone inside `do concurrent` is fine; `transfer` is the
trigger.

Usage: python3 tools/dc_transfer_lint.py [paths...]   (default: src)
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

DC_OPEN = re.compile(r"^\s*do\s+concurrent\b", re.I)
DO_ANY = re.compile(r"^\s*do\b", re.I)
END_DO = re.compile(r"^\s*end\s+do\b", re.I)
TRANSFER = re.compile(r"\btransfer\s*\(", re.I)
WAIVER = "dc-transfer-ok"


def strip_comment(line: str) -> str:
    out, in_str, q = [], False, ""
    for ch in line:
        if in_str:
            out.append(ch)
            if ch == q:
                in_str = False
        elif ch in "'\"":
            in_str, q = True, ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def scan(path: Path) -> list[tuple[int, str]]:
    lines = path.read_text(errors="replace").splitlines()
    hits, depth = [], 0
    for n, raw in enumerate(lines, 1):
        code = strip_comment(raw)
        if depth:
            if TRANSFER.search(code) and WAIVER not in raw:
                hits.append((n, raw.strip()[:88]))
            if END_DO.match(code):
                depth -= 1
                continue
            if DO_ANY.match(code):
                depth += 1
            continue
        if DC_OPEN.match(code):
            depth = 1
    return hits


def main() -> int:
    roots = [a for a in sys.argv[1:] if not a.startswith("-")] or ["src"]
    files: list[Path] = []
    for r in roots:
        p = Path(r)
        files += sorted(p.rglob("*.F90")) if p.is_dir() else [p]
    bad = 0
    for f in files:
        for n, text in scan(f):
            bad += 1
            print(f"ERROR  {f}:{n}  [dc-transfer]  `transfer()` inside a "
                  f"`do concurrent` -- nvfortran miscompiles this on -stdpar=gpu")
            print(f"       {text}")
    if bad:
        print(f"\n{bad} occurrence(s).  FIX: keep the loop as `!$acc parallel loop`, "
              f"or move the bit manipulation host-side.")
        print(f"WAIVER: add `{WAIVER}` in a comment on the offending line if you "
              f"have verified your toolchain.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
