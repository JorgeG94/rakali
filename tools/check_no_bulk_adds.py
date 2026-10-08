#!/usr/bin/env python3
"""Refuse commits that add a suspiciously large number of NEW files.

Why this exists
---------------
This worktree carries dozens of untracked files at any time -- scratch runs,
draft plans, in-progress examples.  `git add -A` sweeps all of them into the
index, and the resulting commit silently carries unrelated work.  That has
happened repeatedly; the standing "use `git add -u`" rule is not a control,
because it depends on remembering it at exactly the wrong moment.

`.gitignore` removes the unambiguous scratch from the addable set, but the
pending work (docs/ plans, examples/, validation_examples/) deliberately stays
visible -- it might still get committed on purpose.  This hook guards THAT
remainder.

The heuristic
-------------
A deliberate commit adds a handful of new files: a module and its test, a doc,
a validation case.  An accidental `git add -A` adds dozens spanning unrelated
top-level directories.  So: fail when a commit adds more than MAX_ADDS new
files, and report them so the author can see what got swept in.

This is intentionally a blunt count rather than a cleverer "are these files
related" test -- a blunt rule that fires reliably beats a subtle one that can
be argued with.

Override
--------
Legitimate bulk additions exist (vendoring, a new validation suite, the
coastal/ocean split).  Set RAKALI_ALLOW_BULK_ADD=1 for that commit:

    RAKALI_ALLOW_BULK_ADD=1 git commit ...

The override is deliberately explicit and per-commit: it forces the decision to
be conscious rather than habitual.
"""
from __future__ import annotations

import os
import subprocess
import sys

MAX_ADDS = 8


def staged_added() -> list[str]:
    """Paths staged as NEW files (status A), via the index rather than argv.

    pre-commit passes the changed filenames, but not their status, so ask git
    directly -- otherwise a large edit-only commit would trip the check.
    """
    out = subprocess.run(
        ["git", "diff", "--cached", "--name-status", "--diff-filter=A"],
        capture_output=True, text=True, check=True,
    ).stdout
    return [ln.split("\t", 1)[1] for ln in out.splitlines() if "\t" in ln]


def main() -> int:
    if os.environ.get("RAKALI_ALLOW_BULK_ADD") == "1":
        return 0

    added = staged_added()
    if len(added) <= MAX_ADDS:
        return 0

    top = sorted({p.split("/", 1)[0] for p in added})
    print("check_no_bulk_adds: FAIL")
    print(f"  this commit adds {len(added)} new files "
          f"across {len(top)} top-level paths: {', '.join(top)}")
    print()
    for p in sorted(added)[:20]:
        print(f"    A  {p}")
    if len(added) > 20:
        print(f"    ... and {len(added) - 20} more")
    print()
    print("  If this was `git add -A`, it almost certainly swept in unrelated")
    print("  untracked files.  Prefer:")
    print("      git reset && git add -u && git add <the new files you meant>")
    print()
    print("  If the bulk addition IS intended, say so explicitly:")
    print("      RAKALI_ALLOW_BULK_ADD=1 git commit ...")
    return 1


if __name__ == "__main__":
    sys.exit(main())
