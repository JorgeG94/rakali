#!/usr/bin/env python3
"""openmp_portability_lint.py — flag !$acc directives that break under OpenMP.

Driven by tools/openmp_portability_rules.yaml. Scans Fortran source files for
directive shapes known to fail on specific OpenMP compilers (today: ifx on
Aurora). Used as a PR gate on `main` so portability-hostile directives are
caught before they propagate into the auto-regenerated `auto/dc-openmp`
variant.

Usage:
    python tools/openmp_portability_lint.py PATH [PATH ...]
    python tools/openmp_portability_lint.py src/ app/ --fail-on-warn
    python tools/openmp_portability_lint.py --diff origin/main src/ app/

Options:
    --fail-on-warn        Exit non-zero when any rule trips at warn severity.
                          Error-severity hits always exit non-zero.
    --diff REF            Only consider lines that differ vs REF (for PRs).
                          Without this flag the entire file is scanned —
                          useful for backfilling rules across the codebase.
    --rules FILE          Override the default rules file location.
    --list-rules          Print loaded rules and exit.

YAML is loaded with PyYAML if available, otherwise a tiny inline parser is
used (sufficient for the simple flat-list schema we have today).
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

RULES_DEFAULT = Path(__file__).parent / "openmp_portability_rules.yaml"
SOURCE_SUFFIXES = (".F90", ".f90")


@dataclass
class Rule:
    id: str
    match: re.Pattern[str]
    compiler: list[str]
    severity: str  # "error" or "warn"
    message: str
    fix: str
    evidence: str


def load_rules(path: Path) -> list[Rule]:
    text = path.read_text()
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
    except ImportError:
        data = _minimal_yaml_load(text)
    rules: list[Rule] = []
    for r in data.get("rules", []):
        rules.append(
            Rule(
                id=r["id"],
                match=re.compile(r["match"]),
                compiler=list(r.get("compiler", [])),
                severity=r.get("severity", "error"),
                message=" ".join(r.get("message", "").split()),
                fix=" ".join(r.get("fix", "").split()),
                evidence=" ".join(r.get("evidence", "").split()),
            )
        )
    return rules


def _minimal_yaml_load(text: str) -> dict:
    """Tiny YAML reader for our flat list-of-dicts-with-scalars schema.

    Handles top-level `rules:`, `- id: ...` list items, `key: value` lines,
    `key: [a, b]` inline lists, and `key: >` folded blocks. Not a general YAML
    parser — purposely minimal so the lint runs without PyYAML.
    """
    out: dict = {"rules": []}
    cur: dict | None = None
    pending_key: str | None = None
    pending_block: list[str] = []

    def flush_block():
        nonlocal pending_key, pending_block, cur
        if pending_key is not None and cur is not None:
            cur[pending_key] = " ".join(s.strip() for s in pending_block).strip()
        pending_key = None
        pending_block = []

    for raw in text.splitlines():
        line = raw.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue
        if pending_key is not None:
            # Folded-block (`>`) continuation lines are indented DEEPER than the
            # 4-space keys (`fix:`/`evidence:` sit at 4 spaces too); keying off a
            # bare 4-space prefix swallows the next key into the block, so a rule
            # loses `fix`/`evidence` -> KeyError. Require indent > 4.
            indent = len(line) - len(line.lstrip())
            if indent > 4 or line.startswith("\t"):
                pending_block.append(line.strip())
                continue
            flush_block()
        if line.startswith("rules:"):
            continue
        if line.startswith("  - id:"):
            if cur is not None:
                out["rules"].append(cur)
            cur = {"id": line.split("id:", 1)[1].strip()}
            continue
        if line.startswith("    ") and cur is not None and ":" in line:
            k, _, v = line.strip().partition(":")
            v = v.strip()
            if v == ">":
                pending_key = k
                pending_block = []
            elif v.startswith("[") and v.endswith("]"):
                cur[k] = [s.strip() for s in v[1:-1].split(",") if s.strip()]
            else:
                cur[k] = v.strip("'\"")
    flush_block()
    if cur is not None:
        out["rules"].append(cur)
    return out


def changed_lines(ref: str) -> dict[Path, set[int]]:
    """Lines added vs `ref` in each file. {path -> {line numbers}}."""
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
                for i in range(start, start + n):
                    result[current].add(i)
    return result


def walk_files(paths: list[Path]) -> Iterable[Path]:
    for p in paths:
        if p.is_file() and p.suffix in SOURCE_SUFFIXES:
            yield p
        elif p.is_dir():
            for f in p.rglob("*"):
                if f.suffix in SOURCE_SUFFIXES:
                    yield f


def lint(paths: list[Path], rules: list[Rule],
         diff_lines: dict[Path, set[int]] | None) -> list[tuple[str, Path, int, Rule, str]]:
    """Returns list of (severity, path, lineno, rule, snippet)."""
    hits: list[tuple[str, Path, int, Rule, str]] = []
    for f in walk_files(paths):
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        line_filter = None
        if diff_lines is not None:
            rel = f.resolve().relative_to(Path.cwd().resolve()) if f.is_absolute() else f
            line_filter = diff_lines.get(rel)
            if not line_filter:
                # --diff was passed AND this file has no added lines vs REF
                # (file not in diff, or only deletions). Skip entirely.
                continue
        for n, line in enumerate(text.splitlines(), 1):
            if line_filter is not None and n not in line_filter:
                continue
            for r in rules:
                if r.match.search(line):
                    hits.append((r.severity, f, n, r, line.strip()))
    return hits


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="*", default=["src", "app"], type=Path)
    ap.add_argument("--fail-on-warn", action="store_true")
    ap.add_argument("--diff", metavar="REF",
                    help="only scan lines that differ vs REF")
    ap.add_argument("--rules", type=Path, default=RULES_DEFAULT)
    ap.add_argument("--list-rules", action="store_true")
    args = ap.parse_args()

    rules = load_rules(args.rules)
    if args.list_rules:
        for r in rules:
            print(f"{r.severity:5s}  {r.id:30s}  [{','.join(r.compiler)}]")
        return 0

    diff_lines = changed_lines(args.diff) if args.diff else None
    hits = lint(args.paths, rules, diff_lines)

    n_err = sum(1 for s, *_ in hits if s == "error")
    n_warn = sum(1 for s, *_ in hits if s == "warn")

    for sev, path, ln, rule, snippet in hits:
        prefix = "ERROR" if sev == "error" else "WARN "
        print(f"{prefix}  {path}:{ln}  [{rule.id}]")
        print(f"        {snippet}")
        print(f"        {rule.message}")
        print(f"        FIX: {rule.fix}")
        if rule.evidence:
            print(f"        EVIDENCE: {rule.evidence}")
        print()

    if not hits:
        print("openmp_portability_lint: clean")
        return 0
    print(f"openmp_portability_lint: {n_err} error(s), {n_warn} warning(s)")
    if n_err > 0:
        return 1
    if n_warn > 0 and args.fail_on_warn:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
