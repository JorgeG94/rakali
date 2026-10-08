#!/usr/bin/env bash
# Pre-commit wrapper around the OpenMP portability lint.
#
# The lint is diff-aware (only NEW offenders fail), which needs a base ref to
# diff against.  Locally that is origin/main; but some checkouts don't have it
# resolvable — the generic `pre-commit run --all-files` CI job uses a shallow
# clone, and fresh worktrees may not track origin/main.  In that case
# `git diff origin/main` aborts with exit 128 and the lint crashes.
#
# Resolve the first base ref that exists and diff against it.  If none is
# available, skip cleanly: the dedicated openmp-portability-lint.yml workflow
# (full-history checkout) is the authoritative gate for PRs to main.
set -euo pipefail

base_ref=""
for cand in origin/main main; do
    if git rev-parse --verify --quiet "${cand}^{commit}" >/dev/null 2>&1; then
        base_ref="${cand}"
        break
    fi
done

if [ -z "${base_ref}" ]; then
    echo "openmp-portability: no base ref (origin/main) in this checkout — skipping."
    echo "  (the openmp-portability-lint.yml CI workflow is the authoritative gate)"
    exit 0
fi

exec python3 tools/openmp_portability_lint.py src/ app/ --diff "${base_ref}" --fail-on-warn
