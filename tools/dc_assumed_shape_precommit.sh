#!/usr/bin/env bash
# Pre-commit wrapper around the assumed-shape lint.
#
# The lint is diff-aware (only NEW offenders fail), which needs a base ref.
# Locally that is origin/main; some checkouts don't have it (shallow CI
# clones, fresh worktrees).  In that case resolve the first available base
# ref.  If none is available, skip cleanly — the dedicated CI workflow is
# the authoritative gate.
#
# Mirror of tools/openmp_portability_precommit.sh.
set -euo pipefail

base_ref=""
for cand in origin/main main; do
    if git rev-parse --verify --quiet "${cand}^{commit}" >/dev/null 2>&1; then
        base_ref="${cand}"
        break
    fi
done

if [ -z "${base_ref}" ]; then
    echo "dc-assumed-shape: no base ref (origin/main) in this checkout — skipping."
    echo "  (the dc-assumed-shape CI workflow is the authoritative gate)"
    exit 0
fi

exec python3 tools/dc_assumed_shape_lint.py src/ app/ --diff "${base_ref}"
