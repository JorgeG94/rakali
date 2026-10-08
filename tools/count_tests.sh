#!/usr/bin/env bash
# count_tests.sh - Count unit tests per testsuite in tests/unit/

UNIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../tests/" && pwd)"

printf "%-45s %s\n" "Testsuite" "Tests"
printf "%-45s %s\n" "-----------------------------------------" "-----"

total_tests=0
total_suites=0

for f in "$UNIT_DIR"/test_*.F90; do
    filename="$(basename "$f")"

    # Count new_unittest calls inside the testsuite array assignment
    count=$(grep -c 'new_unittest(' "$f" 2>/dev/null || echo 0)

    # Extract the collect_* subroutine name as the suite label
    suite=$(grep -oE 'collect_[a-zA-Z0-9_]+' "$f" | head -1)
    suite="${suite:-${filename%.F90}}"

    printf "%-45s %d\n" "$suite" "$count"

    total_tests=$((total_tests + count))
    total_suites=$((total_suites + 1))
done

printf "%-45s %s\n" "-----------------------------------------" "-----"
printf "%-45s %d\n" "TOTAL ($total_suites suites)" "$total_tests"
