#!/bin/bash
# Run the OFFICIAL GIS Cup 2026 evaluator against our submission.
#
# Clones (or reuses) the organizers' evaluator, installs its pinned deps, and
# evaluates results/solutions.txt with their ArcGIS-backed geometry engine.
# Intended for a Codespace / CI runner so the ~400 MB of node deps stay off the
# local machine.
#
# Usage: scripts/run_official_eval.sh [dataset.geojson] [solutions.txt]
set -euo pipefail
cd "$(dirname "$0")/.."
REPO_ROOT="$PWD"

DATASET="${1:-$REPO_ROOT/data/GIS-cup-sample-dataset.geojson}"
SUBMISSION="${2:-$REPO_ROOT/results/solutions.txt}"
EVAL_DIR="${EVAL_DIR:-$HOME/gis-cup-evaluator}"

[ -f "$DATASET" ] || { echo "missing dataset: $DATASET" >&2; exit 1; }
[ -f "$SUBMISSION" ] || { echo "missing submission: $SUBMISSION" >&2; exit 1; }

# The evaluator runs from its own checkout, so relative inputs would resolve
# against the wrong directory once we cd there.
abspath() { (cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")"); }
DATASET="$(abspath "$DATASET")"
SUBMISSION="$(abspath "$SUBMISSION")"

if [ ! -d "$EVAL_DIR/.git" ]; then
    echo "== cloning official evaluator =="
    git clone --depth 1 https://github.com/alowe/gis-cup-2026-evaluator "$EVAL_DIR"
fi

cd "$EVAL_DIR"
echo "== installing pinned evaluator deps =="
if command -v pnpm >/dev/null 2>&1; then
    pnpm install --frozen-lockfile
else
    npx --yes pnpm@11.16.0 install --frozen-lockfile
fi

cp "$REPO_ROOT/scripts/official_eval.test.ts" benchmarks/official_eval.test.ts

# One vitest run per sub-problem. Evaluating all nine in a single test blows
# past vitest's per-test timeout, and a late failure would discard the earlier
# results; per-sub-problem runs finish well inside the limit and each one's
# verdict is durable.
COUNT=$(( $(grep -c '' "$SUBMISSION") / 3 ))
echo "== evaluating $COUNT sub-problems =="
LOG="$REPO_ROOT/results/official-eval.log"
: > "$LOG"
FAILED=0
for i in $(seq 1 "$COUNT"); do
    echo "-- sub-problem $i/$COUNT --" | tee -a "$LOG"
    if ! DATASET="$DATASET" SUBMISSION="$SUBMISSION" SUBPROBLEM="$i" \
        REPORT="$REPO_ROOT/results/official-eval-report-$i.json" \
        npx vitest run --config vitest.benchmark.config.ts \
        --testTimeout=3600000 \
        benchmarks/official_eval.test.ts --reporter=verbose 2>&1 | tee -a "$LOG"
    then
        FAILED=1
    fi
done

echo "== summary ==" | tee -a "$LOG"
grep -oE '"tau":[0-9.]+,"k":[0-9]+.*"lost":[0-9]+' "$LOG" \
    | sed 's/"antennas[A-Za-z]*":[0-9]*,//g' | tee -a "$LOG"
if [ "$FAILED" -ne 0 ]; then
    echo "SOME SUB-PROBLEMS FAILED — do not submit until resolved" | tee -a "$LOG"
    exit 1
fi
echo "all sub-problems verified; reports in $REPO_ROOT/results/" | tee -a "$LOG"
