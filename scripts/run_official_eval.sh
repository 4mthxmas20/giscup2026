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

echo "== evaluating =="
DATASET="$DATASET" SUBMISSION="$SUBMISSION" \
    REPORT="$REPO_ROOT/results/official-eval-report.json" \
    npx vitest run --config vitest.benchmark.config.ts \
    benchmarks/official_eval.test.ts --reporter=verbose 2>&1 | tee "$REPO_ROOT/results/official-eval.log"

echo "== report: $REPO_ROOT/results/official-eval-report.json =="
