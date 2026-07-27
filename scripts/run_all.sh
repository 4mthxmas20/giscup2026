#!/bin/bash
# End-to-end GISCUP 2026 pipeline: geojson -> solutions -> verify -> zip.
#
# Usage: scripts/run_all.sh data/input.geojson [lstime_seconds] [R] [spacing]
set -euo pipefail
cd "$(dirname "$0")/.."

GEOJSON="${1:?usage: run_all.sh input.geojson [lstime] [R] [spacing]}"
LSTIME="${2:-120}"
R="${3:-1200}"
SPACING="${4:-6}"
PY=.venv/bin/python
[ -x "$PY" ] || PY=python3

echo "== build =="
make -C cpp

echo "== prepare =="
mkdir -p results cache
$PY python/prepare.py "$GEOJSON" results/buildings.txt

echo "== solve (R=$R spacing=$SPACING lstime=${LSTIME}s per config) =="
cpp/solver results/buildings.txt --R "$R" --spacing "$SPACING" \
    --lstime "$LSTIME" --witness --out results

echo "== verify =="
$PY python/verify.py "$GEOJSON" results --samples 20

echo "== package =="
$PY python/make_submission.py results results/submission.zip

echo "== done =="
grep -c . results/solutions.txt >/dev/null && echo "submission at results/submission.zip"
