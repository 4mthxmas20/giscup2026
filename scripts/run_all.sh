#!/bin/bash
# End-to-end GISCUP 2026 pipeline: geojson -> solutions -> verify -> zip.
#
# Usage: scripts/run_all.sh data/input.geojson [lstime_seconds] [R] [spacing]
#
# The competition ships its own parameter file, so tau/k are overridable:
#   TAUS=0.3,0.6 KS=100,1000 scripts/run_all.sh comp.geojson 600
set -euo pipefail
cd "$(dirname "$0")/.."

GEOJSON="${1:?usage: run_all.sh input.geojson [lstime] [R] [spacing]}"
LSTIME="${2:-120}"
R="${3:-1200}"
SPACING="${4:-6}"
TAUS="${TAUS:-0.25,0.5,0.75}"
KS="${KS:-50,500,1000}"
PY=.venv/bin/python
[ -x "$PY" ] || PY=python3

echo "== build =="
make -C cpp

echo "== prepare =="
mkdir -p results cache
$PY python/prepare.py "$GEOJSON" results/buildings.txt

echo "== solve (R=$R spacing=$SPACING lstime=${LSTIME}s, taus=$TAUS ks=$KS) =="
rm -f results/sol_t*_k*.txt results/wit_t*_k*.txt
cpp/solver results/buildings.txt --R "$R" --spacing "$SPACING" \
    --lstime "$LSTIME" --witness --taus "$TAUS" --ks "$KS" --out results

echo "== verify =="
$PY python/verify.py "$GEOJSON" results --samples 20

echo "== package =="
$PY python/make_submission.py results results/submission.zip

echo "== done =="
grep -c . results/solutions.txt >/dev/null && echo "submission at results/submission.zip"
