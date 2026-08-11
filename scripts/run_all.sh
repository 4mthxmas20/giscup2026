#!/bin/bash
# End-to-end GISCUP 2026 pipeline: geojson -> solutions -> verify -> zip.
#
# Usage: scripts/run_all.sh data/input.geojson [lstime_seconds] [R] [spacing]
#
# The competition ships its own parameter file, so tau/k are overridable:
#   TAUS=0.3,0.6 KS=100,1000 scripts/run_all.sh comp.geojson 600
set -euo pipefail
cd "$(dirname "$0")/.."

GEOJSON="${1:?usage: run_all.sh input.geojson [budget_s] [R] [spacing]}"
# Per-config search budget, split between 1-opt and LNS. 1-opt saturates after
# a couple of minutes while LNS keeps paying, so cap the 1-opt share.
BUDGET="${2:-240}"
LSTIME=$(( BUDGET / 3 )); [ "$LSTIME" -gt 180 ] && LSTIME=180
LNSTIME=$(( BUDGET - LSTIME ))
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

echo "== solve (R=$R spacing=$SPACING ls=${LSTIME}s lns=${LNSTIME}s, taus=$TAUS ks=$KS) =="
rm -f results/sol_t*_k*.txt results/wit_t*_k*.txt
cpp/solver results/buildings.txt --R "$R" --spacing "$SPACING" \
    --lstime "$LSTIME" --lnstime "$LNSTIME" --witness \
    --taus "$TAUS" --ks "$KS" --out results

echo "== verify =="
$PY python/verify.py "$GEOJSON" results --samples 20

echo "== package =="
$PY python/make_submission.py results results/submission.zip

echo "== done =="
grep -c . results/solutions.txt >/dev/null && echo "submission at results/submission.zip"
