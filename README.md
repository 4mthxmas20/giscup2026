# GISCUP 2026 — Antenna Placement Solver

Solution for the [ACM SIGSPATIAL GIS Cup 2026](https://sigspatial2026.sigspatial.org/giscup.html):
place `k` antennas on building perimeters to maximize the number of buildings whose
boundary is at least a fraction `τ` visible from the antenna set.

## Method

1. **Candidate generation** — every polygon vertex plus extra samples on long
   edges (spacing configurable), each offset 1e-6 m outward for robust
   visibility queries; candidates that would fall inside a touching building
   are discarded.
2. **Exact visibility** — per candidate, a rotational sweep over all boundary
   edges within radius `R` yields the exactly visible sub-segments of every
   building boundary as arc-length intervals. Conservative by construction:
   obstacles are every edge intersecting the disk, and credit is clipped to
   the disk, so credited intervals are always truly visible.
3. **Optimization** — per (τ, k): lazy greedy over a portfolio of surrogate
   objectives (newly-serviced bonus + γ-weighted capped coverage progress for
   a grid of γ values), then swap local search (evict the antenna with the
   smallest exact service loss, re-insert the globally best replacement, 15%
   random evictions to escape plateaus).
4. **Independent verification** — a Python/shapely checker recomputes antenna
   positions, coverage-interval unions, and sampled sight lines from the
   original GeoJSON before anything is packaged.

## Requirements

- C++17 compiler (clang or gcc; tested with Apple clang on macOS)
- Python 3.10+ with `shapely`, `numpy`, `matplotlib` (only for verification,
  visualization, and packaging — the solver itself has zero dependencies)

```
python3 -m venv .venv && .venv/bin/pip install shapely numpy matplotlib
```

## Run (end to end)

```
./scripts/run_all.sh <input.geojson> [lstime_per_config_s] [R_m] [spacing_m]
```

e.g. `./scripts/run_all.sh data/GIS-cup-sample-dataset.geojson 120 1200 6`.
This builds the solver, converts the GeoJSON, solves all nine (τ, k)
sub-problems, verifies every claim independently, and writes
`results/solutions.txt` + `results/submission.zip`.

## Manual steps

```
make -C cpp                                  # build
python3 python/prepare.py in.geojson results/buildings.txt
cpp/solver results/buildings.txt --R 1200 --spacing 6 \
    --lstime 120 --witness --out results     # solve 9 sub-problems
python3 python/verify.py in.geojson results  # independent check
python3 python/make_submission.py results results/submission.zip
python3 python/visualize.py in.geojson results/sol_t0.5_k500.txt out.png
```

Useful solver flags: `--probe N` (time a candidate sample and extrapolate the
visibility pass before committing to R/spacing), `--viscache f.bin` (reuse
visibility across runs; the cache is validated against the prepared dataset
and geometry parameters), `--lstime` / `--lnstime` (per-config search budgets),
`--maxruin`, `--taus`, `--ks`, `--gammas`, `--threads`, `--minivlen`.

## Checking against the official evaluator

`.github/workflows/official-eval.yml` runs the organizers' own evaluator
(ArcGIS geometry, 1 mm spatial tolerance) on GitHub runners, one job per
sub-problem, and fails if any claimed building is rejected:

```bash
gh workflow run official-eval.yml
```

See [RUNBOOK.md](RUNBOOK.md) for the competition-day procedure.

## Layout

- `data/` — input GeoJSON datasets
- `cpp/` — C++17 solver core (geometry, grid index, sweep, optimizer)
- `python/` — prepare / verify / visualize / package + sweep cross-check test
  (`test_visibility.py`, run with a seed argument)
- `scripts/run_all.sh` — one-command pipeline
- `results/` — solutions, witness certificates, submission zip
