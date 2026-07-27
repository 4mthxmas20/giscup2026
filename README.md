# GISCUP 2026 — Antenna Placement Solver

Solution for the [ACM SIGSPATIAL GIS Cup 2026](https://sigspatial2026.sigspatial.org/giscup.html):
place `k` antennas on building perimeters to maximize the number of buildings whose
boundary is at least a fraction `τ` visible from the antenna set.

## Layout

- `data/` — input GeoJSON datasets
- `cpp/` — C++17 solver core (visibility sweep + greedy optimizer), no external deps except OpenMP
- `python/` — data prep, independent verification (shapely), visualization, submission packaging
- `results/` — solver outputs per (τ, k) sub-problem
- `scripts/` — end-to-end pipeline drivers

## Build

```
cd cpp && make        # requires clang/gcc with C++17; OpenMP optional but recommended
```

## Run (end to end)

```
./scripts/run_all.sh data/GIS-cup-sample-dataset.geojson
```

This prepares inputs, solves all nine (τ, k) sub-problems, verifies claimed
buildings with an independent Python/shapely checker, and writes the submission
file under `results/`.
