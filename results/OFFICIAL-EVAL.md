# Official evaluator verdict — sample dataset

Run with the organizers' evaluator (`alowe/gis-cup-2026-evaluator`, ArcGIS
geometry, 0.001 m spatial tolerance) against `results/solutions.txt`, one
vitest run per sub-problem, in the codespace that produced the solution.

Dataset: `data/GIS-cup-sample-dataset.geojson` (12,860 buildings, 78,727 edges;
the corrected release with building 9448's inner ring removed).

| τ | k | claimed | verified | lost | eval seconds |
| --- | --- | --- | --- | --- | --- |
| 0.25 | 50 | 2388 | 2388 | **0** | 56.2 |
| 0.25 | 500 | 10624 | 10624 | **0** | 428.9 |
| 0.25 | 1000 | 12860 | 12860 | **0** | 772.3 |
| 0.5 | 50 | 755 | 755 | **0** | 33.2 |
| 0.5 | 500 | 5939 | 5939 | **0** | 429.6 |
| 0.5 | 1000 | 10165 | 10165 | **0** | 798.3 |
| 0.75 | 50 | 329 | 329 | **0** | 27.8 |
| 0.75 | 500 | 3086 | 3086 | **0** | 331.8 |
| 0.75 | 1000 | 5895 | 5895 | **0** | 724.4 |

Every claimed building survives. No antenna was rejected as off-boundary and
none collided after snapping (`antennasUnique == k` for all nine). The only
warnings were `ANTENNA_SNAPPED` on a handful of entries, where a candidate on
a segment interior differs from the evaluator's recomputed nearest point by
~1e-10 m — well inside the 1 mm tolerance and harmless.

τ=0.25 / k=1000 reaches every building in the dataset (12860/12860).

## How this was produced

```bash
scripts/run_official_eval.sh data/GIS-cup-sample-dataset.geojson results/solutions.txt
```

Three independent implementations agree on these numbers: the C++ solver's own
bookkeeping (checked by recomputing the score from the final antenna set),
`python/verify.py` (shapely, 1.04 M sampled sight lines), and the evaluator
above.
