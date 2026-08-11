# Competition-day runbook

**Data drops 2026-08-15 16:00 UTC (09:00 PDT). Submission closes 2026-08-16
16:00 UTC.** Submit through EasyChair (`giscup2026`) as a single archive
containing the results file and a source directory with build instructions.

Sources: `https://github.com/alowe/gis-cup-2026-evaluator` →
`datasets/GIS-cup-competition-dataset.geojson` and
`datasets/competition-parameters.txt`.

## Where to run

All compute runs on a Codespace, not the local machine.

```bash
gh codespace create -R 4mthxmas20/giscup2026 -b main \
  -m standardLinux32gb --idle-timeout 240m --retention-period 72h
gh codespace ssh -c <name> -- "cd /workspaces/giscup2026 && <command>"
```

`standardLinux32gb` (4 cores, 16 GB) is the largest machine this account can
use. **`.devcontainer/devcontainer.json` must not request more than that** —
`hostRequirements` above the ceiling filters every machine out and creation
fails with "no available machine types", which looks like a quota problem but
is not.

Four cores means roughly double the wall time of an 8-core laptop; scale the
budget estimates below accordingly. Start long runs detached so the SSH session
is not load-bearing:

```bash
gh codespace ssh -c <name> -- "cd /workspaces/giscup2026 && nohup bash -c '<pipeline>' > run.log 2>&1 & echo LAUNCHED"
```

Codespaces bill by the core-hour and the machine keeps running until the idle
timeout, so stop it when the work is done: `gh codespace stop -c <name>`.

## 0. Get the data (T+0)

```bash
cd ~/Desktop/GIScup
curl -L -o data/competition.geojson \
  https://raw.githubusercontent.com/alowe/gis-cup-2026-evaluator/main/datasets/GIS-cup-competition-dataset.geojson
curl -L -o data/competition-parameters.txt \
  https://raw.githubusercontent.com/alowe/gis-cup-2026-evaluator/main/datasets/competition-parameters.txt
cat data/competition-parameters.txt
```

**Read the parameter file before anything else.** τ and k need not match the
sample's 3×3 grid; everything downstream takes them as arguments.

## 1. Size it up (T+5 min)

```bash
python3 python/prepare.py data/competition.geojson results/buildings.txt
make -C cpp
cpp/solver results/buildings.txt --probe 2000 --R 1200 --spacing 6
```

`--probe` times a candidate sample and extrapolates. Pick R and spacing so the
estimated visibility pass lands **under ~1 hour**:

| Dial | Effect |
| --- | --- |
| `--R` | cost grows roughly with R²; R=800 scores ~0.6% below R=1200 for a quarter of the time |
| `--spacing` | candidate count scales inversely; 6→12 roughly halves both cost and candidate density |

Re-probe after changing dials. If the dataset is far larger than the 12,860-building
sample, start at `--R 600 --spacing 10` and only tighten if the estimate is small.

## 2. Solve

```bash
TAUS=<from parameter file> KS=<from parameter file> \
  scripts/run_all.sh data/competition.geojson 600 <R> <spacing>
```

`run_all.sh` prepares, solves, verifies, and packages. The second argument is
the per-config **total** search budget in seconds, split automatically between
1-opt and LNS (1-opt capped at a third, max 180 s — it saturates early while
LNS keeps paying). Nine configs × that budget is the bulk of the wall time —
budget backwards from the deadline and leave **≥3 hours of slack**.

Measured on the sample at τ=0.75, k=500 with a 240 s budget: 1-opt alone
reached 2991, adding LNS reached 3050.

To split the budget explicitly, call the solver directly:

```bash
cpp/solver results/buildings.txt --R <R> --spacing <spacing> \
  --taus <taus> --ks <ks> --lstime 300 --lnstime 600 \
  --viscache cache/comp.bin --witness --out results
```

`--viscache` writes the visibility table once so later re-runs with more search
time skip the expensive pass. **Never reuse a cache across datasets.**

### Spend the budget where it can move the ranking

Scoring is relative — each config contributes `mine / best_submitted`, capped at
1. A config everyone nearly saturates cannot differentiate anyone, so extra
search there buys close to nothing.

On the sample, τ=0.25/k=1000 reaches 12841 of 12860 buildings (99.9%): every
serious entry will score ~1.0 on it. τ=0.75/k=50 reaches 302 of 12860 (2.3%) —
that is where solutions actually spread out.

So do **not** split the budget evenly. Give the low-τ / high-k configs a short
pass and pour the remaining hours into high-τ / low-k, which stay far from
saturation and keep improving under LNS:

```bash
# cheap pass on the near-saturated configs
cpp/solver results/buildings.txt --viscache cache/comp.bin \
  --taus 0.25 --ks 50,500,1000 --lstime 60 --lnstime 120 --witness --out results

# everything else goes to the hard end
cpp/solver results/buildings.txt --viscache cache/comp.bin \
  --taus 0.75 --ks 50 --lstime 180 --lnstime 5400 --witness --out results
```

Check saturation before deciding: if a config already claims >95% of buildings,
it is done.

### Strategy: bank a full result, then intensify

Do a complete pass with a modest budget first so a submittable answer exists
early, then spend the remaining hours improving whichever configs look weakest
(highest τ, lowest k are usually furthest from saturation). `--warmstart`
resumes from a previous solution instead of restarting from greedy:

```bash
cpp/solver results/buildings.txt --viscache cache/comp.bin \
  --R <R> --spacing <spacing> --taus 0.75 --ks 500 \
  --warmstart results/sol_t0.75_k500.txt --lstime 0 --lnstime 3600 \
  --witness --out results
```

The warm-start file must come from a run with the **same dataset, R and
spacing** — the solver checks each coordinate against the candidate set and
refuses to start if they do not match. Re-run `verify.py` afterwards.

## 3. Verify (mandatory before submitting)

```bash
python3 python/verify.py data/competition.geojson results --samples 20
```

Every config must print `OK`. Anything else means claims are being
over-reported — fix before packaging, never submit a `FAIL`.

Run this **on the machine that produced the solution**. The sweep is not
bit-reproducible across platforms — the same dataset at R=1200 yielded
6,011,070 intervals on macOS/clang and 6,010,413 in the Linux codespace,
because `atan2` differs by an ulp at the margins. The gap is 0.01% and both
sides stay conservative, but verifying elsewhere compares against subtly
different geometry for no reason.

Then get the organizers' own verdict (runs on GitHub, not this machine):

```bash
git add -A && git commit -m "competition run" && git push
gh workflow run official-eval.yml \
  -f dataset=data/competition.geojson -f submission=results/solutions.txt
gh run watch
```

Each of the nine jobs must report `lost=0`. Download reports with
`gh run download <run-id>`.

## 4. Package and submit

```bash
python3 python/make_submission.py results results/submission.zip
unzip -l results/submission.zip
```

The archive holds `solutions.txt` plus `source/` (C++ core, Python tooling,
README with build and run instructions). Upload to EasyChair.

## Pre-flight checklist

- [ ] τ/k taken from `competition-parameters.txt`, not assumed
- [ ] `verify.py` prints `OK` for every config
- [ ] Official evaluator workflow: `lost=0` on all nine
- [ ] Antenna count equals k exactly, with no duplicate coordinates
      (both asserted by `make_submission.py`)
- [ ] Archive contains results **and** source with instructions
- [ ] Submitted before 16:00 UTC on 2026-08-16

## Known failure modes

| Symptom | Cause | Fix |
| --- | --- | --- |
| `prepare.py` warns about holes | multi-ring polygon | holes are ignored; the official loader rejects such datasets outright, so this should not occur |
| Visibility pass runs hours | R or spacing too aggressive for the dataset | re-probe, lower R first |
| `verify.py` reports `cov_fail` | solver/verifier threshold drift | do not submit; investigate the witness file for that config |
| Evaluator reports `ANTENNA_OFF_BOUNDARY` | antenna >1 mm from any boundary | should be impossible; candidates are exact boundary points |
