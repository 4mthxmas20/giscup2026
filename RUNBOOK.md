# Competition-day runbook

**Data drops 2026-08-15 16:00 UTC (09:00 PDT). Submission closes 2026-08-16
16:00 UTC.** Submit through EasyChair (`giscup2026`) as a single archive
containing the results file and a source directory with build instructions.

Sources: `https://github.com/alowe/gis-cup-2026-evaluator` →
`datasets/GIS-cup-competition-dataset.geojson` and
`datasets/competition-parameters.txt`.

## Where to run

| | Role |
| --- | --- |
| **GitHub Actions** | the compute: parallel search and batch verification |
| **Codespace** | the control room: observation, tuning, debugging, manual fallback |
| **Local machine** | never a compute environment — viewing and submitting only |

### The main path: one dispatch

```bash
gh workflow run competition.yml \
  -f dataset=data/competition.geojson \
  -f taus=<from parameter file> -f ks=<from parameter file> \
  -f R=<from probe> -f spacing=<from probe> \
  -f budget_easy=900 -f budget_hard=7200 -f seeds_hard=3
gh run watch
```

`prepare` builds the visibility table once into the Actions cache; `search`
fans out one runner per `(tau, k, seed)`, each verifying its own configuration;
`aggregate` picks the best seed per configuration and packages; `official-eval`
re-checks the packaged `solutions.txt` with the organizers' evaluator. Download
`submission` when it is green.

Budgets are per configuration and run **concurrently**, so wall time is roughly
one budget, not nine. Give the saturating configurations a short pass and the
hard ones the hours (see *Spend the budget where it can move the ranking*).

Two independent guards keep the shared inputs honest: a sha256 `manifest.txt`
checked by every downstream job (fails loud), and the solver's own cache
fingerprint over the dataset, R, spacing, eps and min-interval (fails safe by
rebuilding). The cache key additionally hashes the geometry sources, so
changing the sweep cannot silently reuse a table built by the old logic.

### The fallback: Codespace, serial

Use when Actions is unavailable or a run needs hands-on debugging.

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

`--probe` times a candidate sample and extrapolates.

Measured on a 4-core codespace by tiling the sample (`python/synth_scale.py`),
all at `--spacing 6`:

| Buildings | Candidates | R=600 | R=800 | R=1200 | Intervals @R=1200 |
| --- | --- | --- | --- | --- | --- |
| 12,860 (sample) | 183 K | — | — | 13 min | 72 MB |
| 51,440 (2×2) | 734 K | 11.1 min | 21.6 min | 54.4 min | 309 MB |
| 115,740 (3×3) | 1,651 K | 25.7 min | 47.6 min | 125.2 min | 711 MB |

Visibility cost is **linear in building count** (2.25× the buildings cost
2.20–2.31× the time) and grows about R²–R^2.4. Memory is a non-issue.

**Raise `--spacing` before lowering `--R`.** Lower R only shortens the
visibility pass. Candidate count scales exactly with building count, and the
LNS rebuild rescans every candidate on each iteration — at 1.65 M candidates
the search runs ~9× fewer iterations per second than on the sample, and search
iterations are what actually buy score. Spacing is the only dial that cuts
both: it shrinks the visibility pass *and* the per-iteration scan.

So on a large dataset, budget backwards from search, not from visibility: a
125-minute visibility pass is affordable inside 24 hours, but a search that
crawls is not recoverable.

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
time skip the expensive pass. The solver now validates the cache against the
prepared dataset plus R/spacing/epsilon/min-interval settings and rebuilds if
anything changed.

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

Then get the organizers' own verdict. Run it in the same codespace, so the
evaluator sees exactly the geometry that produced the solution:

```bash
scripts/run_official_eval.sh data/competition.geojson results/solutions.txt
```

This clones the evaluator, installs its pinned dependencies (~400 MB, one
time), and evaluates every sub-problem. Each must report `lost=0`; the report
lands in `results/official-eval-report.json`.

The `.github/workflows/official-eval.yml` workflow does the same thing on
GitHub runners with the nine sub-problems in parallel. It is the faster option
when the codespace is busy, but it needs the results committed and pushed
first.

## 4. Package and submit

```bash
python3 python/make_submission.py results results/submission.zip
unzip -l results/submission.zip
```

The archive holds `solutions.txt` plus `source/` (C++ core, Python tooling,
README with build and run instructions). Upload to EasyChair.

## Pre-flight checklist

- [ ] τ/k taken from `competition-parameters.txt`, not assumed
- [ ] `prepare` saved the visibility cache (no `cache-hit` warning downstream)
- [ ] every `search` job logged `loaded vis cache` — a job that recomputes it
      is burning minutes for nothing
- [ ] `verify.py` printed `OK` inside each search job
- [ ] `aggregate` produced both `solutions.txt` and `subproblems.json`
- [ ] Official evaluator: `lost=0` on every index in `subproblems.json`
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
| `results/` suddenly holds only inputs | `git stash -u` swept the untracked outputs away | `git stash pop`; never run `stash -u`, `clean`, or `checkout .` in the workspace holding a finished run |
| Evaluator dies instantly with `ENOENT` | relative dataset/submission path, resolved against the evaluator checkout | fixed in `run_official_eval.sh`; pass absolute paths if calling vitest by hand |
| `git pull` in the codespace asks for a username | the SSH-exec context has no git credentials | push from the local clone and use `gh codespace cp` to copy individual files in |
