#!/usr/bin/env python3
"""Assemble the final submission: solutions.txt (one 3-line block per
sub-problem, in the required format) + an archive with the source code.

The (tau, k) list is discovered from the solver outputs in `soldir` unless
given explicitly -- the competition publishes its own parameter file and its
values need not match the sample's.

Usage: make_submission.py soldir out.zip [--configs 0.25:50 0.25:500 ...]
"""
import argparse
import glob
import json
import os
import re
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def discover_configs(soldir):
    """(tau, k) pairs from sol_t<tau>_k<k>.txt, ordered by tau then k."""
    found = []
    for path in glob.glob(os.path.join(soldir, "sol_t*_k*.txt")):
        m = re.fullmatch(r"sol_t([0-9.]+)_k([0-9]+)\.txt", os.path.basename(path))
        if m:
            found.append((m.group(1), m.group(2)))
    return sorted(found, key=lambda c: (float(c[0]), int(c[1])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("soldir")
    ap.add_argument("out")
    ap.add_argument("--ids", default=None,
                    help="index -> ID map from prepare.py "
                         "(default: <soldir>/buildings.txt.ids.json)")
    ap.add_argument("--configs", nargs="*", default=None,
                    help="ordered tau:k pairs; default discovers them in soldir")
    args = ap.parse_args()

    ids_path = args.ids or os.path.join(args.soldir, "buildings.txt.ids.json")
    with open(ids_path) as f:
        id_map = json.load(f)

    if args.configs:
        configs = [tuple(c.split(":")) for c in args.configs]
    else:
        configs = discover_configs(args.soldir)
    if not configs:
        raise SystemExit(f"no sol_t*_k*.txt found in {args.soldir}")
    print(f"configs: {', '.join(f'({t}, {k})' for t, k in configs)}")

    # Solver output stores 0-based indices, so pairing it with a newer ID map
    # would silently emit the wrong building IDs. Refuse stale files outright.
    ids_mtime = os.path.getmtime(ids_path)
    stale = [f"sol_t{t}_k{k}.txt" for t, k in configs
             if os.path.getmtime(os.path.join(args.soldir, f"sol_t{t}_k{k}.txt"))
             < ids_mtime]
    if stale:
        raise SystemExit(
            f"stale solver output (older than {os.path.basename(ids_path)}): "
            f"{', '.join(stale)}\nre-run the solver before packaging")

    blocks = []
    for tau, k in configs:
        path = os.path.join(args.soldir, f"sol_t{tau}_k{k}.txt")
        with open(path) as f:
            _tau, _k = f.readline().split()
            coords = [tuple(map(float, c.split()))
                      for c in f.readline().split(";")]
            claimed = [int(x) for x in f.readline().split(",") if x.strip()]
        assert len(coords) == int(k), f"{path}: {len(coords)} != k={k}"
        assert len(set(coords)) == len(coords), f"{path}: duplicate antennas"
        line1 = f"({tau}, {k})"
        line2 = ", ".join(f"({x!r}, {y!r})" for x, y in coords)
        line3 = ", ".join(id_map[i] for i in claimed)
        blocks.append(f"{line1}\n{line2}\n{line3}\n")

    sol_txt = os.path.join(args.soldir, "solutions.txt")
    with open(sol_txt, "w") as f:
        f.write("".join(blocks))
    print(f"wrote {sol_txt}")

    with zipfile.ZipFile(args.out, "w", zipfile.ZIP_DEFLATED) as z:
        z.write(sol_txt, "solutions.txt")
        for rel in ["README.md", "cpp/Makefile", "cpp/solver.cpp",
                    "cpp/geometry.hpp", "cpp/grid.hpp", "cpp/visibility.hpp",
                    "python/prepare.py", "python/verify.py",
                    "python/make_submission.py", "scripts/run_all.sh"]:
            z.write(os.path.join(ROOT, rel), f"source/{rel}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
