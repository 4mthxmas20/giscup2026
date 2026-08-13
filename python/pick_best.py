#!/usr/bin/env python3
"""Collapse per-seed solver outputs into one best solution per (tau, k).

The competition matrix may run a hard configuration under several seeds; each
job self-verifies, so the claimed count is trustworthy and the largest wins.

Usage: pick_best.py <indir> <outdir>

`indir` is searched recursively for `sol_t<tau>_k<k>.txt` (typically one per
downloaded artifact directory). The winner of each (tau, k) is copied to
`outdir`, along with its witness file when present.
"""
import argparse
import os
import re
import shutil
import sys


def claimed_count(path):
    with open(path) as f:
        f.readline()  # (tau, k)
        f.readline()  # coordinates
        ids = f.readline().strip()
    return len([x for x in ids.split(",") if x.strip()]) if ids else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("outdir")
    args = ap.parse_args()

    best = {}  # (tau, k) -> (count, path)
    pattern = re.compile(r"sol_t([0-9.]+)_k([0-9]+)\.txt$")
    for root, _dirs, files in os.walk(args.indir):
        for name in files:
            m = pattern.match(name)
            if not m:
                continue
            path = os.path.join(root, name)
            key = (m.group(1), m.group(2))
            n = claimed_count(path)
            if key not in best or n > best[key][0]:
                best[key] = (n, path)

    if not best:
        sys.exit(f"no sol_t*_k*.txt found under {args.indir}")

    os.makedirs(args.outdir, exist_ok=True)
    for (tau, k), (n, path) in sorted(best.items(), key=lambda x: (float(x[0][0]),
                                                                  int(x[0][1]))):
        dest = os.path.join(args.outdir, f"sol_t{tau}_k{k}.txt")
        shutil.copy2(path, dest)
        wit = os.path.join(os.path.dirname(path), f"wit_t{tau}_k{k}.txt")
        if os.path.exists(wit):
            shutil.copy2(wit, os.path.join(args.outdir, f"wit_t{tau}_k{k}.txt"))
        print(f"(tau={tau}, k={k}): {n} claimed  <- {os.path.relpath(path, args.indir)}")


if __name__ == "__main__":
    main()
