#!/usr/bin/env python3
"""Assemble the final submission: solutions.txt (9 sub-problem blocks in the
required 3-line format) + a zip with the source code and instructions.

Usage: make_submission.py soldir out.zip [--configs ...]
"""
import argparse
import os
import zipfile

TAUS = ["0.25", "0.5", "0.75"]
KS = ["50", "500", "1000"]

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("soldir")
    ap.add_argument("out")
    args = ap.parse_args()

    blocks = []
    for tau in TAUS:
        for k in KS:
            path = os.path.join(args.soldir, f"sol_t{tau}_k{k}.txt")
            with open(path) as f:
                _tau, _k = f.readline().split()
                coords = [tuple(map(float, c.split()))
                          for c in f.readline().split(";")]
                ids = f.readline().strip()
            assert len(coords) == int(k)
            line1 = f"({tau}, {k})"
            line2 = ", ".join(f"({x!r}, {y!r})" for x, y in coords)
            blocks.append(f"{line1}\n{line2}\n{ids}\n")

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
