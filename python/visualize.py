#!/usr/bin/env python3
"""Render a solution: buildings colored by serviced status + antenna positions.

Usage: visualize.py geojson sol_file out.png
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


def main():
    gj_path, sol_path, out_path = sys.argv[1:4]
    with open(gj_path) as f:
        gj = json.load(f)
    with open(sol_path) as f:
        tau, k = f.readline().split()
        coords = [tuple(map(float, c.split())) for c in f.readline().split(";")]
        claimed = set(int(x) for x in f.readline().split(","))

    served_polys, unserved_polys = [], []
    for feat in gj["features"]:
        if feat["geometry"]["type"] != "Polygon":
            continue
        ring = feat["geometry"]["coordinates"][0]
        (served_polys if feat["properties"]["id"] in claimed
         else unserved_polys).append(ring)

    fig, ax = plt.subplots(figsize=(14, 12))
    ax.add_collection(PolyCollection(unserved_polys, facecolors="#d0d0d0",
                                     edgecolors="#999999", linewidths=0.2))
    ax.add_collection(PolyCollection(served_polys, facecolors="#7fc97f",
                                     edgecolors="#4a874a", linewidths=0.2))
    xs = [x for x, _ in coords]
    ys = [y for _, y in coords]
    ax.scatter(xs, ys, s=14, c="red", marker="^", zorder=5, label="antenna")
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set_title(f"tau={tau} k={k}: {len(claimed)} serviced")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
