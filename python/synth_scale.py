#!/usr/bin/env python3
"""Tile a building dataset into a larger synthetic one, to measure how the
solver scales before competition day (the evaluation dataset's size is
unknown).

Tiles are laid out on an n x n grid with a gap, so footprints stay disjoint
and each copy keeps the original's local density.

Usage: synth_scale.py in.geojson out.geojson --tiles 2 [--gap 200]
"""
import argparse
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--tiles", type=int, default=2,
                    help="grid side; total size multiplies by tiles^2")
    ap.add_argument("--gap", type=float, default=200.0,
                    help="metres of empty space between tiles")
    args = ap.parse_args()

    with open(args.infile) as f:
        data = json.load(f)
    feats = data["features"]

    xs = [p[0] for ft in feats for p in ft["geometry"]["coordinates"][0]]
    ys = [p[1] for ft in feats for p in ft["geometry"]["coordinates"][0]]
    width = max(xs) - min(xs) + args.gap
    height = max(ys) - min(ys) + args.gap

    out = []
    for row in range(args.tiles):
        for col in range(args.tiles):
            dx, dy = col * width, row * height
            for ft in feats:
                ring = [[p[0] + dx, p[1] + dy]
                        for p in ft["geometry"]["coordinates"][0]]
                out.append({
                    "type": "Feature",
                    "properties": {"id": f"{ft['properties']['id']}-r{row}c{col}"},
                    "geometry": {"type": "Polygon", "coordinates": [ring]},
                })

    data["features"] = out
    with open(args.outfile, "w") as f:
        json.dump(data, f)
    print(f"wrote {args.outfile}: {len(out)} buildings "
          f"({args.tiles}x{args.tiles} tiles of {len(feats)}), "
          f"extent {width * args.tiles:.0f} x {height * args.tiles:.0f} m")


if __name__ == "__main__":
    main()
