#!/usr/bin/env python3
"""Convert a GISCUP GeoJSON building-footprint file into the plain-text
format consumed by the C++ solver, and print dataset statistics.

Output format (text, full double precision):
    n_buildings
    id n_vertices
    x y            (n_vertices lines, exterior ring, closing vertex dropped)
    ...
"""
import argparse
import json
import math
import sys


def ring_area(coords):
    """Signed area of a ring (positive = CCW)."""
    s = 0.0
    n = len(coords)
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return 0.5 * s


def ring_perimeter(coords):
    s = 0.0
    n = len(coords)
    for i in range(n):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % n]
        s += math.hypot(x2 - x1, y2 - y1)
    return s


def clean_ring(coords):
    """Drop the closing vertex and any consecutive duplicate points."""
    if len(coords) >= 2 and coords[0] == coords[-1]:
        coords = coords[:-1]
    out = []
    for p in coords:
        if not out or (p[0] != out[-1][0] or p[1] != out[-1][1]):
            out.append(p)
    # trailing point equal to first
    while len(out) >= 2 and out[0] == out[-1]:
        out.pop()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geojson")
    ap.add_argument("out")
    args = ap.parse_args()

    with open(args.geojson) as f:
        data = json.load(f)

    feats = data["features"]
    buildings = []  # (id, ring)
    n_holes = 0
    n_multi = 0
    for feat in feats:
        geom = feat["geometry"]
        bid = feat["properties"]["id"]
        if geom["type"] == "Polygon":
            rings = [geom["coordinates"]]
        elif geom["type"] == "MultiPolygon":
            rings = geom["coordinates"]
            n_multi += 1
        else:
            print(f"warning: skipping geometry type {geom['type']} (id={bid})",
                  file=sys.stderr)
            continue
        if len(rings) != 1:
            print(f"warning: id={bid} has {len(rings)} polygons; using first",
                  file=sys.stderr)
        poly = rings[0]
        if len(poly) > 1:
            n_holes += 1
            print(f"warning: id={bid} has {len(poly) - 1} hole(s); ignored",
                  file=sys.stderr)
        ring = clean_ring(poly[0])
        if len(ring) < 3:
            print(f"warning: id={bid} degenerate ring; skipped", file=sys.stderr)
            continue
        if ring_area(ring) < 0:  # normalize to CCW
            ring.reverse()
        buildings.append((bid, ring))

    with open(args.out, "w") as f:
        f.write(f"{len(buildings)}\n")
        for bid, ring in buildings:
            f.write(f"{bid} {len(ring)}\n")
            for x, y in ring:
                f.write(f"{x!r} {y!r}\n")

    # stats
    nv = [len(r) for _, r in buildings]
    pers = [ring_perimeter(r) for _, r in buildings]
    areas = [ring_area(r) for _, r in buildings]
    xs = [p[0] for _, r in buildings for p in r]
    ys = [p[1] for _, r in buildings for p in r]
    nv_sorted = sorted(nv)
    print(f"buildings: {len(buildings)}  (multi: {n_multi}, with holes: {n_holes})")
    print(f"vertices: total {sum(nv)}, per-building min {min(nv)} "
          f"median {nv_sorted[len(nv)//2]} max {max(nv)}")
    print(f"perimeter: min {min(pers):.1f} median {sorted(pers)[len(pers)//2]:.1f} "
          f"max {max(pers):.1f} (m)")
    print(f"area: min {min(areas):.1f} median {sorted(areas)[len(areas)//2]:.1f} "
          f"max {max(areas):.1f} (m^2)")
    w, h = max(xs) - min(xs), max(ys) - min(ys)
    print(f"bbox: {w:.0f} x {h:.0f} m "
          f"({min(xs):.0f},{min(ys):.0f})..({max(xs):.0f},{max(ys):.0f})")
    print(f"density: {len(buildings) / (w * h / 1e6):.0f} buildings/km^2")


if __name__ == "__main__":
    main()
