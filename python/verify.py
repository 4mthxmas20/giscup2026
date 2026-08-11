#!/usr/bin/env python3
"""Independent verification of solver solutions against the original GeoJSON.

For each claimed serviced building the verifier checks, with shapely and
without reusing any C++ geometry code:
  1. every antenna lies on some building boundary (distance < 1e-6);
  2. the witness intervals, unioned and clipped to [0, perimeter], have total
     length >= tau * perimeter (computed from the GeoJSON, not the solver);
  3. sampled points inside witness intervals are actually visible from their
     witness antenna (segment penetrates no building interior beyond 1e-6 m).

Usage: verify.py geojson soldir [--samples N] [--configs t0.25_k50 ...]
"""
import argparse
import glob
import json
import math
import os
import random
import sys

from shapely.geometry import LineString, Point, Polygon
from shapely.strtree import STRtree


def load_buildings(path):
    with open(path) as f:
        data = json.load(f)
    blds = {}
    for feat in data["features"]:
        geom = feat["geometry"]
        if geom["type"] != "Polygon":
            continue
        bid = str(feat["properties"]["id"]).strip()
        blds[bid] = Polygon(geom["coordinates"][0])
    return blds


def arc_point(ring, cum, t):
    """Point at arc length t along ring (list of coords, closing point dropped)."""
    n = len(ring)
    t = t % cum[-1]
    # binary search
    lo, hi = 0, n
    while lo < hi:
        mid = (lo + hi) // 2
        if cum[mid + 1] < t:
            lo = mid + 1
        else:
            hi = mid
    j = lo
    seg = cum[j + 1] - cum[j]
    u = (t - cum[j]) / seg if seg > 0 else 0
    a, b = ring[j], ring[(j + 1) % n]
    return (a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u)


def penetration_blocked(q, p, tree, shrunk):
    """Blocked iff segment qp penetrates a building interior by > 1e-6 m.

    Uses polygons shrunk by 1e-7 m so that sight lines grazing along a wall
    (collinear within float rounding) are not misclassified as blocked."""
    seg = LineString([q, p])
    if seg.length == 0:
        return False
    for idx in tree.query(seg):
        s = shrunk[idx]
        if s.is_empty:
            continue
        if seg.intersection(s).length > 1e-6:
            return True
    return False


def union_len(ivs, per):
    """Union length of intervals clipped to [0, per] (t values may exceed per
    slightly due to float rounding at the wrap point)."""
    cl = []
    for t0, t1 in ivs:
        t0, t1 = max(0.0, t0), min(per, t1)
        if t1 > t0:
            cl.append((t0, t1))
    cl.sort()
    total, ce = 0.0, -1.0
    for t0, t1 in cl:
        if t0 > ce:
            total += t1 - t0
            ce = t1
        elif t1 > ce:
            total += t1 - ce
            ce = t1
    return total


def verify_config(gj_blds, tree, polys, shrunk, soldir, tag, nsamples, rng, ids):
    """Solver files identify buildings by 0-based index; `ids` maps those back
    to the dataset's own IDs."""
    sol = os.path.join(soldir, f"sol_{tag}.txt")
    wit = os.path.join(soldir, f"wit_{tag}.txt")
    with open(sol) as f:
        tau, k = f.readline().split()
        tau, k = float(tau), int(k)
        coords = [tuple(map(float, c.split())) for c in f.readline().split(";")]
        claimed = [int(x) for x in f.readline().split(",")]
    assert len(coords) == k, f"{tag}: {len(coords)} coords != k={k}"
    assert len(set(coords)) == k, f"{tag}: duplicate antenna coordinates"

    # 1. antennas on boundaries
    bad_pos = 0
    for q in coords:
        pt = Point(q)
        ok = any(polys[i].exterior.distance(pt) < 1e-6
                 for i in tree.query(pt.buffer(1e-5)))
        if not ok:
            bad_pos += 1
    # 2 & 3. witness check
    witness = {}
    with open(wit) as f:
        for line in f:
            parts = line.split()
            bid = int(parts[0])
            ivs = []
            for w in parts[3:]:
                ai, t0, t1 = w.split(":")
                ivs.append((int(ai), float(t0), float(t1)))
            witness[bid] = ivs

    n_cov_fail = 0
    n_vis_fail = 0
    n_checked_pts = 0
    valid = []
    for bid in claimed:
        poly = gj_blds[ids[bid]]
        ring = list(poly.exterior.coords)[:-1]
        # match the solver's CCW normalization so arc params line up
        area = sum(ring[j][0] * ring[(j + 1) % len(ring)][1] -
                   ring[(j + 1) % len(ring)][0] * ring[j][1]
                   for j in range(len(ring)))
        if area < 0:
            ring = ring[::-1]
        cum = [0.0]
        for j in range(len(ring)):
            cum.append(cum[-1] + math.dist(ring[j], ring[(j + 1) % len(ring)]))
        per = cum[-1]
        ivs = witness.get(bid, [])
        cov = union_len([(t0, t1) for _, t0, t1 in ivs], per)
        if cov < tau * per:
            n_cov_fail += 1
            continue
        ok = True
        for _ in range(nsamples):
            ai, t0, t1 = ivs[rng.randrange(len(ivs))]
            t = rng.uniform(t0, min(t1, per))
            p = arc_point(ring, cum, t)
            q = coords[ai]
            if math.dist(q, p) < 1e-12:
                continue
            n_checked_pts += 1
            if penetration_blocked(q, p, tree, shrunk):
                n_vis_fail += 1
                ok = False
                break
        if ok:
            valid.append(bid)

    status = "OK " if (bad_pos == 0 and n_cov_fail == 0 and n_vis_fail == 0) \
        else "FAIL"
    print(f"{status} {tag}: claimed {len(claimed)}, valid {len(valid)}, "
          f"cov_fail {n_cov_fail}, vis_fail {n_vis_fail}, "
          f"bad_antenna_pos {bad_pos}, sampled {n_checked_pts} pts")
    return bad_pos == 0 and n_cov_fail == 0 and n_vis_fail == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geojson")
    ap.add_argument("soldir")
    ap.add_argument("--samples", type=int, default=20,
                    help="visibility samples per claimed building")
    ap.add_argument("--configs", nargs="*", default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ids", default=None,
                    help="index -> ID map from prepare.py "
                         "(default: <soldir>/buildings.txt.ids.json)")
    args = ap.parse_args()

    ids_path = args.ids or os.path.join(args.soldir, "buildings.txt.ids.json")
    with open(ids_path) as f:
        ids = json.load(f)

    gj_blds = load_buildings(args.geojson)
    poly_ids = sorted(gj_blds)
    polys = [gj_blds[i] for i in poly_ids]
    tree = STRtree(polys)
    shrunk = [p.buffer(-1e-7) for p in polys]
    rng = random.Random(args.seed)

    tags = args.configs
    if not tags:
        tags = sorted(
            os.path.basename(p)[4:-4]
            for p in glob.glob(os.path.join(args.soldir, "sol_*.txt")))
    all_ok = True
    for tag in tags:
        all_ok &= verify_config(gj_blds, tree, polys, shrunk, args.soldir,
                                tag, args.samples, rng, ids)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
