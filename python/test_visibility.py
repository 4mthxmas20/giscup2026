#!/usr/bin/env python3
"""Cross-check the C++ visibility sweep against a shapely brute force.

Generates random scenes of disjoint boxes, queries visibility from random
outside points via `solver --testvis`, then verifies:
  1. (hard) every point the sweep claims visible IS visible and within R;
  2. (soft) points that are visible & within R are claimed, except within
     a tiny tolerance of interval endpoints.
"""
import math
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from shapely.geometry import LineString, Point, Polygon

ROOT = Path(__file__).resolve().parent.parent
SOLVER = ROOT / "cpp" / "solver"


def rand_scene(rng, nbox, extent=100.0):
    polys = []
    tries = 0
    while len(polys) < nbox and tries < 1000:
        tries += 1
        w, h = rng.uniform(3, 18), rng.uniform(3, 18)
        cx, cy = rng.uniform(0, extent), rng.uniform(0, extent)
        ang = rng.uniform(0, math.pi)
        ca, sa = math.cos(ang), math.sin(ang)
        pts = []
        for dx, dy in [(-w/2, -h/2), (w/2, -h/2), (w/2, h/2), (-w/2, h/2)]:
            pts.append((cx + dx*ca - dy*sa, cy + dx*sa + dy*ca))
        p = Polygon(pts)
        if all(p.distance(o) > 0.5 for o in polys):
            polys.append(p)
    return polys


def write_scene(polys, path):
    rings = []
    with open(path, "w") as f:
        f.write(f"{len(polys)}\n")
        for i, p in enumerate(polys):
            ring = list(p.exterior.coords)[:-1]
            # ensure CCW
            area = sum(ring[j][0]*ring[(j+1) % len(ring)][1] -
                       ring[(j+1) % len(ring)][0]*ring[j][1]
                       for j in range(len(ring)))
            if area < 0:
                ring.reverse()
            f.write(f"{i+1} {len(ring)}\n")
            for x, y in ring:
                f.write(f"{x!r} {y!r}\n")
            rings.append(ring)
    return rings


def arc_point(ring, t):
    """Point at arc length t along ring."""
    for j in range(len(ring)):
        a, b = ring[j], ring[(j+1) % len(ring)]
        L = math.dist(a, b)
        if t <= L:
            u = t / L if L > 0 else 0
            return (a[0] + (b[0]-a[0])*u, a[1] + (b[1]-a[1])*u)
        t -= L
    return ring[0]


def visible_brute(q, p, polys):
    """True if segment qp does not pass through any building interior.

    Tolerance-based: endpoint coordinates carry ~1e-15 rounding error, so an
    exact interior test misfires on points meant to lie on a boundary. We
    measure penetration length instead (running along a wall is allowed)."""
    seg = LineString([q, p])
    if seg.length == 0:
        return True
    for poly in polys:
        inter = seg.intersection(poly)
        pen = inter.length - seg.intersection(poly.exterior).length
        if pen > 1e-6:
            return False
    return True


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    rng = random.Random(seed)
    hard_fail = 0
    soft_fail = 0
    total_checked = 0
    for trial in range(20):
        polys = rand_scene(rng, rng.randint(3, 12))
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            scene_path = tf.name
        rings = write_scene(polys, scene_path)
        perims = [sum(math.dist(r[j], r[(j+1) % len(r)]) for j in range(len(r)))
                  for r in rings]

        for _ in range(8):
            # random q outside; occasionally epsilon-close to a boundary
            while True:
                if rng.random() < 0.3:
                    poly = rng.choice(polys)
                    t = rng.uniform(0, poly.exterior.length)
                    bp = poly.exterior.interpolate(t)
                    c = poly.centroid
                    d = math.dist((bp.x, bp.y), (c.x, c.y))
                    q = (bp.x + (bp.x-c.x)/d*1e-6, bp.y + (bp.y-c.y)/d*1e-6)
                else:
                    q = (rng.uniform(-10, 110), rng.uniform(-10, 110))
                if not any(pl.contains(Point(q)) for pl in polys):
                    break
            R = rng.uniform(20, 120)
            res = subprocess.run(
                [str(SOLVER), scene_path, "--testvis", repr(q[0]), repr(q[1]),
                 repr(R)],
                capture_output=True, text=True, check=True)
            ivs = {}
            for line in res.stdout.splitlines():
                b, t0, t1 = line.split()
                ivs.setdefault(int(b), []).append((float(t0), float(t1)))

            tol = 1e-5
            for bi, ring in enumerate(rings):
                L = perims[bi]
                biv = ivs.get(bi, [])
                # sample boundary points
                ts = [rng.uniform(0, L) for _ in range(25)]
                for t0, t1 in biv:
                    for e in (t0, t1):
                        ts += [max(0, e - 10*tol), min(L, e + 10*tol),
                               min(L, (t0+t1)/2)]
                for t in ts:
                    p = arc_point(ring, t)
                    claimed = any(a - 1e-9 <= t <= b + 1e-9 for a, b in biv)
                    near_edge = any(abs(t-a) < tol or abs(t-b) < tol
                                    for a, b in biv)
                    dist = math.dist(q, p)
                    vis = visible_brute(q, p, polys) and dist <= R + 1e-9
                    total_checked += 1
                    if claimed and not near_edge and not vis:
                        if hard_fail == 0:
                            import json
                            with open("/tmp/failcase.json", "w") as ff:
                                json.dump({"rings": rings, "q": q, "R": R,
                                           "bld": bi, "t": t}, ff)
                        hard_fail += 1
                        print(f"HARD FAIL trial={trial} q={q} R={R} bld={bi} "
                              f"t={t} p={p} dist={dist}")
                    if (not claimed) and (not near_edge) and vis and dist < R - tol:
                        soft_fail += 1
                        if soft_fail <= 10:
                            print(f"soft miss trial={trial} q={q} R={R} bld={bi} "
                                  f"t={t} p={p} dist={dist}")
    print(f"checked {total_checked} samples: hard_fail={hard_fail} "
          f"soft_miss={soft_fail}")
    sys.exit(1 if hard_fail else 0)


if __name__ == "__main__":
    main()
