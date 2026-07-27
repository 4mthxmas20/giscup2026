#!/usr/bin/env python3
"""Minimal repro harness: reuses test_visibility's scene generator to dump a
given (seed, trial, query#) case, prints claimed intervals vs brute truth per
sampled point."""
import math
import random
import subprocess
import sys
from pathlib import Path

from shapely.geometry import Point

import test_visibility as tv

ROOT = Path(__file__).resolve().parent.parent


def main():
    seed, trial_want, q_want = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(seed)
    for trial in range(trial_want + 1):
        polys = tv.rand_scene(rng, rng.randint(3, 12))
        scene_path = "/tmp/debug_scene.txt"
        rings = tv.write_scene(polys, scene_path)
        perims = [sum(math.dist(r[j], r[(j+1) % len(r)]) for j in range(len(r)))
                  for r in rings]
        for qi in range(8):
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
            if trial != trial_want or qi != q_want:
                continue
            print(f"scene at {scene_path}, q={q}, R={R}")
            print(f"polys: {[list(r) for r in rings]}")
            res = subprocess.run(
                [str(ROOT / "cpp" / "solver"), scene_path, "--testvis",
                 repr(q[0]), repr(q[1]), repr(R)],
                capture_output=True, text=True, check=True)
            ivs = {}
            for line in res.stdout.splitlines():
                b, t0, t1 = line.split()
                ivs.setdefault(int(b), []).append((float(t0), float(t1)))
            for bi, ring in enumerate(rings):
                biv = ivs.get(bi, [])
                print(f"bld {bi} perim={perims[bi]:.2f} claimed={biv}")
                for t0, t1 in biv:
                    tm = (t0 + t1) / 2
                    p = tv.arc_point(ring, tm)
                    ok = tv.visible_brute(q, p, polys)
                    d = math.dist(q, p)
                    print(f"   mid t={tm:.3f} p=({p[0]:.3f},{p[1]:.3f}) "
                          f"dist={d:.2f} brute_visible={ok}")
            return


if __name__ == "__main__":
    main()
