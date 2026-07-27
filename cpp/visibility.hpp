#pragma once
// Rotational-sweep visibility: from a query point q, compute the exactly
// visible sub-segments of all building boundary edges within radius R.
//
// Conservative-by-construction: obstacles considered are all edges
// intersecting disk(q, R), and visibility is only credited within that disk,
// so every credited interval is truly visible (never overcounts).
#include <algorithm>
#include <cmath>
#include <vector>

#include "geometry.hpp"
#include "grid.hpp"

// A visible piece of a building boundary, as an arc-length interval
// [t0, t1] along the (CCW) boundary parameterization, in meters.
// Never wraps: intervals crossing 0 are emitted split.
struct VisIv {
    uint32_t bld;
    double t0, t1;
};

struct SweepScratch {
    std::vector<uint32_t> edgeIds;
    std::vector<uint32_t> stamp;
    uint32_t stampCur = 0;
    struct Ev {
        double ang;
        uint32_t seg;   // index into segs
        bool add;
    };
    struct Seg {
        Pt a, b;        // oriented so that sweep angle increases a -> b
        bool flip;      // true if (a,b) is reversed vs boundary direction
        uint32_t edge;  // dataset edge index
    };
    std::vector<Ev> evs;
    std::vector<Seg> segs;
    std::vector<uint32_t> active;   // indices into segs
    std::vector<int32_t> activePos; // segs -> position in active, or -1
    std::vector<VisIv> raw;
};

// Ray from q with direction d against segment ab: returns distance r along d
// and sets s in [0,1] (clamped) along a->b. Returns +inf when parallel.
static inline double raySeg(Pt q, Pt d, Pt a, Pt b, double& s) {
    Pt ab = b - a;
    Pt aq = a - q;
    double den = cross(d, ab);
    if (std::fabs(den) < 1e-15) {
        s = 0;
        return HUGE_VAL;
    }
    s = -cross(d, aq) / den;
    if (s < 0) s = 0;
    else if (s > 1) s = 1;
    // recompute r from the clamped point for robustness
    Pt p = a + ab * s;
    double r = dot(p - q, d);
    return r;
}

// Compute visible intervals from q within radius R. q must lie strictly
// outside every building (callers offset boundary candidates outward).
inline void sweepVisibility(Pt q, double R, const Dataset& ds, const EdgeGrid& grid,
                            SweepScratch& sc, std::vector<VisIv>& out,
                            double minIvLen = 0.0) {
    out.clear();
    grid.collect(q, R, ds, sc.edgeIds, sc.stamp, sc.stampCur);
    if (sc.edgeIds.empty()) return;

    sc.evs.clear();
    sc.segs.clear();
    sc.active.clear();
    sc.raw.clear();

    const double R2 = R * R;
    // Build oriented segments + events
    for (uint32_t ei : sc.edgeIds) {
        const Edge& e = ds.edges[ei];
        Pt ra = e.a - q, rb = e.b - q;
        double cr = cross(ra, rb);
        // collinear with q: radial wall, does not block and gets no credit
        if (cr == 0) continue;
        SweepScratch::Seg s;
        s.edge = ei;
        if (cr > 0) { s.a = e.a; s.b = e.b; s.flip = false; }
        else        { s.a = e.b; s.b = e.a; s.flip = true; }
        double angA = std::atan2(s.a.y - q.y, s.a.x - q.x);
        double angB = std::atan2(s.b.y - q.y, s.b.x - q.x);
        uint32_t si = (uint32_t)sc.segs.size();
        sc.segs.push_back(s);
        if (angA <= angB) {
            sc.evs.push_back({angA, si, true});
            sc.evs.push_back({angB, si, false});
        } else {
            // wraps across the -pi/pi cut: active from start
            sc.active.push_back(si);
            sc.evs.push_back({angB, si, false});
            sc.evs.push_back({angA, si, true});
        }
    }
    if (sc.segs.empty()) return;

    sc.activePos.assign(sc.segs.size(), -1);
    for (size_t i = 0; i < sc.active.size(); i++)
        sc.activePos[sc.active[i]] = (int32_t)i;

    std::sort(sc.evs.begin(), sc.evs.end(),
              [](const SweepScratch::Ev& x, const SweepScratch::Ev& y) {
                  if (x.ang != y.ang) return x.ang < y.ang;
                  return x.add < y.add;  // removes before adds at equal angle
              });

    auto addActive = [&](uint32_t si) {
        if (sc.activePos[si] >= 0) return;
        sc.activePos[si] = (int32_t)sc.active.size();
        sc.active.push_back(si);
    };
    auto delActive = [&](uint32_t si) {
        int32_t p = sc.activePos[si];
        if (p < 0) return;
        uint32_t last = sc.active.back();
        sc.active[p] = last;
        sc.activePos[last] = p;
        sc.active.pop_back();
        sc.activePos[si] = -1;
    };

    // Emit the visible portion of the nearest segment(s) in window [t0,t1].
    auto processWindow = [&](double w0, double w1) {
        if (w1 - w0 < 1e-12 || sc.active.empty()) return;
        double mid = 0.5 * (w0 + w1);
        Pt dm{std::cos(mid), std::sin(mid)};
        double rmin = HUGE_VAL;
        for (uint32_t si : sc.active) {
            const auto& s = sc.segs[si];
            double sp;
            double r = raySeg(q, dm, s.a, s.b, sp);
            if (r > 1e-12 && r < rmin) rmin = r;
        }
        if (rmin == HUGE_VAL || rmin > R) return;
        Pt d0{std::cos(w0), std::sin(w0)}, d1{std::cos(w1), std::sin(w1)};
        // credit every segment tied for nearest (handles duplicated shared walls)
        for (uint32_t si : sc.active) {
            const auto& s = sc.segs[si];
            double sp;
            double r = raySeg(q, dm, s.a, s.b, sp);
            if (!(r > 1e-12) || r > rmin + 1e-9) continue;
            double s0, s1;
            raySeg(q, d0, s.a, s.b, s0);
            raySeg(q, d1, s.a, s.b, s1);
            if (s0 > s1) std::swap(s0, s1);
            // clip to distance <= R (|a-q + s*ab|^2 <= R^2)
            Pt aq = s.a - q, ab = s.b - s.a;
            double A = norm2(ab), B = 2 * dot(aq, ab), C = norm2(aq) - R2;
            if (A > 0) {
                double disc = B * B - 4 * A * C;
                if (disc <= 0) continue;  // entirely outside the disk
                double sq = std::sqrt(disc);
                double lo = (-B - sq) / (2 * A), hi = (-B + sq) / (2 * A);
                s0 = std::max(s0, lo);
                s1 = std::min(s1, hi);
            }
            if (s1 <= s0) continue;
            const Edge& e = ds.edges[s.edge];
            double u0 = s.flip ? 1.0 - s1 : s0;   // param along boundary dir
            double u1 = s.flip ? 1.0 - s0 : s1;
            sc.raw.push_back({e.bld, e.arc0 + u0 * e.len, e.arc0 + u1 * e.len});
        }
    };

    double prev = -M_PI;
    size_t i = 0;
    while (i < sc.evs.size()) {
        double ang = sc.evs[i].ang;
        processWindow(prev, ang);
        while (i < sc.evs.size() && sc.evs[i].ang == ang) {
            if (sc.evs[i].add) addActive(sc.evs[i].seg);
            else delActive(sc.evs[i].seg);
            i++;
        }
        prev = ang;
    }
    processWindow(prev, M_PI);

    // merge per building: sort by (bld, t0), coalesce touching intervals
    std::sort(sc.raw.begin(), sc.raw.end(), [](const VisIv& x, const VisIv& y) {
        if (x.bld != y.bld) return x.bld < y.bld;
        return x.t0 < y.t0;
    });
    for (const VisIv& iv : sc.raw) {
        if (!out.empty() && out.back().bld == iv.bld &&
            iv.t0 <= out.back().t1 + 1e-9) {
            if (iv.t1 > out.back().t1) out.back().t1 = iv.t1;
        } else {
            if (!out.empty() && out.back().t1 - out.back().t0 < minIvLen)
                out.pop_back();
            out.push_back(iv);
        }
    }
    if (!out.empty() && out.back().t1 - out.back().t0 < minIvLen) out.pop_back();
}
