#pragma once
// Uniform grid over boundary edges for disk range queries.
#include <algorithm>
#include <cmath>
#include <vector>

#include "geometry.hpp"

class EdgeGrid {
  public:
    void build(const Dataset& ds, double cell) {
        cell_ = cell;
        x0_ = ds.minx - cell;
        y0_ = ds.miny - cell;
        nx_ = (int)std::ceil((ds.maxx - x0_ + cell) / cell) + 1;
        ny_ = (int)std::ceil((ds.maxy - y0_ + cell) / cell) + 1;
        cells_.assign((size_t)nx_ * ny_, {});
        for (size_t ei = 0; ei < ds.edges.size(); ei++) {
            const Edge& e = ds.edges[ei];
            // rasterize edge bbox (edges are short; bbox overestimate is fine)
            int cx0 = clampx(ix(std::min(e.a.x, e.b.x)));
            int cx1 = clampx(ix(std::max(e.a.x, e.b.x)));
            int cy0 = clampy(iy(std::min(e.a.y, e.b.y)));
            int cy1 = clampy(iy(std::max(e.a.y, e.b.y)));
            for (int cy = cy0; cy <= cy1; cy++)
                for (int cx = cx0; cx <= cx1; cx++) {
                    // skip cells the segment doesn't actually come near
                    double bx0 = x0_ + cx * cell_, by0 = y0_ + cy * cell_;
                    if (segNearBox(e.a, e.b, bx0, by0, bx0 + cell_, by0 + cell_))
                        cells_[(size_t)cy * nx_ + cx].push_back((uint32_t)ei);
                }
        }
        stamp_.assign(ds.edges.size(), 0);
        cur_ = 0;
    }

    // Collect indices of edges whose distance to q is <= R (superset is OK,
    // exact filtering happens in the caller). NOT thread-safe on stamp_; each
    // thread must own its own EdgeGrid copy or use collectInto with own stamp.
    void collect(Pt q, double R, const Dataset& ds, std::vector<uint32_t>& out,
                 std::vector<uint32_t>& stamp, uint32_t& cur) const {
        out.clear();
        if (stamp.size() != ds.edges.size()) stamp.assign(ds.edges.size(), 0);
        ++cur;
        int cx0 = clampx(ix(q.x - R)), cx1 = clampx(ix(q.x + R));
        int cy0 = clampy(iy(q.y - R)), cy1 = clampy(iy(q.y + R));
        double R2 = R * R;
        for (int cy = cy0; cy <= cy1; cy++) {
            // prune whole rows of cells outside the disk
            double celly0 = y0_ + cy * cell_, celly1 = celly0 + cell_;
            double dy = q.y < celly0 ? celly0 - q.y : (q.y > celly1 ? q.y - celly1 : 0);
            for (int cx = cx0; cx <= cx1; cx++) {
                double cellx0 = x0_ + cx * cell_, cellx1 = cellx0 + cell_;
                double dx = q.x < cellx0 ? cellx0 - q.x
                                         : (q.x > cellx1 ? q.x - cellx1 : 0);
                if (dx * dx + dy * dy > R2) continue;
                for (uint32_t ei : cells_[(size_t)cy * nx_ + cx]) {
                    if (stamp[ei] == cur) continue;
                    stamp[ei] = cur;
                    const Edge& e = ds.edges[ei];
                    if (distSeg2(q, e.a, e.b) <= R2) out.push_back(ei);
                }
            }
        }
    }

  private:
    static bool segNearBox(Pt a, Pt b, double x0, double y0, double x1, double y1) {
        // conservative: does segment ab intersect the (slightly inflated) box?
        double e = 1e-9;
        x0 -= e; y0 -= e; x1 += e; y1 += e;
        // quick accept: endpoint inside
        if ((a.x >= x0 && a.x <= x1 && a.y >= y0 && a.y <= y1) ||
            (b.x >= x0 && b.x <= x1 && b.y >= y0 && b.y <= y1))
            return true;
        // separating axis on box sides
        if (std::max(a.x, b.x) < x0 || std::min(a.x, b.x) > x1 ||
            std::max(a.y, b.y) < y0 || std::min(a.y, b.y) > y1)
            return false;
        // segment line vs box corners
        Pt d = b - a;
        double c1 = cross(d, Pt{x0, y0} - a), c2 = cross(d, Pt{x1, y0} - a);
        double c3 = cross(d, Pt{x0, y1} - a), c4 = cross(d, Pt{x1, y1} - a);
        double mx = std::max(std::max(c1, c2), std::max(c3, c4));
        double mn = std::min(std::min(c1, c2), std::min(c3, c4));
        return mx >= 0 && mn <= 0;
    }
    int ix(double x) const { return (int)std::floor((x - x0_) / cell_); }
    int iy(double y) const { return (int)std::floor((y - y0_) / cell_); }
    int clampx(int c) const { return c < 0 ? 0 : (c >= nx_ ? nx_ - 1 : c); }
    int clampy(int c) const { return c < 0 ? 0 : (c >= ny_ ? ny_ - 1 : c); }

    double cell_ = 50, x0_ = 0, y0_ = 0;
    int nx_ = 1, ny_ = 1;
    std::vector<std::vector<uint32_t>> cells_;
    std::vector<uint32_t> stamp_;
    uint32_t cur_ = 0;
};
