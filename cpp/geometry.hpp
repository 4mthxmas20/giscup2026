#pragma once
// Basic 2D geometry types and the building dataset.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

struct Pt {
    double x = 0, y = 0;
};
static inline Pt operator-(Pt a, Pt b) { return {a.x - b.x, a.y - b.y}; }
static inline Pt operator+(Pt a, Pt b) { return {a.x + b.x, a.y + b.y}; }
static inline Pt operator*(Pt a, double s) { return {a.x * s, a.y * s}; }
static inline double cross(Pt a, Pt b) { return a.x * b.y - a.y * b.x; }
static inline double dot(Pt a, Pt b) { return a.x * b.x + a.y * b.y; }
static inline double norm2(Pt a) { return a.x * a.x + a.y * a.y; }
static inline double norm(Pt a) { return std::sqrt(norm2(a)); }

struct Building {
    int64_t id = 0;
    std::vector<Pt> pts;        // CCW exterior ring, no closing vertex
    std::vector<double> cum;    // cum[i] = arc length from pts[0] to pts[i]; size n+1
    double perimeter = 0;
};

// One directed boundary edge pts[e] -> pts[e+1] of some building.
struct Edge {
    Pt a, b;
    uint32_t bld;       // building index
    double arc0;        // arc-length of point a along the building boundary
    double len;         // edge length
};

struct Dataset {
    std::vector<Building> blds;
    std::vector<Edge> edges;
    double minx = 0, miny = 0, maxx = 0, maxy = 0;

    void load(const std::string& path) {
        FILE* f = std::fopen(path.c_str(), "r");
        if (!f) throw std::runtime_error("cannot open " + path);
        size_t n = 0;
        if (std::fscanf(f, "%zu", &n) != 1) throw std::runtime_error("bad header");
        blds.resize(n);
        minx = miny = 1e300;
        maxx = maxy = -1e300;
        for (size_t i = 0; i < n; i++) {
            Building& b = blds[i];
            size_t nv = 0;
            if (std::fscanf(f, "%lld %zu", (long long*)&b.id, &nv) != 2)
                throw std::runtime_error("bad building header");
            b.pts.resize(nv);
            for (size_t j = 0; j < nv; j++) {
                if (std::fscanf(f, "%lf %lf", &b.pts[j].x, &b.pts[j].y) != 2)
                    throw std::runtime_error("bad vertex");
                minx = std::min(minx, b.pts[j].x);
                maxx = std::max(maxx, b.pts[j].x);
                miny = std::min(miny, b.pts[j].y);
                maxy = std::max(maxy, b.pts[j].y);
            }
            b.cum.resize(nv + 1);
            b.cum[0] = 0;
            for (size_t j = 0; j < nv; j++)
                b.cum[j + 1] = b.cum[j] + norm(b.pts[(j + 1) % nv] - b.pts[j]);
            b.perimeter = b.cum[nv];
        }
        std::fclose(f);
        for (size_t i = 0; i < blds.size(); i++) {
            const Building& b = blds[i];
            size_t nv = b.pts.size();
            for (size_t j = 0; j < nv; j++) {
                Edge e;
                e.a = b.pts[j];
                e.b = b.pts[(j + 1) % nv];
                e.bld = (uint32_t)i;
                e.arc0 = b.cum[j];
                e.len = b.cum[j + 1] - b.cum[j];
                if (e.len > 0) edges.push_back(e);
            }
        }
    }
};

// Squared distance from point p to segment ab.
static inline double distSeg2(Pt p, Pt a, Pt b) {
    Pt d = b - a;
    double L2 = norm2(d);
    double t = L2 > 0 ? dot(p - a, d) / L2 : 0;
    t = t < 0 ? 0 : (t > 1 ? 1 : t);
    return norm2(p - (a + d * t));
}
