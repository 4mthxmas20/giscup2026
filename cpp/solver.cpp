// GISCUP 2026 antenna-placement solver.
//
// Pipeline: load buildings -> generate boundary candidates -> batch
// visibility sweep (shared across all sub-problems) -> lazy greedy per tau
// with solution snapshots at each requested k.
#include <atomic>
#include <chrono>
#include <cstring>
#include <queue>
#include <thread>

#include "geometry.hpp"
#include "grid.hpp"
#include "visibility.hpp"

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
}

// ---------------------------------------------------------------- buildings
// Coarse grid of building indices (by bbox) for point-in-building tests.
class BuildingGrid {
  public:
    void build(const Dataset& ds, double cell) {
        cell_ = cell;
        x0_ = ds.minx - cell;
        y0_ = ds.miny - cell;
        nx_ = (int)std::ceil((ds.maxx - x0_ + cell) / cell) + 1;
        ny_ = (int)std::ceil((ds.maxy - y0_ + cell) / cell) + 1;
        cells_.assign((size_t)nx_ * ny_, {});
        for (size_t bi = 0; bi < ds.blds.size(); bi++) {
            double bx0 = 1e300, by0 = 1e300, bx1 = -1e300, by1 = -1e300;
            for (Pt p : ds.blds[bi].pts) {
                bx0 = std::min(bx0, p.x); bx1 = std::max(bx1, p.x);
                by0 = std::min(by0, p.y); by1 = std::max(by1, p.y);
            }
            int cx0 = cl(ix(bx0), nx_), cx1 = cl(ix(bx1), nx_);
            int cy0 = cl(iy(by0), ny_), cy1 = cl(iy(by1), ny_);
            for (int cy = cy0; cy <= cy1; cy++)
                for (int cx = cx0; cx <= cx1; cx++)
                    cells_[(size_t)cy * nx_ + cx].push_back((uint32_t)bi);
        }
    }
    // strictly-inside test (boundary counts as outside)
    bool insideAny(Pt p, const Dataset& ds) const {
        int cx = cl(ix(p.x), nx_), cy = cl(iy(p.y), ny_);
        for (uint32_t bi : cells_[(size_t)cy * nx_ + cx]) {
            const auto& pts = ds.blds[bi].pts;
            size_t n = pts.size();
            bool in = false;
            for (size_t i = 0, j = n - 1; i < n; j = i++) {
                if (distSeg2(p, pts[j], pts[i]) < 1e-18) { in = false; break; }
                if ((pts[i].y > p.y) != (pts[j].y > p.y)) {
                    double xx = pts[i].x + (pts[j].x - pts[i].x) * (p.y - pts[i].y) /
                                               (pts[j].y - pts[i].y);
                    if (p.x < xx) in = !in;
                }
            }
            if (in) return true;
        }
        return false;
    }

  private:
    int ix(double x) const { return (int)std::floor((x - x0_) / cell_); }
    int iy(double y) const { return (int)std::floor((y - y0_) / cell_); }
    static int cl(int c, int n) { return c < 0 ? 0 : (c >= n ? n - 1 : c); }
    double cell_ = 100, x0_ = 0, y0_ = 0;
    int nx_ = 1, ny_ = 1;
    std::vector<std::vector<uint32_t>> cells_;
};

// ---------------------------------------------------------------- candidates
struct Candidate {
    Pt p;        // exact point on the boundary (output coordinate)
    Pt q;        // epsilon-offset query point, strictly outside all buildings
    uint32_t bld;
};

static void genCandidates(const Dataset& ds, const BuildingGrid& bg,
                          double spacing, double eps,
                          std::vector<Candidate>& out) {
    size_t skipped = 0;
    for (size_t bi = 0; bi < ds.blds.size(); bi++) {
        const auto& pts = ds.blds[bi].pts;
        size_t n = pts.size();
        auto outwardEdge = [&](size_t j) {  // outward normal of edge j->j+1 (CCW ring)
            Pt d = pts[(j + 1) % n] - pts[j];
            double L = norm(d);
            return L > 0 ? Pt{d.y / L, -d.x / L} : Pt{0, 0};
        };
        for (size_t j = 0; j < n; j++) {
            // vertex candidate: offset along bisector of adjacent outward normals
            Pt n1 = outwardEdge((j + n - 1) % n), n2 = outwardEdge(j);
            Pt bis = n1 + n2;
            double L = norm(bis);
            bis = L > 1e-12 ? bis * (1.0 / L) : n2;
            Candidate c{pts[j], pts[j] + bis * eps, (uint32_t)bi};
            if (!bg.insideAny(c.q, ds)) out.push_back(c);
            else skipped++;
            // extra samples on long edges
            Pt a = pts[j], b = pts[(j + 1) % n];
            double elen = norm(b - a);
            int extra = (int)std::floor(elen / spacing);
            Pt nrm = outwardEdge(j);
            for (int t = 1; t <= extra; t++) {
                double u = (double)t / (extra + 1);
                Pt p = a + (b - a) * u;
                Candidate cc{p, p + nrm * eps, (uint32_t)bi};
                if (!bg.insideAny(cc.q, ds)) out.push_back(cc);
                else skipped++;
            }
        }
    }
    std::fprintf(stderr, "candidates: %zu (skipped %zu inside a neighbor)\n",
                 out.size(), skipped);
}

// ---------------------------------------------------------------- CSR store
// Packed visible-interval record.
struct IvRec {
    uint32_t bld;
    float t0, t1;
};

struct VisTable {
    std::vector<uint64_t> off;   // size ncand+1
    std::vector<IvRec> data;
};

static void batchVisibility(const Dataset& ds, const EdgeGrid& grid,
                            const std::vector<Candidate>& cands, double R,
                            double minIvLen, int nthreads, VisTable& vt) {
    size_t n = cands.size();
    std::vector<std::vector<IvRec>> results(n);
    std::atomic<size_t> next{0};
    auto worker = [&]() {
        SweepScratch sc;
        std::vector<VisIv> out;
        for (;;) {
            size_t i = next.fetch_add(64);
            if (i >= n) break;
            size_t hi = std::min(i + 64, n);
            for (; i < hi; i++) {
                sweepVisibility(cands[i].q, R, ds, grid, sc, out, minIvLen);
                auto& r = results[i];
                r.reserve(out.size());
                for (const VisIv& iv : out)
                    r.push_back({iv.bld, (float)iv.t0, (float)iv.t1});
            }
        }
    };
    std::vector<std::thread> ths;
    for (int t = 0; t < nthreads; t++) ths.emplace_back(worker);
    for (auto& t : ths) t.join();

    vt.off.assign(n + 1, 0);
    for (size_t i = 0; i < n; i++) vt.off[i + 1] = vt.off[i] + results[i].size();
    vt.data.resize(vt.off[n]);
    for (size_t i = 0; i < n; i++) {
        std::copy(results[i].begin(), results[i].end(), vt.data.begin() + vt.off[i]);
        results[i].clear();
        results[i].shrink_to_fit();
    }
}

// ---------------------------------------------------------------- greedy
// Sorted disjoint interval union with length tracking.
struct IvUnion {
    std::vector<std::pair<double, double>> ivs;
    double len = 0;

    // length that inserting [t0,t1] would add (no mutation)
    double addedLen(double t0, double t1) const {
        double add = t1 - t0;
        for (const auto& iv : ivs) {
            if (iv.second <= t0) continue;
            if (iv.first >= t1) break;
            add -= std::min(iv.second, t1) - std::max(iv.first, t0);
        }
        return add > 0 ? add : 0;
    }
    void insert(double t0, double t1) {
        size_t i = 0;
        while (i < ivs.size() && ivs[i].second < t0) i++;
        double lo = t0, hi = t1;
        size_t j = i;
        while (j < ivs.size() && ivs[j].first <= t1) {
            lo = std::min(lo, ivs[j].first);
            hi = std::max(hi, ivs[j].second);
            len -= ivs[j].second - ivs[j].first;
            j++;
        }
        ivs.erase(ivs.begin() + i, ivs.begin() + j);
        ivs.insert(ivs.begin() + i, {lo, hi});
        len += hi - lo;
    }
};

struct GreedyResult {
    std::vector<uint32_t> picks;                    // candidate indices, in order
    // snapshots[k] = serviced building indices right after pick #k
    std::vector<std::pair<int, std::vector<uint32_t>>> snapshots;
};

static void runGreedy(const Dataset& ds, const std::vector<Candidate>& cands,
                      const VisTable& vt, double tau, const std::vector<int>& ks,
                      double gamma, int nthreads, GreedyResult& res) {
    size_t nb = ds.blds.size(), nc = cands.size();
    std::vector<double> need(nb);
    for (size_t b = 0; b < nb; b++) need[b] = tau * ds.blds[b].perimeter;
    std::vector<IvUnion> cov(nb);
    std::vector<uint8_t> served(nb, 0);

    auto gainOf = [&](uint32_t c) {
        double g = 0;
        uint64_t lo = vt.off[c], hi = vt.off[c + 1];
        uint64_t i = lo;
        while (i < hi) {
            uint32_t b = vt.data[i].bld;
            if (served[b]) {
                while (i < hi && vt.data[i].bld == b) i++;
                continue;
            }
            double add = 0;
            while (i < hi && vt.data[i].bld == b) {
                add += cov[b].addedLen(vt.data[i].t0, vt.data[i].t1);
                i++;
            }
            if (add <= 0) continue;
            double c0 = cov[b].len;
            double eff = std::min(c0 + add, need[b]) - std::min(c0, need[b]);
            g += gamma * eff / need[b];
            if (c0 < need[b] && c0 + add >= need[b]) g += 1.0;
        }
        return g;
    };

    // initial gains in parallel
    std::vector<double> gain0(nc);
    {
        std::atomic<size_t> nextIdx{0};
        auto w = [&]() {
            for (;;) {
                size_t i = nextIdx.fetch_add(256);
                if (i >= nc) break;
                size_t hi2 = std::min(i + 256, nc);
                for (; i < hi2; i++) gain0[i] = gainOf((uint32_t)i);
            }
        };
        std::vector<std::thread> ths;
        for (int t = 0; t < nthreads; t++) ths.emplace_back(w);
        for (auto& t : ths) t.join();
    }

    struct HeapEnt {
        double g;
        uint32_t c;
        uint32_t round;
        bool operator<(const HeapEnt& o) const { return g < o.g; }
    };
    std::priority_queue<HeapEnt> heap;
    std::vector<uint8_t> chosen(nc, 0);
    for (size_t i = 0; i < nc; i++)
        heap.push({gain0[i], (uint32_t)i, 0});

    int kmax = 0;
    for (int k : ks) kmax = std::max(kmax, k);
    uint32_t round = 0;
    std::vector<int> ksSorted(ks);
    std::sort(ksSorted.begin(), ksSorted.end());
    size_t ksNext = 0;

    while ((int)res.picks.size() < kmax && !heap.empty()) {
        HeapEnt e = heap.top();
        heap.pop();
        if (chosen[e.c]) continue;
        if (e.round != round) {
            e.g = gainOf(e.c);
            e.round = round;
            heap.push(e);
            continue;
        }
        // commit
        chosen[e.c] = 1;
        uint64_t lo = vt.off[e.c], hi = vt.off[e.c + 1];
        for (uint64_t i = lo; i < hi; i++) {
            uint32_t b = vt.data[i].bld;
            cov[b].insert(vt.data[i].t0, vt.data[i].t1);
            if (!served[b] && cov[b].len >= need[b]) served[b] = 1;
        }
        res.picks.push_back(e.c);
        round++;
        while (ksNext < ksSorted.size() &&
               (int)res.picks.size() == ksSorted[ksNext]) {
            std::vector<uint32_t> sv;
            for (size_t b = 0; b < nb; b++)
                if (served[b]) sv.push_back((uint32_t)b);
            res.snapshots.push_back({ksSorted[ksNext], std::move(sv)});
            ksNext++;
        }
    }
}

// ---------------------------------------------------------------- main
static std::vector<double> parseList(const char* s) {
    std::vector<double> v;
    for (const char* p = s; *p;) {
        v.push_back(std::strtod(p, (char**)&p));
        if (*p == ',') p++;
    }
    return v;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
                     "usage: solver buildings.txt [--R m] [--spacing m] [--gamma g]\n"
                     "       [--taus 0.25,0.5,0.75] [--ks 50,500,1000] [--out dir]\n"
                     "       [--threads n] [--minivlen m]\n"
                     "       solver buildings.txt --testvis qx qy R\n");
        return 1;
    }
    double R = 400, spacing = 12, eps = 1e-6, gamma = 0.4, minIvLen = 0.02;
    std::vector<double> taus = {0.25, 0.5, 0.75};
    std::vector<int> ks = {50, 500, 1000};
    std::string outDir = "results";
    int nthreads = (int)std::thread::hardware_concurrency();
    bool testvis = false;
    double tvx = 0, tvy = 0, tvR = 0;

    for (int i = 2; i < argc; i++) {
        auto is = [&](const char* o) { return std::strcmp(argv[i], o) == 0; };
        if (is("--R")) R = std::atof(argv[++i]);
        else if (is("--spacing")) spacing = std::atof(argv[++i]);
        else if (is("--gamma")) gamma = std::atof(argv[++i]);
        else if (is("--minivlen")) minIvLen = std::atof(argv[++i]);
        else if (is("--threads")) nthreads = std::atoi(argv[++i]);
        else if (is("--out")) outDir = argv[++i];
        else if (is("--taus")) taus = parseList(argv[++i]);
        else if (is("--ks")) {
            ks.clear();
            for (double v : parseList(argv[++i])) ks.push_back((int)v);
        } else if (is("--testvis")) {
            testvis = true;
            tvx = std::atof(argv[++i]);
            tvy = std::atof(argv[++i]);
            tvR = std::atof(argv[++i]);
        }
    }

    Dataset ds;
    ds.load(argv[1]);
    std::fprintf(stderr, "loaded %zu buildings, %zu edges\n", ds.blds.size(),
                 ds.edges.size());
    EdgeGrid grid;
    grid.build(ds, 40.0);

    if (testvis) {
        SweepScratch sc;
        std::vector<VisIv> out;
        sweepVisibility({tvx, tvy}, tvR, ds, grid, sc, out, 0.0);
        for (const VisIv& iv : out)
            std::printf("%u %.17g %.17g\n", iv.bld, iv.t0, iv.t1);
        return 0;
    }

    auto t0 = Clock::now();
    BuildingGrid bg;
    bg.build(ds, 100.0);
    std::vector<Candidate> cands;
    genCandidates(ds, bg, spacing, eps, cands);

    VisTable vt;
    batchVisibility(ds, grid, cands, R, minIvLen, nthreads, vt);
    auto t1 = Clock::now();
    std::fprintf(stderr, "visibility: %.1fs, %zu candidates, %zu intervals (%.0f MB)\n",
                 secs(t0, t1), cands.size(), vt.data.size(),
                 vt.data.size() * sizeof(IvRec) / 1e6);

    for (double tau : taus) {
        auto tg0 = Clock::now();
        GreedyResult res;
        runGreedy(ds, cands, vt, tau, ks, gamma, nthreads, res);
        auto tg1 = Clock::now();
        for (auto& [k, servedIdx] : res.snapshots) {
            char path[512];
            std::snprintf(path, sizeof path, "%s/sol_t%g_k%d.txt", outDir.c_str(),
                          tau, k);
            FILE* f = std::fopen(path, "w");
            if (!f) { std::fprintf(stderr, "cannot write %s\n", path); return 1; }
            std::fprintf(f, "%g %d\n", tau, k);
            for (int i = 0; i < k; i++) {
                const Candidate& c = cands[res.picks[i]];
                std::fprintf(f, "%s%.17g %.17g", i ? ";" : "", c.p.x, c.p.y);
            }
            std::fprintf(f, "\n");
            for (size_t i = 0; i < servedIdx.size(); i++)
                std::fprintf(f, "%s%lld", i ? "," : "",
                             (long long)ds.blds[servedIdx[i]].id);
            std::fprintf(f, "\n");
            std::fclose(f);
            std::fprintf(stderr, "tau=%.2f k=%d: served %zu / %zu  (%.1fs)\n", tau,
                         k, servedIdx.size(), ds.blds.size(), secs(tg0, tg1));
        }
    }
    return 0;
}
