// GISCUP 2026 antenna-placement solver.
//
// Pipeline: load buildings -> generate boundary candidates -> batch
// visibility sweep (shared across all sub-problems) -> lazy greedy per tau
// with solution snapshots at each requested k.
#include <atomic>
#include <chrono>
#include <cstring>
#include <map>
#include <queue>
#include <set>
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

// `only` (when non-null) restricts generation to the flagged buildings, which
// is how the refinement pass densifies just the region that still matters.
static void genCandidates(const Dataset& ds, const BuildingGrid& bg,
                          double spacing, double eps,
                          std::vector<Candidate>& out,
                          const std::vector<uint8_t>* only = nullptr,
                          bool quiet = false) {
    size_t skipped = 0;
    for (size_t bi = 0; bi < ds.blds.size(); bi++) {
        if (only && !(*only)[bi]) continue;
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
    if (!quiet)
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

// Mutable optimization state for one tau. Supports add/remove of antennas
// with exact per-building covered-interval unions.
struct OptState {
    const Dataset* ds = nullptr;
    const VisTable* vt = nullptr;
    double tau = 0, gamma = 0;
    std::vector<double> need;
    std::vector<IvUnion> cov;
    std::vector<uint8_t> served;
    std::vector<std::vector<uint32_t>> touch;  // building -> selected cands
    std::vector<uint8_t> isSel;
    std::vector<uint32_t> selected;
    int score = 0;

    void init(const Dataset& d, const VisTable& v, double tau_, double gamma_) {
        ds = &d;
        vt = &v;
        tau = tau_;
        gamma = gamma_;
        size_t nb = d.blds.size();
        need.resize(nb);
        // +1mm slack: guards against float32 interval rounding so that every
        // claimed building clears the exact threshold with margin
        for (size_t b = 0; b < nb; b++)
            need[b] = tau * d.blds[b].perimeter + 1e-3;
        cov.assign(nb, {});
        served.assign(nb, 0);
        touch.assign(nb, {});
        isSel.assign(v.off.size() - 1, 0);
        selected.clear();
        score = 0;
    }

    // surrogate gain: 1.0 per newly-served building + gamma * capped progress
    double gainOf(uint32_t c) const {
        double g = 0;
        uint64_t i = vt->off[c], hi = vt->off[c + 1];
        while (i < hi) {
            uint32_t b = vt->data[i].bld;
            if (served[b]) {
                while (i < hi && vt->data[i].bld == b) i++;
                continue;
            }
            double add = 0;
            while (i < hi && vt->data[i].bld == b) {
                add += cov[b].addedLen(vt->data[i].t0, vt->data[i].t1);
                i++;
            }
            if (add <= 0) continue;
            double c0 = cov[b].len;
            double eff = std::min(c0 + add, need[b]) - std::min(c0, need[b]);
            g += gamma * eff / need[b];
            if (c0 < need[b] && c0 + add >= need[b]) g += 1.0;
        }
        return g;
    }

    void add(uint32_t c) {
        isSel[c] = 1;
        selected.push_back(c);
        for (uint64_t i = vt->off[c]; i < vt->off[c + 1]; i++) {
            uint32_t b = vt->data[i].bld;
            if (touch[b].empty() || touch[b].back() != c) touch[b].push_back(c);
            cov[b].insert(vt->data[i].t0, vt->data[i].t1);
            if (!served[b] && cov[b].len >= need[b]) {
                served[b] = 1;
                score++;
            }
        }
    }

    // intervals of candidate c on building b (vt columns are sorted by bld)
    std::pair<uint64_t, uint64_t> colOn(uint32_t c, uint32_t b) const {
        uint64_t lo = vt->off[c], hi = vt->off[c + 1];
        const IvRec* base = vt->data.data();
        const IvRec* s = std::lower_bound(base + lo, base + hi, b,
            [](const IvRec& r, uint32_t bb) { return r.bld < bb; });
        const IvRec* e = s;
        while (e < base + hi && e->bld == b) e++;
        return {(uint64_t)(s - base), (uint64_t)(e - base)};
    }

    void remove(uint32_t c) {
        isSel[c] = 0;
        selected.erase(std::find(selected.begin(), selected.end(), c));
        uint64_t i = vt->off[c], hi = vt->off[c + 1];
        while (i < hi) {
            uint32_t b = vt->data[i].bld;
            while (i < hi && vt->data[i].bld == b) i++;
            auto& tb = touch[b];
            tb.erase(std::find(tb.begin(), tb.end(), c));
            // rebuild b's union from remaining antennas
            cov[b].ivs.clear();
            cov[b].len = 0;
            for (uint32_t a : tb) {
                auto [s, e] = colOn(a, b);
                for (uint64_t j = s; j < e; j++)
                    cov[b].insert(vt->data[j].t0, vt->data[j].t1);
            }
            bool sv = cov[b].len >= need[b];
            if (served[b] && !sv) { served[b] = 0; score--; }
            else if (!served[b] && sv) { served[b] = 1; score++; }
        }
    }

    // exact drop in served count if c were removed
    int removalLoss(uint32_t c) const {
        int loss = 0;
        uint64_t i = vt->off[c], hi = vt->off[c + 1];
        while (i < hi) {
            uint32_t b = vt->data[i].bld;
            while (i < hi && vt->data[i].bld == b) i++;
            if (!served[b]) continue;
            IvUnion u;
            for (uint32_t a : touch[b]) {
                if (a == c) continue;
                auto [s, e] = colOn(a, b);
                for (uint64_t j = s; j < e; j++)
                    u.insert(vt->data[j].t0, vt->data[j].t1);
                if (u.len >= need[b]) break;
            }
            if (u.len < need[b]) loss++;
        }
        return loss;
    }

    std::vector<uint32_t> servedIdx() const {
        std::vector<uint32_t> sv;
        for (size_t b = 0; b < served.size(); b++)
            if (served[b]) sv.push_back((uint32_t)b);
        return sv;
    }

    // Rebuild the whole state from an antenna set. Cheaper than it looks:
    // dominated by the adds, which cost the same as one greedy step each.
    void setSelection(const std::vector<uint32_t>& sel) {
        std::fill(cov.begin(), cov.end(), IvUnion{});
        std::fill(served.begin(), served.end(), (uint8_t)0);
        for (auto& t : touch) t.clear();
        std::fill(isSel.begin(), isSel.end(), (uint8_t)0);
        selected.clear();
        score = 0;
        for (uint32_t c : sel) add(c);
    }
};

// Parallel argmax of state.gainOf over all non-selected candidates.
static uint32_t bestCandidate(const OptState& st, int nthreads,
                              const std::vector<uint8_t>* exclude = nullptr) {
    size_t nc = st.vt->off.size() - 1;
    std::vector<std::pair<double, uint32_t>> best(nthreads, {-1.0, UINT32_MAX});
    std::atomic<size_t> next{0};
    auto w = [&](int tid) {
        std::pair<double, uint32_t> loc{-1.0, UINT32_MAX};
        for (;;) {
            size_t i = next.fetch_add(1024);
            if (i >= nc) break;
            size_t hi = std::min(i + 1024, nc);
            for (; i < hi; i++) {
                if (st.isSel[i] || (exclude && (*exclude)[i])) continue;
                double g = st.gainOf((uint32_t)i);
                if (g > loc.first) loc = {g, (uint32_t)i};
            }
        }
        best[tid] = loc;
    };
    std::vector<std::thread> ths;
    for (int t = 0; t < nthreads; t++) ths.emplace_back(w, t);
    for (auto& t : ths) t.join();
    std::pair<double, uint32_t> r{-1.0, UINT32_MAX};
    for (auto& b : best)
        if (b.first > r.first) r = b;
    return r.second;
}

struct GreedyResult {
    std::vector<uint32_t> picks;
    std::vector<std::pair<int, int>> scores;  // (k, served) at snapshots
};

static void runGreedy(OptState& st, const std::vector<int>& ks, int nthreads,
                      GreedyResult& res) {
    size_t nc = st.vt->off.size() - 1;
    std::vector<double> gain0(nc);
    {
        std::atomic<size_t> nextIdx{0};
        auto w = [&]() {
            for (;;) {
                size_t i = nextIdx.fetch_add(256);
                if (i >= nc) break;
                size_t hi2 = std::min(i + 256, nc);
                for (; i < hi2; i++) gain0[i] = st.gainOf((uint32_t)i);
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
    for (size_t i = 0; i < nc; i++) heap.push({gain0[i], (uint32_t)i, 0});

    int kmax = 0;
    for (int k : ks) kmax = std::max(kmax, k);
    std::vector<int> ksSorted(ks);
    std::sort(ksSorted.begin(), ksSorted.end());
    size_t ksNext = 0;
    uint32_t round = 0;

    while ((int)res.picks.size() < kmax && !heap.empty()) {
        HeapEnt e = heap.top();
        heap.pop();
        if (st.isSel[e.c]) continue;
        if (e.round != round) {
            e.g = st.gainOf(e.c);
            e.round = round;
            heap.push(e);
            continue;
        }
        st.add(e.c);
        res.picks.push_back(e.c);
        round++;
        while (ksNext < ksSorted.size() &&
               (int)res.picks.size() == ksSorted[ksNext]) {
            res.scores.push_back({ksSorted[ksNext], st.score});
            ksNext++;
        }
    }
}

// ---------------------------------------------------------------- refinement
// Densify candidates only where they can still buy buildings: around footprints
// whose coverage sits just short of tau. Halving `spacing` globally would
// multiply the candidate set and slow every LNS iteration; aiming the extra
// resolution at the margin keeps the scan cheap. Self-limiting by construction
// — a saturated configuration has no near-miss buildings and refines nothing.
static size_t refineCandidates(const Dataset& ds, const BuildingGrid& bg,
                               const EdgeGrid& grid, const OptState& st,
                               double band, double radius, double refSpacing,
                               double eps, double R, double minIvLen,
                               int nthreads, std::vector<Candidate>& cands,
                               VisTable& vt) {
    size_t nb = ds.blds.size();
    struct Box { double x0, y0, x1, y1; };
    std::vector<Box> box(nb);
    for (size_t b = 0; b < nb; b++) {
        Box& e = box[b];
        e.x0 = e.y0 = 1e300;
        e.x1 = e.y1 = -1e300;
        for (Pt p : ds.blds[b].pts) {
            e.x0 = std::min(e.x0, p.x); e.x1 = std::max(e.x1, p.x);
            e.y0 = std::min(e.y0, p.y); e.y1 = std::max(e.y1, p.y);
        }
    }

    std::vector<uint32_t> nearMiss;
    for (size_t b = 0; b < nb; b++) {
        double have = st.cov[b].len, want = st.need[b];
        if (have < want && have >= (1.0 - band) * want)
            nearMiss.push_back((uint32_t)b);
    }
    if (nearMiss.empty()) return 0;

    // Bucket the near-miss footprints so the neighbourhood test stays linear.
    double cell = std::max(radius, 1.0);
    std::map<std::pair<long, long>, std::vector<uint32_t>> buckets;
    auto cellOf = [&](double x, double y) {
        return std::make_pair((long)std::floor(x / cell), (long)std::floor(y / cell));
    };
    for (uint32_t b : nearMiss) {
        const Box& e = box[b];
        buckets[cellOf((e.x0 + e.x1) / 2, (e.y0 + e.y1) / 2)].push_back(b);
    }

    std::vector<uint8_t> refine(nb, 0);
    size_t nRegion = 0;
    for (size_t b = 0; b < nb; b++) {
        const Box& e = box[b];
        auto [cx, cy] = cellOf((e.x0 + e.x1) / 2, (e.y0 + e.y1) / 2);
        bool hit = false;
        for (long dx = -1; dx <= 1 && !hit; dx++) {
            for (long dy = -1; dy <= 1 && !hit; dy++) {
                auto it = buckets.find({cx + dx, cy + dy});
                if (it == buckets.end()) continue;
                for (uint32_t m : it->second) {
                    const Box& o = box[m];
                    double ddx = std::max({o.x0 - e.x1, e.x0 - o.x1, 0.0});
                    double ddy = std::max({o.y0 - e.y1, e.y0 - o.y1, 0.0});
                    if (ddx * ddx + ddy * ddy <= radius * radius) { hit = true; break; }
                }
            }
        }
        if (hit) { refine[b] = 1; nRegion++; }
    }

    std::vector<Candidate> fresh;
    genCandidates(ds, bg, refSpacing, eps, fresh, &refine, true);

    // Vertices are regenerated unconditionally, so without this the refined set
    // would duplicate existing points — and two antennas on one coordinate are
    // counted once by the evaluator while still consuming k.
    std::set<std::pair<double, double>> known;
    for (const Candidate& c : cands) known.insert({c.p.x, c.p.y});
    std::vector<Candidate> added;
    for (const Candidate& c : fresh)
        if (known.insert({c.p.x, c.p.y}).second) added.push_back(c);
    if (added.empty()) return 0;

    VisTable nvt;
    batchVisibility(ds, grid, added, R, minIvLen, nthreads, nvt);

    uint64_t base = vt.off.back();
    for (size_t i = 1; i < nvt.off.size(); i++) vt.off.push_back(base + nvt.off[i]);
    vt.data.insert(vt.data.end(), nvt.data.begin(), nvt.data.end());
    cands.insert(cands.end(), added.begin(), added.end());

    std::fprintf(stderr,
                 "  refine: %zu near-miss buildings, %zu in region, "
                 "+%zu candidates (total %zu)\n",
                 nearMiss.size(), nRegion, added.size(), cands.size());
    return added.size();
}

// Large-neighborhood search: ruin a handful of antennas and rebuild greedily.
// A 1-opt swap can only cross a barrier one antenna at a time, which stalls at
// high tau where a building needs several antennas together; removing m at once
// lets the rebuild discover those combinations.
static void lnsSearch(OptState& st, double seconds, int nthreads, uint64_t seed,
                      int maxRuin, const std::vector<Candidate>& cands) {
    auto tstart = Clock::now();
    uint64_t rng = seed ? seed : 0x9e3779b97f4a7c15ull;
    auto rnd = [&]() {
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        return rng;
    };
    std::vector<uint32_t> best = st.selected;
    int bestScore = st.score;
    std::vector<uint8_t> excl(st.isSel.size(), 0);
    std::vector<uint32_t> victims, added;
    std::vector<std::pair<double, uint32_t>> nearest;
    int iters = 0, improves = 0, undos = 0;

    while (secs(tstart, Clock::now()) < seconds) {
        if (st.score == (int)st.served.size()) break;  // every building served
        iters++;
        int scoreBefore = st.score;
        int m = 1 + (int)(rnd() % (uint64_t)std::max(1, maxRuin));
        m = std::min(m, (int)st.selected.size());

        victims.clear();
        added.clear();
        if (st.selected.empty()) break;
        uint32_t seedAnt = st.selected[rnd() % st.selected.size()];
        if (rnd() % 100 < 60) {
            // Related removal: take the seed's nearest selected neighbours.
            // Antennas serving the same block are the ones worth recombining;
            // a scattered random set rarely has a better joint arrangement.
            nearest.clear();
            Pt sp = cands[seedAnt].p;
            for (uint32_t s : st.selected)
                nearest.push_back({norm2(cands[s].p - sp), s});
            size_t take = std::min((size_t)m, nearest.size());
            std::partial_sort(nearest.begin(), nearest.begin() + take, nearest.end());
            for (size_t i = 0; i < take; i++) victims.push_back(nearest[i].second);
        } else {
            victims.push_back(seedAnt);
            for (int i = 1; i < m && (int)victims.size() < m; i++) {
                uint32_t v = st.selected[rnd() % st.selected.size()];
                if (std::find(victims.begin(), victims.end(), v) == victims.end())
                    victims.push_back(v);
            }
        }
        for (uint32_t v : victims) st.remove(v);
        if (victims.empty()) break;

        // Forbid the victims on the first rebuild pick so each iteration
        // genuinely moves; afterwards they may come back if they are best.
        for (uint32_t v : victims) excl[v] = 1;
        for (size_t i = 0; i < victims.size(); i++) {
            uint32_t c = bestCandidate(st, nthreads, &excl);
            if (i == 0) {
                for (uint32_t v : victims) excl[v] = 0;
            }
            if (c == UINT32_MAX) break;
            st.add(c);
            added.push_back(c);
        }
        for (uint32_t v : victims) excl[v] = 0;

        // A short rebuild would silently shrink the antenna set below k, so
        // treat it as a failed move even when the score held.
        if (st.score < scoreBefore || added.size() != victims.size()) {
            // Undo incrementally: m removes + m adds, not a full k-add rebuild.
            for (uint32_t c : added) st.remove(c);
            for (uint32_t v : victims) st.add(v);
            undos++;
        } else if (st.score > bestScore) {
            bestScore = st.score;
            best = st.selected;
            improves++;
        }
        // equal score: keep the new configuration and keep walking the plateau
    }
    if (st.score < bestScore) st.setSelection(best);
    std::fprintf(stderr, "  lns: %d iters (%d improves, %d undos), %.1fs\n",
                 iters, improves, undos, secs(tstart, Clock::now()));
}

// Swap-based local search: repeatedly evict the antenna with the smallest
// exact removal loss and re-add the globally best replacement.
static void localSearch(OptState& st, double seconds, int nthreads,
                        uint64_t seed) {
    auto tstart = Clock::now();
    uint64_t rng = seed ? seed : 0x9e3779b97f4a7c15ull;
    auto rnd = [&]() {
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        return rng;
    };
    std::vector<uint8_t> excl(st.isSel.size(), 0);
    int sinceImprove = 0;
    int iters = 0;
    while (secs(tstart, Clock::now()) < seconds && sinceImprove < 5000) {
        if (st.score == (int)st.served.size()) break;  // every building served
        iters++;
        // victim = random antenna among those with minimal removal loss;
        // occasionally a fully random one to escape plateaus
        int minLoss = INT32_MAX;
        std::vector<uint32_t> cands;
        if (rnd() % 100 < 15) {
            cands.push_back(st.selected[rnd() % st.selected.size()]);
            minLoss = 0;
        } else {
            std::vector<int> loss(st.selected.size());
            std::atomic<size_t> next{0};
            auto w = [&]() {
                for (;;) {
                    size_t i = next.fetch_add(8);
                    if (i >= st.selected.size()) break;
                    size_t hi = std::min(i + 8, st.selected.size());
                    for (; i < hi; i++) loss[i] = st.removalLoss(st.selected[i]);
                }
            };
            std::vector<std::thread> ths;
            for (int t = 0; t < nthreads; t++) ths.emplace_back(w);
            for (auto& t : ths) t.join();
            for (size_t i = 0; i < st.selected.size(); i++)
                minLoss = std::min(minLoss, loss[i]);
            for (size_t i = 0; i < st.selected.size(); i++)
                if (loss[i] == minLoss) cands.push_back(st.selected[i]);
        }
        uint32_t victim = cands[rnd() % cands.size()];
        int before = st.score;
        st.remove(victim);
        excl[victim] = 1;
        uint32_t repl = bestCandidate(st, nthreads, &excl);
        excl[victim] = 0;
        if (repl == UINT32_MAX) { st.add(victim); break; }
        st.add(repl);
        if (st.score < before) {  // revert
            st.remove(repl);
            st.add(victim);
            sinceImprove++;
        } else if (st.score == before) {
            sinceImprove++;
        } else {
            sinceImprove = 0;
        }
    }
    std::fprintf(stderr, "  localsearch: %d iters, %.1fs\n", iters,
                 secs(tstart, Clock::now()));
}

// Load a previous solution's antennas back into candidate indices, so a run
// can pick up where an earlier one stopped instead of re-deriving from greedy.
// Requires an identical candidate set (same dataset, R and spacing).
static bool loadWarmStart(const std::string& path,
                          const std::vector<Candidate>& cands,
                          std::vector<uint32_t>& picks) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return false;
    char line[64];
    if (!std::fgets(line, sizeof line, f)) { std::fclose(f); return false; }
    std::string coords;
    int ch;
    while ((ch = std::fgetc(f)) != EOF && ch != '\n') coords.push_back((char)ch);
    std::fclose(f);

    std::map<std::pair<double, double>, uint32_t> byPoint;
    for (size_t i = 0; i < cands.size(); i++)
        byPoint.emplace(std::make_pair(cands[i].p.x, cands[i].p.y), (uint32_t)i);

    picks.clear();
    size_t pos = 0;
    while (pos < coords.size()) {
        size_t sep = coords.find(';', pos);
        std::string tok = coords.substr(pos, sep == std::string::npos
                                                 ? std::string::npos : sep - pos);
        double x = 0, y = 0;
        if (std::sscanf(tok.c_str(), "%lf %lf", &x, &y) == 2) {
            auto it = byPoint.find({x, y});
            if (it == byPoint.end()) {
                std::fprintf(stderr,
                             "warmstart: (%.17g, %.17g) is not in the candidate "
                             "set; R/spacing must match the original run\n", x, y);
                return false;
            }
            picks.push_back(it->second);
        }
        if (sep == std::string::npos) break;
        pos = sep + 1;
    }
    return !picks.empty();
}

// ---------------------------------------------------------------- vis cache
static uint64_t hashBytes(uint64_t h, const void* data, size_t n) {
    const unsigned char* p = (const unsigned char*)data;
    for (size_t i = 0; i < n; i++) {
        h ^= (uint64_t)p[i];
        h *= 1099511628211ull;
    }
    return h;
}
static uint64_t hashU64(uint64_t h, uint64_t v) {
    return hashBytes(h, &v, sizeof v);
}
static uint64_t hashDouble(uint64_t h, double v) {
    return hashBytes(h, &v, sizeof v);
}
static uint64_t datasetFingerprint(const Dataset& ds) {
    uint64_t h = 1469598103934665603ull;
    h = hashU64(h, (uint64_t)ds.blds.size());
    for (const Building& b : ds.blds) {
        h = hashU64(h, (uint64_t)b.id);
        h = hashU64(h, (uint64_t)b.pts.size());
        for (Pt p : b.pts) {
            h = hashDouble(h, p.x);
            h = hashDouble(h, p.y);
        }
    }
    return h;
}

// Cache format v2:
// [magic 8][dataset fingerprint u64][R double][spacing double][eps double]
// [minIvLen double][ncand u64][off (ncand+1) u64][ndata u64][IvRec...]
// [cand p/q/bld]. Mismatched caches are ignored and rebuilt.
static bool saveVis(const std::string& path, const VisTable& vt,
                    const std::vector<Candidate>& cands, uint64_t fingerprint,
                    double R, double spacing, double eps, double minIvLen) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    const char magic[8] = {'G', 'I', 'S', 'V', 'I', 'S', '2', '\0'};
    uint64_t nc = cands.size(), nd = vt.data.size();
    std::fwrite(magic, 1, sizeof magic, f);
    std::fwrite(&fingerprint, 8, 1, f);
    std::fwrite(&R, sizeof R, 1, f);
    std::fwrite(&spacing, sizeof spacing, 1, f);
    std::fwrite(&eps, sizeof eps, 1, f);
    std::fwrite(&minIvLen, sizeof minIvLen, 1, f);
    std::fwrite(&nc, 8, 1, f);
    std::fwrite(vt.off.data(), 8, vt.off.size(), f);
    std::fwrite(&nd, 8, 1, f);
    std::fwrite(vt.data.data(), sizeof(IvRec), nd, f);
    std::fwrite(cands.data(), sizeof(Candidate), nc, f);
    std::fclose(f);
    return true;
}
static bool loadVis(const std::string& path, VisTable& vt,
                    std::vector<Candidate>& cands, uint64_t fingerprint,
                    double R, double spacing, double eps, double minIvLen) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return false;
    char magic[8];
    const char want[8] = {'G', 'I', 'S', 'V', 'I', 'S', '2', '\0'};
    if (std::fread(magic, 1, sizeof magic, f) != sizeof magic ||
        std::memcmp(magic, want, sizeof magic) != 0) {
        std::fprintf(stderr, "ignoring legacy or invalid vis cache: %s\n",
                     path.c_str());
        std::fclose(f);
        return false;
    }
    uint64_t gotFingerprint = 0;
    double gotR = 0, gotSpacing = 0, gotEps = 0, gotMinIvLen = 0;
    if (std::fread(&gotFingerprint, 8, 1, f) != 1 ||
        std::fread(&gotR, sizeof gotR, 1, f) != 1 ||
        std::fread(&gotSpacing, sizeof gotSpacing, 1, f) != 1 ||
        std::fread(&gotEps, sizeof gotEps, 1, f) != 1 ||
        std::fread(&gotMinIvLen, sizeof gotMinIvLen, 1, f) != 1) {
        std::fclose(f);
        return false;
    }
    if (gotFingerprint != fingerprint || gotR != R || gotSpacing != spacing ||
        gotEps != eps || gotMinIvLen != minIvLen) {
        std::fprintf(stderr,
                     "ignoring mismatched vis cache: %s "
                     "(dataset/R/spacing/eps/minivlen changed)\n",
                     path.c_str());
        std::fclose(f);
        return false;
    }
    uint64_t nc = 0, nd = 0;
    if (std::fread(&nc, 8, 1, f) != 1) { std::fclose(f); return false; }
    vt.off.resize(nc + 1);
    if (std::fread(vt.off.data(), 8, nc + 1, f) != nc + 1) { std::fclose(f); return false; }
    if (std::fread(&nd, 8, 1, f) != 1) { std::fclose(f); return false; }
    vt.data.resize(nd);
    if (std::fread(vt.data.data(), sizeof(IvRec), nd, f) != nd) { std::fclose(f); return false; }
    cands.resize(nc);
    if (std::fread(cands.data(), sizeof(Candidate), nc, f) != nc) { std::fclose(f); return false; }
    std::fclose(f);
    return true;
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
    double R = 400, spacing = 12, eps = 1e-6, minIvLen = 0.02, lsTime = 0;
    double lnsTime = 0;
    // Tuned on the sample at tau=0.75/k=500: ruin 3 beat 6 and 12, because
    // cheaper iterations buy far more of them and the greedy rebuild cannot
    // reassemble a large ruin well.
    int maxRuin = 3;
    std::vector<double> taus = {0.25, 0.5, 0.75};
    std::vector<double> gammas = {0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5};
    std::vector<int> ks = {50, 500, 1000};
    std::string outDir = "results";
    std::string visCache, warmStart;
    int nthreads = (int)std::thread::hardware_concurrency();
    bool testvis = false, witness = false;
    int probeN = 0;
    // Refinement defaults: a building within 35% of its threshold is close
    // enough to be worth chasing, and 150 m covers the sight lines that can
    // realistically close such a gap.
    bool refine = false;
    double refineBand = 0.35, refineRadius = 150.0, refineSpacing = 2.0;
    bool visOnly = false;      // build the cache and exit, for a shared CI stage
    uint64_t seedBase = 0;     // vary to run independent searches in parallel
    double tvx = 0, tvy = 0, tvR = 0;

    for (int i = 2; i < argc; i++) {
        auto is = [&](const char* o) { return std::strcmp(argv[i], o) == 0; };
        if (is("--R")) R = std::atof(argv[++i]);
        else if (is("--spacing")) spacing = std::atof(argv[++i]);
        else if (is("--gammas")) gammas = parseList(argv[++i]);
        else if (is("--lstime")) lsTime = std::atof(argv[++i]);
        else if (is("--lnstime")) lnsTime = std::atof(argv[++i]);
        else if (is("--maxruin")) maxRuin = std::atoi(argv[++i]);
        else if (is("--minivlen")) minIvLen = std::atof(argv[++i]);
        else if (is("--threads")) nthreads = std::atoi(argv[++i]);
        else if (is("--out")) outDir = argv[++i];
        else if (is("--viscache")) visCache = argv[++i];
        else if (is("--witness")) witness = true;
        else if (is("--probe")) probeN = std::atoi(argv[++i]);
        else if (is("--warmstart")) warmStart = argv[++i];
        else if (is("--refine")) refine = true;
        else if (is("--visonly")) visOnly = true;
        else if (is("--seed")) seedBase = (uint64_t)std::strtoull(argv[++i], nullptr, 10);
        else if (is("--refineband")) refineBand = std::atof(argv[++i]);
        else if (is("--refineradius")) refineRadius = std::atof(argv[++i]);
        else if (is("--refinespacing")) refineSpacing = std::atof(argv[++i]);
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
    uint64_t visFingerprint = datasetFingerprint(ds);
    EdgeGrid grid;
    grid.build(ds, 40.0);

    if (probeN > 0) {
        // Time a random sample of candidates so the operator can pick R and
        // spacing that fit the 24-hour window on an unseen dataset.
        BuildingGrid bg;
        bg.build(ds, 100.0);
        std::vector<Candidate> cands;
        genCandidates(ds, bg, spacing, eps, cands);
        size_t n = std::min((size_t)probeN, cands.size());
        std::vector<Candidate> sample;
        uint64_t rs = 12345;
        for (size_t i = 0; i < n; i++) {
            rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
            sample.push_back(cands[rs % cands.size()]);
        }
        VisTable svt;
        auto p0 = Clock::now();
        batchVisibility(ds, grid, sample, R, minIvLen, nthreads, svt);
        double el = secs(p0, Clock::now());
        double perCand = el / (double)n;
        double ivPerCand = (double)svt.data.size() / (double)n;
        std::printf("probe R=%g spacing=%g threads=%d\n", R, spacing, nthreads);
        std::printf("  candidates total : %zu\n", cands.size());
        std::printf("  sampled          : %zu in %.2fs (%.3f ms/candidate)\n",
                    n, el, perCand * 1000);
        std::printf("  est. visibility  : %.0fs (%.1f min)\n",
                    perCand * cands.size(), perCand * cands.size() / 60);
        std::printf("  est. intervals   : %.0f (%.0f MB)\n",
                    ivPerCand * cands.size(),
                    ivPerCand * cands.size() * sizeof(IvRec) / 1e6);
        return 0;
    }

    if (testvis) {
        SweepScratch sc;
        std::vector<VisIv> out;
        sweepVisibility({tvx, tvy}, tvR, ds, grid, sc, out, 0.0);
        for (const VisIv& iv : out)
            std::printf("%u %.17g %.17g\n", iv.bld, iv.t0, iv.t1);
        return 0;
    }

    auto t0 = Clock::now();
    std::vector<Candidate> cands;
    VisTable vt;
    BuildingGrid bg;  // also needed by refinement, so build it either way
    bg.build(ds, 100.0);
    if (!visCache.empty() &&
        loadVis(visCache, vt, cands, visFingerprint, R, spacing, eps, minIvLen)) {
        std::fprintf(stderr, "loaded vis cache: %zu candidates, %zu intervals\n",
                     cands.size(), vt.data.size());
    } else {
        genCandidates(ds, bg, spacing, eps, cands);
        batchVisibility(ds, grid, cands, R, minIvLen, nthreads, vt);
        auto t1 = Clock::now();
        std::fprintf(stderr,
                     "visibility: %.1fs, %zu candidates, %zu intervals (%.0f MB)\n",
                     secs(t0, t1), cands.size(), vt.data.size(),
                     vt.data.size() * sizeof(IvRec) / 1e6);
        if (!visCache.empty())
            saveVis(visCache, vt, cands, visFingerprint, R, spacing, eps,
                    minIvLen);
    }

    if (visOnly) {
        std::fprintf(stderr, "visonly: cache ready (%zu candidates, %zu intervals)\n",
                     cands.size(), vt.data.size());
        return 0;
    }

    for (double tau : taus) {
        // greedy portfolio over gamma values; keep best picks per k
        struct Best {
            int score = -1;
            double gamma = 0;
            std::vector<uint32_t> picks;
        };
        std::map<int, Best> best;
        std::vector<uint32_t> warmPicks;
        bool warm = !warmStart.empty() && loadWarmStart(warmStart, cands, warmPicks);
        if (warm) {
            std::fprintf(stderr, "warmstart: %zu antennas from %s\n",
                         warmPicks.size(), warmStart.c_str());
            bool matched = false;
            for (int k : ks)
                if ((int)warmPicks.size() == k) { best[k] = {0, 0, warmPicks}; matched = true; }
            if (!matched) {
                std::fprintf(stderr,
                             "FATAL: warmstart has %zu antennas, which matches no "
                             "requested k; pass --ks %zu\n",
                             warmPicks.size(), warmPicks.size());
                return 2;
            }
        }
        for (double g : gammas) {
            if (warm) break;  // the warm solution already beats a fresh greedy
            OptState st;
            st.init(ds, vt, tau, g);
            GreedyResult res;
            runGreedy(st, ks, nthreads, res);
            for (auto& [k, sc] : res.scores) {
                if (sc > best[k].score) {
                    best[k] = {sc, g,
                               std::vector<uint32_t>(res.picks.begin(),
                                                     res.picks.begin() + k)};
                }
            }
        }
        for (int k : ks) {
            Best& b = best[k];
            OptState st;
            st.init(ds, vt, tau, 0.4);  // small progress term as LS tiebreak
            for (uint32_t c : b.picks) st.add(c);
            int g0 = st.score;
            if (lsTime > 0)
                localSearch(st, lsTime, nthreads, seedBase * 1000003 + 12345 + k);
            int g1 = st.score;
            // With refinement on, the LNS budget is split: search once on the
            // coarse set, densify around what nearly made it, then search again
            // with the finer positions available.
            double lns1 = refine ? lnsTime * 0.5 : lnsTime;
            if (lns1 > 0)
                lnsSearch(st, lns1, nthreads, seedBase * 7919 + 999 + k, maxRuin,
                          cands);
            int g2 = st.score;
            if (refine) {
                std::vector<uint32_t> keep = st.selected;
                size_t grew = refineCandidates(ds, bg, grid, st, refineBand,
                                               refineRadius, refineSpacing, eps,
                                               R, minIvLen, nthreads, cands, vt);
                if (grew) {
                    st.init(ds, vt, tau, 0.4);  // isSel must cover the new size
                    for (uint32_t c : keep) st.add(c);
                }
                if (lnsTime - lns1 > 0)
                    lnsSearch(st, lnsTime - lns1, nthreads,
                              seedBase * 6151 + 777 + k, maxRuin, cands);
            }

            // Incremental add/remove bookkeeping drives the whole search; a
            // drift here would silently optimize a phantom objective, so
            // recompute the score from the final antenna set and compare.
            {
                int incremental = st.score;
                std::vector<uint32_t> finalSel = st.selected;
                st.setSelection(finalSel);
                if (st.score != incremental) {
                    std::fprintf(stderr,
                                 "FATAL: score drift tau=%g k=%d incremental=%d "
                                 "recomputed=%d\n", tau, k, incremental, st.score);
                    return 2;
                }
            }
            if (refine)
                std::fprintf(stderr,
                             "tau=%.2f k=%d: greedy %d (gamma=%.2f) -> ls %d "
                             "-> lns %d -> refined %d / %zu\n",
                             tau, k, g0, b.gamma, g1, g2, st.score, ds.blds.size());
            else
                std::fprintf(stderr,
                             "tau=%.2f k=%d: greedy %d (gamma=%.2f) -> ls %d -> lns %d / %zu\n",
                             tau, k, g0, b.gamma, g1, st.score, ds.blds.size());

            char path[512];
            std::snprintf(path, sizeof path, "%s/sol_t%g_k%d.txt", outDir.c_str(),
                          tau, k);
            FILE* f = std::fopen(path, "w");
            if (!f) { std::fprintf(stderr, "cannot write %s\n", path); return 1; }
            std::fprintf(f, "%g %d\n", tau, k);
            for (size_t i = 0; i < st.selected.size(); i++) {
                const Candidate& c = cands[st.selected[i]];
                std::fprintf(f, "%s%.17g %.17g", i ? ";" : "", c.p.x, c.p.y);
            }
            std::fprintf(f, "\n");
            auto sv = st.servedIdx();
            for (size_t i = 0; i < sv.size(); i++)
                std::fprintf(f, "%s%lld", i ? "," : "",
                             (long long)ds.blds[sv[i]].id);
            std::fprintf(f, "\n");
            std::fclose(f);

            if (witness) {
                // certificate: for each served building, the antenna intervals
                // (antenna index = position in the coordinate list) covering it
                std::snprintf(path, sizeof path, "%s/wit_t%g_k%d.txt",
                              outDir.c_str(), tau, k);
                FILE* wf = std::fopen(path, "w");
                std::map<uint32_t, size_t> selPos;
                for (size_t i = 0; i < st.selected.size(); i++)
                    selPos[st.selected[i]] = i;
                for (uint32_t b : sv) {
                    std::fprintf(wf, "%lld %.17g %.17g", (long long)ds.blds[b].id,
                                 ds.blds[b].perimeter, st.need[b]);
                    for (uint32_t a : st.touch[b]) {
                        auto [s, e] = st.colOn(a, b);
                        for (uint64_t j = s; j < e; j++)
                            std::fprintf(wf, " %zu:%.9g:%.9g", selPos[a],
                                         (double)vt.data[j].t0,
                                         (double)vt.data[j].t1);
                    }
                    std::fprintf(wf, "\n");
                }
                std::fclose(wf);
            }
        }
    }
    return 0;
}
