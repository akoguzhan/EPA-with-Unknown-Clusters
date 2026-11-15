// Walk of one Panel Kmeans initialization along X(phi) = U + phi V (see R/partition_walk.R).
// Returns the segments [lo, hi] on which the Lloyd path is constant, with the final partition
// and validity of each segment (adjacent segments with the same final partition are merged).
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <map>
#include <algorithm>
using namespace Rcpp;

static void lloyd_continue(const std::vector<double>& x, int N, int P, int K,
                           std::vector<int> g, int nmax,
                           std::vector< std::vector<int> >& out, bool& valid) {
  valid = true;
  std::vector<double> ctr(K * P);
  std::vector<int> cnt(K);
  for (int m = 0; m < nmax; ++m) {
    std::fill(ctr.begin(), ctr.end(), 0.0);
    std::fill(cnt.begin(), cnt.end(), 0);
    for (int i = 0; i < N; ++i) {
      int k = g[i];
      if (k <= 0) continue;
      cnt[k - 1]++;
      for (int p = 0; p < P; ++p) ctr[(k - 1) * P + p] += x[i * P + p];
    }
    for (int k = 0; k < K; ++k) for (int p = 0; p < P; ++p) ctr[k * P + p] /= cnt[k];
    std::vector<int> gnew(N);
    std::vector<int> cnew(K, 0);
    for (int i = 0; i < N; ++i) {
      double best = std::numeric_limits<double>::infinity(); int bk = 1;
      for (int k = 0; k < K; ++k) {
        double d = 0.0;
        for (int p = 0; p < P; ++p) { double z = x[i * P + p] - ctr[k * P + p]; d += z * z; }
        if (d < best) { best = d; bk = k + 1; }
      }
      gnew[i] = bk; cnew[bk - 1]++;
    }
    out.push_back(gnew);
    bool empty = false;
    for (int k = 0; k < K; ++k) if (cnew[k] == 0) empty = true;
    if (empty) { valid = false; return; }
    bool same = true;
    for (int i = 0; i < N; ++i) if (gnew[i] != g[i]) { same = false; break; }
    if (same) return;
    g = gnew;
  }
}

// positive roots of the assignment constraints of one step
static std::vector<double> step_roots(const NumericMatrix& U, const NumericMatrix& V, int K,
                                      const std::vector<int>& gprev, const std::vector<int>& anew,
                                      double P1, double eps) {
  int N = U.nrow(), P = U.ncol();
  std::vector<double> Ub(K * P, 0.0), Vb(K * P, 0.0);
  std::vector<int> cnt(K, 0);
  for (int i = 0; i < N; ++i) {
    int k = gprev[i];
    if (k <= 0) continue;
    cnt[k - 1]++;
    for (int p = 0; p < P; ++p) { Ub[(k - 1) * P + p] += U(i, p); Vb[(k - 1) * P + p] += V(i, p); }
  }
  for (int k = 0; k < K; ++k) for (int p = 0; p < P; ++p) { Ub[k * P + p] /= cnt[k]; Vb[k * P + p] /= cnt[k]; }
  std::vector<double> cb(K, 0.0), ca(K, 0.0), cab(K, 0.0);
  for (int k = 0; k < K; ++k) for (int p = 0; p < P; ++p) {
    cb[k] += Vb[k * P + p] * Vb[k * P + p];
    ca[k] += Ub[k * P + p] * Ub[k * P + p];
    cab[k] += Ub[k * P + p] * Vb[k * P + p];
  }
  std::vector<double> r;
  std::vector<double> fA(K), fB(K), fC(K);
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      double vv = 0, uu = 0, uv = 0, vu = 0;
      for (int p = 0; p < P; ++p) {
        vv += V(i, p) * Vb[k * P + p]; uu += U(i, p) * Ub[k * P + p];
        uv += U(i, p) * Vb[k * P + p]; vu += V(i, p) * Ub[k * P + p];
      }
      fA[k] = cb[k] - 2 * vv; fC[k] = ca[k] - 2 * uu; fB[k] = cab[k] - uv - vu;
    }
    int a = anew[i] - 1;
    for (int g = 0; g < K; ++g) {
      if (g == a) continue;
      double A = fA[a] - fA[g], B = 2 * (fB[a] - fB[g]), C = fC[a] - fC[g];
      double sc = std::fabs(A) * P1 * P1 + std::fabs(B) * P1 + std::fabs(C);
      if (sc == 0) sc = 1;
      C -= eps * sc;
      if (std::fabs(A) * P1 * P1 <= 1e-12 * sc) {
        if (std::fabs(B) * P1 > 1e-12 * sc) { double z = -C / B; if (z > 0) r.push_back(z); }
      } else {
        double disc = B * B - 4 * A * C;
        if (disc >= 0) {
          double sq = std::sqrt(disc);
          double qq = -0.5 * (B + (B >= 0 ? 1.0 : -1.0) * sq);
          double z1 = qq / A, z2 = (qq == 0) ? 0.0 : C / qq;
          if (z1 > 0) r.push_back(z1);
          if (z2 > 0) r.push_back(z2);
        }
      }
    }
  }
  return r;
}

struct Walk {
  std::vector<double> lo, hi;
  std::vector<int> val;
  std::vector< std::vector<int> > parts;
};

// Segments of one initialization over [phimin, phimax]
static Walk walk_one(const NumericMatrix& U, const NumericMatrix& V, int K, const std::vector<int>& init,
                     int itermax, double phimin, double phimax, double eps, double phiref) {
  int N = U.nrow(), P = U.ncol();
  double P1 = std::max(phiref, 1.0);   // scale of the numerical tolerance
  std::vector< std::vector<int> > path;
  std::vector< std::vector<double> > roots;
  bool valid = true;
  Walk w;
  double phi = phimin;
  int mstar = 1, nstep = 0;
  std::vector<double> x(N * P);
  while (phi < phimax && nstep < 100000) {
    ++nstep;
    for (int i = 0; i < N; ++i) for (int p = 0; p < P; ++p) x[i * P + p] = U(i, p) + phi * V(i, p);
    std::vector<int> gstart = (mstar == 1) ? init : path[mstar - 2];
    path.resize(mstar - 1);
    roots.resize(mstar - 1);
    std::vector< std::vector<int> > cont;
    lloyd_continue(x, N, P, K, gstart, itermax - (mstar - 1), cont, valid);
    for (size_t m = 0; m < cont.size(); ++m) path.push_back(cont[m]);
    for (size_t m = mstar - 1; m < path.size(); ++m) {
      const std::vector<int>& prev = (m == 0) ? init : path[m - 1];
      roots.push_back(step_roots(U, V, K, prev, path[m], P1, eps));
    }
    double tl = phi + 1e-10 * std::max(1.0, phi);
    double nb = std::numeric_limits<double>::infinity(); int ms = 1;
    for (size_t m = 0; m < roots.size(); ++m) {
      for (size_t j = 0; j < roots[m].size(); ++j) {
        double z = roots[m][j];
        if (z > tl && z < nb) { nb = z; ms = (int)m + 1; }
      }
    }
    mstar = ms;
    const std::vector<int>& fin = path.back();
    int n = w.lo.size();
    bool merge = false;
    if (n > 0 && w.val[n - 1] == (valid ? 1 : 0)) {
      merge = true;
      if (valid) { for (int i = 0; i < N; ++i) if (w.parts[n - 1][i] != fin[i]) { merge = false; break; } }
    }
    if (merge) w.hi[n - 1] = nb;
    else { w.lo.push_back(phi); w.hi.push_back(nb); w.val.push_back(valid ? 1 : 0); w.parts.push_back(fin); }
    if (!std::isfinite(nb)) break;
    phi = nb + 1e-8 * std::max(1.0, nb);
  }
  return w;
}

// [[Rcpp::export]]
List run_segments_cpp(NumericMatrix U, NumericMatrix V, int K, IntegerVector g0, int itermax,
                      double phimin, double phimax, double eps) {
  int N = U.nrow();
  std::vector<int> init(g0.begin(), g0.end());
  Walk w = walk_one(U, V, K, init, itermax, phimin, phimax, eps, phimax);
  int ns = w.lo.size();
  IntegerMatrix pm(N, ns);
  for (int s = 0; s < ns; ++s) for (int i = 0; i < N; ++i) pm(i, s) = w.parts[s][i];
  return List::create(Named("lo") = wrap(w.lo), Named("hi") = wrap(w.hi), Named("valid") = wrap(w.val),
                      Named("parts") = pm);
}

// ---- retained partition of one K, piecewise in phi (steps 1-2 of .good_block) ----------------

// relabel by order of first appearance (as .partition_key)
static std::vector<int> canonical(const std::vector<int>& g) {
  std::vector<int> map, out(g.size());
  std::vector<int> seen;
  for (size_t i = 0; i < g.size(); ++i) {
    int lab = -1;
    for (size_t j = 0; j < seen.size(); ++j) if (seen[j] == g[i]) { lab = (int)j + 1; break; }
    if (lab < 0) { seen.push_back(g[i]); lab = seen.size(); }
    out[i] = lab;
  }
  return out;
}

// WSS(phi; g) = c1 phi^2 + c2 phi + c3
static void wss_coefs(const NumericMatrix& U, const NumericMatrix& V, const std::vector<int>& g, int K,
                      double* c) {
  int N = U.nrow(), P = U.ncol();
  std::vector<double> Ub(K * P, 0.0), Vb(K * P, 0.0);
  std::vector<int> cnt(K, 0);
  for (int i = 0; i < N; ++i) {
    int k = g[i] - 1; cnt[k]++;
    for (int p = 0; p < P; ++p) { Ub[k * P + p] += U(i, p); Vb[k * P + p] += V(i, p); }
  }
  for (int k = 0; k < K; ++k) for (int p = 0; p < P; ++p) { Ub[k * P + p] /= cnt[k]; Vb[k * P + p] /= cnt[k]; }
  long double s1 = 0, s2 = 0, s3 = 0;
  for (int p = 0; p < P; ++p) for (int i = 0; i < N; ++i) {
    int k = g[i] - 1;
    double ru = U(i, p) - Ub[k * P + p], rv = V(i, p) - Vb[k * P + p];
    s1 += (long double)rv * rv; s2 += (long double)ru * rv; s3 += (long double)ru * ru;
  }
  c[0] = (double)s1; c[1] = 2.0 * (double)s2; c[2] = (double)s3;
}

// positive roots of A x^2 + B x + C (as .quad_roots with eps = 0)
static void quad_roots0(double A, double B, double C, double phi_ref, std::vector<double>& r) {
  double P1 = std::max(phi_ref, 1.0);
  double sc = std::fabs(A) * P1 * P1 + std::fabs(B) * P1 + std::fabs(C);
  if (sc == 0) sc = 1;
  if (std::fabs(A) * P1 * P1 <= 1e-12 * sc) {
    if (std::fabs(B) * P1 > 1e-12 * sc) { double z = -C / B; if (std::isfinite(z) && z > 0) r.push_back(z); }
  } else {
    double disc = B * B - 4 * A * C;
    if (disc >= 0) {
      double sq = std::sqrt(disc);
      double qq = -0.5 * (B + (B >= 0 ? 1.0 : -1.0) * sq);
      double z1 = qq / A, z2 = (qq == 0) ? 0.0 : C / qq;
      if (std::isfinite(z1) && z1 > 0) r.push_back(z1);
      if (std::isfinite(z2) && z2 > 0) r.push_back(z2);
    }
  }
}

// [[Rcpp::export]]
List winners_cpp(NumericMatrix U, NumericMatrix V, int K, List inits, int nmain, int itermax,
                 double lo_b, double hi_b, double eps, double phiref) {
  int N = U.nrow();
  int S = inits.size();
  std::map< std::vector<int>, int > pid_of;
  std::vector< std::vector<int> > plist;
  std::vector< std::vector<double> > tlo(S), thi(S);
  std::vector< std::vector<int> > tid(S);
  int nseg = 0;
  for (int s = 0; s < S; ++s) {
    IntegerVector gi = inits[s];
    std::vector<int> init(gi.begin(), gi.end());
    Walk w = walk_one(U, V, K, init, itermax, lo_b, hi_b, eps, phiref);
    for (size_t j = 0; j < w.lo.size(); ++j) {
      int id = -1;
      if (w.val[j] == 1) {
        std::vector<int> c = canonical(w.parts[j]);
        std::map< std::vector<int>, int >::iterator it = pid_of.find(c);
        if (it == pid_of.end()) { id = plist.size(); pid_of[c] = id; plist.push_back(c); }
        else id = it->second;
      }
      size_t n = tlo[s].size();
      if (n > 0 && tid[s][n - 1] == id) thi[s][n - 1] = std::max(thi[s][n - 1], w.hi[j]);
      else { tlo[s].push_back(w.lo[j]); thi[s].push_back(w.hi[j]); tid[s].push_back(id); }
    }
    nseg += tlo[s].size();
  }
  // elementary pieces
  std::vector<double> br;
  br.push_back(lo_b); br.push_back(hi_b);
  for (int s = 0; s < S; ++s) { br.insert(br.end(), tlo[s].begin(), tlo[s].end()); br.insert(br.end(), thi[s].begin(), thi[s].end()); }
  std::sort(br.begin(), br.end());
  br.erase(std::unique(br.begin(), br.end()), br.end());
  std::vector<double> b2;
  for (size_t j = 0; j < br.size(); ++j) if (br[j] >= lo_b && br[j] <= hi_b) b2.push_back(br[j]);
  br.swap(b2);
  int nm = (int)br.size() - 1;
  // candidate sets on each piece (-2 = not covered, -1 = invalid)
  std::vector< std::vector<int> > cands(std::max(nm, 0));
  for (int m = 0; m < nm; ++m) {
    double mid = (br[m + 1] + br[m]) / 2;
    std::vector<int> ks(S);
    for (int s = 0; s < S; ++s) {
      const std::vector<double>& lo = tlo[s];
      int j = std::upper_bound(lo.begin(), lo.end(), mid) - lo.begin();   // number of lo <= mid
      int k = -2;
      if (j > 0 && !(mid > thi[s][j - 1] + 1e-7 * std::max(1.0, mid))) k = tid[s][j - 1];
      ks[s] = k;
    }
    std::vector<int> c;
    for (int s = 0; s < std::min(nmain, S); ++s) if (ks[s] >= 0 && std::find(c.begin(), c.end(), ks[s]) == c.end()) c.push_back(ks[s]);
    if (c.empty() && S > nmain) for (int s = nmain; s < S; ++s) if (ks[s] >= 0) { c.push_back(ks[s]); break; }
    cands[m] = c;
  }
  std::vector<double> olo, ohi;
  std::vector<int> okey;
  std::map<int, std::vector<double> > wcf;
  int m = 0;
  while (m < nm) {
    int e = m;
    while (e + 1 < nm && cands[e + 1] == cands[m]) ++e;
    double a = br[m], b = br[e + 1];
    const std::vector<int>& uk = cands[m];
    m = e + 1;
    if (b <= a) continue;
    if (uk.empty()) { olo.push_back(a); ohi.push_back(b); okey.push_back(-1); continue; }
    if (uk.size() == 1) { olo.push_back(a); ohi.push_back(b); okey.push_back(uk[0]); continue; }
    int nk = uk.size();
    std::vector<double> cf(3 * nk);
    for (int k = 0; k < nk; ++k) {
      if (!wcf.count(uk[k])) { std::vector<double> c(3); wss_coefs(U, V, plist[uk[k]], K, &c[0]); wcf[uk[k]] = c; }
      const std::vector<double>& c = wcf[uk[k]];
      cf[3 * k] = c[0]; cf[3 * k + 1] = c[1]; cf[3 * k + 2] = c[2];
    }
    std::vector<double> rr;
    for (int i = 0; i < nk - 1; ++i) for (int j = i + 1; j < nk; ++j)
      quad_roots0(cf[3 * i] - cf[3 * j], cf[3 * i + 1] - cf[3 * j + 1], cf[3 * i + 2] - cf[3 * j + 2], std::max(b, 1.0), rr);
    std::vector<double> cuts;
    cuts.push_back(a); cuts.push_back(b);
    for (size_t j = 0; j < rr.size(); ++j) if (rr[j] > a && rr[j] < b) cuts.push_back(rr[j]);
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    for (size_t j = 0; j + 1 < cuts.size(); ++j) {
      double mm = (cuts[j + 1] + cuts[j]) / 2;
      int win = 0; double best = std::numeric_limits<double>::infinity();
      for (int k = 0; k < nk; ++k) {
        double v = mm * mm * cf[3 * k] + mm * cf[3 * k + 1] + cf[3 * k + 2];
        if (v < best) { best = v; win = k; }
      }
      olo.push_back(cuts[j]); ohi.push_back(cuts[j + 1]); okey.push_back(uk[win]);
    }
  }
  int np = plist.size();
  IntegerMatrix pm(N, np);
  for (int s = 0; s < np; ++s) for (int i = 0; i < N; ++i) pm(i, s) = plist[s][i];
  return List::create(Named("lo") = wrap(olo), Named("hi") = wrap(ohi), Named("key") = wrap(okey),
                      Named("parts") = pm, Named("nseg") = nseg);
}
