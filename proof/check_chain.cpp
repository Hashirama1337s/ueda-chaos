// Checker A (chain version): rigorous CAPD verification of covering relations for Ueda's map
//   P = time-2*pi map of  x'' + k x' + x^3 = B cos t,  k = 1/20, B = 15/2, state (x, v = x'), phase 0.
// Integrated as the autonomous system (tau = t / (2*pi), W = 2*pi as an interval)
//   x' = W v,  v' = W(-k v - x^3 + B c),  c' = -W s,  s' = W c,   (c, s)(0) = (1, 0) exactly,
// so P = (x, v)-part of the exact time-1 map (c, s return to (1, 0)).
//
// Charts: "curved"  X(u,s) = Q + (u + p(s)) U + s S   (p: fixed polynomial, double coefficients)
//         "affine"  X(u,s) = Q + u U + s S
// Inverse (rigorous): (a, s) = [U S]^(-1) (X - Q) in interval arithmetic; u = a - p(s) (curved) or a (affine).
// h-set N = [c - w, c + w] x [-b, b] in its chart; exit direction u, entry direction s.
// For every transition N => M (map P, orientation o = +1 or -1) the program proves, with rigorous enclosures:
//   (E-) P(u-edge of N at c - w) lies in M-chart {u < cM - wM} (o = +1) or {u > cM + wM} (o = -1),
//   (E+) P(u-edge of N at c + w) lies on the opposite side,
//   (S)  P(N) lies in M-chart {|s| < bM}.
// These imply the covering relation N =P=> M for h-sets with one exit direction (Zgliczynski & Gidea, J. Differential
// Equations 202 (2004), 32-58): the straight-line homotopy (in M's chart) from P to the linear map (u,s) -> (o*L*u, 0),
// L large, never moves exit-edge images into M (their u stays on the same side) and never meets M's entry set, since
// (S) holds for the whole image and the homotopy only shrinks |s|.
// The h-sets are also proved pairwise disjoint. With the design's transition graph this gives a compact P-invariant set
// semiconjugate to the subshift of finite type of the graph, hence h_top(P) >= log(spectral radius) > 0.
// usage: check_chain <design.txt> [initial pieces along s] [max depth]
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <memory>
#include "capd/capdlib.h"
using namespace capd; using namespace std;

struct Chart { string name; bool curved; double Q[2], U[2], S[2]; IVector Q4; IMatrix Minv4; };
struct HSet { string name; double c, w, b; int chart; };
struct Trans { int src, dst, o; };
vector<double> pc; vector<Chart> charts; vector<HSet> sets; vector<Trans> trans;

interval P_(interval s) { interval v(0.); for (double q : pc) v = v * s + interval(q); return v; }
interval dP_(interval s) { interval v(0.); int n = pc.size(); for (int k = 0; k < n - 1; ++k) v = v * s + interval(pc[k]) * interval(n - 1 - k); return v; }
interval ddP_(interval s) { interval v(0.); int n = pc.size(); for (int k = 0; k < n - 2; ++k) v = v * s + interval(pc[k]) * interval(n - 1 - k) * interval(n - 2 - k); return v; }
double Pd(double s) { double v = 0; for (double q : pc) v = v * s + q; return v; }

int findSet(const string& n) { for (size_t i = 0; i < sets.size(); ++i) if (sets[i].name == n) return i; cerr << "no set " << n << endl; exit(2); }

void readDesign(const char* fn) {
  ifstream in(fn); string tok;
  if (!in) { cerr << "FATAL: cannot open design file " << fn << endl; exit(2); }
  while (in >> tok) {
    if (tok[0] == '#') { getline(in, tok); continue; }
    if (tok == "m") { int m; in >> m; if (m != 1) { cerr << "only m = 1 supported" << endl; exit(2); } }
    else if (tok == "p") { int n; in >> n; pc.resize(n); for (auto& q : pc) in >> q; }
    else if (tok == "chart") { Chart c; string ty; in >> c.name >> ty >> c.Q[0] >> c.Q[1] >> c.U[0] >> c.U[1] >> c.S[0] >> c.S[1]; c.curved = ty == "curved"; charts.push_back(c); }
    else if (tok == "set") { HSet h; in >> h.name >> h.c >> h.w >> h.b; h.chart = -1;
      for (size_t i = 0; i < charts.size(); ++i) if (charts[i].name == h.name) h.chart = i;
      if (h.chart < 0) { cerr << "no chart for set " << h.name << endl; exit(2); } sets.push_back(h); }
    else if (tok == "trans") { string a, b; int o; in >> a >> b >> o; if (o != 1 && o != -1) { cerr << "FATAL: orientation must be +1 or -1" << endl; exit(2); } trans.push_back({findSet(a), findSet(b), o}); }
    else { cerr << "FATAL: unknown token '" << tok << "' in design file" << endl; exit(2); }
    if (in.fail()) { cerr << "FATAL: parse error after token '" << tok << "'" << endl; exit(2); }
  }
  if (!in.eof()) { cerr << "FATAL: design file not fully read" << endl; exit(2); }
  if (sets.empty() || trans.empty()) { cerr << "FATAL: design has no h-sets or no transitions" << endl; exit(2); }
  for (size_t i = 0; i < sets.size(); ++i) { bool out = false; for (auto& t : trans) out = out || t.src == (int)i;
    if (!out) { cerr << "FATAL: h-set " << sets[i].name << " has no outgoing transition" << endl; exit(2); } }
  for (auto& c : charts) {
    IMatrix E(2, 2); E[0][0] = interval(c.U[0]); E[1][0] = interval(c.U[1]); E[0][1] = interval(c.S[0]); E[1][1] = interval(c.S[1]);
    interval det = E[0][0] * E[1][1] - E[0][1] * E[1][0];
    c.Minv4 = IMatrix::Identity(4);
    c.Minv4[0][0] = E[1][1] / det; c.Minv4[0][1] = -E[0][1] / det; c.Minv4[1][0] = -E[1][0] / det; c.Minv4[1][1] = E[0][0] / det;
    c.Q4 = IVector({interval(c.Q[0]), interval(c.Q[1]), interval(1.), interval(0.)});
  }
}

IMap* field; IOdeSolver* solver; ITimeMap* tmap;
DMap* dfield; DOdeSolver* dsolver; DTimeMap* dtmap;

// chart coordinates of a 4-D enclosure given as an affine set
pair<interval, interval> toChart(C0TripletonSet& set, const Chart& ch) {
  IVector as = set.affineTransformation(ch.Minv4, ch.Q4);
  interval s = as[1], u = ch.curved ? as[0] - P_(as[1]) : as[0];
  return {u, s};
}

// rigorous image of the chart box (u in Uu, s in Ss) of chart ch under P
C0TripletonSet* imageSet(const Chart& ch, interval Uu, interval Ss) {  // caller owns the result
  interval s0(Ss.mid().leftBound()), u0(Uu.mid().leftBound()), dS = Ss - s0, dU = Uu - u0;
  interval pa(0.), pd(0.), R(0.);
  if (ch.curved) { pa = P_(s0); pd = dP_(s0); R = interval(0.5) * ddP_(Ss) * sqr(dS); }   // Lagrange remainder on the hull Ss
  interval a0 = u0 + pa;
  IVector X0({interval(ch.Q[0]) + a0 * interval(ch.U[0]) + s0 * interval(ch.S[0]),
              interval(ch.Q[1]) + a0 * interval(ch.U[1]) + s0 * interval(ch.S[1]), interval(1.), interval(0.)});
  IMatrix C = IMatrix::Identity(4);
  for (int i = 0; i < 2; ++i) { C[i][0] = interval(ch.U[i]); C[i][1] = pd * interval(ch.U[i]) + interval(ch.S[i]); }
  IVector r({dU + R, dS, interval(0.), interval(0.)});
  C0TripletonSet* set = new C0TripletonSet(X0, C, r);
  (*tmap)(interval(1.), *set);
  return set;
}

// non-rigorous point image in a target chart (used only to detect a definite failure early)
pair<double, double> dimage(const Chart& ch, double u, double s, const Chart& tg) {
  double a = u + (ch.curved ? Pd(s) : 0.);
  DVector X({ch.Q[0] + a * ch.U[0] + s * ch.S[0], ch.Q[1] + a * ch.U[1] + s * ch.S[1]});
  DVector Y = (*dtmap)(1.0, X);
  double det = tg.U[0] * tg.S[1] - tg.U[1] * tg.S[0], dx = Y[0] - tg.Q[0], dv = Y[1] - tg.Q[1];
  double aa = (dx * tg.S[1] - dv * tg.S[0]) / det, ss = (-dx * tg.U[1] + dv * tg.U[0]) / det;
  return {aa - (tg.curved ? Pd(ss) : 0.), ss};
}

enum Kind { LOW, HIGH, SBAND };
struct Cond { Kind k; double val; int tgt; };
long nPieces = 0, nFailInt = 0; int maxDepth = 40; bool fatal = false;

bool decided(const pair<interval, interval>& us, const Cond& c) {
  if (c.k == LOW) return us.first < interval(c.val);
  if (c.k == HIGH) return us.first > interval(c.val);
  return us.second > interval(-c.val) && us.second < interval(c.val);
}
bool pointOK(const Chart& ch, double u, double s, const Cond& c) {
  auto q = dimage(ch, u, s, charts[sets[c.tgt].chart]);
  if (c.k == LOW) return q.first < c.val; if (c.k == HIGH) return q.first > c.val; return fabs(q.second) < c.val;
}

bool verify(const Chart& ch, interval Uu, interval Ss, const vector<Cond>& conds, bool splitU, int depth) {
  if (fatal) return false;
  vector<Cond> rest; bool ok = true; double uwidth = 0;
  try {
    std::unique_ptr<C0TripletonSet> set(imageSet(ch, Uu, Ss)); ++nPieces;
    for (auto& c : conds) {
      auto us = toChart(*set, charts[sets[c.tgt].chart]);
      uwidth = std::max(uwidth, (double)(us.first.rightBound() - us.first.leftBound()));
      if (!decided(us, c)) rest.push_back(c);
    }
  } catch (exception& e) { ok = false; ++nFailInt; rest = conds; }
  if (rest.empty()) return true;
  for (auto& c : rest) if (!pointOK(ch, Uu.mid().leftBound(), Ss.mid().leftBound(), c)) {
    cerr << "DEFINITE FAIL: chart " << ch.name << " u=" << Uu.mid() << " s=" << Ss.mid() << " kind " << c.k << " val " << c.val
         << " target " << sets[c.tgt].name << endl;
    fatal = true; return false;
  }
  if (depth >= maxDepth) { cerr << "max depth: chart " << ch.name << " U=" << Uu << " S=" << Ss << endl; fatal = true; return false; }
  bool cutU = splitU && (depth % 2 == 1) && (!ok || uwidth > 1e-3);
  if (cutU) { double m = Uu.mid().leftBound();
    return verify(ch, interval(Uu.leftBound(), m), Ss, rest, splitU, depth + 1) && verify(ch, interval(m, Uu.rightBound()), Ss, rest, splitU, depth + 1); }
  double m = Ss.mid().leftBound();
  return verify(ch, Uu, interval(Ss.leftBound(), m), rest, splitU, depth + 1) && verify(ch, Uu, interval(m, Ss.rightBound()), rest, splitU, depth + 1);
}

// rigorous test that the (x, v)-box B misses h-set h (B mapped into h's chart has u or s range disjoint from h's)
bool missesSet(const IVector& B, const HSet& h) {
  const Chart& ch = charts[h.chart];
  interval dx = B[0] - interval(ch.Q[0]), dv = B[1] - interval(ch.Q[1]);
  interval a = ch.Minv4[0][0] * dx + ch.Minv4[0][1] * dv, s = ch.Minv4[1][0] * dx + ch.Minv4[1][1] * dv;
  interval u = ch.curved ? a - P_(s) : a;
  return u.rightBound() < h.c - h.w || u.leftBound() > h.c + h.w || s.rightBound() < -h.b || s.leftBound() > h.b;
}
bool disjoint(const HSet& a, const HSet& b, int depth = 0) {
  // subdivide a until each piece's (x, v) hull misses b
  vector<pair<interval, interval>> st{{interval(a.c - a.w, a.c + a.w), interval(-a.b, a.b)}};
  const Chart& ch = charts[a.chart];
  int guard = 0;
  while (!st.empty()) {
    auto pr = st.back(); st.pop_back();
    interval A = pr.first + (ch.curved ? P_(pr.second) : interval(0.));
    IVector B({interval(ch.Q[0]) + A * interval(ch.U[0]) + pr.second * interval(ch.S[0]), interval(ch.Q[1]) + A * interval(ch.U[1]) + pr.second * interval(ch.S[1])});
    if (missesSet(B, b)) continue;
    if (++guard > 200000) return false;
    double m = pr.second.mid().leftBound(), mu = pr.first.mid().leftBound();
    if (pr.second.rightBound() - pr.second.leftBound() > pr.first.rightBound() - pr.first.leftBound()) {
      st.push_back({pr.first, interval(pr.second.leftBound(), m)}); st.push_back({pr.first, interval(m, pr.second.rightBound())}); }
    else { st.push_back({interval(pr.first.leftBound(), mu), pr.second}); st.push_back({interval(mu, pr.first.rightBound()), pr.second}); }
  }
  return true;
}

int main(int argc, char** argv) {
  cout.precision(10);
  if (argc < 2) { cerr << "usage: check_chain design.txt [pieces] [maxdepth]" << endl; return 2; }
  readDesign(argv[1]);
  int n0 = argc > 2 ? atoi(argv[2]) : 64; if (argc > 3) maxDepth = atoi(argv[3]);
  if (n0 < 1 || maxDepth < 1) { cerr << "FATAL: need pieces >= 1 and max depth >= 1" << endl; return 2; }
  field = new IMap("par:W,K,B;var:x,v,c,s;fun:W*v,W*(-K*v-x^3+B*c),-W*s,W*c;");
  field->setParameter("W", 2 * interval::pi()); field->setParameter("K", interval(1) / interval(20)); field->setParameter("B", interval(15) / interval(2));
#ifdef UEDA_NEGATIVE_CONTROL_BUILD
  // negative-control binary only (compiled with -DUEDA_NEGATIVE_CONTROL_BUILD); never prints a verification
  field->setParameter("B", interval(atof(getenv("UEDA_NEGCTL_B") ? getenv("UEDA_NEGCTL_B") : "7.4")));
  cout << "NEGATIVE CONTROL BUILD: B = " << (getenv("UEDA_NEGCTL_B") ? getenv("UEDA_NEGCTL_B") : "7.4") << endl;
#endif
  int ord = getenv("UEDA_ORDER") ? atoi(getenv("UEDA_ORDER")) : 12;
  solver = new IOdeSolver(*field, ord); tmap = new ITimeMap(*solver);
  dfield = new DMap("time:t;par:W,K,B;var:x,v;fun:W*v,W*(-K*v-x^3+B*cos(W*t));");
  dfield->setParameter("W", 2 * M_PI); dfield->setParameter("K", 0.05); dfield->setParameter("B", 7.5);
  dsolver = new DOdeSolver(*dfield, 20); dtmap = new DTimeMap(*dsolver);
  cout << "Ueda map P (k = 1/20, B = 15/2), C0 Lohner tripleton, Taylor order " << ord << "; " << sets.size() << " h-sets, "
       << trans.size() << " transitions" << endl;
  for (auto& c : charts) cout << "chart " << c.name << (c.curved ? " curved" : " affine") << " Q=(" << c.Q[0] << ", " << c.Q[1] << ") U=(" << c.U[0] << ", " << c.U[1] << ") S=(" << c.S[0] << ", " << c.S[1] << ")" << endl;
  for (auto& h : sets) cout << "h-set " << h.name << ": |u - " << h.c << "| <= " << h.w << ", |s| <= " << h.b << endl;
  for (auto& t : trans) cout << "transition " << sets[t.src].name << " => " << sets[t.dst].name << " (orientation " << t.o << ")" << endl;
  bool all = true;
  for (size_t i = 0; i < sets.size(); ++i) for (size_t j = i + 1; j < sets.size(); ++j) {
    bool d = disjoint(sets[i], sets[j]) || disjoint(sets[j], sets[i]);
    cout << "disjoint " << sets[i].name << " / " << sets[j].name << ": " << (d ? "VERIFIED" : "NOT VERIFIED") << endl;
    all = all && d;
  }
  for (size_t si = 0; si < sets.size(); ++si) {
    const HSet& h = sets[si]; const Chart& ch = charts[h.chart];
    vector<Cond> lo, hi, whole;
    for (auto& t : trans) if (t.src == (int)si) {
      const HSet& M = sets[t.dst];
      Cond L{LOW, M.c - M.w, t.dst}, H{HIGH, M.c + M.w, t.dst};
      if (t.o > 0) { lo.push_back(L); hi.push_back(H); } else { lo.push_back(H); hi.push_back(L); }
      whole.push_back(Cond{SBAND, M.b, t.dst});
    }
    if (whole.empty()) continue;
    const char* nm[3] = {"edge u=c-w", "edge u=c+w", "whole set (S)"};
    for (int part = 0; part < 3 && all; ++part) {
      long before = nPieces, fb = nFailInt; bool ok = true;
      for (int k = 0; k < n0 && ok; ++k) {
        double sLo = (k == 0) ? -h.b : -h.b + 2 * h.b * k / n0, sHi = (k == n0 - 1) ? h.b : -h.b + 2 * h.b * (k + 1) / n0;
        interval Sk(sLo, sHi);  // consecutive pieces share endpoints; first/last are exactly -b and b
        if (part == 0) ok = verify(ch, interval(h.c - h.w), Sk, lo, false, 0);
        else if (part == 1) ok = verify(ch, interval(h.c + h.w), Sk, hi, false, 0);
        else ok = verify(ch, interval(h.c - h.w, h.c + h.w), Sk, whole, true, 0);
      }
      cout << "set " << h.name << " " << nm[part] << ": " << (ok ? "VERIFIED" : "FAILED") << " (" << nPieces - before
           << " enclosures, " << nFailInt - fb << " integration failures resolved by subdivision)" << endl;
      if (nPieces - before < 1) { cout << "FAIL: no enclosure computed for this part" << endl; ok = false; }
      all = all && ok;
    }
  }
  cout << (all ? "ALL COVERING RELATIONS AND DISJOINTNESS VERIFIED" : "NOT VERIFIED") << endl;
  return all ? 0 : 1;
}
