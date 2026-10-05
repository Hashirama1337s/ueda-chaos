// Independent verifier B for a chain of h-set covering relations of Ueda's stroboscopic map.
// Interval arithmetic and the C0 Lohner integrator are CAPD/filib. No other verifier is used.
//
// Rigour notes baked into the checks:
// - Decimal data are enclosed by outward string conversion, then treated as those exact decimals.
// - k = 1/20 and the nominal B = 15/2 are exact rationals. 2*pi is an interval enclosure.
// - A box proves a strict inequality only when its outward image bound lies strictly on the
//   inner side of the threshold interval. That implies the inequality for every point of the box.
// - An enclosure lying entirely on the wrong side refutes the universal statement.
// - Parameter boxes are closed interval products whose union contains the h-set (or its edge).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "capd/capdlib.h"

using namespace std;
using namespace capd;

static const int kTaylorOrder = 12;
static const double kInitialParamStep = 0.01;
static const long kMaxIntegrationsPerGoal = 100000;

static interval enclose_decimal(const string& s) {
  return interval(s, s);
}

static double mid_inside(const interval& z) {
  double a = z.leftBound();
  double b = z.rightBound();
  double m = a + (b - a) * 0.5;
  if (m < a) m = a;
  if (m > b) m = b;
  return m;
}

static interval force_contains_zero(interval r) {
  if (r.leftBound() > 0.0) r = interval(0.0, r.rightBound());
  if (r.rightBound() < 0.0) r = interval(r.leftBound(), 0.0);
  return r;
}

static bool split_interval(const interval& z, interval& L, interval& R) {
  double a = z.leftBound();
  double b = z.rightBound();
  if (!(b > a)) return false;
  double m = a + (b - a) * 0.5;
  if (!(m > a && m < b)) {
    m = nextafter(a, b);
    if (!(m > a && m < b)) return false;
  }
  L = interval(a, m);
  R = interval(m, b);
  return true;
}

static interval hull2(const interval& a, const interval& b) {
  return interval(std::min(a.leftBound(), b.leftBound()), std::max(a.rightBound(), b.rightBound()));
}

static bool separated_intervals(const interval& a, const interval& b) {
  return a.rightBound() < b.leftBound() || b.rightBound() < a.leftBound();
}

struct Poly {
  vector<interval> c, d1, d2;
  explicit Poly(vector<interval> coeffs) : c(std::move(coeffs)) {
    int deg = (int)c.size() - 1;
    for (int i = 0; i < deg; ++i) d1.push_back(interval(double(deg - i)) * c[i]);
    int deg1 = (int)d1.size() - 1;
    for (int i = 0; i < deg1; ++i) d2.push_back(interval(double(deg1 - i)) * d1[i]);
  }
  interval eval(const vector<interval>& coeff, interval s) const {
    interval acc(0.0);
    for (size_t i = 0; i < coeff.size(); ++i) acc = acc * s + coeff[i];
    return acc;
  }
  interval p(interval s) const { return eval(c, s); }
  interval dp(interval s) const { return eval(d1, s); }
  interval ddp(interval s) const { return eval(d2, s); }

  // Outward enclosure of p on s, tightened by chopping so Horner dependency stays small.
  interval range(interval s) const {
    double a = s.leftBound();
    double b = s.rightBound();
    double width = b - a;
    if (!(width > 0.0)) return p(s);
    int n = (int)ceil(width / 0.01);
    if (n < 1) n = 1;
    if (n > 10000) n = 10000;
    interval acc;
    bool any = false;
    double prev = a;
    for (int i = 1; i <= n; ++i) {
      double r = (i == n) ? b : (a + (b - a) * double(i) / double(n));
      if (r < prev) r = prev;
      if (r > b) r = b;
      interval piece = p(interval(prev, r));
      acc = any ? hull2(acc, piece) : piece;
      any = true;
      prev = r;
    }
    return acc;
  }
};

struct Chart {
  string name;
  bool curved = false;
  IVector Q, U, S;
  IMatrix Minv;
  interval det;
  interval dotUS;
  Chart() : Q(2), U(2), S(2), Minv(2, 2) {}
};

struct HSet {
  string name;
  int chart = -1;
  interval c, w, b;
};

struct Trans {
  int src = -1;
  int dst = -1;
  int orient = 0;
};

struct Box {
  interval u, s;
};

struct Image {
  interval u, s;
};

struct GoalResult {
  bool ok = false;
  bool refuted = false;
  long boxes_proved = 0;
  long integrations = 0;
  long integrator_failures = 0;
  double min_slack = 1e300;
  interval worst_u, worst_s;
  string reason;
};

static vector<Box> initial_grid(interval u, interval s, double step) {
  double uL = u.leftBound(), uR = u.rightBound();
  double sL = s.leftBound(), sR = s.rightBound();
  auto parts = [&](double a, double b) {
    vector<interval> out;
    if (!(b > a)) {
      out.push_back(interval(a, a));
      return out;
    }
    int n = (int)ceil((b - a) / step);
    if (n < 1) n = 1;
    double prev = a;
    for (int i = 1; i <= n; ++i) {
      double r = (i == n) ? b : (a + (b - a) * double(i) / double(n));
      if (r < prev) r = prev;
      if (r > b) r = b;
      out.push_back(interval(prev, r));
      prev = r;
    }
    return out;
  };
  vector<interval> us = parts(uL, uR);
  vector<interval> ss = parts(sL, sR);
  vector<Box> boxes;
  boxes.reserve(us.size() * ss.size());
  for (const auto& uu : us)
    for (const auto& ss_i : ss) boxes.push_back(Box{uu, ss_i});
  return boxes;
}

class Verifier {
 public:
  Poly poly;
  vector<Chart> charts;
  vector<HSet> sets;
  vector<Trans> trans;
  int m_header = 1;
  interval B;
  IMap vectorField;
  IOdeSolver solver;
  ITimeMap timeMap;
  interval period;

  Verifier(Poly p, interval Bparam)
      : poly(std::move(p)),
        B(Bparam),
        vectorField("par:k,B;time:t;var:x,v;fun:v,-k*v-x*x*x+B*cos(t);"),
        solver(vectorField, kTaylorOrder),
        timeMap(solver),
        period(2.0 * interval::pi()) {
    vectorField.setParameter("k", interval(1) / interval(20));
    vectorField.setParameter("B", B);
    solver.setMaxStep(interval(1) / interval(5));
  }

  int chart_index(const string& name) const {
    for (int i = 0; i < (int)charts.size(); ++i)
      if (charts[i].name == name) return i;
    return -1;
  }
  int set_index(const string& name) const {
    for (int i = 0; i < (int)sets.size(); ++i)
      if (sets[i].name == name) return i;
    return -1;
  }

  C0Rect2Set make_set(const Chart& ch, const interval& u, const interval& s) const {
    double um = mid_inside(u);
    double sm = mid_inside(s);
    interval du = force_contains_zero(u - interval(um));
    interval ds = force_contains_zero(s - interval(sm));
    IVector xc(2);
    IMatrix C(2, 2);
    if (!ch.curved) {
      xc[0] = ch.Q[0] + interval(um) * ch.U[0] + interval(sm) * ch.S[0];
      xc[1] = ch.Q[1] + interval(um) * ch.U[1] + interval(sm) * ch.S[1];
      C[0][0] = ch.U[0];
      C[1][0] = ch.U[1];
      C[0][1] = ch.S[0];
      C[1][1] = ch.S[1];
    } else {
      interval rem = interval(0.5) * poly.ddp(s) * sqr(ds);
      interval ac = interval(um) + poly.p(interval(sm)) + rem;
      xc[0] = ch.Q[0] + ac * ch.U[0] + interval(sm) * ch.S[0];
      xc[1] = ch.Q[1] + ac * ch.U[1] + interval(sm) * ch.S[1];
      interval dp = poly.dp(interval(sm));
      C[0][0] = ch.U[0];
      C[1][0] = ch.U[1];
      C[0][1] = ch.S[0] + dp * ch.U[0];
      C[1][1] = ch.S[1] + dp * ch.U[1];
    }
    IVector r0(2);
    r0[0] = du;
    r0[1] = ds;
    IVector xMid(2), xRad(2);
    split(xc, xMid, xRad);
    xRad[0] = force_contains_zero(xRad[0]);
    xRad[1] = force_contains_zero(xRad[1]);
    return C0Rect2Set(xMid, C, r0, IMatrix::Identity(2), xRad);
  }

  // Returns false on integrator failure. On success, img encloses P(chart(u,s)) in the target chart.
  bool integrate_box(const Chart& src, const Chart& dst, const interval& u, const interval& s, Image& img, string& err) {
    solver.setMaxStep(interval(1) / interval(5));
    try {
      C0Rect2Set set = make_set(src, u, s);
      timeMap(period, set);
      if (!subset(period, set.getCurrentTime())) {
        err = "integrator time interval does not contain 2*pi";
        return false;
      }
      IVector as = set.affineTransformation(dst.Minv, dst.Q);
      img.s = as[1];
      if (dst.curved) img.u = as[0] - poly.p(img.s);
      else img.u = as[0];
      return true;
    } catch (const exception& e) {
      solver.clearCoefficients();
      err = e.what();
      return false;
    }
  }

  enum class Goal { LeftExit, RightExit, Strip };

  // orient +1: left exit requires u < c-w, right exit requires u > c+w.
  // orient -1: the two inequalities are swapped.
  GoalResult prove(Goal goal, int orient, const Chart& src, const HSet& dstSet, const Chart& dst,
                   interval uDom, interval sDom) {
    GoalResult R;
    const interval thrL = dstSet.c - dstSet.w;
    const interval thrR = dstSet.c + dstSet.w;
    const interval bM = dstSet.b;
    const bool left_wants_low = (goal == Goal::LeftExit && orient == +1) || (goal == Goal::RightExit && orient == -1);
    const bool wants_high = (goal == Goal::LeftExit && orient == -1) || (goal == Goal::RightExit && orient == +1);

    vector<Box> stack = initial_grid(uDom, sDom, kInitialParamStep);
    while (!stack.empty()) {
      if (R.integrations >= kMaxIntegrationsPerGoal) {
        R.ok = false;
        R.reason = "integration budget exhausted";
        return R;
      }
      Box box = stack.back();
      stack.pop_back();
      Image img;
      string err;
      ++R.integrations;
      bool ok_int = integrate_box(src, dst, box.u, box.s, img, err);
      if (!ok_int) {
        ++R.integrator_failures;
        interval uL, uR, sL, sR;
        bool su = split_interval(box.u, uL, uR);
        bool ss = split_interval(box.s, sL, sR);
        double wu = diam(box.u).rightBound();
        double ws = diam(box.s).rightBound();
        if (wu >= ws && su) {
          stack.push_back(Box{uL, box.s});
          stack.push_back(Box{uR, box.s});
          continue;
        }
        if (ss) {
          stack.push_back(Box{box.u, sL});
          stack.push_back(Box{box.u, sR});
          continue;
        }
        if (su) {
          stack.push_back(Box{uL, box.s});
          stack.push_back(Box{uR, box.s});
          continue;
        }
        R.ok = false;
        R.reason = "integrator failed on an unsplittable box: " + err;
        R.worst_u = box.u;
        R.worst_s = box.s;
        return R;
      }

      bool proved = false;
      bool refuted = false;
      double slack = 0.0;
      if (goal == Goal::Strip) {
        // |s| < b. Sufficient: s < b.left and s > -b.left.
        double slack_hi = bM.leftBound() - img.s.rightBound();
        double slack_lo = img.s.leftBound() - (-bM.leftBound());
        slack = std::min(slack_hi, slack_lo);
        proved = slack > 0.0;
        // Entire image has |s| >= b_exact.
        refuted = img.s.leftBound() >= bM.rightBound() || img.s.rightBound() <= -bM.rightBound();
      } else if (left_wants_low || (goal != Goal::Strip && !wants_high)) {
        // u < c - w. Sufficient: u.right < (c-w).left.
        slack = thrL.leftBound() - img.u.rightBound();
        proved = slack > 0.0;
        refuted = img.u.leftBound() > thrL.rightBound();
      } else {
        // u > c + w. Sufficient: u.left > (c+w).right.
        slack = img.u.leftBound() - thrR.rightBound();
        proved = slack > 0.0;
        refuted = img.u.rightBound() < thrR.leftBound();
      }

      if (proved) {
        ++R.boxes_proved;
        if (slack < R.min_slack) {
          R.min_slack = slack;
          R.worst_u = img.u;
          R.worst_s = img.s;
        }
        continue;
      }
      if (refuted) {
        R.ok = false;
        R.refuted = true;
        R.reason = "enclosure lies entirely outside the required side";
        R.min_slack = slack;
        R.worst_u = img.u;
        R.worst_s = img.s;
        return R;
      }

      interval uL, uR, sL, sR;
      bool su = split_interval(box.u, uL, uR);
      bool ss = split_interval(box.s, sL, sR);
      double wu = diam(box.u).rightBound();
      double ws = diam(box.s).rightBound();
      bool too_small = wu < 1e-14 && ws < 1e-12;
      if (too_small || (!su && !ss)) {
        R.ok = false;
        R.reason = "enclosure still straddles the threshold at minimum box size";
        R.min_slack = slack;
        R.worst_u = img.u;
        R.worst_s = img.s;
        return R;
      }
      if (wu >= ws && su) {
        stack.push_back(Box{uL, box.s});
        stack.push_back(Box{uR, box.s});
      } else if (ss) {
        stack.push_back(Box{box.u, sL});
        stack.push_back(Box{box.u, sR});
      } else {
        stack.push_back(Box{uL, box.s});
        stack.push_back(Box{uR, box.s});
      }
    }
    if (R.boxes_proved <= 0) {
      R.ok = false;
      R.reason = "empty cover";
      return R;
    }
    R.ok = true;
    return R;
  }
};

static void load_chain(const string& path, Verifier& V) {
  ifstream in(path);
  if (!in) throw runtime_error("cannot open " + path);
  string line;
  int lineno = 0;
  while (getline(in, line)) {
    ++lineno;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    istringstream ls(line);
    string tag;
    if (!(ls >> tag)) continue;
    if (tag[0] == '#') continue;
    auto need = [&](string& tok) {
      if (!(ls >> tok)) throw runtime_error("line " + to_string(lineno) + " truncated");
    };
    if (tag == "m") {
      string tok;
      need(tok);
      V.m_header = stoi(tok);
    } else if (tag == "p") {
      string tok;
      need(tok);
      int n = stoi(tok);
      vector<interval> coeffs;
      for (int i = 0; i < n; ++i) {
        need(tok);
        coeffs.push_back(enclose_decimal(tok));
      }
      if (n < 1) throw runtime_error("polynomial needs at least one coefficient");
      V.poly = Poly(std::move(coeffs));
    } else if (tag == "chart") {
      string name, kind, tok;
      need(name);
      need(kind);
      Chart ch;
      ch.name = name;
      if (kind == "curved") ch.curved = true;
      else if (kind == "affine") ch.curved = false;
      else throw runtime_error("unknown chart kind " + kind);
      need(tok); ch.Q[0] = enclose_decimal(tok);
      need(tok); ch.Q[1] = enclose_decimal(tok);
      need(tok); ch.U[0] = enclose_decimal(tok);
      need(tok); ch.U[1] = enclose_decimal(tok);
      need(tok); ch.S[0] = enclose_decimal(tok);
      need(tok); ch.S[1] = enclose_decimal(tok);
      ch.det = ch.U[0] * ch.S[1] - ch.U[1] * ch.S[0];
      ch.dotUS = ch.U[0] * ch.S[0] + ch.U[1] * ch.S[1];
      if (subset(interval(0.0), ch.det)) throw runtime_error("chart " + name + " has det containing 0");
      ch.Minv[0][0] = ch.S[1] / ch.det;
      ch.Minv[0][1] = -ch.S[0] / ch.det;
      ch.Minv[1][0] = -ch.U[1] / ch.det;
      ch.Minv[1][1] = ch.U[0] / ch.det;
      if (V.chart_index(name) >= 0) throw runtime_error("duplicate chart " + name);
      V.charts.push_back(ch);
    } else if (tag == "set") {
      string name, tok;
      need(name);
      HSet h;
      h.name = name;
      h.chart = V.chart_index(name);
      if (h.chart < 0) throw runtime_error("set " + name + " has no chart");
      need(tok); h.c = enclose_decimal(tok);
      need(tok); h.w = enclose_decimal(tok);
      need(tok); h.b = enclose_decimal(tok);
      if (!(h.w.leftBound() > 0.0) || !(h.b.leftBound() > 0.0))
        throw runtime_error("set " + name + " widths must be positive");
      if (V.set_index(name) >= 0) throw runtime_error("duplicate set " + name);
      V.sets.push_back(h);
    } else if (tag == "trans") {
      string a, b, tok;
      need(a);
      need(b);
      need(tok);
      Trans t;
      t.src = V.set_index(a);
      t.dst = V.set_index(b);
      if (t.src < 0 || t.dst < 0) throw runtime_error("trans references unknown set on line " + to_string(lineno));
      t.orient = stoi(tok);
      if (t.orient != 1 && t.orient != -1) throw runtime_error("orientation must be +1 or -1");
      V.trans.push_back(t);
    } else {
      throw runtime_error("unknown tag " + tag + " on line " + to_string(lineno));
    }
  }
}

struct Zono {
  IVector Q, U, S;
  interval a, s;
  Zono() : Q(2), U(2), S(2) {}
  interval xbox() const { return Q[0] + a * U[0] + s * S[0]; }
  interval vbox() const { return Q[1] + a * U[1] + s * S[1]; }
};

static interval project(const Zono& Z, double dx, double dy) {
  interval dirx(dx), diry(dy);
  return Z.Q[0] * dirx + Z.Q[1] * diry
      + Z.a * (Z.U[0] * dirx + Z.U[1] * diry)
      + Z.s * (Z.S[0] * dirx + Z.S[1] * diry);
}

static bool zono_separated(const Zono& A, const Zono& B, double& gap, string& how) {
  if (separated_intervals(A.xbox(), B.xbox())) {
    gap = std::max(A.xbox().leftBound() - B.xbox().rightBound(), B.xbox().leftBound() - A.xbox().rightBound());
    how = "bounding boxes (x)";
    return true;
  }
  if (separated_intervals(A.vbox(), B.vbox())) {
    gap = std::max(A.vbox().leftBound() - B.vbox().rightBound(), B.vbox().leftBound() - A.vbox().rightBound());
    how = "bounding boxes (v)";
    return true;
  }
  auto consider = [&](double dx, double dy, const string& name) -> bool {
    double n2 = dx * dx + dy * dy;
    if (!(n2 > 0.0)) return false;
    interval pA = project(A, dx, dy);
    interval pB = project(B, dx, dy);
    if (separated_intervals(pA, pB)) {
      gap = std::max(pA.leftBound() - pB.rightBound(), pB.leftBound() - pA.rightBound());
      how = name;
      return true;
    }
    return false;
  };
  auto perp = [](const IVector& V, double& dx, double& dy) {
    double mx = mid_inside(V[0]);
    double my = mid_inside(V[1]);
    dx = -my;
    dy = mx;
  };
  double dx, dy;
  perp(A.U, dx, dy);
  if (consider(dx, dy, "normal to first U")) return true;
  perp(A.S, dx, dy);
  if (consider(dx, dy, "normal to first S")) return true;
  perp(B.U, dx, dy);
  if (consider(dx, dy, "normal to second U")) return true;
  perp(B.S, dx, dy);
  if (consider(dx, dy, "normal to second S")) return true;
  return false;
}

static Zono affine_zono(const Chart& ch, const HSet& h) {
  Zono Z;
  Z.Q = ch.Q;
  Z.U = ch.U;
  Z.S = ch.S;
  Z.a = interval((h.c - h.w).leftBound(), (h.c + h.w).rightBound());
  Z.s = interval(-h.b.rightBound(), h.b.rightBound());
  return Z;
}

// Fat parallelogram in (a,s) containing the curved h-set slice s in sRange.
static Zono curved_slice_zono(const Chart& ch, const HSet& h, const Poly& poly, interval sRange) {
  Zono Z;
  Z.Q = ch.Q;
  Z.U = ch.U;
  Z.S = ch.S;
  interval uRange((h.c - h.w).leftBound(), (h.c + h.w).rightBound());
  Z.a = uRange + poly.range(sRange);
  Z.s = sRange;
  return Z;
}

static bool curved_vs_zono(const Chart& ch, const HSet& h, const Poly& poly, const Zono& other,
                           interval sRange, int depth, double& gap, string& how) {
  Zono mine = curved_slice_zono(ch, h, poly, sRange);
  if (zono_separated(mine, other, gap, how)) return true;
  if (depth >= 60 || !(diam(sRange).rightBound() > 1e-14)) return false;
  interval L, R;
  if (!split_interval(sRange, L, R)) return false;
  double g1 = 0, g2 = 0;
  string h1, h2;
  if (!curved_vs_zono(ch, h, poly, other, L, depth + 1, g1, h1)) return false;
  if (!curved_vs_zono(ch, h, poly, other, R, depth + 1, g2, h2)) return false;
  gap = std::min(g1, g2);
  how = "s-subdivision, " + h1;
  return true;
}

static void print_goal(const string& name, const GoalResult& R) {
  cout << "  " << name << ": ";
  if (R.ok) {
    cout << "VERIFIED  boxes=" << R.boxes_proved
         << "  integrations=" << R.integrations
         << "  integrator_retries=" << R.integrator_failures
         << "  min_slack=" << R.min_slack << "\n";
    cout << "    worst image u=" << R.worst_u << "  s=" << R.worst_s << "\n";
  } else {
    cout << "NOT VERIFIED";
    if (R.refuted) cout << " (refuted)";
    cout << "  boxes_proved=" << R.boxes_proved
         << "  integrations=" << R.integrations
         << "  integrator_retries=" << R.integrator_failures << "\n";
    cout << "    reason: " << R.reason << "\n";
    cout << "    witness u=" << R.worst_u << "  s=" << R.worst_s << "\n";
  }
}

static bool taylor_self_test(const Poly& poly) {
  interval samples[] = {interval(-0.15, 0.15), interval(0.0, 0.01), interval(-0.12, -0.11), interval(0.14, 0.15)};
  for (interval srange : samples) {
    double sm = mid_inside(srange);
    interval ds = srange - interval(sm);
    interval rem = interval(0.5) * poly.ddp(srange) * sqr(ds);
    interval form = poly.p(interval(sm)) + poly.dp(interval(sm)) * ds + rem;
    for (int i = 0; i <= 8; ++i) {
      double t = srange.leftBound() + (srange.rightBound() - srange.leftBound()) * double(i) / 8.0;
      interval pt = poly.p(interval(t, t));
      if (!subset(pt, form)) {
        cout << "polynomial Taylor self-test FAILED at t=" << t << " pt=" << pt << " form=" << form << "\n";
        return false;
      }
    }
  }
  if (!subset(interval(0.0), poly.p(interval(0.0)))) {
    cout << "polynomial Taylor self-test FAILED: p(0) does not contain 0, p(0)=" << poly.p(interval(0.0)) << "\n";
    return false;
  }
  return true;
}

int main(int argc, char** argv) {
  cout.precision(17);
  cout << unitbuf;
  string path = "design_chain.txt";
  bool flip = false;
  bool customB = false;
  string Btext;
  for (int i = 1; i < argc; ++i) {
    string a = argv[i];
    if (a == "--flip-last") flip = true;
    else if (a == "--B") {
      if (i + 1 >= argc) {
        cerr << "missing value after --B\n";
        return 2;
      }
      Btext = argv[++i];
      customB = true;
    } else {
      path = a;
    }
  }

  interval B = customB ? enclose_decimal(Btext) : (interval(15) / interval(2));
  // Placeholder polynomial; load_chain overwrites coefficients and derivatives.
  Poly placeholder({interval(0.0)});
  try {
    Verifier V(placeholder, B);
    load_chain(path, V);
    if (V.m_header != 1) {
      cout << "NOT VERIFIED\n";
      cout << "m=" << V.m_header << " but this verifier only implements a single period of P (spec defines P, and the file's m is unspecified except that the supplied chain has m=1).\n";
      return 2;
    }
    if (flip) {
      if (V.trans.empty()) throw runtime_error("no transitions to flip");
      V.trans.back().orient = -V.trans.back().orient;
    }

    cout << "checker B\n";
    cout << "chain=" << path << "\n";
    cout << "B=" << V.B << (customB ? " (override)\n" : " (15/2)\n");
    cout << "k=" << (interval(1) / interval(20)) << "\n";
    cout << "period=" << V.period << "\n";
    cout << "taylor_order=" << kTaylorOrder << "  max_step<=1/5\n";
    cout << "flip_last=" << (flip ? "yes" : "no") << "\n";
    cout << "m=" << V.m_header << " (single application of P per link)\n";
    if (!taylor_self_test(V.poly)) {
      cout << "NOT VERIFIED\n";
      return 2;
    }
    cout << "polynomial Taylor self-test: passed\n";

    cout << "charts:\n";
    for (const auto& ch : V.charts) {
      cout << "  " << ch.name << (ch.curved ? " curved" : " affine")
           << " det=" << ch.det << " U·S=" << ch.dotUS << "\n";
    }
    cout << "transitions:\n";
    for (const auto& t : V.trans) {
      cout << "  " << V.sets[t.src].name << " -> " << V.sets[t.dst].name
           << " orient " << t.orient << "\n";
    }

    // Point self-test: the chart-A origin should come back near itself (periodic point to the printed precision).
    {
      int ia = V.chart_index("A");
      if (ia < 0) throw runtime_error("no chart A");
      Image img;
      string err;
      if (!V.integrate_box(V.charts[ia], V.charts[ia], interval(0.0), interval(0.0), img, err)) {
        cout << "point self-test FAILED: " << err << "\nNOT VERIFIED\n";
        return 2;
      }
      cout << "point self-test P(A(0,0)) in chart A: u=" << img.u << " s=" << img.s << "\n";
    }

    auto started = chrono::steady_clock::now();
    bool all = true;
    long total_integrations = 1;  // the point self-test
    for (int ti = 0; ti < (int)V.trans.size(); ++ti) {
      const Trans& t = V.trans[ti];
      const HSet& N = V.sets[t.src];
      const HSet& M = V.sets[t.dst];
      const Chart& src = V.charts[N.chart];
      const Chart& dst = V.charts[M.chart];
      auto t0 = chrono::steady_clock::now();
      cout << "\ntransition " << (ti + 1) << ": " << N.name << " -> " << M.name
           << " orient " << t.orient << "\n";

      interval uL = N.c - N.w;
      interval uR = N.c + N.w;
      interval sDom(-N.b.rightBound(), N.b.rightBound());
      interval uDom((N.c - N.w).leftBound(), (N.c + N.w).rightBound());

      GoalResult eL = V.prove(Verifier::Goal::LeftExit, t.orient, src, M, dst, uL, sDom);
      print_goal("(E-) left edge u=c-w", eL);
      GoalResult eR = V.prove(Verifier::Goal::RightExit, t.orient, src, M, dst, uR, sDom);
      print_goal("(E+) right edge u=c+w", eR);
      GoalResult sG = V.prove(Verifier::Goal::Strip, t.orient, src, M, dst, uDom, sDom);
      print_goal("(S)  |s| < b_dst on the whole h-set", sG);

      double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
      cout << "  transition_seconds=" << sec << "\n";
      total_integrations += eL.integrations + eR.integrations + sG.integrations;
      if (!(eL.ok && eR.ok && sG.ok)) all = false;
    }

    cout << "\n(D) pairwise disjointness\n";
    bool disjoint_ok = true;
    for (int i = 0; i < (int)V.sets.size(); ++i) {
      for (int j = i + 1; j < (int)V.sets.size(); ++j) {
        const HSet& A = V.sets[i];
        const HSet& B = V.sets[j];
        const Chart& cA = V.charts[A.chart];
        const Chart& cB = V.charts[B.chart];
        double gap = 0;
        string how;
        bool sep = false;
        if (!cA.curved && !cB.curved) {
          sep = zono_separated(affine_zono(cA, A), affine_zono(cB, B), gap, how);
        } else if (cA.curved && !cB.curved) {
          interval sRange(-A.b.rightBound(), A.b.rightBound());
          sep = curved_vs_zono(cA, A, V.poly, affine_zono(cB, B), sRange, 0, gap, how);
        } else if (!cA.curved && cB.curved) {
          interval sRange(-B.b.rightBound(), B.b.rightBound());
          sep = curved_vs_zono(cB, B, V.poly, affine_zono(cA, A), sRange, 0, gap, how);
        } else {
          // Both curved: slice the first and test against fat slices of the second by recursion on the first only,
          // using a fat parallelogram of the whole second set (sufficient).
          interval sRange(-A.b.rightBound(), A.b.rightBound());
          interval sB(-B.b.rightBound(), B.b.rightBound());
          Zono fatB = curved_slice_zono(cB, B, V.poly, sB);
          sep = curved_vs_zono(cA, A, V.poly, fatB, sRange, 0, gap, how);
        }
        cout << "  " << A.name << " vs " << B.name << ": ";
        if (sep) cout << "VERIFIED  by " << how << "  gap=" << gap << "\n";
        else {
          cout << "NOT VERIFIED  (no separating direction found)\n";
          disjoint_ok = false;
        }
      }
    }
    if (!disjoint_ok) all = false;

    double total_sec = chrono::duration<double>(chrono::steady_clock::now() - started).count();
    cout << "\nintegrations=" << total_integrations << "\n";
    cout << "seconds=" << total_sec << "\n";
    cout << (all ? "ALL VERIFIED" : "NOT VERIFIED") << "\n";
    return all ? 0 : 2;
  } catch (const exception& e) {
    cout << "fatal: " << e.what() << "\nNOT VERIFIED\n";
    return 2;
  }
}
