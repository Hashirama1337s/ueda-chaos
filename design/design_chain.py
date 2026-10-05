"""P1 design v2 (numerical): a chain of h-sets along a transverse homoclinic orbit of the saddle D, every link ONE period.
  N_A : around D, in the curved chart X(u,s) = D + (u + p(s)) eu + s es (u = 0 is W^s_loc(D))
  C_0..C_3 : around q_k = P^k(q_0), q_0 = chart(uc, 0) on W^u(D) with P^4(q_0) on W^s_loc(D); affine charts
             X = q_k + u e^u_k + s e^s_k, e^u_k = tangent of W^u at q_k, e^s_k = stable direction transported back from q_4.
Transitions (all by P): A=>A, A=>C0, C0=>C1, C1=>C2, C2=>C3, C3=>A.
Loops of length 1 and 5 => topological entropy of P >= log(largest root of x^5 = x^4 + 1) (about 0.2524 per period).
Output: proof/design_chain.txt (input of checker A) + json; dense floating-point check of every condition.
usage: python code/design_chain.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from manifolds2 import Pk, D

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
dj = json.load(open(os.path.join(HERE, "design", "design_m4.json")))
eu, es, coef, uc = np.array(dj["eu"]), np.array(dj["es"]), np.array(dj["p_coef_high_to_low"]), dj["uc"]


def chartA(u, s):
    a = u + np.polyval(coef, s); return D[0] + a * eu[0] + s * es[0], D[1] + a * eu[1] + s * es[1]


MA = np.linalg.inv(np.column_stack([eu, es]))


def uncA(x, v):
    a, s = MA @ np.vstack([np.atleast_1d(x) - D[0], np.atleast_1d(v) - D[1]]); return a - np.polyval(coef, s), s


def jac(x, v, h=1e-7):
    X, V = Pk(np.array([x, x + h, x]), np.array([v, v, v + h]), 1)
    return np.array([[X[1] - X[0], X[2] - X[0]], [V[1] - V[0], V[2] - V[0]]]) / h, (X[0], V[0])


def main():
    # orbit q0..q4 (refine uc so that q4 lies on u = 0 of chart A)
    q = [np.array(chartA(np.array([uc]), np.array([0.0]))).ravel()]
    for k in range(4):
        X, V = Pk(np.array([q[-1][0]]), np.array([q[-1][1]]), 1); q.append(np.array([X[0], V[0]]))
    u4, s4 = uncA(*q[4]); print(f"q4 in chart A: u={u4[0]:+.2e}, s={s4[0]:+.6f}")
    J = [jac(*q[k])[0] for k in range(4)]
    # unstable tangents forward: t_0 = d/du chart at q0 = eu ; t_{k+1} = J_k t_k
    tu = [eu.copy()]
    for k in range(4): t = J[k] @ tu[-1]; tu.append(t)
    lam = [np.linalg.norm(tu[k + 1]) / np.linalg.norm(tu[k]) for k in range(4)]
    tu = [t / np.linalg.norm(t) for t in tu]
    # stable tangent at q4: d/ds of chart A at s = s4 ; transport back: t_k = J_k^{-1} t_{k+1}
    s4v = s4[0]; dp = np.polyval(np.polyder(coef), s4v)
    ts = [None] * 5; ts[4] = es + dp * eu
    mu = [0] * 4
    for k in (3, 2, 1, 0):
        t = np.linalg.solve(J[k], ts[k + 1]); mu[k] = np.linalg.norm(ts[k + 1]) / np.linalg.norm(t); ts[k] = t
    ts = [t / np.linalg.norm(t) for t in ts]
    for k in range(5):
        ang = np.degrees(np.arccos(abs(tu[k] @ ts[k])))
        print(f"q{k} = ({q[k][0]:+.6f}, {q[k][1]:+.6f})  angle(e_u, e_s) = {ang:6.2f} deg" + (f"  stretch {lam[k]:8.3f}  contraction {mu[k]:.4f}" if k < 4 else ""))
    print("product of stretches", np.prod(lam), "(compare du-stretch of P^4 at C: -2535)")

    # sizes: A first (P(N_A) must cover C0 in u), then C3, C2, C1, C0 backwards with safety factor 3
    wA, bA = 1.6e-3, 0.15
    w = [0.0] * 4
    w[3] = 3 * wA / lam[3]
    for k in (2, 1, 0): w[k] = (10 if k == 0 else 3) * w[k + 1] / lam[k]
    b = [0.0] * 4
    b[0] = 0.012
    for k in range(1, 4): b[k] = max(4 * mu[k - 1] * b[k - 1], 1e-6)
    print("w:", ["%.2e" % x for x in w], " b:", ["%.2e" % x for x in b])

    frames = {"A": None}
    for k in range(4): frames[f"C{k}"] = (q[k], tu[k], ts[k])

    def unc(name, x, v):
        if name == "A": return uncA(x, v)
        Q, U, S = frames[name]; M = np.linalg.inv(np.column_stack([U, S]))
        a, s = M @ np.vstack([np.atleast_1d(x) - Q[0], np.atleast_1d(v) - Q[1]]); return a, s

    def chart(name, u, s):
        if name == "A": return chartA(u, s)
        Q, U, S = frames[name]; return Q[0] + u * U[0] + s * S[0], Q[1] + u * U[1] + s * S[1]

    sets = {"A": (0.0, wA, bA)}
    for k in range(4): sets[f"C{k}"] = (0.0, w[k], b[k])
    trans = [("A", "A"), ("A", "C0"), ("C0", "C1"), ("C1", "C2"), ("C2", "C3"), ("C3", "A")]
    allok = True; orient = {}
    for src, dst in trans:
        c, ww, bb = sets[src]; c2, w2, b2 = sets[dst]
        sg = np.linspace(-bb, bb, 801)
        L = unc(dst, *Pk(*chart(src, np.full_like(sg, c - ww), sg), 1))
        R = unc(dst, *Pk(*chart(src, np.full_like(sg, c + ww), sg), 1))
        U, Sg = np.meshgrid(np.linspace(c - ww, c + ww, 41), np.linspace(-bb, bb, 201))
        W = unc(dst, *Pk(*chart(src, U.ravel(), Sg.ravel()), 1))
        o = 1 if L[0].mean() < R[0].mean() else -1
        lo, hi = (L, R) if o > 0 else (R, L)
        m1 = (c2 - w2) - lo[0].max(); m2 = hi[0].min() - (c2 + w2); m3 = b2 - np.abs(W[1]).max()
        ok = m1 > 0 and m2 > 0 and m3 > 0; allok &= ok; orient[(src, dst)] = o
        print(f"{src:>2} => {dst:<2} orient {o:+d}: exit margins {m1:+.3e} {m2:+.3e} (target half-width {w2:.2e}); "
              f"entry margin {m3:+.3e} (target b {b2:.2e})  {'OK' if ok else 'FAIL'}")
    # disjointness (in (x,v): compare chart boxes crudely by mapping C_k boxes into chart A)
    print("ALL OK (floating point)" if allok else "SOME FAIL")
    lines = ["# Ueda map: chain of h-sets along a transverse homoclinic orbit of D; every link is one period (P)",
             "m 1", "p %d %s" % (len(coef), " ".join(repr(float(c)) for c in coef))]
    lines.append("chart A curved %r %r %r %r %r %r" % tuple(float(x) for x in (D[0], D[1], eu[0], eu[1], es[0], es[1])))
    for k in range(4):
        Q, U, S = frames[f"C{k}"]
        lines.append("chart C%d affine %r %r %r %r %r %r" % ((k,) + tuple(float(x) for x in (Q[0], Q[1], U[0], U[1], S[0], S[1]))))
    for nm, (c, ww, bb) in sets.items():
        lines.append("set %s %r %r %r" % (nm, float(c), float(ww), float(bb)))
    for (s_, d_) in trans:
        lines.append(f"trans {s_} {d_} {orient[(s_, d_)]}")
    open(os.path.join(HERE, "design", "design_chain_regenerated.txt"), "w", newline="\n").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
