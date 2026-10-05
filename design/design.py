"""P1 design (numerical): a 2-symbol horseshoe for P^m near the saddle D of Ueda's map.
Chart (global homeomorphism of R^2):  X(u, s) = D + (u + p(s)) eu + s es,
p = polynomial fit of the local stable manifold W^s_loc(D) (so u = 0 is W^s up to fit error).
h-sets (exit direction u, entry direction s):
  N_A = [-wA, wA] x [-b, b]          (around D)
  N_C = [uc - wC, uc + wC] x [-b, b] (around the preimage of a transverse homoclinic point, uc < 0 on W^u)
Covering conditions for P^m (checked here by dense floating-point sampling; checker A redoes them rigorously):
  (i)   |s| < b on P^m(N_A) and P^m(N_C)                                  (images stay inside in the entry direction)
  (ii)  each X in {A, C}: P^m(one u-edge of X) has u < uc - wC  and  P^m(the other u-edge) has u > wA
=> N_A => N_A, N_A => N_C, N_C => N_A, N_C => N_C : full 2-shift for P^m, h_top(P) >= ln2 / m.
usage: python code/design.py
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from manifolds2 import Pk, D, eigvecs

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def Pinv(x, v, k=1, nst=2500):
    """k-fold inverse map: integrate from t=0 backwards to t=-2pi (field is 2pi-periodic)."""
    K, B = 0.05, 7.5
    h = -2 * np.pi / nst
    for _ in range(k):
        t = 0.0
        for _ in range(nst):
            a1x, a1v = v, -K * v - x ** 3 + B * np.cos(t)
            x2, v2 = x + h / 2 * a1x, v + h / 2 * a1v; a2x, a2v = v2, -K * v2 - x2 ** 3 + B * np.cos(t + h / 2)
            x3, v3 = x + h / 2 * a2x, v + h / 2 * a2v; a3x, a3v = v3, -K * v3 - x3 ** 3 + B * np.cos(t + h / 2)
            x4, v4 = x + h * a3x, v + h * a3v; a4x, a4v = v4, -K * v4 - x4 ** 3 + B * np.cos(t + h)
            x = x + h / 6 * (a1x + 2 * a2x + 2 * a3x + a4x); v = v + h / 6 * (a1v + 2 * a2v + 2 * a3v + a4v)
            t += h
    return x, v


def main():
    eu, es, ev = eigvecs()
    M = np.linalg.inv(np.column_stack([eu, es]))           # (x,v)-D -> (a, s) with a along eu, s along es

    def to_as(x, v):
        a, s = M @ np.vstack([np.atleast_1d(x) - D[0], np.atleast_1d(v) - D[1]]); return a, s
    # --- local stable manifold: backward iterates of a tiny stable segment, then fit a = p(s)
    s0 = np.r_[-np.geomspace(1e-5, 1e-5 / 0.056169, 400)[::-1], np.geomspace(1e-5, 1e-5 / 0.056169, 400)]
    xs, vs = D[0] + s0 * es[0], D[1] + s0 * es[1]
    A, S = [], []
    a, s = to_as(xs, vs); A.append(a); S.append(s)
    for k in range(1, 4):
        X, V = Pinv(xs, vs, k, nst=12500); a, s = to_as(X, V); A.append(a); S.append(s)
    A = np.concatenate(A); S = np.concatenate(S)
    keep = np.abs(S) < 0.3
    deg = 8
    coef = np.polyfit(S[keep], A[keep], deg)
    coef[-1] = 0.0; coef[-2] = 0.0                            # p(0) = p'(0) = 0 (tangency at D)
    resid = np.abs(np.polyval(coef, S[keep]) - A[keep]).max()
    print(f"W^s_loc fit: |s|<0.3, {keep.sum()} pts, deg {deg}, max residual {resid:.2e}; p(+-0.2) = {np.polyval(coef,0.2):+.3e}, {np.polyval(coef,-0.2):+.3e}")

    def chart(u, s):
        a = u + np.polyval(coef, s)
        return D[0] + a * eu[0] + s * es[0], D[1] + a * eu[1] + s * es[1]

    def unchart(x, v):
        a, s = to_as(x, v); return a - np.polyval(coef, s), s

    m = 4
    # --- locate uc: on the line s=0, u<0, find u with u-coordinate of P^m(chart(u,0)) == 0 (closest to the k=4 return)
    us = -np.linspace(0.010, 0.020, 4001)
    X, V = Pk(*chart(us, np.zeros_like(us)), m); uu, ss = unchart(X, V)
    sgn = np.nonzero(np.sign(uu[:-1]) != np.sign(uu[1:]))[0]
    print("zeros of u(P^4) on s=0, u in [-0.02,-0.01]:", [(round(us[i], 7), round(float(ss[i]), 4)) for i in sgn])
    if len(sgn) == 0:
        print("no return found"); return
    i = sgn[0]; uc = us[i] - uu[i] * (us[i + 1] - us[i]) / (uu[i + 1] - uu[i])
    X, V = Pk(*chart(np.array([uc]), np.array([0.0])), m); u1, s1 = unchart(X, V)
    # local stretch factor at C along u
    du = 1e-9; X2, V2 = Pk(*chart(np.array([uc + du]), np.array([0.0])), m); u2, _ = unchart(X2, V2)
    lamC = (u2[0] - u1[0]) / du
    print(f"uc = {uc:.9f}: P^4 lands at (u, s) = ({u1[0]:+.2e}, {s1[0]:+.5f}); du-stretch at C = {lamC:+.1f}; at D = {13.003744**4:.1f}")

    def edges(lo, hi, bb, ns):
        sgrid = np.linspace(-bb, bb, ns); uu_ = np.linspace(lo, hi, 41)
        U, Sg = np.meshgrid(uu_, sgrid); X, V = Pk(*chart(U.ravel(), Sg.ravel()), m); iu, isg = unchart(X, V)
        iu = iu.reshape(U.shape); isg = isg.reshape(U.shape)
        left, right = iu[:, 0], iu[:, -1]
        if left.mean() < right.mean(): return np.abs(isg).max(), left.max(), right.min(), "+"
        return np.abs(isg).max(), right.max(), left.min(), "flip"

    def check(wA, wC, bA, bC, ns=401):
        """golden-mean design: A=>A, A=>C (P^4(N_A) thin, u-edges beyond both strips); C=>A (P^4(N_C) inside |s|<bA,
        u-edges beyond N_A on both sides). Margins > 0 = satisfied (floating-point sampling)."""
        sA, loA, hiA, oA = edges(-wA, wA, bA, ns)
        sC, loC, hiC, oC = edges(uc - wC, uc + wC, bC, ns)
        return {"A_s": float(min(bA, bC) - sA), "A_u": float(min((uc - wC) - loA, hiA - wA)), "A_or": oA,
                "C_s": float(bA - sC), "C_u": float(min(-wA - loC, hiC - wA)), "C_or": oC}

    best = None
    for bA in (0.15,):
        for bC in (5e-4, 1e-3, 2e-3):
            for wA in (7e-7, 1e-6, 1.4e-6):
                for wC in (1e-5, 2e-5, 4e-5):
                    r = check(wA, wC, bA, bC, ns=81)
                    mn = min(r["A_s"] / bC, r["C_s"] / bA, r["A_u"] / 0.016, r["C_u"] / (2535 * wC))
                    print(f"bA={bA} bC={bC:.0e} wA={wA:.0e} wC={wC:.0e}: {r}  norm-min {mn:+.3f}")
                    if best is None or mn > best[0]: best = (mn, bA, bC, wA, wC)
    mn, bA, bC, wA, wC = best
    print(f"BEST: bA={bA} bC={bC} wA={wA} wC={wC} normalised min margin {mn:+.3f}")
    r2 = check(wA, wC, bA, bC, ns=2001)
    print("dense re-check (2001 s-samples):", r2)
    json.dump({"m": m, "D": D.tolist(), "eu": eu.tolist(), "es": es.tolist(), "p_coef_high_to_low": coef.tolist(),
               "uc": float(uc), "bA": bA, "bC": bC, "wA": wA, "wC": wC, "check": r2, "fit_residual": float(resid),
               "transitions": "A=>A, A=>C, C=>A (golden mean), map P^4"},
              open(os.path.join(HERE, "design", "design_m4.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
