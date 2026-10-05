"""P0(c) v2: adaptive growth of W^u(D) for Ueda's map, and its returns to the LOCAL stable segment of D.
A return of P^m(fundamental domain of W^u_loc) across the local stable segment = a numerical transverse homoclinic point,
i.e. where a horseshoe for P^m lives (Smale-Birkhoff). Numerical only (vectorised RK4); rigorous work is CAPD's job.
usage: python code/manifolds2.py [max_iter]
"""
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

K, B = 0.05, 7.5
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NST = 2500


def Pk(x, v, k):
    """k-fold time-2pi map, vectorised RK4."""
    h = 2 * np.pi / NST
    for _ in range(k):
        t = 0.0
        for _ in range(NST):
            a1x, a1v = v, -K * v - x ** 3 + B * np.cos(t)
            x2, v2 = x + h / 2 * a1x, v + h / 2 * a1v; a2x, a2v = v2, -K * v2 - x2 ** 3 + B * np.cos(t + h / 2)
            x3, v3 = x + h / 2 * a2x, v + h / 2 * a2v; a3x, a3v = v3, -K * v3 - x3 ** 3 + B * np.cos(t + h / 2)
            x4, v4 = x + h * a3x, v + h * a3v; a4x, a4v = v4, -K * v4 - x4 ** 3 + B * np.cos(t + h)
            x = x + h / 6 * (a1x + 2 * a2x + 2 * a3x + a4x); v = v + h / 6 * (a1v + 2 * a2v + 2 * a3v + a4v)
            t += h
    return x, v


D = np.array([2.8913632426, 0.2469647907])
LAM_U, LAM_S = 13.003744, 0.056169


def eigvecs():
    eps = 1e-6
    x = np.array([D[0], D[0] + eps, D[0]]); v = np.array([D[1], D[1], D[1] + eps])
    X, V = Pk(x, v, 1)
    J = np.array([[X[1] - X[0], X[2] - X[0]], [V[1] - V[0], V[2] - V[0]]]) / eps
    ev, E = np.linalg.eig(J)
    iu = int(np.argmax(abs(ev)))
    eu, es = E[:, iu].real, E[:, 1 - iu].real
    return eu / np.linalg.norm(eu), es / np.linalg.norm(es), ev


def grow(eu, k, sign, s0=2e-3, delta=0.004, maxpts=400000):
    """images under P^k of the fundamental segment {D + s eu : s in [s0, LAM_U s0]} (one branch), adaptively refined."""
    s = s0 * LAM_U ** np.linspace(0, 1, 2001)
    for rnd in range(40):
        X, V = Pk(D[0] + sign * s * eu[0], D[1] + sign * s * eu[1], k)
        gap = np.hypot(np.diff(X), np.diff(V))
        bad = np.nonzero(gap > delta)[0]
        if len(bad) == 0 or len(s) > maxpts:
            return s, X, V, len(bad)
        s = np.sort(np.r_[s, 0.5 * (s[bad] + s[bad + 1])])
    return s, X, V, len(bad)


def cross_stable(X, V, es, eu, L=0.3):
    """crossings of the polyline (X,V) with the local stable segment D + t es, |t| <= L. Returns t, angle, index."""
    # coordinates relative to D in the (eu, es) frame: crossing when the eu-component changes sign
    M = np.linalg.inv(np.column_stack([eu, es]))
    a, b = M @ np.vstack([X - D[0], V - D[1]])  # a: unstable coordinate, b: stable coordinate
    out = []
    idx = np.nonzero(np.sign(a[:-1]) * np.sign(a[1:]) < 0)[0]
    for i in idx:
        t = a[i] / (a[i] - a[i + 1]); bb = b[i] + t * (b[i + 1] - b[i])
        if abs(bb) <= L:
            d = np.array([X[i + 1] - X[i], V[i + 1] - V[i]]); d /= np.linalg.norm(d)
            ang = np.degrees(np.arcsin(min(1, abs(d[0] * es[1] - d[1] * es[0]))))
            out.append((float(bb), float(ang), int(i)))
    return out


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    eu, es, ev = eigvecs()
    print(f"eigen: {ev}, eu={eu}, es={es}, angle(eu,es)={np.degrees(np.arccos(abs(eu @ es))):.2f} deg")
    res = {"D": D.tolist(), "eu": eu.tolist(), "es": es.tolist(), "returns": []}
    fig, ax = plt.subplots(figsize=(7, 7))
    cols = plt.cm.viridis(np.linspace(0, 1, kmax + 1))
    for sign in (+1, -1):
        for k in range(1, kmax + 1):
            s, X, V, nbad = grow(eu, k, sign)
            cr = cross_stable(X, V, es, eu)
            print(f"branch {'+' if sign > 0 else '-'} P^{k}: {len(s)} pts (unresolved gaps {nbad}); "
                  f"x in [{X.min():.2f},{X.max():.2f}] v in [{V.min():.2f},{V.max():.2f}]; crossings of local W^s: {len(cr)}"
                  + (f"; first: t={cr[0][0]:+.4f}, angle {cr[0][1]:.1f} deg" if cr else ""))
            for c in cr:
                res["returns"].append({"branch": sign, "k": k, "t_stable": c[0], "angle_deg": c[1],
                                       "s": float(s[c[2]])})
            ax.plot(X, V, "-", lw=0.4, color=cols[k], alpha=0.9)
    ax.plot(*D, "r*", ms=10)
    L = 0.3; ax.plot([D[0] - L * es[0], D[0] + L * es[0]], [D[1] - L * es[1], D[1] + L * es[1]], "r-", lw=1.5)
    ax.set_xlabel("x"); ax.set_ylabel("dx/dt"); ax.set_title("W^u(D) of Ueda's map (colour = iterate), red = local W^s(D)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "unstable_manifold.png"), dpi=170)
    json.dump(res, open(os.path.join(HERE, "design", "manifolds2.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
