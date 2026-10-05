"""P0 numerics for Ueda's map P: (x, v) at t=0 -> (x, v) at t=2*pi for x'' + K x' + x^3 = B cos t.
Positive control: reproduce Ueda's attractor picture; Lyapunov exponent; fixed/period-2 points and their stability.
Vectorised RK4 (many initial conditions at once) for exploration only - NOT rigorous.
usage: python code/explore.py
"""
import os, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

K, B = 0.05, 7.5
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NSTEP = 4000  # RK4 steps per period (h = 2pi/4000 ~ 1.6e-3); checked against 8000 below


def f(t, x, v):
    return v, -K * v - x ** 3 + B * np.cos(t)


def P(x, v, n=NSTEP, var=False):
    """time-2pi map of arrays x, v; optionally the 2x2 Jacobian (variational equations)."""
    h = 2 * np.pi / n
    x = np.array(x, float); v = np.array(v, float)
    if var:
        J = np.zeros(x.shape + (2, 2)); J[..., 0, 0] = 1; J[..., 1, 1] = 1
    t = 0.0
    for _ in range(n):
        if not var:
            k1 = f(t, x, v)
            k2 = f(t + h / 2, x + h / 2 * k1[0], v + h / 2 * k1[1])
            k3 = f(t + h / 2, x + h / 2 * k2[0], v + h / 2 * k2[1])
            k4 = f(t + h, x + h * k3[0], v + h * k3[1])
            x = x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            v = v + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        else:
            def g(t, x, v, J):
                dx, dv = f(t, x, v)
                A = np.zeros(x.shape + (2, 2)); A[..., 0, 1] = 1; A[..., 1, 0] = -3 * x ** 2; A[..., 1, 1] = -K
                return dx, dv, A @ J
            a1 = g(t, x, v, J)
            a2 = g(t + h / 2, x + h / 2 * a1[0], v + h / 2 * a1[1], J + h / 2 * a1[2])
            a3 = g(t + h / 2, x + h / 2 * a2[0], v + h / 2 * a2[1], J + h / 2 * a2[2])
            a4 = g(t + h, x + h * a3[0], v + h * a3[1], J + h * a3[2])
            x = x + h / 6 * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0])
            v = v + h / 6 * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
            J = J + h / 6 * (a1[2] + 2 * a2[2] + 2 * a3[2] + a4[2])
        t += h
    return (x, v, J) if var else (x, v)


def attractor(npts=3000, burn=300):
    x, v = np.array([2.5]), np.array([0.0])
    pts = []
    for i in range(burn + npts):
        x, v = P(x, v)
        if i >= burn:
            pts.append((x[0], v[0]))
    return np.array(pts)


def lyapunov(x0, v0, n=2000):
    x, v = np.array([x0]), np.array([v0]); w = np.array([1.0, 0.0]); s = 0.0
    for _ in range(n):
        x, v, J = P(x, v, var=True)
        w = J[0] @ w; nw = np.linalg.norm(w); s += np.log(nw); w /= nw
    return s / (n * 2 * np.pi), s / n  # per unit time, per period


def periodic_points(k, grid=60, box=(-1.0, 4.5, -7.0, 8.0)):
    """Newton on P^k from a grid; return unique points with eigenvalues."""
    xs = np.linspace(box[0], box[1], grid); vs = np.linspace(box[2], box[3], grid)
    X, V = np.meshgrid(xs, vs); x = X.ravel(); v = V.ravel()
    for it in range(25):
        xx, vv = x.copy(), v.copy(); J = np.broadcast_to(np.eye(2), x.shape + (2, 2)).copy()
        for _ in range(k):
            xx, vv, Jk = P(xx, vv, var=True); J = Jk @ J
        F = np.stack([xx - x, vv - v], -1)
        M = J - np.eye(2)
        try:
            d = np.linalg.solve(M, F[..., None])[..., 0]
        except np.linalg.LinAlgError:
            d = np.zeros_like(F)
        step = np.clip(d, -0.5, 0.5)
        x = x - step[:, 0]; v = v - step[:, 1]
        ok = np.isfinite(x) & np.isfinite(v) & (np.abs(x) < 10) & (np.abs(v) < 20)
        x, v = x[ok], v[ok]
    xx, vv = x.copy(), v.copy(); J = np.broadcast_to(np.eye(2), x.shape + (2, 2)).copy()
    for _ in range(k):
        xx, vv, Jk = P(xx, vv, var=True); J = Jk @ J
    res = np.hypot(xx - x, vv - v)
    good = res < 1e-9
    pts = []
    for xi, vi, Ji in zip(x[good], v[good], J[good]):
        if all(np.hypot(xi - p[0], vi - p[1]) > 1e-5 for p in pts):
            ev = np.linalg.eigvals(Ji)
            pts.append((float(xi), float(vi), [complex(e) for e in ev]))
    return pts


def main():
    # accuracy check of the RK4 map: compare 4000 vs 8000 steps on a few points
    x0 = np.array([2.5, 1.0, 3.0]); v0 = np.array([0.0, 2.0, -3.0])
    a = P(x0, v0, NSTEP); b = P(x0, v0, 2 * NSTEP)
    print("RK4 step check |P4000-P8000| =", float(np.max(np.abs(np.r_[a[0] - b[0], a[1] - b[1]]))))
    pts = attractor()
    print(f"attractor: x in [{pts[:,0].min():.3f}, {pts[:,0].max():.3f}], v in [{pts[:,1].min():.3f}, {pts[:,1].max():.3f}]")
    lam_t, lam_p = lyapunov(pts[-1, 0], pts[-1, 1], 1500)
    print(f"largest Lyapunov exponent: {lam_t:.4f} per unit time ({lam_p:.4f} per forcing period)")
    out = {"attractor_box": [float(pts[:, 0].min()), float(pts[:, 0].max()), float(pts[:, 1].min()), float(pts[:, 1].max())],
           "lyap_per_time": lam_t, "lyap_per_period": lam_p}
    for k in (1, 2):
        pp = periodic_points(k)
        out[f"period{k}"] = [{"x": p[0], "v": p[1], "eig": [[e.real, e.imag] for e in p[2]]} for p in pp]
        print(f"\nperiod-{k} points of P (k={k}): {len(pp)}")
        for p in pp:
            ev = p[2]; kind = "saddle" if (abs(ev[0]) > 1) != (abs(ev[1]) > 1) else ("sink" if max(map(abs, ev)) < 1 else "source")
            print(f"  x={p[0]:+.6f} v={p[1]:+.6f}  eig={[f'{e.real:+.4g}{e.imag:+.3g}j' for e in ev]}  {kind}")
    json.dump(out, open(os.path.join(HERE, "design", "explore.json"), "w"), indent=1)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(pts[:, 0], pts[:, 1], ",", color="k", alpha=0.6)
    for k, mk in ((1, "o"), (2, "s")):
        for p in out[f"period{k}"]:
            ax.plot(p["x"], p["v"], mk, ms=6, mfc="none", color="tab:red" if k == 1 else "tab:blue")
    ax.set_xlabel("x"); ax.set_ylabel("dx/dt"); ax.set_title("Ueda's map, x'' + 0.05x' + x^3 = 7.5 cos t (stroboscopic, t = 2πn)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "explore_attractor.png"), dpi=160)


if __name__ == "__main__":
    main()
