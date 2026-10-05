"""Figure: Ueda's attractor with the verified chain of h-sets (A around the saddle D, C0..C3 along the homoclinic loop)."""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from manifolds2 import Pk
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = [l.split() for l in open(os.path.join(HERE, "proof", "design_chain.txt")) if l.strip() and not l.startswith("#")]
p = [float(x) for x in [l for l in L if l[0] == "p"][0][2:]]
charts = {l[1]: (l[2], np.array(list(map(float, l[3:9])))) for l in L if l[0] == "chart"}
sets = {l[1]: tuple(map(float, l[2:5])) for l in L if l[0] == "set"}
def X(name, u, s):
    ty, c = charts[name]; a = u + (np.polyval(p, s) if ty == "curved" else 0)
    return c[0] + a * c[2] + s * c[4], c[1] + a * c[3] + s * c[5]
rng = np.random.default_rng(2); x = rng.uniform(1.6, 3.3, 3000); v = rng.uniform(-4, 5, 3000)
for _ in range(40): x, v = Pk(x, v, 1)
pts = []
for _ in range(20): x, v = Pk(x, v, 1); pts.append(np.c_[x, v])
pts = np.vstack(pts)
fig, axs = plt.subplots(1, 2, figsize=(12, 6))
for ax, zoom in zip(axs, (False, True)):
    ax.plot(pts[:, 0], pts[:, 1], ",", color="0.55", alpha=0.6)
    cols = {"A": "tab:red", "C0": "tab:orange", "C1": "tab:green", "C2": "tab:blue", "C3": "tab:purple"}
    for nm, (c, w, b) in sets.items():
        s = np.linspace(-b, b, 200); u = np.linspace(c - w, c + w, 50)
        bx = np.r_[X(nm, np.full_like(s, c - w), s)[0], X(nm, u, np.full_like(u, b))[0], X(nm, np.full_like(s, c + w), s[::-1])[0], X(nm, u[::-1], np.full_like(u, -b))[0]]
        by = np.r_[X(nm, np.full_like(s, c - w), s)[1], X(nm, u, np.full_like(u, b))[1], X(nm, np.full_like(s, c + w), s[::-1])[1], X(nm, u[::-1], np.full_like(u, -b))[1]]
        ax.fill(bx, by, color=cols[nm], alpha=0.5, lw=1.5, ec=cols[nm], label=nm)
        cx, cy = X(nm, np.array([c]), np.array([0.0]))
        if zoom and nm in ("A", "C0"):
            ax.annotate(nm, (cx[0], cy[0]), xytext=(8, 6), textcoords="offset points", color=cols[nm], fontsize=11, weight="bold")
        elif not zoom and nm in ("C1", "C2", "C3"):
            ax.annotate(nm, (cx[0], cy[0]), xytext=(8, 6), textcoords="offset points", color=cols[nm], fontsize=11, weight="bold")
        elif not zoom and nm == "A":
            ax.annotate("A, C0 (at D)", (cx[0], cy[0]), xytext=(-95, -22), textcoords="offset points", color="tab:red", fontsize=10, weight="bold")
    q = [np.array(X("C0", np.array([0.]), np.array([0.]))).ravel()]
    for k in range(4):
        a, bb = Pk(np.array([q[-1][0]]), np.array([q[-1][1]]), 1); q.append(np.array([a[0], bb[0]]))
    q = np.array(q)
    ax.plot(q[:, 0], q[:, 1], "k--", lw=0.8, alpha=0.7)
    D = charts["A"][1][:2]; ax.plot(*D, "k*", ms=12)
    if zoom:
        ax.set_xlim(2.86, 2.93); ax.set_ylim(0.05, 0.48); ax.set_title("zoom near the saddle D: A (red) and C0 (orange)", fontsize=10)
    else:
        ax.set_title("Ueda's attractor (grey) and the verified chain A -> C0 -> C1 -> C2 -> C3 -> A\n(dashed: joins the 5 successive points of one numerically computed homoclinic loop of the saddle D, star)", fontsize=10)
    ax.set_xlabel("x"); ax.set_ylabel("dx/dt")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "proof_chain.png"), dpi=170)
print("ok")
