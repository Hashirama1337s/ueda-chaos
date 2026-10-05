"""Positive control: Ueda's attractor picture (vectorised RK4, many orbits at once) + Lyapunov estimate. Numerical only."""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from manifolds2 import Pk
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(1)
x = rng.uniform(-3, 4, 3000); v = rng.uniform(-6, 6, 3000)
for _ in range(60): x, v = Pk(x, v, 1)
pts = []
for _ in range(25):
    x, v = Pk(x, v, 1); pts.append(np.c_[x, v])
pts = np.vstack(pts)
print("attractor box x [%.3f, %.3f] v [%.3f, %.3f], n=%d" % (pts[:,0].min(), pts[:,0].max(), pts[:,1].min(), pts[:,1].max(), len(pts)))
# Lyapunov per period from two nearby orbits (renormalised), 10 orbits
x0, v0 = x[:10].copy(), v[:10].copy(); d0 = 1e-8; xs, vs = x0 + d0, v0.copy(); acc = np.zeros(10); n = 300
for _ in range(n):
    x0, v0 = Pk(x0, v0, 1); xs, vs = Pk(xs, vs, 1)
    d = np.hypot(xs - x0, vs - v0); acc += np.log(d / d0)
    xs, vs = x0 + (xs - x0) * d0 / d, v0 + (vs - v0) * d0 / d
lam = acc / n
print("Lyapunov per forcing period: mean %.4f (min %.4f max %.4f); per unit time %.4f" % (lam.mean(), lam.min(), lam.max(), lam.mean() / (2*np.pi)))
fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.plot(pts[:,0], pts[:,1], ",", color="k", alpha=0.5)
ax.plot(2.8913632426, 0.2469647907, "r*", ms=9)
ax.set_xlabel("x"); ax.set_ylabel("dx/dt"); ax.set_title("Ueda's attractor: x'' + 0.05x' + x^3 = 7.5 cos t, sampled at t = 2πn\n(red star: saddle D)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "attractor.png"), dpi=170)
