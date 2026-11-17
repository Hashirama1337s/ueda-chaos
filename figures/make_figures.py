"""Figures for the ueda-chaos repository.

Usage:  python figures/make_figures.py <repo root>
Writes figures/*.svg next to this script. Inputs, all shipped in the repository:
  proof/design_chain.txt    the five proof boxes (charts, the curve polynomial p, sizes) and the six links
  independent/SPEC.md       the equation's parameters k = 1/20, B = 15/2 (parsed, not typed)
The attractor dots are a float illustration recomputed here from the equation (fixed-step RK4, fixed seed);
the boxes are drawn exactly as defined in the design file. Nothing here is part of the proof.
Checked copy (figure referee, 2026-11-17): the hero subtitle marks the grey dots as simulation, not proved.
"""
import os
import re
import sys
from fractions import Fraction

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon

BG, BORDER, GRID, AXC, TXT, ACC = "#1c2026", "#2f353d", "#22272e", "#9aa4ae", "#e8edf2", "#39ff14"
CAT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c", "#9085e9", "#e66767"]
NEUTRAL = "#5b636c"


def style():
    plt.rcParams.update({
        "svg.fonttype": "none", "svg.hashsalt": "ueda-chaos-figures",
        "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14, "xtick.labelsize": 12.5, "ytick.labelsize": 12.5,
        "legend.fontsize": 12.5, "figure.facecolor": "none", "axes.facecolor": BG, "axes.edgecolor": AXC,
        "axes.labelcolor": TXT, "axes.titlecolor": TXT, "xtick.color": AXC, "ytick.color": AXC, "text.color": TXT,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.spines.top": False,
        "axes.spines.right": False, "legend.facecolor": BG, "legend.edgecolor": BORDER, "legend.labelcolor": TXT,
        "legend.framealpha": 1.0, "lines.linewidth": 2.0,
    })


def frame(fig):
    w, h = fig.get_size_inches()
    fig.patches.append(FancyBboxPatch((0.03, 0.03), w - 0.06, h - 0.06, boxstyle="round,pad=0,rounding_size=0.18",
                                      transform=fig.dpi_scale_trans, fc=BG, ec=BORDER, lw=1.5, zorder=-100))


def save(fig, name):
    here = os.path.dirname(os.path.abspath(__file__))
    out = here if os.path.basename(here) == "figures" else os.path.join(here, "figures")   # repo: figures/make_figures.py
    os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, name), format="svg", dpi=200, metadata={"Date": None, "Creator": None})
    plt.close(fig)
    print("wrote figures/" + name)


def read_design(root):
    L = [l.split() for l in open(os.path.join(root, "proof", "design_chain.txt"), encoding="utf-8")
         if l.strip() and not l.startswith("#")]
    p = [float(x) for x in [t for t in L if t[0] == "p"][0][2:]]
    charts = {t[1]: (t[2], np.array(list(map(float, t[3:9])))) for t in L if t[0] == "chart"}
    sets = {t[1]: tuple(map(float, t[2:5])) for t in L if t[0] == "set"}
    trans = [(t[1], t[2]) for t in L if t[0] == "trans"]
    return p, charts, sets, trans


def read_params(root):
    s = open(os.path.join(root, "independent", "SPEC.md"), encoding="utf-8").read()
    g = re.search(r"k = (\d+/\d+) and B = (\d+/\d+)", s)
    return float(Fraction(g.group(1))), float(Fraction(g.group(2)))


def box_outline(p, chart, st, n=120):
    """Outline of an h-set (as in design/fig_proof.py): affine X = Q + uU + sS, curved X = Q + (u + p(s))U + sS."""
    ty, c = chart
    cc, w, b = st

    def X(u, s):
        a = u + (np.polyval(p, s) if ty == "curved" else 0.0)
        return np.c_[c[0] + a * c[2] + s * c[4], c[1] + a * c[3] + s * c[5]]

    s = np.linspace(-b, b, n)
    u = np.linspace(cc - w, cc + w, n)
    return np.vstack([X(np.full(n, cc - w), s), X(u, np.full(n, b)), X(np.full(n, cc + w), s[::-1]),
                      X(u[::-1], np.full(n, -b))]), X(np.array([cc]), np.array([0.0]))[0]


def strobe(x, v, k, B, periods, nst=600):
    """time-2pi map of x'' + k x' + x^3 = B cos t, vectorised fixed-step RK4 (float illustration only)."""
    h = 2 * np.pi / nst
    for _ in range(periods):
        t = 0.0
        for _ in range(nst):
            f = lambda t, x, v: (v, -k * v - x ** 3 + B * np.cos(t))
            a1 = f(t, x, v)
            a2 = f(t + h / 2, x + h / 2 * a1[0], v + h / 2 * a1[1])
            a3 = f(t + h / 2, x + h / 2 * a2[0], v + h / 2 * a2[1])
            a4 = f(t + h, x + h * a3[0], v + h * a3[1])
            x = x + h / 6 * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0])
            v = v + h / 6 * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
            t += h
    return x, v


def fig_hero(root, p, charts, sets):
    k, B = read_params(root)
    rng = np.random.default_rng(1)
    x, v = rng.uniform(1.6, 3.3, 4000), rng.uniform(-4.0, 5.0, 4000)
    x, v = strobe(x, v, k, B, 50)
    pts = []
    for _ in range(12):
        x, v = strobe(x, v, k, B, 1)
        pts.append(np.c_[x, v])
    pts = np.vstack(pts)
    names = ["A", "C0", "C1", "C2", "C3"]
    col = {n: CAT[i] for i, n in enumerate(names)}
    mk = {"A": "D", "C0": "s", "C1": "s", "C2": "s", "C3": "s"}
    outl = {n: box_outline(p, charts[n], sets[n]) for n in names}

    fig = plt.figure(figsize=(12.5, 7.0))
    frame(fig)
    ax = fig.add_axes([0.075, 0.10, 0.50, 0.70])
    az = fig.add_axes([0.665, 0.10, 0.31, 0.70])
    fig.text(0.03, 0.94, "Five tiny boxes, drawn over Ueda's attractor, carry a proof of chaos", fontsize=17, weight="bold")
    fig.text(0.03, 0.845, "One grey dot per push of the forcing (strobe view), x'' + %g x' + x^3 = %g cos t: simulation, "
             "not proved.\nColoured marks: the proof boxes. Proved route: A -> A (stay), or A -> C0 -> C1 -> C2 -> C3 -> A."
             % (k, B),
             fontsize=13, color=AXC, linespacing=1.4)

    ax.plot(pts[:, 0], pts[:, 1], ".", ms=1.2, color=AXC, alpha=0.35, rasterized=True, zorder=1)
    for n in names:
        ax.plot(*outl[n][1], mk[n], ms=9, color=col[n], mec=BG, mew=1.0, zorder=4)
    lab = {"A": (-14, -4, "right", "A and C0 (zoom, right)"), "C1": (12, 0, "left", "C1"), "C2": (12, 0, "left", "C2"),
           "C3": (12, 0, "left", "C3")}
    for n, (dx, dy, ha, s) in lab.items():
        ax.annotate(s, outl[n][1], xytext=(dx, dy), textcoords="offset points", ha=ha, va="center",
                    color=col[n] if n != "A" else TXT, fontsize=13, weight="bold", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.2", fc=BG, ec="none", alpha=0.85))
    zx, zy = (2.868, 2.914), (0.08, 0.42)
    ax.add_patch(Rectangle((zx[0], zy[0]), zx[1] - zx[0], zy[1] - zy[0], fill=False, ec=TXT, lw=1.2, zorder=5))
    ax.set_xlabel("position x")
    ax.set_ylabel("speed dx/dt")
    ax.set_xlim(1.45, 3.5)
    ax.set_ylim(-4.6, 5.8)

    az.plot(pts[:, 0], pts[:, 1], ".", ms=2.0, color=AXC, alpha=0.5, rasterized=True, zorder=1)
    for n in ("A", "C0"):
        P = outl[n][0]
        az.add_patch(Polygon(P, closed=True, fc=col[n], ec=col[n], alpha=0.35, lw=0, zorder=3))
        az.add_patch(Polygon(P, closed=True, fill=False, ec=col[n], lw=1.4, zorder=4))
    az.annotate("A", outl["A"][0][40], xytext=(-10, 0), textcoords="offset points", ha="right", va="center",
                color=col["A"], fontsize=13, weight="bold")
    az.annotate("C0", outl["C0"][1], xytext=(10, 6), textcoords="offset points", ha="left", va="bottom",
                color=col["C0"], fontsize=13, weight="bold")
    az.set_xlim(*zx)
    az.set_ylim(*zy)
    az.set_title("zoom: boxes at true size", fontsize=13.5, loc="left")
    az.set_xlabel("position x")
    az.set_ylabel("speed dx/dt")
    az.locator_params(axis="x", nbins=4)
    for s in az.spines.values():
        s.set_visible(True)
        s.set_color(TXT)
    save(fig, "ueda-proof-boxes-on-attractor.svg")


def closed_routes(nodes, trans, nmax):
    idx = {n: i for i, n in enumerate(nodes)}
    G = [[0] * len(nodes) for _ in nodes]
    for a, b in trans:
        G[idx[a]][idx[b]] += 1
    M = [[int(i == j) for j in range(len(nodes))] for i in range(len(nodes))]
    out = []
    for _ in range(nmax):
        M = [[sum(M[i][q] * G[q][j] for q in range(len(nodes))) for j in range(len(nodes))] for i in range(len(nodes))]
        out.append(sum(M[i][i] for i in range(len(nodes))))
    return out, np.array(G, float)


def fig_count(sets, trans):
    N = 40
    cnt, G = closed_routes(list(sets), trans, N)
    r = max(abs(np.linalg.eigvals(G)))
    assert abs(r ** 3 - r - 1) < 1e-9, "growth rate is not the plastic number"
    n = np.arange(1, N + 1)
    fig = plt.figure(figsize=(12.5, 5.6))
    frame(fig)
    ax = fig.add_axes([0.10, 0.14, 0.85, 0.66])
    fig.text(0.03, 0.92, "Infinitely many repeating motions, and the count keeps growing", fontsize=17, weight="bold")
    fig.text(0.03, 0.865, "Different box-routes that close up after n pushes, counted exactly from the route map "
             "(log scale)", fontsize=13, color=AXC)
    ax.semilogy(n, cnt, "o", ms=8, color=ACC, zorder=3)
    for q in (1, 5, 10, 20, 40):
        ax.annotate("{:,}".format(cnt[q - 1]), (q, cnt[q - 1]), xytext=(0, 10), textcoords="offset points",
                    ha="center", fontsize=12.5, color=TXT)
    ax.text(21, 3, "long-run growth: a factor of about %.4f per push\n(the plastic number, the real root of x^3 = x + 1)"
            % r, fontsize=12.5, color=AXC, ha="left", va="center")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: "{:,.0f}".format(v)))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("n (pushes until the motion repeats)")
    ax.set_ylabel("number of routes")
    ax.set_xlim(0, N + 1.5)
    ax.set_ylim(0.6, max(cnt) * 3)
    save(fig, "ueda-repeating-routes.svg")
    return cnt, r


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
    style()
    p, charts, sets, trans = read_design(root)
    assert list(sets) == ["A", "C0", "C1", "C2", "C3"] and len(trans) == 6
    fig_hero(root, p, charts, sets)
    cnt, r = fig_count(sets, trans)
    print("routes n=1,5,10,20,40:", cnt[0], cnt[4], cnt[9], cnt[19], cnt[39], " r=%.6f ln r=%.6f" % (r, np.log(r)))


if __name__ == "__main__":
    main()
