# ueda-chaos: figure captions and sources

Both figures are drawn by `figures/make_figures.py` (run: `python figures/make_figures.py .` from the repo root;
Python 3.11, numpy, matplotlib 3.11, about 15 s). No number is typed into a plot.

## 1. `figures/ueda-proof-boxes-on-attractor.svg` (hero)

**Point.** The whole proof rests on five tiny boxes placed on the famous picture of Ueda's chaos.

**Caption.** Ueda's oscillator x'' + 0.05 x' + x^3 = 7.5 cos t, seen once per push of the forcing (position across,
speed up). Grey: the familiar picture of Ueda's attractor, recomputed here by simulation (an illustration, not part of
the proof). Coloured marks: where the five proof boxes sit; the proved route is A to A (stay one push), or A to C0, C1,
C2, C3 and back to A (a five-push loop). The right panel zooms in on A and C0, drawn at true size exactly as defined in
`proof/design_chain.txt`: A is a long thin strip, C0 a short sliver. The theorem does not say that the attractor
contains these boxes or the chaotic set.

**Alt text.** The swirling, layered crescent of Ueda's attractor in grey, with five small coloured markers labelled A and
C0 (together, near the middle), C1, C2 (upper) and C3 (lower). A zoom panel shows box A as a long blue strip and C0 as a
short orange sliver beside it.

**Sources.**
| element | file | field |
|---|---|---|
| k = 1/20 = 0.05, B = 15/2 = 7.5 | `independent/SPEC.md` | "k = 1/20 and B = 15/2" (parsed) |
| box positions and outlines | `proof/design_chain.txt` | `p` (curve polynomial), `chart` (type, Q, U, S), `set` (c, w, b); outline as in `design/fig_proof.py` |
| route A ⇒ A, A ⇒ C0 ⇒ C1 ⇒ C2 ⇒ C3 ⇒ A | `proof/design_chain.txt` | `trans` lines (script asserts 5 sets, 6 links) |
| grey attractor | computed | 4,000 starts (seed 1), 50 transient pushes, 12 recorded; fixed-step RK4, 600 steps per push |
| "not claimed: attractor contains the set" | README "What we do not claim" | text |

## 2. `figures/ueda-repeating-routes.svg`

**Point.** The proof gives infinitely many repeating motions, and their number grows by the plastic number per push.

**Caption.** How many different box-routes close up after n pushes, counted exactly from the proof's route map (log
scale). For every such route the theorem gives an exact solution that repeats after n pushes, and because the five
boxes do not overlap, different routes give different starting states: 1 for n = 1 (staying in A), 6 for n = 5, 16
for n = 10, 276 for n = 20 and 76,724 for n = 40. In the long run the count grows by a factor of about 1.3247 per
push, the plastic number (the real root of x^3 = x + 1); its logarithm, about 0.2812, is the proved lower bound on
the topological entropy.

**Alt text.** Dot chart on a log scale: the number of closed routes is 1 for n = 1 to 4, jumps to 6 at n = 5 and then
rises almost along a straight line to 76,724 at n = 40.

**Sources.**
| element | file | field |
|---|---|---|
| route map | `proof/design_chain.txt` | `trans` lines |
| counts 1, 6, 16, 276, 76,724 | computed | trace of G^n, exact integers (script prints them) |
| 1.3247, 0.2812 | computed; matches README "How fast the patterns multiply" | spectral radius of G (script asserts r^3 = r + 1); ln r = 0.281200 |

Dropped: a margin chart of the 18 proved margins (`independent/run_B.log`): an audit table, already in TECHNICAL.md.
Dropped: drawing the guide loop as lines between strobe points (the old `proof_chain.png` does this); straight segments
between strobe points are not paths of the motion and could be read as such.
