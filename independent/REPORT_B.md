# Checker B report

**In plain words.** The six covering relations in `design_chain.txt`, and the pairwise disjointness of the five h-sets, all hold under interval arithmetic and a validated integrator. The run is `ALL VERIFIED`. Both required negative controls come back `NOT VERIFIED`.

## What was proved

Ueda's equation \(x'' + k x' + x^3 = B \cos t\) with \(k = 1/20\) and \(B = 15/2\). The map \(P\) is the time-\(2\pi\) stroboscopic map from phase 0. For each `trans SRC DST o` the program proves (E−), (E+), and (S) from `SPEC.md`, and it proves (D) for every pair of h-sets.

Every accepted inequality is strict and is witnessed by an outward enclosure of the image sitting strictly on the inner side of the threshold. A threshold such as \(c - w\) is an interval around the exact decimal; the test uses the inward endpoint, so the comparison implies the inequality for the exact threshold.

## Method and choices

The verifier is `verify_b.cpp`. It uses CAPD with filib, Taylor order 12, and the C0 doubleton (`C0Rect2Set`, Lohner / QR). Time dependence is CAPD's `time:t` parameter in the vector field, with the state kept in \((x, v)\). The step controller is automatic, with each step at most \(1/5\).

That combination is deliberate. The specification warns that an extra time state plus a high Taylor order makes the step controller blow up on this cubic, and that a Lagrange remainder has to be taken on the whole step hull. CAPD's C0 step builds the remainder from the rough enclosure of the step (`enclosure`, then `computeRemainderCoefficients` on that hull). Order 12 sits in the range the specification calls reliable. No threads are used.

A parameter rectangle in a chart is sent to phase space as a doubleton. An affine chart is exactly the parallelogram \(Q + u U + s S\). The curved chart uses the Taylor expansion of \(p\) in \(s\), with the Lagrange remainder in \(p''\) carried as an interval and added into the center, so the doubleton contains the curved piece. After integration, \((a, s) = [U\ S]^{-1}(X - Q)\) is evaluated on the doubleton (center, \(C r_0\), and \(B r\) together). On the curved target, \(u = a - p(s)\).

The integrator is asked for the time interval \(2 \cdot \mathtt{interval::pi()}\), which contains the real \(2\pi\). A box is accepted only when the set's time interval contains that enclosure. The returned set encloses the solution throughout that time interval, hence at time exactly \(2\pi\).

Subdivision starts at parameter width \(0.01\) and splits any box whose enclosure still crosses the threshold, or whose integration throws. A box whose enclosure lies entirely on the wrong side refutes the statement and stops that goal. Nothing in the successful run needed a retry: `integrator_retries=0` on every goal. Decimal inputs are parsed by filib's outward string constructor, so each coordinate interval contains the exact printed decimal.

Disjointness does not integrate. An affine h-set is a parallelogram. A curved h-set is covered by a parallelogram in the \((a, s)\) coordinates with \(a \in [c-w, c+w] + p([-b, b])\), which contains the curved strip; \(s\) is chopped until the pieces separate if the fat parallelogram does not. Two parallelograms are declared disjoint when their projections onto a tested direction have disjoint outward enclosures. Directions tried: the coordinate axes, and the edge normals.

## Results

Log: `run_b.log`. 205 integrations (one of them the point check below). Timed section: **1.553 seconds**. Final line: **ALL VERIFIED**.

Point check, not part of the covering statement: \(P\) of the chart-A origin returns
\(u \in [1.493 \times 10^{-9},\ 1.506 \times 10^{-9}]\),
\(s \in [-8.45 \times 10^{-11},\ -8.41 \times 10^{-11}]\)
in chart A. The anchor is a periodic point up to the printed truncation of \(Q\).

Minimum slack on each goal (outward image bound versus inward threshold). The smallest slack in the whole chain is the C1 → C2 strip, \(2.87 \times 10^{-4}\).

| link | orient | (E−) slack | (E+) slack | (S) slack | boxes (E− / E+ / S) |
|---|---:|---:|---:|---:|---|
| A → A | +1 | 1.79e-2 | 1.80e-2 | 1.41e-1 | 31 / 31 / 31 |
| A → C0 | +1 | 3.92e-3 | 3.48e-2 | 3.20e-3 | 31 / 31 / 31 |
| C0 → C1 | +1 | 1.85e-3 | 1.35e-3 | 2.00e-3 | 3 / 3 / 3 |
| C1 → C2 | +1 | 1.58e-3 | 1.58e-3 | 2.87e-4 | 1 / 1 / 1 |
| C2 → C3 | +1 | 4.79e-4 | 4.79e-4 | 1.39e-3 | 1 / 1 / 1 |
| C3 → A | −1 | 3.29e-3 | 3.28e-3 | 3.28e-2 | 1 / 1 / 1 |

Per-transition times in the log: 0.643 s, 0.638 s, 0.087 s, 0.050 s, 0.048 s, 0.049 s.

(D), all ten pairs:

| pair | separated by | gap |
|---|---|---:|
| A vs C0 | normal to A's \(S\) | 2.68e-3 |
| A vs C1 | \(v\)-bounding boxes | 5.01e-2 |
| A vs C2 | \(x\)-bounding boxes | 5.92e-3 |
| A vs C3 | \(x\)-bounding boxes | 5.36e-2 |
| C0 vs C1 | \(x\)-bounding boxes | 1.09e-2 |
| C0 vs C2 | \(x\)-bounding boxes | 2.47e-2 |
| C0 vs C3 | \(x\)-bounding boxes | 6.95e-2 |
| C1 vs C2 | \(x\)-bounding boxes | 3.86e-2 |
| C1 vs C3 | \(x\)-bounding boxes | 5.78e-2 |
| C2 vs C3 | \(x\)-bounding boxes | 9.72e-2 |

The A–C0 gap is the projection onto \((-S_y, S_x)\), a vector of length 1. In the chart-A coordinate \(a\) the same separation is about \(0.013\), because \(U \cdot n \approx |\det| \approx 0.199\).

## Negative controls

Both end with `NOT VERIFIED`.

1. Last orientation flipped (`trans C3 A 1`). Log `run_b_flip.log`. **1.807 seconds**, 205 integrations, exit status 2. Transitions 1–5 still verify. On C3 → A with the wrong sign, (E−) and (E+) are refuted by a single box each: the left-edge image has \(u \in [0.00489,\ 0.00500]\), which is entirely above \(c + w\), so it fails \(u < c - w\). The right edge is entirely below \(c - w\). The strip (S) still holds; orientation does not enter (S). Disjointness is unchanged.

2. \(B = 7.4\) instead of \(7.5\), enclosed as the exact decimal \([7.3999999999999995,\ 7.4000000000000012]\). Log `run_b_B74.log`. **0.516 seconds**, 81 integrations, exit status 2. The chart-A origin is sent to \(u \approx -0.867\), \(s \approx 0.319\), far from the h-set. Every transition has at least one refuted goal (right exit and strip on the first five; left exit and strip on C3 → A). Disjointness still verifies, because (D) does not use \(B\).

## What looked suspicious

- The file line `m 1` is not defined in `SPEC.md`. The header of the chain says each link is one period of \(P\). This checker treats `m` as the number of periods and refuses any value other than 1.
- On every chart, \(U\) and \(S\) are unit vectors with dot product about \(\pm 0.98\) and \(|\det|\) about \(0.20\) to \(0.37\). The frames are invertible (0 is outside every determinant interval) and strongly skewed. The sets are long thin regions in the plane. This is a property of the design data. It does not break the charts, and the disjointness proof accounts for it. A and C0 are the only pair whose axis-aligned boxes meet; they separate on the normal to A's \(S\).
- The anchor \(Q\) of chart A is printed with fewer digits than the other anchors. The validated return map of that point misses the origin by about \(1.5 \times 10^{-9}\) in \(u\). That is consistent with a truncated printing of a periodic point, and it is far inside the h-set half-width \(1.6 \times 10^{-3}\).
- The proved slacks are at least \(2.87 \times 10^{-4}\). The data intervals are about one unit in the last place wide. The proof does not sit on a rounding boundary.
- No integration in the successful run hit the step-size blow-up described in the specification. The initial parameter width \(0.01\) was already fine enough for every goal.

## Limits

The dynamical part trusts CAPD's C0 Taylor integrator at order 12, which is the tool the specification allows. The geometric reductions (chart inverse, polynomial Taylor remainder, separating projections) are in this file and were checked against the same log: the point return is near the origin, the C3 → A exit under the flipped orientation is refuted on the nose, and \(B = 7.4\) moves the images cleanly off the sets.
