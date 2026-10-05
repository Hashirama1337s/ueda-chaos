# Checker B specification (blind): verify a chain of covering relations for Ueda's map

You are writing an INDEPENDENT rigorous verifier. Do not look for, open or copy any other verifier code (in particular
nothing under ../proof/). Work only from this specification and design_chain.txt in this folder.

## 1. The map
- Ueda's equation: x'' + k x' + x^3 = B cos t, with k = 1/20 and B = 15/2 (exactly; use interval enclosures of 1/20 and 2*pi).
- State (x, v) with v = x'.
- P: R^2 -> R^2 is the time-2*pi map started at t = 0 (stroboscopic map at phase 0).

## 2. Charts (all coordinates below are doubles read from design_chain.txt; treat them as exact numbers)
- Polynomial p(s) = sum_j p_j s^(n-1-j), coefficients listed high-to-low on the "p" line. Used only by the curved chart.
- `chart NAME curved Qx Qv Ux Uv Sx Sv`: X(u, s) = Q + (u + p(s)) U + s S.
- `chart NAME affine Qx Qv Ux Uv Sx Sv`: X(u, s) = Q + u U + s S.
- Inverse of a chart:
  - (a, s) = [U S]^(-1) (X - Q), where [U S] is the 2x2 matrix with columns U and S;
  - u = a - p(s) (curved) or u = a (affine).
  - Each chart is a homeomorphism of R^2.

## 3. h-sets
- `set NAME c w b`: N = {X(u, s) : |u - c| <= w, |s| <= b} in chart NAME.
- Exit direction u, entry direction s.
- Left exit edge u = c - w; right exit edge u = c + w (each for |s| <= b).

## 4. What to prove
For every line `trans SRC DST o` (o = +1 or -1), with N = SRC and M = DST in their own charts:
- **(E-)** for every point of N's left edge, P(point) in M's chart has u < cM - wM if o = +1, or u > cM + wM if o = -1.
- **(E+)** for every point of N's right edge, the opposite inequality (u > cM + wM if o = +1, u < cM - wM if o = -1).
- **(S)** for every point of N, P(point) in M's chart has |s| < bM.

Also prove:
- **(D)** all h-sets are pairwise disjoint.

Every inequality must be established with rigorous interval arithmetic (outward rounding) and a validated ODE
integrator, for ALL points of the edges/sets (cover them by finitely many boxes and enclose each box's image).

## 5. Tools
You may use CAPD (installed in WSL at ~/capd_src, built library in ~/capd_src/build).
- Example build line: g++ -I$HOME/capd_src/capdDynSys/include -I$HOME/capd_src/capdAlg/include -I$HOME/capd_src/capdAux/include
  -I$HOME/capd_src/capdExt/include -I$HOME/capd_src/capdExt/filibsrc -std=c++17 -O2 -frounding-math -D__USE_FILIB__
  -DFILIB_EXTENDED -DFILIB_HAVE_SSE file.cpp -L$HOME/capd_src/build -L$HOME/capd_src/build/capdExt/filibsrc -lcapd -lfilib -o out
- Or any other rigorous method you implement yourself (e.g. your own interval Taylor integrator).

**Independence matters more than speed.** Make your own choices for:
- how you handle the time dependence;
- the set representation;
- subdivision strategy;
- the inverse-chart computation.

**Known CAPD pitfalls:**
- Do NOT use multithreading with CAPD (not thread-safe).
- With the explicit time variable and high Taylor orders, automatic step control can blow up erratically on this stiff
  cubic. Moderate orders (about 10-14) are reliable.
- A Lagrange remainder must be evaluated on the whole interval (hull), not at a point.

## 6. Deliverables (write in this folder only)
- **verify_b.cpp** (or .py): the verifier.
- **run_b.log:** the full output of a run on design_chain.txt, with per-transition and per-edge results and a final line
  ALL VERIFIED / NOT VERIFIED.
- **Negative controls, each must report NOT VERIFIED:**
  1. design_chain.txt with the last transition's orientation flipped (`trans C3 A 1`);
  2. B = 7.4 instead of 7.5.
- **REPORT_B.md:** method, choices, results, run times, anything that looked suspicious in the design or this spec.

Stop after these deliverables; do not modify anything outside this folder.
