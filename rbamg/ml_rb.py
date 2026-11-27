"""Idea 13: multilevel N, level 0 = red-black, k=2, w_F=1, w_C=0.2, single-entry lumped S~ (exact on 1);
coarse levels = greedy robust N.  Constant flow, cell Pe 6."""
import sys, time, numpy as np
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); cfg = sys.argv[2]; methods = sys.argv[3].split(",") if len(sys.argv) > 3 else ["N", "twin", "markovJ"]
pe = float(sys.argv[4]) if len(sys.argv) > 4 else 6.0
A, xy = cd(m, "const", pe); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
I, J = np.divmod(np.arange(n), m); isF = (I + J) % 2 == 1
coarse = {"robust": dict(k1=8, tauP_deep=1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=0),
          "mid": dict(k1=4, tauP_deep=1e-3, tau_filter=3e-3, tau_deep=3e-3, refine=0),
          "cheap": dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)}[cfg]
t0 = time.time()
ml = ML(A, omega=1.0, omega_c=0.2, wc_levels=1, k0=2, split0=isF, delta=0.3, mode="lin", xy=xy, nkeep=1, lin_levels=1, **coarse)
ts = time.time() - t0
levnnz = " ".join(f"{L.A.nnz/A.nnz:.2f}" for L in ml.levels)
for wh in methods:
    t0 = time.time(); hist = ml.solve(b, wh, tol=1e-8, maxcyc=60); tsol = time.time() - t0
    nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {m} | {n} | {cfg} | {len(ml.levels)} | {wh} | {f:.3g}{'' if conv else ' (no conv., ' + format(hist[-1], '.0e') + ')'} | "
          f"{nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} | {ml.opcx()[0]:.2f} | {ts:.0f} | {tsol/nc:.2f} | {levnnz} |"); sys.stdout.flush()
