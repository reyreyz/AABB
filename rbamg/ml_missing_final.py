"""O(n) multilevel N where twin and Markov-Jacobi diverge.
Level 0: Jacobi w_F=1, w_C=0.8, k=1, S~ = 2 entries/row exact on 1,x,y.  Levels >= 1: robust N
(Jacobi w=1, k=8, S at 1e-5, filter 1e-3 / 3e-4 from level 8, sparse LU below 16 sqrt(n)).  R = P~^T."""
import sys, time, numpy as np
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); refine = int(sys.argv[2]) if len(sys.argv) > 2 else 0
methods = sys.argv[3].split(",") if len(sys.argv) > 3 else ["N", "twin", "markovJ"]
A = recirc(m, 10.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
t0 = time.time()
ml = ML(A, omega=1.0, omega_c=0.8, wc_levels=1, k0=1, k1=8, delta=0.3, mode="lin", xy=xy, nkeep=2,
        lin_levels=1, tauP_deep=1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=refine)
ts = time.time() - t0
levnnz = " ".join(f"{L.A.nnz/A.nnz:.2f}" for L in ml.levels)
for wh in methods:
    t0 = time.time(); hist = ml.solve(b, wh, tol=1e-8, maxcyc=60); tsol = time.time() - t0
    nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {m} | {n} | {refine} | {len(ml.levels)} | {wh} | {f:.3g}{'' if conv else ' (no conv., ' + format(hist[-1], '.0e') + ')'} | "
          f"{nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} | {ml.opcx()[0]:.2f} | {ts:.0f} | {tsol/nc:.2f} | {levnnz} |"); sys.stdout.flush()
