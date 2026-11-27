"""Multilevel: finest level = Claim 10 (k=1, w_C=0.8, 2-entry linear-exact S~), coarse levels = robust N."""
import sys, time, numpy as np
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
t0 = time.time()
ml = ML(A, omega=1.0, omega_c=0.8, wc_levels=1, k0=1, k1=8, delta=0.3, mode="lin", xy=xy, nkeep=2,
        lin_levels=1, tauP_deep=1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=1)
ts = time.time() - t0
oc = ml.opcx()
for wh in ("N", "twin", "markovJ"):
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=40)
    nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {m} | {n} | {len(ml.levels)} | {wh} | {f:.3g} | {nc if conv else 'no conv. (' + format(hist[-1], '.1e') + ')'} | "
          f"{ml.wpc:.0f} | {ml.wpc * nc if conv else float('inf'):.0f} | {oc[0]:.2f} | {ts:.0f} |"); sys.stdout.flush()
