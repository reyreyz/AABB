"""Claim 4 setting (recirc Pe 20, Jacobi w=1, robust hierarchy): N vs MJF (black-box sweeps + algebraic F-Jacobi)."""
import sys, numpy as np
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 20.0); b = np.random.default_rng(0).standard_normal(A.shape[0])
ml = ML(A, omega=1.0, k0=4, k1=8, delta=0.3, tauP=1e-5, trunc_from=1, mode="drop", tau_filter=1e-3, tau_deep=3e-4, refine=1)
for wh in ("N", "markovJF", "markovJ", "twin"):
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=40); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {m} | {wh} | {f:.3g}{'' if conv else ' (no conv.)'} | {nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} |"); sys.stdout.flush()
