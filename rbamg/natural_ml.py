"""Multilevel natural examples.  A: constant flow Pe 8, damped Jacobi w=0.1 (inside its stability range), k=16,
greedy split, one-point P~ on level 0; coarse levels robust N (w=1, k=2, S at 1e-2, filter 1e-2).
B: constant flow Pe 6, Jacobi w=0.8 on all levels, greedy, k=2, S truncated at 0.3 on all levels, filter 1e-2."""
import sys, time, numpy as np
from ml_twins import ML
from flows import cd
m, case = int(sys.argv[1]), sys.argv[2]; methods = sys.argv[3].split(",") if len(sys.argv) > 3 else ["N", "twin", "markovJ"]
if case == "A":
    A, xy = cd(m, "const", 8.0)
    kw = dict(omega=0.1, omega_coarse=1.0, k0=16, k1=2, delta=0.3, mode="lin", xy=xy, nkeep=1, lin_levels=1,
              tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)
elif case == "A2":
    A, xy = cd(m, "const", 8.0)
    kw = dict(omega=0.1, k0=16, k1=16, delta=0.3, mode="lin", xy=xy, nkeep=1, lin_levels=1,
              tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)
elif case == "B2":
    A, xy = cd(m, "const", 6.0)
    kw = dict(omega=0.8, k0=2, k1=2, delta=0.3, mode="drop", tauP=0.1, tau_filter=1e-3, tau_deep=1e-3, refine=0)
else:
    A, xy = cd(m, "const", 6.0)
    kw = dict(omega=0.8, k0=2, k1=2, delta=0.3, mode="drop", tauP=0.3, tau_filter=1e-2, tau_deep=1e-2, refine=0)
b = np.random.default_rng(0).standard_normal(A.shape[0]); t0 = time.time(); ml = ML(A, **kw); ts = time.time() - t0
for wh in methods:
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=80); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {case} | {m} | {A.shape[0]} | {len(ml.levels)} | {wh} | {f:.3g}{'' if conv else ' (no conv.)'} | {nc if conv else '-'} | "
          f"{ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} | {ml.opcx()[0]:.2f} | {ts:.0f} |"); sys.stdout.flush()
