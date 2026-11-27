"""Multilevel N, constant flow cell Pe 12: level 0 red-black, k, w_F=1, w_C, blended S~ (theta); cheap coarse levels."""
import sys, time, numpy as np
from ml_twins import ML
from flows import cd
m, pe, th, k, wc = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
cfg = sys.argv[6] if len(sys.argv) > 6 else "cheap"; methods = sys.argv[7].split(",") if len(sys.argv) > 7 else ["N", "twin", "markovJ"]
A, xy = cd(m, "const", pe); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
coarse = {"robust": dict(k1=8, tauP_deep=1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=0),
          "mid": dict(k1=4, tauP_deep=1e-3, tau_filter=3e-3, tau_deep=3e-3, refine=0),
          "cheap": dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0),
          "cheap-r": dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=1)}[cfg]
t0 = time.time()
ml = ML(A, omega=1.0, omega_c=wc, wc_levels=1, k0=k, split0=rb, mode="lin", xy=xy, theta=th, lin_levels=1, **coarse)
ts = time.time() - t0
for wh in methods:
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=80); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {pe:g} | {m} | {cfg} | {len(ml.levels)} | {wh} | {f:.3g}{'' if conv else ' (no conv.)'} | {nc if conv else '-'} | {ml.wpc:.0f} | "
          f"{ml.wpc*nc if conv else float('inf'):.0f} | {ml.opcx()[0]:.2f} | {ts:.0f} |"); sys.stdout.flush()
