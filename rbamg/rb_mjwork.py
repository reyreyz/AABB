import sys, numpy as np
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); A, xy = cd(m, "const", 6.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
I, J = np.divmod(np.arange(n), m); isF = (I + J) % 2 == 1
coarse = dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)
for lab, wf, wc, k, wh in (("N (wF 1, wC 0.2, k 2)", 1.0, 0.2, 2, "N"), ("MJ wF 0.8 wC 0.16 k 8", 0.8, 0.16, 8, "markovJ"),
                           ("MJ wF 0.6 wC 0.12 k 8", 0.6, 0.12, 8, "markovJ"), ("MJ wF 0.8 wC 0.16 k 16", 0.8, 0.16, 16, "markovJ"),
                           ("MJ wF 0.6 wC 0.12 k 16", 0.6, 0.12, 16, "markovJ")):
    ml = ML(A, omega=wf, omega_c=wc, wc_levels=1, k0=k, split0=isF, mode="lin", xy=xy, nkeep=1, lin_levels=1, **coarse)
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=200); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {m} | {lab} | {f:.3g} | {nc if conv else 'no conv.'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} |"); sys.stdout.flush()
