import sys, numpy as np
from ml_twins import ML
from flows import cd
m, pe, th = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
A, xy = cd(m, "const", pe); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
coarse = dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)
for wf, wc, k in ((0.15, 0.03, 32), (1.0, 0.03, 16), (1.0, 0.03, 32), (1.0, 0.02, 32), (0.25, 0.05, 32), (1.0, 0.05, 16), (1.0, 0.05, 32)):
    ml = ML(A, omega=wf, omega_c=wc, wc_levels=1, k0=k, split0=rb, mode="lin", xy=xy, theta=th, lin_levels=1, **coarse)
    hist = ml.solve(b, "markovJ", tol=1e-8, maxcyc=120); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {pe:g} | {m} | MJ wF {wf} wC {wc} k {k} | {f:.3g} | {nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} |"); sys.stdout.flush()
