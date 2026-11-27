"""Strong test: multilevel N vs every damped Markovian cycle on the same hierarchy (WU to 1e-8)."""
import sys, itertools, numpy as np
from ml_twins import ML
from flows import cd
m, pe, th, k, wc = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
A, xy = cd(m, "const", pe); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
coarse = dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)
def solve(ml, wh, maxcyc=150):
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=maxcyc); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    return f, nc, conv, ml.wpc * nc if conv else float("inf"), ml.wpc
ml = ML(A, omega=1.0, omega_c=wc, wc_levels=1, k0=k, split0=rb, mode="lin", xy=xy, theta=th, lin_levels=1, **coarse)
fN = solve(ml, "N")
print(f"m={m} Pe={pe} theta={th}: N (k={k}, wC={wc}): factor {fN[0]:.3g}, cycles {fN[1]}, WU {fN[3]:.0f} ({fN[4]:.0f}/cycle), op.cx {ml.opcx()[0]:.2f}")
res = []
for wf, r, kk in itertools.product((0.03, 0.06, 0.1, 0.15, 0.25, 0.4, 0.7, 1.0), (0.0, 0.2, 0.5, 1.0), (4, 8, 16, 32)):
    mm = ML(A, omega=wf, omega_c=r * wf, wc_levels=1, k0=kk, split0=rb, mode="lin", xy=xy, theta=th, lin_levels=1, **coarse)
    f = solve(mm, "markovJ", 80); res.append((f[3], f[0], "MJ", wf, r * wf, kk))
    if r == 0.0:
        f = solve(mm, "twin", 80); res.append((f[3], f[0], "twin", wf, "-", kk))
for wcc, kk in itertools.product((0.01, 0.02, 0.03, 0.05, 0.1), (8, 16, 32)):
    mm = ML(A, omega=1.0, omega_c=wcc, wc_levels=1, k0=kk, split0=rb, mode="lin", xy=xy, theta=th, lin_levels=1, **coarse)
    f = solve(mm, "markovJ", 80); res.append((f[3], f[0], "MJ", 1.0, wcc, kk))
res.sort(key=lambda x: (x[0], x[1]))
print(f"best Markovian cycles on the same hierarchy (WU total, factor, method, w_F, w_C, k), of {len(res)}:")
for x in res[:6]: print("  ", (round(x[0]) if np.isfinite(x[0]) else "inf", round(float(x[1]), 3)) + x[2:])
