import sys, itertools, numpy as np
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
print(f"m={m}\n| k1 | tauP_deep | filter | w_C coarse | levels | op.cx | N factor | cycles | WU/cycle | WU total |\n|---|---|---|---|---|---|---|---|---|---|")
for k1, tp, tf, wcc in itertools.product((2, 4, 8), (1e-5, 1e-3, 1e-2), (1e-3, 1e-2), (1.0, 0.8)):
    ml = ML(A, omega=1.0, omega_c=0.8, wc_levels=1 if wcc == 1.0 else 99, k0=1, k1=k1, delta=0.3, mode="lin", xy=xy,
            nkeep=2, lin_levels=1, tauP_deep=tp, tau_filter=tf, tau_deep=min(tf, 3e-4) if tf == 1e-3 else tf, refine=1)
    hist = ml.solve(b, "N", tol=1e-8, maxcyc=60); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else float("inf")
    print(f"| {k1} | {tp:g} | {tf:g} | {wcc} | {len(ml.levels)} | {ml.opcx()[0]:.2f} | {f:.3g} | {nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} |"); sys.stdout.flush()
