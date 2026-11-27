import sys, itertools, numpy as np
from ml_twins import ML
from robust_n import recirc
ms = [int(x) for x in sys.argv[1].split(",")]
cfgs = [(int(c.split(":")[0]), float(c.split(":")[1]), float(c.split(":")[2]), int(c.split(":")[3]), c.split(":")[4]) for c in sys.argv[2].split(",")]
print("| m | n | k1 | tauP_deep | filter | refine | coarse S | levels | op.cx | N factor | cycles | WU/cycle | WU total |\n|" + "---|" * 13)
for m in ms:
    A = recirc(m, 10.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
    h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
    for k1, tp, tf, rf, cs in cfgs:
        ml = ML(A, omega=1.0, omega_c=0.8, wc_levels=1 if cs == "drop" else 99, k0=1, k1=k1, delta=0.3, mode="lin", xy=xy,
                nkeep=2, lin_levels=1 if cs == "drop" else 99, tauP_deep=tp, tau_filter=tf, tau_deep=tf, refine=rf)
        hist = ml.solve(b, "N", tol=1e-8, maxcyc=80); nc = len(hist) - 1; conv = hist[-1] < 1e-8
        f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
        print(f"| {m} | {n} | {k1} | {tp:g} | {tf:g} | {rf} | {cs} | {len(ml.levels)} | {ml.opcx()[0]:.2f} | {f:.3g} | {nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} |"); sys.stdout.flush()
