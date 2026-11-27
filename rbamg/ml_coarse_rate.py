"""Rate of the level-1 cycle alone on A_1 = R A P~ (exact, matrix-free), for cheap and robust coarse settings."""
import sys, numpy as np
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
for lab, k1, tp, tf, td, rf in (("cheap k1=2, 1e-2, 1e-2", 2, 1e-2, 1e-2, 1e-2, 0),
                                ("mid k1=4, 1e-3, 3e-3", 4, 1e-3, 3e-3, 3e-3, 0),
                                ("robust k1=8, 1e-5, 1e-3/3e-4", 8, 1e-5, 1e-3, 3e-4, 0)):
    ml = ML(A, omega=1.0, omega_c=0.8, wc_levels=1, k0=1, k1=k1, delta=0.3, mode="lin", xy=xy, nkeep=2,
            lin_levels=1, tauP_deep=tp, tau_filter=tf, tau_deep=td, refine=rf)
    L0 = ml.levels[0]; nc = len(L0.C); rc = np.random.default_rng(1).standard_normal(nc)
    e = np.zeros(nc); hist = [1.0]; ml.work = 0.0
    for it in range(12):
        e = e + ml.cycle(1, rc - ml.RAP(L0, e), np.zeros(nc), "N")
        hist.append(np.linalg.norm(rc - ml.RAP(L0, e)) / np.linalg.norm(rc))
    f = (hist[-1] / hist[2]) ** (1 / 10)
    print(f"| {m} | {lab} | {len(ml.levels)} | {f:.3g} | {ml.work/12:.1f} |"); sys.stdout.flush()
