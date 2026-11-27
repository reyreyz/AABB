"""N (truncated, multilevel) == stationary Markovian cycle with M^ = M_twin + Q(S - S~)P^# sum_j T^j M."""
import sys, numpy as np
from ml_twins import ML
from robust_n import recirc
for m in (63, 127):
    A = recirc(m, 10.0); n = A.shape[0]; rng = np.random.default_rng(3)
    h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
    for lab, kw in (("Claim 11 hierarchy (w_C=0.6)", dict(omega=1.0, omega_c=0.6, wc_levels=1, k0=1, k1=8, delta=0.3, mode="lin", xy=xy,
                     nkeep=2, lin_levels=1, tauP_deep=1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=1)),
                    ("Pe 10, lump 0.5 all levels, k=4", dict(omega=1.0, k0=4, k1=4, tauP=0.5, mode="lump", refine=1))):
        ml = ML(A, **kw); b = rng.standard_normal(n)
        xN = xH = np.zeros(n); diffs = []
        for c in range(6):
            xN = ml.cycle(0, b, xN, "N"); xH = ml.cycle(0, b, xH, "hat")
            diffs.append(np.linalg.norm(xN - xH) / max(np.linalg.norm(xN), 1e-300))
        print(f"| {m} | {lab} | {len(ml.levels)} | {max(diffs):.1e} | {np.linalg.norm(b - A @ xN)/np.linalg.norm(b):.1e} |")
