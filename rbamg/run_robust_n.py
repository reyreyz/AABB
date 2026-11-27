"""Robust non-Markovian multilevel cycle (Algorithm 7.3 family) on the section 7.3
problem: centered recirculating flow, cell Peclet 20 (Jacobi and the Markovian
cycle diverge).  Greedy F-points (delta = 0.3), Jacobi omega = 1, P0, injection,
exact P' = P0 + Q0 sum_{l<k} K^l H, Petrov-Galerkin coarse operators filtered
(tau = 1e-3), one refinement step on the finest levels.
Sweep of k (finest / coarse) and refinement depth; total work to 1e-8.
"""

import sys
import time
import numpy as np
import ml_z73 as mzz
from run_base_nz import problem
from run_ml_nz import rate

B = dict(k=4, k_coarse=8, refine_levels=1, tau_deep=3e-5)
T = dict(B, tauP=1e-5)
T3 = dict(T, tau_deep=3e-4)
T3 = dict(T, tau_deep=3e-4)
T3 = dict(T, tau_deep=3e-4)
T3 = dict(T, tau_deep=3e-4)
T3 = dict(T, tau_deep=3e-4)
CONFIGS = [
    ("P'^T, exact solve at 16 sqrt(n)", dict(T3, restrict="ppt", coarsest=-16)),
    ("P'^T, recursive, refine 0-1", dict(T3, restrict="ppt")),
    ("P'^T, recursive, 1 extra cycle at 16 sqrt(n) level", dict(T3, restrict="ppt", refine_small=-16, n_refine_small=1)),
    ("P'^T, recursive, 3 extra cycles at 16 sqrt(n) level", dict(T3, restrict="ppt", refine_small=-16, n_refine_small=3)),
    ("P'^T, recursive, 7 extra cycles at 16 sqrt(n) level", dict(T3, restrict="ppt", refine_small=-16, n_refine_small=7)),
]


def main(ms=(63, 127, 255), only=None):
    print("# Robust N: section 7.3 recirculating flow, Pe 20\n")
    print("| m | n | configuration | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] | solve [s] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = problem(m)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for lab, kw in CONFIGS:
            if only and lab not in only:
                continue
            t = time.time()
            kw = dict(kw)
            if kw.get('coarsest', 0) and kw['coarsest'] < 0:
                kw['coarsest'] = int(-kw['coarsest'] * np.sqrt(A.shape[0]))
            if kw.get('refine_small', 0) and kw['refine_small'] < 0:
                kw['refine_small'] = int(-kw['refine_small'] * np.sqrt(A.shape[0]))
            c = mzz.Alg73Z(A, zmode="Q0", delta=0.3, **kw)
            ts = time.time() - t
            t = time.time()
            h, w = c.solve(b, tol=1e-8, maxcyc=30)
            tsol = time.time() - t
            conv = h[-1] < 1e-8
            print(f"| {m} | {A.shape[0]} | {lab} | {len(c.levels)} | {rate(h):.3g} | {len(h)-1}{'' if conv else ' (no conv.)'} | "
                  f"{w[-1]/max(len(h)-1,1):.0f} | {w[-1] if conv else float('inf'):.0f} | {c.op_complexity():.2f} | "
                  f"{c.stored_complexity():.1f} | {ts:.1f} | {tsol:.1f} |")
            sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (63, 127, 255)
    only = sys.argv[2].split(";") if len(sys.argv) > 2 else None
    main(ms, only)
