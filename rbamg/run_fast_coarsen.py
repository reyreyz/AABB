"""Faster coarsening on the coarse levels, with N_Z restoring the CR there.

Section 7.3 problem (recirculating flow, Pe 20).  Finest level: greedy split
delta = 0.3 (A_ff diagonal, K = 0).  Coarse levels: delta_coarse (larger = fewer
C-points).  Bounded settings: k = 4 / 8, S truncated at 1e-5, filter 1e-3 on
levels < 8 and 3e-4 below.  N_Z: bounded Z (<= 4 per column) on the coarse
levels, safeguarded by theta (CR rate of the level pair).
Refinement on levels 0-1 ('r01') or on every level ('rall', W-like).
"""

import sys
import time
import numpy as np
import ml_z73 as mzz
from run_base_nz import problem
from run_ml_nz import rate

B = dict(k=4, k_coarse=8, tauP=1e-5, tau_deep=3e-4, delta=0.3)
CONFIGS = [
    ("N dc=0.3 r01 (bounded base)", dict(B, zmode="Q0", refine_levels=1)),
    ("N dc=0.45 r01", dict(B, zmode="Q0", delta_coarse=0.45, refine_levels=1)),
    ("N_Z dc=0.45 r01", dict(B, q=4, z_from=1, theta="auto", delta_coarse=0.45, refine_levels=1)),
    ("N dc=0.45 rall", dict(B, zmode="Q0", delta_coarse=0.45, refine_levels=99)),
    ("N_Z dc=0.45 rall", dict(B, q=4, z_from=1, theta="auto", delta_coarse=0.45, refine_levels=99)),
    ("N dc=0.6 r01", dict(B, zmode="Q0", delta_coarse=0.6, refine_levels=1)),
    ("N_Z dc=0.6 r01", dict(B, q=4, z_from=1, theta="auto", delta_coarse=0.6, refine_levels=1)),
    ("N dc=0.6 rall", dict(B, zmode="Q0", delta_coarse=0.6, refine_levels=99)),
    ("N_Z dc=0.6 rall", dict(B, q=4, z_from=1, theta="auto", delta_coarse=0.6, refine_levels=99)),
]


def main(m, only=None):
    A = problem(m)
    b = np.random.default_rng(0).standard_normal(A.shape[0])
    for lab, kw in CONFIGS:
        if only and lab not in only:
            continue
        t = time.time()
        c = mzz.Alg73Z(A, **kw)
        ts = time.time() - t
        t = time.time()
        h, w = c.solve(b, tol=1e-8, maxcyc=30)
        tsol = time.time() - t
        conv = h[-1] < 1e-8
        ratios = [len(L.C) / L.A.shape[0] for L in c.levels[1:-1]]
        rc = np.mean(ratios) if ratios else 0.0
        th = [getattr(L, "theta", None) for L in c.levels[1:-1]]
        thm = np.mean([t_ for t_ in th if t_ is not None]) if any(t_ is not None for t_ in th) else float("nan")
        print(f"| {m} | {lab} | {len(c.levels)} | {rc:.2f} | {thm:.2f} | {rate(h):.3g} | {len(h)-1}{'' if conv else ' (no conv.)'} | "
              f"{w[-1]/max(len(h)-1,1):.0f} | {w[-1] if conv else float('inf'):.0f} | {c.op_complexity():.2f} | "
              f"{c.stored_complexity():.1f} | {ts:.1f} | {tsol:.1f} |")
        sys.stdout.flush()


if __name__ == "__main__":
    m = int(sys.argv[1])
    only = sys.argv[2].split(";") if len(sys.argv) > 2 else None
    if not only or only == [""]:
        print("| m | configuration | levels | mean nc/n (coarse) | mean theta | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] | solve [s] |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    main(m, only if only and only != [""] else None)
