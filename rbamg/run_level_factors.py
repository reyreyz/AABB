"""Per-level two-level factors of the robust N hierarchy as m grows.

Hierarchy: bounded base (greedy delta = 0.3, k = 4 / 8, S truncated at 1e-5,
filter 1e-3 / 3e-4, refinement 0-1).  For each of the first levels, the level
operator A_l is treated as a two-level problem with the same construction
(k = 4 on level 0, 8 below; same split rule and truncation) and an exact,
unfiltered coarse solve; its convergence factor is measured from 30 cycles.
"""

import sys
import numpy as np
import ml_z73 as mzz
import ml_nz as mz
from run_base_nz import problem
from run_ml_nz import rate

B = dict(zmode="Q0", k=4, k_coarse=8, tauP=1e-5, tau_deep=3e-4, delta=0.3, refine_levels=1)


def main(ms=(63, 127, 255, 511), nlev=6):
    print("# Per-level two-level factors (exact coarse solve) of the robust N hierarchy\n")
    print("| m | level | n_l | nnz/row | nc/n | rho(K) | two-level factor |")
    print("|---|---|---|---|---|---|---|")
    for m in ms:
        A = problem(m)
        c = mzz.Alg73Z(A, **B)
        for lev, L in enumerate(c.levels[:min(nlev, len(c.levels) - 1)]):
            Al = L.A
            k = 4 if lev == 0 else 8
            t = mzz.Alg73Z(Al, zmode="Q0", k=k, delta=0.3, tauA=0.0, tauP=1e-5, tauP_from=0 if lev else 1,
                           max_levels=2)
            bl = np.random.default_rng(1).standard_normal(Al.shape[0])
            h, _ = t.solve(bl, tol=1e-12, maxcyc=30)
            hh = np.asarray(h)
            nc_ = len(hh) - 1
            fac = (hh[-1] / hh[1]) ** (1.0 / max(nc_ - 1, 1)) if nc_ > 1 else hh[-1]
            L0 = t.levels[0]
            print(f"| {m} | {lev} | {Al.shape[0]} | {Al.nnz/Al.shape[0]:.1f} | {len(L0.C)/Al.shape[0]:.2f} | "
                  f"{mz.cr_rate(L0.KZ):.3f} | {fac:.3g} ({nc_} cycles to 1e-12) |")
            sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (63, 127, 255, 511)
    main(ms)
