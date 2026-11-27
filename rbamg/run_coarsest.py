"""Robust N with an early exact coarse solve: stop the hierarchy at
n_coarsest = c * sqrt(n) and factor that level (sparse LU).  Work counts the LU
triangular solves (nnz(L) + nnz(U)) at every visit; storage includes the factors.
Bounded base: k = 4 / 8, refinement on levels 0-1, S truncated at 1e-5,
coarse filter 1e-3 (levels < 8) / 3e-4 (below), greedy delta = 0.3.
"""

import sys
import time
import numpy as np
import ml_z73 as mzz
from run_base_nz import problem
from run_ml_nz import rate

B = dict(zmode="Q0", k=4, k_coarse=8, tauP=1e-5, tau_deep=3e-4, delta=0.3, refine_levels=1)


def main(m, cs):
    A = problem(m)
    n = A.shape[0]
    b = np.random.default_rng(0).standard_normal(n)
    for c in cs:
        nco = 150 if c == 0 else int(c * np.sqrt(n))
        t = time.time()
        cy = mzz.Alg73Z(A, coarsest=nco, **B)
        ts = time.time() - t
        t = time.time()
        h, w = cy.solve(b, tol=1e-8, maxcyc=30)
        tsol = time.time() - t
        conv = h[-1] < 1e-8
        Lc = cy.levels[-1]
        lab = "150 (base)" if c == 0 else f"{c} sqrt(n) = {nco}"
        print(f"| {m} | {n} | {lab} | {len(cy.levels)} | {Lc.A.shape[0]} | {Lc.lu_nnz/A.nnz:.2f} | {rate(h):.3g} | "
              f"{len(h)-1}{'' if conv else ' (no conv.)'} | {w[-1]/max(len(h)-1,1):.0f} | {w[-1] if conv else float('inf'):.0f} | "
              f"{cy.op_complexity():.2f} | {cy.stored_complexity():.1f} | {ts:.1f} | {tsol:.1f} |")
        sys.stdout.flush()


if __name__ == "__main__":
    m = int(sys.argv[1])
    cs = [float(a) for a in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0, 4, 8, 16]
    main(m, cs)
