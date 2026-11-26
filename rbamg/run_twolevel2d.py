"""O(N) 2D example: centered CD (wind along x, cell Peclet 8), semi-coarsening by 3,
stencil decomposition (W, Z) with 48 fixed numbers, two-level non-Markovian cycle
with coarse solves by PyAMG AIR (O(N)); sparse LU as reference on small grids.
M and N (Q0, P0) use the same split, relaxation and coarse solver.
"""

import sys
import time
import numpy as np
from stencil2d import cd2d_x, layout2d, fit2d
from twolevel_wz import TwoLevelWZ
from run_ml_nz import rate


def main(pe=8, k=4, sizes=(47, 95, 191, 383), coarse="air", coarse_cycles=2):
    t = time.time()
    th = fit2d(pe, k, wy=1)
    print(f"# O(N) two-level example: centered CD, wind along x, Pe = {pe}, k = {k}\n")
    print(f"stencil (W, Z): {len(th)} numbers, fitted once on a 23 x 12 grid in {time.time()-t:.0f} s; "
          f"coarse solver: {coarse} ({coarse_cycles} cycles)\n")
    print("| m | n | cycle | factor | cycles to 1e-10 | WU/cycle | stored nnz / nnz(A) | setup [µs/unknown] |")
    print("|---|---|---|---|---|---|---|---|")
    for m in sizes:
        A = cd2d_x(m, m, pe)
        C, F, e, nslot = layout2d(m, m, 1, 2, 1)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for lab, mode in (("M", "M"), ("N (Q0, P0)", "Q0"), ("N_WZ", "WZ")):
            tl = TwoLevelWZ(A, C, F, e, nslot, th, k, mode=mode, coarse=coarse, coarse_cycles=coarse_cycles)
            h, w = tl.solve(b, maxcyc=40)
            conv = h[-1] < 1e-10
            print(f"| {m} | {A.shape[0]} | {lab} | {rate(h):.3g} | {len(h)-1}{'' if conv else ' (no conv.)'} | "
                  f"{w[-1]/max(len(h)-1,1):.1f} | {tl.complexity():.2f} | {1e6*tl.setup_time/A.shape[0]:.1f} |")
            sys.stdout.flush()


if __name__ == "__main__":
    coarse = sys.argv[1] if len(sys.argv) > 1 else "air"
    main(coarse=coarse)
