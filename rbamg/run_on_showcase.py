"""O(N) showcase: multilevel N_Z V-cycles, centered CD at fixed cell Peclet number.

Everything sparse and bounded:
  * Z on a fixed graph-distance pattern (bounded nnz per column),
    Z = argmin ||K_Z^k||_F over the pattern: L-BFGS with sparse gradients and a
    fixed iteration cap (O(N) per level);
  * K_Z, H, P' = P + Q sum_{l<k} K_Z^l H sparse with fill bounded by k and the pattern;
  * coarse operators P^T A P' filtered (row-max, lumped), coarsest < 150 factored;
  * optional damping Z <- theta Z by the CR rate of the level pair.
Reported: factor, cycles to 1e-10, WU/cycle (nnz(A0) units), operator complexity,
max nnz(Z)/n_f over levels, setup time per unknown.
"""

import sys
import time
import numpy as np
import ml_nz as mz
from run_ml_nz import rate

VARIANTS = [
    ("M (Markovian)", dict(zmode="Q0", markov=True)),
    ("N (Q0)", dict(zmode="Q0")),
    ("N_Z Kk dist3", dict(zmode="Kk", pattern="dist3")),
    ("N_Z Kk dist3, theta auto", dict(zmode="Kk", pattern="dist3", theta="auto")),
    ("N_Z Kk dist5, theta auto", dict(zmode="Kk", pattern="dist5", theta="auto")),
]


def main(ms=(31, 63, 127, 255), pe=6.0, ks=(5,), tauA=1e-3, gamma=1):
    print(f"# O(N) N_Z showcase: centered CD, cell Peclet {pe}, bilinear hierarchy, Galerkin, "
          f"Jacobi 2/3, tauA = {tauA}, {'V' if gamma == 1 else 'W'}-cycles\n")
    print("| m | n | k | cycle | theta per level | factor | cycles | WU/cycle | op.cx | max nnz(Z)/nf | setup [µs/unknown] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = mz.cd_centered_sparse(m, pe)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for k in ks:
            for lab, kw in VARIANTS:
                t = time.time()
                c = mz.NZCycle(A, m, k=k, tauA=tauA, gamma=gamma, **kw)
                ts = time.time() - t
                h, w = c.solve(b, maxcyc=60)
                conv = h[-1] < 1e-10
                th = ",".join(f"{getattr(L, 'theta', 1.0):g}" for L in c.levels[:-1]) if kw.get("zmode") != "Q0" else "–"
                nz = max((L.Z.nnz / len(L.F) for L in c.levels[:-1]), default=0.0)
                print(f"| {m} | {A.shape[0]} | {k} | {lab} | {th} | {rate(h):.3g} | {len(h)-1}{'' if conv else ' (no conv.)'} | "
                      f"{w[-1]/(len(h)-1):.1f} | {c.op_complexity():.2f} | {nz:.1f} | {1e6*ts/A.shape[0]:.0f} |")
                sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (31, 63, 127, 255)
    main(ms)
