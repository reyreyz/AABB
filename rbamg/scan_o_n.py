"""Two-level scan for a mesh-independent N_Z example (bounded Z).

Constant-wind centered CD on (-1,1)^2, beta = (cos pi/6, sin pi/6), cell Peclet
|beta| h / eps = Pe (fixed under refinement), greedy F-points with delta,
Jacobi omega = 1, P0, injection, Z with <= q entries per column (ml_z73).
Two levels, exact coarse solve.  Reported: CR rates of K (Q0) and K_Z, and the
convergence factors of M, N (Q0), N_Z.
"""

import sys
import numpy as np
import sparse_problems as spp
import ml_z73 as mzz
import ml_nz as mz
from run_ml_nz import rate

WIND = lambda x, y: (np.cos(np.pi / 6) + 0 * x, np.sin(np.pi / 6) + 0 * x)


def problem(m, pe):
    return spp.cd(m, WIND, spp.eps_for_peclet(m, WIND, pe, "2", False), "centered")


def main(ms=(31, 63, 127), pes=(8, 20), deltas=(0.6, 1.0, 1.5), k=4, q=4):
    print(f"# Two-level scan, constant-wind centered CD, k = {k}, q = {q}\n")
    print("| Pe | delta | m | nc/n | rho(K) Q0 | rho(K_Z) | M | N (Q0) | N_Z |")
    print("|---|---|---|---|---|---|---|---|---|")
    for pe in pes:
        for dl in deltas:
            for m in ms:
                A = problem(m, pe)
                b = np.random.default_rng(0).standard_normal(A.shape[0])
                out = []
                for kw in (dict(markov=True), dict(zmode="Q0"), dict(q=q)):
                    c = mzz.Alg73Z(A, k=k, delta=dl, max_levels=2, **kw)
                    h, _ = c.solve(b, maxcyc=40)
                    out.append((c, rate(h)))
                L0q, L0z = out[1][0].levels[0], out[2][0].levels[0]
                f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                print(f"| {pe} | {dl} | {m} | {len(L0q.C)/A.shape[0]:.2f} | {mz.cr_rate(L0q.KZ):.3f} | "
                      f"{mz.cr_rate(L0z.KZ):.3f} | {f(out[0][1])} | {f(out[1][1])} | {f(out[2][1])} |")
                sys.stdout.flush()


if __name__ == "__main__":
    main()
