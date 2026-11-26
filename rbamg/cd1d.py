"""1D centered convection-diffusion, coarsening by 3: an exactly analysable N_Z example.

    -eps u'' + u' on (0,1), centered differences, cell Peclet Pe = h / eps (fixed).
C-points: every third point; between two C-points sit two F-points, so A_ff is
block diagonal with 2x2 blocks.  Jacobi (omega = 1) gives the local CR blocks
    K_blk = -D^{-1} A_ff,blk,   |mu| = sqrt(Pe^2/4 - 1) / 2   (Pe > 2),
so rho(K) > 1 for Pe > sqrt(20) ~ 4.47: the local complement Q0 cannot converge.
Z: each F-point uses its two neighbouring C-points (2 entries per column).
Two-level, P0 = [I; 0], injection (setting of Algorithm 7.3).
"""

import sys
import numpy as np
import scipy.sparse as sp
import rbamg as rb
import nz_opt as no


def cd1d(n, pe):
    h = 1.0 / (n + 1)
    eps = h / pe
    lo = -eps - 0.5 * h            # coefficient of u_{i-1} (h^2-scaled: -eps - h b/2)
    up = -eps + 0.5 * h
    return sp.diags([lo * np.ones(n - 1), 2 * eps * np.ones(n), up * np.ones(n - 1)], [-1, 0, 1]).tocsr()


def split3(n):
    """C at indices 2, 5, 8, ... (two F-points between C-points, F at both ends)."""
    isC = np.zeros(n, dtype=bool)
    isC[2::3] = True
    return ~isC


def setup(n, pe, k):
    A = cd1d(n, pe).toarray()
    isF = split3(n)
    C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
    perm = np.concatenate([C, F])
    A = A[np.ix_(perm, perm)]
    nc = len(C)
    P = np.vstack([np.eye(nc), np.zeros((len(F), nc))])
    s = rb.Setup(A, np.diag(1 / np.diag(A)), P, P.T.copy(), k)
    return s, C, F


def pattern_nn(s, C, F):
    """two nearest C-points of each F-point."""
    return [np.argsort(np.abs(C - f))[:2] for f in F]


def main(n=191, pes=(3, 5, 8, 20), ks=(2, 3, 4)):
    print(f"# 1D centered CD, coarsening by 3, n = {n}\n")
    print("| Pe | k | rho(T) | rho(K) Q0 | LFA bound | objective | rho(K_Z) | M | N (Q0) | N_Z |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for pe in pes:
        for k in ks:
            s, C, F = setup(n, pe, k)
            lfa = np.sqrt(max(pe ** 2 / 4 - 1, 0)) / 2
            pr = no.Problem(s, pattern_nn(s, C, F))
            z0 = pr.ls_start()
            for ob in ("LS", "Kk", "K_rho"):
                z = z0 if ob == "LS" else pr.optimize(ob, z0)
                r = pr.report(z)
                f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                print(f"| {pe} | {k} | {rb.rho(s.T):.2f} | {rb.rho(s.K):.3f} | {lfa:.3f} | {ob} | {f(r['rhoK'])} | "
                      f"{f(rb.rho(rb.E_markov(s)))} | {f(rb.rho(rb.E_nonmarkov(s)))} | {f(r['rhoE'])} |")
                sys.stdout.flush()


if __name__ == "__main__":
    main()
