"""Does a sparse complement Z allow more aggressive C/F splits than Q0?

Setting of Algorithm 7.3 (two-level, exact coarse solve): centered recirculating
flow on (-1,1)^2, cell Peclet 20, Jacobi omega = 1, P = P0 = [I; 0], injection.
F-points greedy with ||I - D_ff^{-1} A_ff||_inf <= delta (larger delta = fewer
C-points).  CR propagator K = I - D_ff^{-1} A_ff, H = -D_ff^{-1} A_fc.
N_Z with Q = Q0 + P0 Z, Z on the coarse points within graph distance <= d of
each fine point (bounded), Z = argmin ||K_Z^k||_F (JAX + L-BFGS).
"""

import sys
import numpy as np
import scipy.sparse as sp
import rbamg as rb
import sparse_problems as spp
import ml_cycles as mc
import nz_opt as no


def split_points(A, m, kind):
    """'greedyX': greedy with delta = X; 'full': standard coarsening (C at odd-odd,
    n_c/n ~ 1/4); 'semi-x': C on every other column (n_c/n ~ 1/2)."""
    if kind.startswith("greedy"):
        return mc.greedy_f_points(A, float(kind[6:]))
    I, J = np.divmod(np.arange(m * m), m)
    if kind == "full":
        return ~((I % 2 == 1) & (J % 2 == 1))
    return ~(J % 2 == 1)


def setup(m, kind, k):
    A = spp.cd(m, spp.RECIRC, spp.eps_for_peclet(m, spp.RECIRC, 20, "2", False), "centered")
    isF = split_points(A, m, kind)
    C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
    perm = np.concatenate([C, F])
    Ad = A.toarray()[np.ix_(perm, perm)]
    nc = len(C)
    P = np.vstack([np.eye(nc), np.zeros((len(F), nc))])
    R = P.T.copy()                                     # injection
    M = np.diag(1 / np.diag(Ad))                       # Jacobi, omega = 1
    return rb.Setup(Ad, M, P, R, k), A, perm


def pattern_graph(s, Asp, perm, d):
    G = (abs(sp.csr_matrix(Asp)) > 0).astype(float)
    Gd = G.copy()
    for _ in range(d - 1):
        Gd = ((Gd @ G) > 0).astype(float)
    Gd = Gd.toarray()[np.ix_(perm, perm)]
    return [np.nonzero(Gd[s.nc + j, :s.nc])[0] for j in range(s.nf)]


KINDS = ("greedy1.5", "semi-x", "full")
DS = (2, 3, 4)


def main(m=31, k=4):
    print(f"# Aggressive splits with a sparse complement (recirculating flow, Pe 20, m = {m}, k = {k})\n")
    print("| split | nc/n | rho(T) | rho(K) Q0 | M | N (Q0) | d | nnz(Z)/nf | rho(K_Z) | N_Z |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for kind in KINDS:
        s, A, perm = setup(m, kind, k)
        base = (s.nc / s.n, rb.rho(s.T), rb.rho(s.K), rb.rho(rb.E_markov(s)), rb.rho(rb.E_nonmarkov(s)))
        for d in DS:
            pr = no.Problem(s, pattern_graph(s, A, perm, d))
            z = pr.optimize("Kk", pr.ls_start())
            r = pr.report(z)
            f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
            print(f"| {kind} | {base[0]:.2f} | {base[1]:.2f} | {base[2]:.3f} | {f(base[3])} | {f(base[4])} | {d} | "
                  f"{r['nnz']:.1f} | {f(r['rhoK'])} | {f(r['rhoE'])} |")
            sys.stdout.flush()


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 31)
