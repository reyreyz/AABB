"""Search for a problem where M diverges, N with Q0 fails, and N_Z converges.

Exact two-grid analysis, bilinear P (C/F order), Galerkin Rhat = P^T, Jacobi omega.
For each problem:
  rho(T)                     relaxation      (M ~ Pi_R T^k diverges if rho(T) > 1)
  rho(K), rho(K_Z)           CR propagators of Q0 and of Q0 + P Z
  rho(E_M), rho(E_S), rho(E_N(Q0)), rho(E_N(Z))
Z candidates (pattern: interpolatory coarse points): LS (paper), LSk (k-step CR LS).
"""

import sys
import numpy as np
import rbamg as rb
import problems as pb
from reproduce_sm4 import ls_feedback, nonmarkov_with

mpow = np.linalg.matrix_power


def centered_cd(m, pe, theta=np.pi / 6):
    """-eps Lap u + beta.grad u, centered, unit square, cell Peclet |beta|_inf h / eps = pe."""
    h = 1.0 / (m + 1)
    bx, by = np.cos(theta), np.sin(theta)
    eps = max(abs(bx), abs(by)) * h / pe
    I = np.eye(m)
    D0 = 0.5 * (np.eye(m, k=1) - np.eye(m, k=-1))
    return eps * pb.poisson(m) + h * (bx * np.kron(I, D0) + by * np.kron(D0, I))


def lsk_feedback(s):
    """column j: min ||K^{k-1}(K e_j + H z)||, z on the interpolatory pattern."""
    nc, nf = s.nc, s.nf
    W = s.P[nc:]
    L = mpow(s.K, s.k - 1)
    LK, LH = L @ s.K, L @ s.H
    Z = np.zeros((nc, nf))
    for j in range(nf):
        Cj = np.nonzero(W[j])[0]
        z, *_ = np.linalg.lstsq(LH[:, Cj], -LK[:, j], rcond=None)
        Z[Cj, j] = z
    return Z


def analyse(A, m, k, omega=2 / 3):
    perm, P = pb.bilinear_cf(m)
    A = A[np.ix_(perm, perm)]
    M = omega * np.diag(1 / np.diag(A))
    s = rb.Setup(A, M, P, P.T, k)
    out = dict(rhoT=rb.rho(s.T), rhoK=rb.rho(s.K), EM=rb.rho(rb.E_markov(s)),
               ES=rb.rho(rb.E_semi(s)), EN=rb.rho(rb.E_nonmarkov(s)))
    for lab, Z in (("LS", ls_feedback(s)), ("LSk", lsk_feedback(s))):
        rK, nK, rE = nonmarkov_with(s, Z)
        out["rhoK_" + lab], out["EN_" + lab] = rK, rE
    return out


def main(m=31):
    print(f"# Showcase scan (m = {m}, bilinear P, Galerkin, Jacobi 2/3)\n")
    print("| problem | k | rho(T) | rho(K) Q0 | rho(K_Z) LS | rho(K_Z) LSk | M | S | N (Q0) | N_Z LS | N_Z LSk |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    cases = [(f"centered CD Pe={pe}", centered_cd(m, pe)) for pe in (1, 2, 5, 10, 20)]
    lam = np.sort(np.linalg.eigvalsh(pb.poisson(m)))
    for nneg in (3, 6, 10, 15):
        cases.append((f"shifted Laplacian, {nneg} neg. modes",
                      pb.poisson(m) - 0.5 * (lam[nneg - 1] + lam[nneg]) * np.eye(m * m)))
    for lab, A in cases:
        for k in (2, 4):
            r = analyse(A, m, k)
            f = lambda v: f"**{v:.3f}**" if v < 1 else f"{v:.3g}"
            print(f"| {lab} | {k} | {r['rhoT']:.3f} | {r['rhoK']:.3f} | {r['rhoK_LS']:.3f} | {r['rhoK_LSk']:.3f} | "
                  f"{f(r['EM'])} | {f(r['ES'])} | {f(r['EN'])} | {f(r['EN_LS'])} | {f(r['EN_LSk'])} |")
            sys.stdout.flush()


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 31)
