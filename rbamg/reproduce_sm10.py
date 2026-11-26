"""Table SM10: GMRES iterations (relative residual 1e-8), one two-grid cycle
(k = 2) as right preconditioner, Ruge-Stueben interpolation, Galerkin restriction.

Preconditioner B b := one cycle of the algorithm from x0 = 0 (so B = (I - E) A^{-1}).
Algorithms used: M (Petrov-Galerkin), M' (P' in place of P), S (Alg. 4.1),
N (Alg. 4.2, Q0), NZ (Alg. 4.2 with the LS feedback complement Q0 + P Z).
"""

import numpy as np
import rbamg as rb
from reproduce_table1 import rs_setup
from reproduce_sm4 import ls_feedback

PAPER = {("Poisson", 63): (2, 6, 6, 6, 6), ("Rotated aniso.", 63): (20, 11, 14, 11, 12)}


def gmres_right(A, b, prec, tol=1e-8, maxit=200):
    """Right-preconditioned full GMRES; returns iterations to ||r||/||b|| <= tol."""
    n = len(b)
    beta = np.linalg.norm(b)
    V = np.zeros((maxit + 1, n))
    H = np.zeros((maxit + 1, maxit))
    V[0] = b / beta
    for j in range(maxit):
        w = A @ prec(V[j])
        for _ in range(2):
            c = V[: j + 1] @ w
            w -= c @ V[: j + 1]
            H[: j + 1, j] += c
        H[j + 1, j] = np.linalg.norm(w)
        e1 = np.zeros(j + 2)
        e1[0] = beta
        y, *_ = np.linalg.lstsq(H[: j + 2, : j + 1], e1, rcond=None)
        res = np.linalg.norm(H[: j + 2, : j + 1] @ y - e1) / beta
        if res <= tol:
            return j + 1
        if H[j + 1, j] < 1e-14:
            return j + 1
        V[j + 1] = w / H[j + 1, j]
    return maxit


def nz_alg(s, Z):
    """Alg. 4.2 with complement Q = Q0 + P Z: Qsh = Q0sh, Psh = P0sh - Z Q0sh,
    K_Z = K + H Z, H_Z = H."""
    KZ = s.K + s.H @ Z
    Q = s.Q0 + s.P @ Z
    Psh = s.P0sh - Z @ s.Q0sh
    Pp = s.P + Q @ sum(np.linalg.matrix_power(KZ, l) for l in range(s.k)) @ s.H
    RAPp = s.Rhat @ s.A @ Pp

    def alg(X, B):
        xs = rb.relax(s, X, B)
        xk = xs[s.k]
        y = np.zeros((s.nf,) + X.shape[1:])
        for l in range(s.k):
            y = KZ @ y + s.H @ (Psh @ (xk - xs[l]))
        h = Q @ y
        eps = np.linalg.solve(RAPp, s.Rhat @ (B - s.A @ xk) - s.Rhat @ s.A @ h)
        return xk + Pp @ eps + h
    return alg


def main(m=63, k=2):
    print(f"# Table SM10 (m = {m}, k = {k}): GMRES iterations, ours / paper\n")
    print("| problem | M | M' | S | N | NZ |\n|---|---|---|---|---|---|")
    rng = np.random.default_rng(0)
    for name in ("Poisson", "Rotated aniso."):
        s = rs_setup(name, m, k)
        b = rng.standard_normal(s.n)
        z = np.zeros(s.n)
        nzalg = nz_alg(s, ls_feedback(s))
        precs = [lambda v: rb.alg_markov(s, z, v),
                 lambda v: rb.alg_markov(s, z, v, s.Pp),
                 lambda v: rb.alg_semi(s, z, v),
                 lambda v: rb.alg_nonmarkov(s, z, v),
                 lambda v: nzalg(z, v)]
        its = [gmres_right(s.A, b, p) for p in precs]
        print(f"| {name} | " + " | ".join(f"{a} / {p}" for a, p in zip(its, PAPER[(name, m)])) + " |")


if __name__ == "__main__":
    main()
