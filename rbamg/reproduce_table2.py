"""Table 2 (m = 63): best-tuned two-grid configuration of each method on the same
Ruge-Stueben splitting, evaluated exactly (dense error propagators).

Relaxation step "F^nu C^mu, omega": nu Jacobi sweeps on the F-points, then mu on
the C-points, all with weight omega:  S_F = I - omega I_F D^{-1} A,  T = S_C^mu S_F^nu.
  lAIR  : P one-point, R = PyAMG local_air (distance d), E = (I - P(RAP)^{-1}RA) T
  MGR   : P = [I; -D_ff^{-1}A_fc] (Jacobi approximate ideal), injection, E = (I - P(A_C P)^{-1}A_C) T
  N F   : Cor. 4.7, k steps of F^nu: P' = [I; W_id + J^{nu k}(W - W_id)], injection,
          E = (I - P'(A_C P')^{-1} A_C) Q0 J^{nu k} [-W I]
  N F/C : Thm. 4.5 with T = S_C S_F^nu, Q0, injection, P' = P + Q0 sum_{l<k} K^l H,
          E = (I - P'(A_C P')^{-1} A_C) Q0 K^k Q0sh
Work per cycle (WU, units nnz(A)): sweeps (F: nnz(A_F:), C: nnz(A_C:)), residual
(nnz(A), or nnz(A_C:) for injection), restriction nnz(R), interpolation nnz(W).
"""

import sys
import numpy as np
import scipy.sparse as sp
import pyamg
from pyamg.strength import classical_strength_of_connection
from pyamg.classical.split import RS
from pyamg.classical.interpolate import one_point_interpolation, local_air
import sparse_problems as spp
from reproduce_sec63 import recirc

mpow = np.linalg.matrix_power


def cd_unit(m, eps, theta=np.pi / 6):
    """constant wind CD on the unit square (problems.convection_diffusion, sparse)."""
    import problems as pb
    return sp.csr_matrix(pb.convection_diffusion(m, eps, theta))


PROBLEMS = {
    "Recirc. 1e-2": lambda m: recirc(m, 1e-2),
    "Recirc. 1e-3": lambda m: recirc(m, 1e-3),
    "Recirc. 1e-4": lambda m: recirc(m, 1e-4),
    "Recirc. 1e-6": lambda m: recirc(m, 1e-6),
    "CD 1e-3": lambda m: cd_unit(m, 1e-3),
    "CD 1e-5": lambda m: cd_unit(m, 1e-5),
}

# (method, params, paper (rho, WU/cycle, nnz ratio, WU/digit, solves/digit))
CONFIGS = {
    "Recirc. 1e-2": [("lAIR", dict(d=2, w=1.0, nu=1, mu=0), (2.2e-9, 2.10, 0.89, 0.243, 0.116)),
                     ("MGR", dict(w=1.0, nu=1), (4.4e-15, 1.40, 0.89, 0.097, 0.070)),
                     ("N F", dict(W0="P0", w=1.0, nu=1, k=1), (4.4e-15, 1.40, 0.89, 0.097, 0.070)),
                     ("N F/C", dict(W0="P0", nu=1, mu=1, k=1), (4.0e-15, 2.70, 0.89, 0.187, 0.069))],
    "Recirc. 1e-3": [("lAIR", dict(d=2, w=1.0, nu=10, mu=0), (1.1e-9, 6.48, 0.90, 0.722, 0.112)),
                     ("MGR", dict(w=1.1, nu=2), (1.5e-1, 1.88, 0.90, 2.257, 1.203)),
                     ("N F", dict(W0="1pt", w=1.0, nu=4, k=3), (1.4e-11, 6.76, 0.90, 0.623, 0.092)),
                     ("N F/C", dict(W0="P0", nu=4, mu=1, k=3), (2.7e-11, 9.89, 0.90, 0.937, 0.095))],
    "Recirc. 1e-4": [("lAIR", dict(d=2, w=1.0, nu=6, mu=0), (1.3e-8, 4.51, 0.90, 0.571, 0.127)),
                     ("MGR", dict(w=1.1, nu=1), (1.4e-1, 1.39, 0.90, 1.653, 1.192)),
                     ("N F", dict(W0="P0", w=1.0, nu=3, k=3), (1.5e-12, 5.28, 0.90, 0.446, 0.084)),
                     ("N F/C", dict(W0="P0", nu=4, mu=1, k=2), (5.9e-11, 7.00, 0.90, 0.685, 0.098))],
    "Recirc. 1e-6": [("lAIR", dict(d=2, w=1.0, nu=4, mu=0), (1.5e-9, 3.54, 0.90, 0.401, 0.113)),
                     ("MGR", dict(w=1.1, nu=1), (1.2e-1, 1.39, 0.90, 1.492, 1.075)),
                     ("N F", dict(W0="P0", w=1.0, nu=2, k=3), (1.4e-14, 3.82, 0.90, 0.276, 0.072)),
                     ("N F/C", dict(W0="P0", nu=3, mu=1, k=2), (2.4e-14, 6.03, 0.90, 0.443, 0.073))],
    "CD 1e-3": [("lAIR", dict(d=2, w=1.0, nu=1, mu=0), (3.4e-9, 2.10, 0.89, 0.248, 0.118)),
                ("MGR", dict(w=1.0, nu=1), (1.6e-15, 1.40, 0.89, 0.095, 0.068)),
                ("N F", dict(W0="1pt", w=1.0, nu=1, k=1), (1.3e-15, 1.40, 0.89, 0.094, 0.067)),
                ("N F/C", dict(W0="P0", nu=1, mu=1, k=1), (1.5e-15, 2.70, 0.89, 0.182, 0.067))],
    "CD 1e-5": [("lAIR", dict(d=2, w=1.0, nu=1, mu=0), (2.0e-9, 2.10, 0.89, 0.242, 0.115)),
                ("MGR", dict(w=1.0, nu=1), (1.5e-15, 1.40, 0.89, 0.094, 0.068)),
                ("N F", dict(W0="P0", w=1.0, nu=1, k=1), (1.5e-15, 1.40, 0.89, 0.094, 0.068)),
                ("N F/C", dict(W0="P0", nu=1, mu=1, k=1), (1.4e-15, 2.70, 0.89, 0.182, 0.067))],
}


def rho(E):
    return np.max(np.abs(np.linalg.eigvals(E)))


def nnzd(X, tol=1e-14):
    return int(np.count_nonzero(np.abs(X) > tol * max(np.abs(X).max(), 1e-300)))


def evaluate(As, method, p, theta=0.25):
    S = classical_strength_of_connection(As, theta=theta, norm="min")
    split = RS(S, second_pass=True)
    C = np.nonzero(split == 1)[0]
    F = np.nonzero(split == 0)[0]
    perm = np.concatenate([C, F])
    A = As.toarray()[np.ix_(perm, perm)]
    n, nc = A.shape[0], len(C)
    nnzA = As.nnz
    nnzAF, nnzAC = nnzd(A[nc:]), nnzd(A[:nc])
    D = np.diag(A)
    I = np.eye(n)
    Aff, Afc = A[nc:, nc:], A[nc:, :nc]
    Wid = -np.linalg.solve(Aff, Afc)

    def sweeps(w, nu, mu):
        SF = I.copy()
        SF[nc:] -= w * (A[nc:] / D[nc:, None])
        SC = I.copy()
        SC[:nc] -= w * (A[:nc] / D[:nc, None])
        return mpow(SC, mu) @ mpow(SF, nu)

    def coarse_corr(P, R):
        RA = R @ A
        return I - P @ np.linalg.solve(RA @ P, RA), RA @ P

    if method == "lAIR":
        P1 = sp.csr_matrix(one_point_interpolation(As, S, split)).toarray()[perm]
        Rn = sp.csr_matrix(local_air(As, split, theta=0.05, degree=p["d"])).toarray()
        R = Rn[:, perm]
        T = sweeps(p["w"], p["nu"], p["mu"])
        CC, Asig = coarse_corr(P1, R)
        E = CC @ T
        wu = (p["nu"] * nnzAF + p["mu"] * nnzAC + nnzA + nnzd(R) + nnzd(P1[nc:])) / nnzA
    elif method == "MGR":
        W = -Afc / D[nc:, None]
        P = np.vstack([np.eye(nc), W])
        R = np.hstack([np.eye(nc), np.zeros((nc, n - nc))])
        T = sweeps(p["w"], p["nu"], 0)
        CC, Asig = coarse_corr(P, R)
        E = CC @ T
        wu = (p["nu"] * nnzAF + nnzAC + nnzd(W)) / nnzA
    elif method == "N F":
        if p["W0"] == "P0":
            W = np.zeros((n - nc, nc))
        else:
            W = sp.csr_matrix(one_point_interpolation(As, S, split)).toarray()[perm][nc:]
        J = np.eye(n - nc) - p["w"] * Aff / D[nc:, None]
        Jk = mpow(J, p["nu"] * p["k"])
        Wp = Wid + Jk @ (W - Wid)                              # Cor. 4.7
        Pp = np.vstack([np.eye(nc), Wp])
        R = np.hstack([np.eye(nc), np.zeros((nc, n - nc))])
        CC, Asig = coarse_corr(Pp, R)
        Q0 = np.vstack([np.zeros((nc, n - nc)), np.eye(n - nc)])
        E = CC @ Q0 @ Jk @ np.hstack([-W, np.eye(n - nc)])
        wu = (p["nu"] * p["k"] * nnzAF + nnzAC + nnzd(Wp)) / nnzA
    else:  # N F/C
        W = np.zeros((n - nc, nc))
        P = np.vstack([np.eye(nc), W])
        T = sweeps(1.0, p["nu"], p["mu"])
        Q0 = np.vstack([np.zeros((nc, n - nc)), np.eye(n - nc)])
        Q0sh = np.hstack([-W, np.eye(n - nc)])
        K, H = Q0sh @ T @ Q0, Q0sh @ T @ P
        Pp = P + Q0 @ sum(mpow(K, l) for l in range(p["k"])) @ H
        R = np.hstack([np.eye(nc), np.zeros((nc, n - nc))])
        CC, Asig = coarse_corr(Pp, R)
        E = CC @ Q0 @ mpow(K, p["k"]) @ Q0sh
        # k steps of (nu F + mu C) sweeps, C-residual, interpolation with P',
        # history recursion: k applications of K and H
        wu = (p["k"] * (p["nu"] * nnzAF + p["mu"] * nnzAC) + nnzAC + nnzd(Pp[nc:])
              + p["k"] * (nnzd(K) + nnzd(H))) / nnzA
    r = rho(E)
    dig = -np.log10(max(r, 1e-16))
    return r, wu, nnzd(Asig) / nnzA, wu / dig, 1 / dig


def main(m=63):
    print(f"# Table 2 reproduction (m = {m}), paper configurations evaluated exactly\n")
    print("Recirculating flow: stand-in discretization. Each cell: ours / paper.\n")
    print("| problem | method | rho | WU/cycle | nnz(Aσ)/nnz(A) | WU/digit | solves/digit |")
    print("|---|---|---|---|---|---|---|")
    for pname, rows in CONFIGS.items():
        As = PROBLEMS[pname](m)
        for meth, p, pap in rows:
            ours = evaluate(As, meth, p)
            fm = ("{:.1e}", "{:.2f}", "{:.2f}", "{:.3f}", "{:.3f}")
            print(f"| {pname} | {meth} {p} | " + " | ".join(
                f.format(o) + " / " + f.format(q) for f, o, q in zip(fm, ours, pap)) + " |")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
