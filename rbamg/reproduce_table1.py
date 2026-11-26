"""Section 6.1 with Ruge-Stueben coarsening: Table 1, Fig. 1, Table SM11 (rho),
SM6 ideal-P reference, and the decomposition (4.4) of Fig. SM16.

Setup (identified by matching Table 1): PyAMG-style RS splitting (rs.py) with
classical strength theta = 0.25 and the 'min' measure, direct interpolation,
Galerkin Rhat = P^T, damped Jacobi omega = 2/3, local complement Q0.
Rotated anisotropy = PyAMG diffusion_stencil_2d(1e-3, pi/6, 'FD').
Recirculating flow: stand-in (node upwind on (-1,1)^2); not identified.
"""

import sys
import numpy as np
import rbamg as rb
import problems as pb
import rs
from reproduce_sm4 import ls_feedback, nonmarkov_with

PAPER_T1 = {  # (problem, m): (nc/n, M, M', S, N, NZ, rho(K)^k), k = 2
    ("Poisson", 31): (0.50, 0.111, 0.111, 0.111, 0.111, 0.111, 0.111),
    ("Poisson", 63): (0.50, 0.111, 0.111, 0.111, 0.111, 0.111, 0.111),
    ("Rotated aniso.", 31): (0.33, 0.733, 0.507, 0.526, 0.507, 0.505, 0.505),
    ("Rotated aniso.", 63): (0.33, 0.840, 0.511, 0.664, 0.511, 0.511, 0.511),
    ("Recirc. 1e-4 (stand-in)", 31): (0.50, 0.709, 0.403, 0.111, 0.257, 0.114, 0.843),
    ("Recirc. 1e-4 (stand-in)", 63): (0.51, 0.690, 0.403, 0.112, 0.254, 0.114, 0.650),
}
PAPER_SM11_K4 = (0.833, 0.261, 0.441, 0.261, 0.260)   # rotated aniso, m = 63, k = 4
PAPER_SM6_IDEAL = {2: 0.510, 4: 0.261}                 # rotated aniso, m = 63, ideal P


def recirc(m, eps):
    """Stand-in: -eps Lap u + beta . grad u on (-1,1)^2 nodes, beta =
    (2y(1-x^2), -2x(1-y^2)), first-order upwind, h^2-scaled with h of the unit square."""
    h = 1.0 / (m + 1)
    t = 2 * np.arange(1, m + 1) * h - 1
    n = m * m
    A = np.zeros((n, n))
    idx = lambda i, j: i * m + j
    for i in range(m):
        for j in range(m):
            x, y, p = t[j], t[i], idx(i, j)
            A[p, p] += 4 * eps
            for ii, jj in ((i, j - 1), (i, j + 1), (i - 1, j), (i + 1, j)):
                if 0 <= ii < m and 0 <= jj < m:
                    A[p, idx(ii, jj)] -= eps
            for b, lo, hi in ((2 * y * (1 - x * x), (i, j - 1), (i, j + 1)),
                              (-2 * x * (1 - y * y), (i - 1, j), (i + 1, j))):
                nb = lo if b > 0 else hi
                A[p, p] += h * abs(b)
                if 0 <= nb[0] < m and 0 <= nb[1] < m:
                    A[p, idx(*nb)] -= h * abs(b)
    return A


MATS = {
    "Poisson": pb.poisson,
    "Rotated aniso.": lambda m: rs.stencil_grid(rs.diffusion_stencil_2d(1e-3, np.pi / 6, "FD"), m),
    "Recirc. 1e-4 (stand-in)": lambda m: recirc(m, 1e-4),
    "CD eps=1e-5": lambda m: pb.convection_diffusion(m, 1e-5),
}


def rs_setup(name, m, k, A=None):
    A = MATS[name](m) if A is None else A
    perm, P, _ = rs.ruge_stuben(A, 0.25, "min")
    A = A[np.ix_(perm, perm)]
    return rb.Setup(A, pb.jacobi(A), P, P.T, k)


def all_factors(s):
    nz = nonmarkov_with(s, ls_feedback(s))[2]
    return (rb.rho(rb.E_markov(s)), rb.rho(rb.E_markov_pp(s)), rb.rho(rb.E_semi(s)),
            rb.rho(rb.E_nonmarkov(s)), nz, rb.rho(s.K) ** s.k)


def table1(ms=(31, 63)):
    print("\n## Table 1 (k = 2): ours / paper\n")
    print("| problem | m | nc/n | M | M' | S | N | NZ | rho(K)^k |")
    print("|---|---|---|---|---|---|---|---|---|")
    for name in ("Poisson", "Rotated aniso.", "Recirc. 1e-4 (stand-in)"):
        for m in ms:
            s = rs_setup(name, m, 2)
            ours = (s.nc / s.n,) + all_factors(s)
            pap = PAPER_T1[(name, m)]
            print(f"| {name} | {m} | " + " | ".join(
                f"{o:.{2 if i == 0 else 3}f} / {p:.{2 if i == 0 else 3}f}"
                for i, (o, p) in enumerate(zip(ours, pap))) + " |")
            sys.stdout.flush()


def sm11_and_ideal(m=63):
    print(f"\n## Rotated anisotropy, m = {m}: Table SM11 (k = 4) and SM6 ideal P\n")
    s = rs_setup("Rotated aniso.", m, 4)
    ours = all_factors(s)[:5]
    print("| k=4 | M | M' | S | N | NZ |\n|---|---|---|---|---|---|")
    print("| ours / paper | " + " | ".join(f"{o:.3f} / {p:.3f}" for o, p in zip(ours, PAPER_SM11_K4)) + " |")
    print("\n| ideal P, Galerkin | k = 2 | k = 4 |\n|---|---|---|")
    vals = []
    for k in (2, 4):
        s = rs_setup("Rotated aniso.", m, k)
        nc = s.nc
        A = s.A
        Wid = -np.linalg.solve(A[nc:, nc:], A[nc:, :nc])
        P = np.vstack([np.eye(nc), Wid])
        si = rb.Setup(A, s.M, P, P.T, k)
        vals.append(rb.rho(rb.E_markov(si)))
    print(f"| ours / paper | {vals[0]:.3f} / {PAPER_SM6_IDEAL[2]:.3f} | {vals[1]:.3f} / {PAPER_SM6_IDEAL[4]:.3f} |")


def fig1(m=31, ks=range(1, 7)):
    print(f"\n## Fig. 1 data (m = {m}): rho(E) versus k\n")
    for name in ("Poisson", "Rotated aniso.", "CD eps=1e-5", "Recirc. 1e-4 (stand-in)"):
        print(f"\n{name}\n\n| k | M | M' | S | N | NZ | rho(K)^k |\n|---|---|---|---|---|---|---|")
        for k in ks:
            print(f"| {k} | " + " | ".join(f"{v:.2e}" for v in all_factors(rs_setup(name, m, k))) + " |")
            sys.stdout.flush()
    print("\nPaper (text): Poisson M 0.043 vs 0.0014 (memory schemes) at k = 6; rotated "
          "anisotropy M stagnates at 0.72, M', N, NZ reach 0.13; recirc. 1e-4 M 0.20 at k = 6.")


def sm16(m=31, ks=range(1, 7)):
    """(4.4): E^M = sum_l P^(l) Psh T^{k-l} (memory) + Q K^k Qsh (CR noise), canonical Q."""
    print(f"\n## Fig. SM16 data (m = {m}): rho(E^M), ||memory||_2, ||CR noise||_2\n")
    for name in ("Poisson", "Rotated aniso.", "Recirc. 1e-4 (stand-in)"):
        print(f"\n{name}\n\n| k | rho(E^M) | memory | CR noise |\n|---|---|---|---|")
        for k in ks:
            s = rs_setup(name, m, k)
            Q = s.PiR @ s.Q0
            Binv = np.linalg.inv(np.hstack([s.P, Q]))
            Psh, Qsh = Binv[: s.nc], Binv[s.nc:]
            K, H = s.blocks(Q, Qsh)
            mem = sum(Q @ np.linalg.matrix_power(K, l - 1) @ H @ Psh @ np.linalg.matrix_power(s.T, k - l)
                      for l in range(1, k + 1))
            noise = Q @ np.linalg.matrix_power(K, k) @ Qsh
            print(f"| {k} | {rb.rho(rb.E_markov(s)):.3e} | {np.linalg.norm(mem, 2):.3e} | {np.linalg.norm(noise, 2):.3e} |")
            sys.stdout.flush()


if __name__ == "__main__":
    print("# Section 6.1 reproduction (Ruge-Stueben interpolation, Galerkin restriction)")
    table1()
    sm11_and_ideal()
    fig1()
    sm16()
