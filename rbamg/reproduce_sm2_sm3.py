"""Tables SM2 (mesh dependence) and SM3 (cost at m = 63), bilinear P.

SM2: rho(E) of M, M', S, N for CD eps = 1e-3 and random diffusion,
     k = 2, 4, m = 31, 63 (dense eig) and m = 127 (matrix-free Arnoldi, sparse.py).
SM3: operator complexities nnz(P')/nnz(P), nnz(A'_sigma)/nnz(A_sigma), and
     work per cycle in units of nnz(A).

Work rule, derived from the M rows of SM3 (3.9 = 3 + 2 * 0.451 at k = 2):
  one application of a matrix X costs nnz(X)/nnz(A) WU; relaxation costs k+1.
  M : (k+1) + [P] + [Rhat]
  M': (k+1) + [P'] + [Rhat]
"""

import sys
import numpy as np
import rbamg as rb
import problems as pb
from reproduce_sm1 import setup

PAPER_SM2 = {  # (problem, k, m): (M, M', S, N)
    ("CD eps=1e-3", 2, 31): (0.341, 0.218, 0.207, 0.232),
    ("CD eps=1e-3", 2, 63): (0.267, 0.237, 0.216, 0.323),
    ("CD eps=1e-3", 2, 127): (0.275, 0.345, 0.313, 0.368),
    ("CD eps=1e-3", 4, 31): (0.159, 0.095, 0.043, 0.050),
    ("CD eps=1e-3", 4, 63): (0.195, 0.099, 0.047, 0.111),
    ("CD eps=1e-3", 4, 127): (0.186, 0.087, 0.096, 0.146),
    ("Random diffusion", 2, 31): (0.822, 0.793, 0.810, 0.798),
    ("Random diffusion", 2, 63): (0.833, 0.799, 0.816, 0.814),
    ("Random diffusion", 2, 127): (0.861, 0.827, 0.845, 0.848),
    ("Random diffusion", 4, 31): (0.702, 0.630, 0.657, 0.628),
    ("Random diffusion", 4, 63): (0.730, 0.639, 0.666, 0.661),
    ("Random diffusion", 4, 127): (0.768, 0.690, 0.714, 0.722),
}

PAPER_SM3 = {  # (problem, k, scheme): (rho, WU/cycle, nnzP'/nnzP, nnzA'/nnzA)
    ("CD eps=1e-3", 2, "M"): (0.267, 3.9, 1.00, 1.00),
    ("CD eps=1e-3", 2, "M'"): (0.237, 5.8, 5.22, 3.92),
    ("CD eps=1e-3", 4, "M"): (0.195, 5.9, 1.00, 1.00),
    ("CD eps=1e-3", 4, "M'"): (0.099, 12.5, 15.55, 8.19),
    ("Random diffusion", 2, "M"): (0.833, 3.9, 1.00, 1.00),
    ("Random diffusion", 2, "M'"): (0.799, 6.4, 6.47, 3.92),
    ("Random diffusion", 4, "M"): (0.730, 5.9, 1.00, 1.00),
    ("Random diffusion", 4, "M'"): (0.639, 12.8, 16.33, 8.19),
}


def nnz(X, rel=1e-13):
    """Nonzeros above rel * max|X| (exact cancellations count as zeros)."""
    return int(np.count_nonzero(np.abs(X) > rel * np.abs(X).max()))


def sm3(m=63):
    print(f"\n## Table SM3 (m = {m}): ours / paper\n")
    print("| problem | k | scheme | rho | WU/cycle | nnz(P')/nnz(P) | nnz(A'σ)/nnz(Aσ) |")
    print("|---|---|---|---|---|---|---|")
    for name in ("CD eps=1e-3", "Random diffusion"):
        for k in (2, 4):
            s = setup(name, m, k)
            nA, nP, nPp = nnz(s.A), nnz(s.P), nnz(s.Pp)
            As, Asp = s.Rhat @ s.A @ s.P, s.Rhat @ s.A @ s.Pp
            rows = {
                "M": (rb.rho(rb.E_markov(s)), (k + 1) + 2 * nP / nA, 1.0, 1.0),
                "M'": (rb.rho(rb.E_markov_pp(s)), (k + 1) + (nPp + nP) / nA,
                       nPp / nP, nnz(Asp) / nnz(As)),
            }
            for sch, ours in rows.items():
                pap = PAPER_SM3[(name, k, sch)]
                fmt = ("{:.3f}", "{:.1f}", "{:.2f}", "{:.2f}")
                cells = " | ".join(f.format(o) + " / " + f.format(p) for f, o, p in zip(fmt, ours, pap))
                print(f"| {name} | {k} | {sch} | {cells} |")
            sys.stdout.flush()


def sm2(ms=(31, 63)):
    print("\n## Table SM2: ours / paper (M, M', S, N)\n")
    print("| problem | k | m | M | M' | S | N |")
    print("|---|---|---|---|---|---|---|")
    for name in ("CD eps=1e-3", "Random diffusion"):
        for k in (2, 4):
            for m in ms:
                s = setup(name, m, k)
                ours = (rb.rho(rb.E_markov(s)), rb.rho(rb.E_markov_pp(s)),
                        rb.rho(rb.E_semi(s)), rb.rho(rb.E_nonmarkov(s)))
                pap = PAPER_SM2[(name, k, m)]
                print(f"| {name} | {k} | {m} | " +
                      " | ".join(f"{o:.3f} / {p:.3f}" for o, p in zip(ours, pap)) + " |")
                sys.stdout.flush()


if __name__ == "__main__":
    print("# Tables SM2 and SM3 (bilinear P). Random diffusion uses the stand-in field.")
    sm3()
    sm2()
