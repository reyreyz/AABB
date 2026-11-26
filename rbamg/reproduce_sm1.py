"""Table SM1 (m = 31): two-grid factors rho(E) of M, M', S, N and rho(K)^k.

Setup: coarse points at even indices, damped Jacobi omega = 2/3, Galerkin
Rhat = P^T, local complement Q0 for N and rho(K).
P is bilinear, or P1 = [I; sum_{l<3} T_ff^l T_fc]: three steps of projected
relaxation from P0 = [I; 0] (flow step (5.2) with k = 3; identified by
matching all 15 Poisson/P1 entries).
"""

import numpy as np
import rbamg as rb
import problems as pb

PAPER = {  # (problem, P, k): (M, M', S, N, rho(K)^k)
    ("Poisson", "bilinear", 1): (0.662, 0.661, 0.662, 0.662, 0.662),
    ("Poisson", "bilinear", 2): (0.439, 0.437, 0.438, 0.438, 0.438),
    ("Poisson", "bilinear", 4): (0.195, 0.192, 0.192, 0.192, 0.192),
    ("Poisson", "P1", 1): (0.975, 0.950, 0.975, 0.957, 0.664),
    ("Poisson", "P1", 2): (0.967, 0.884, 0.951, 0.891, 0.440),
    ("Poisson", "P1", 4): (0.958, 0.679, 0.904, 0.663, 0.194),
    ("Random diffusion", "bilinear", 1): (0.900, 0.891, 0.900, 0.895, 0.866),
    ("Random diffusion", "bilinear", 2): (0.822, 0.793, 0.810, 0.798, 0.750),
    ("Random diffusion", "bilinear", 4): (0.702, 0.630, 0.657, 0.628, 0.562),
    ("Random diffusion", "P1", 1): (0.991, 0.988, 0.991, 0.989, 0.934),
    ("Random diffusion", "P1", 2): (0.986, 0.974, 0.983, 0.975, 0.873),
    ("Random diffusion", "P1", 4): (0.979, 0.940, 0.965, 0.938, 0.761),
    ("CD eps=1e-2", "bilinear", 1): (0.589, 0.598, 0.589, 0.606, 0.626),
    ("CD eps=1e-2", "bilinear", 2): (0.334, 0.356, 0.347, 0.372, 0.392),
    ("CD eps=1e-2", "bilinear", 4): (0.078, 0.119, 0.120, 0.143, 0.153),
    ("CD eps=1e-2", "P1", 1): (0.589, 0.603, 0.589, 0.611, 0.641),
    ("CD eps=1e-2", "P1", 2): (0.425, 0.359, 0.346, 0.384, 0.412),
    ("CD eps=1e-2", "P1", 4): (0.297, 0.120, 0.120, 0.158, 0.169),
    ("CD eps=1e-3", "bilinear", 1): (0.455, 0.402, 0.455, 0.430, 0.471),
    ("CD eps=1e-3", "bilinear", 2): (0.341, 0.218, 0.207, 0.232, 0.222),
    ("CD eps=1e-3", "bilinear", 4): (0.159, 0.095, 0.043, 0.050, 0.049),
    ("CD eps=1e-3", "P1", 1): (0.459, 0.416, 0.458, 0.420, 0.527),
    ("CD eps=1e-3", "P1", 2): (0.310, 0.210, 0.210, 0.189, 0.278),
    ("CD eps=1e-3", "P1", 4): (0.176, 0.092, 0.045, 0.048, 0.077),
    ("Shifted Laplacian", "bilinear", 1): (0.667, 0.666, 0.667, 0.667, 0.667),
    ("Shifted Laplacian", "bilinear", 2): (0.446, 0.444, 0.445, 0.445, 0.445),
    ("Shifted Laplacian", "bilinear", 4): (0.201, 0.198, 0.198, 0.198, 0.198),
    ("Shifted Laplacian", "P1", 1): (1.068, 1.163, 1.068, 1.140, 0.669),
    ("Shifted Laplacian", "P1", 2): (1.089, 1.534, 1.140, 1.503, 0.447),
    ("Shifted Laplacian", "P1", 4): (1.114, 2.380, 1.299, 2.542, 0.200),
}

MATS = {
    "Poisson": pb.poisson,
    "Random diffusion": pb.random_diffusion,
    "CD eps=1e-2": lambda m: pb.convection_diffusion(m, 1e-2),
    "CD eps=1e-3": lambda m: pb.convection_diffusion(m, 1e-3),
    "Shifted Laplacian": pb.shifted_laplacian,
}


def p1(A, M, nc, sweeps=3):
    """P1 = [I; W],  W = sum_{l<sweeps} T_ff^l T_fc   (5.2)."""
    T = np.eye(len(A)) - M @ A
    Tff, Tfc = T[nc:, nc:], T[nc:, :nc]
    W = sum(np.linalg.matrix_power(Tff, l) @ Tfc for l in range(sweeps))
    return np.vstack([np.eye(nc), W])


def setup(name, m, k, interp="bilinear", A=None):
    perm, Pb = pb.bilinear_cf(m)
    A = (MATS[name](m) if A is None else A)[np.ix_(perm, perm)]
    M = pb.jacobi(A)
    P = Pb if interp == "bilinear" else p1(A, M, Pb.shape[1])
    return rb.Setup(A, M, P, P.T, k)


def factors(s):
    return (rb.rho(rb.E_markov(s)), rb.rho(rb.E_markov_pp(s)), rb.rho(rb.E_semi(s)),
            rb.rho(rb.E_nonmarkov(s)), rb.rho(s.K) ** s.k)


def main(m=31):
    print(f"# Table SM1 reproduction (m = {m})\n")
    print("Each cell: ours / paper. `!` marks |diff| > 0.0015.\n")
    print("| problem | P | k | M | M' | S | N | rho(K)^k |")
    print("|---|---|---|---|---|---|---|---|")
    stats = {}
    for name in MATS:
        for interp in ("bilinear", "P1"):
            for k in (1, 2, 4):
                ours = factors(setup(name, m, k, interp))
                paper = PAPER[(name, interp, k)]
                cells, nbad = [], 0
                for o, p in zip(ours, paper):
                    bad = abs(o - p) > 1.5e-3
                    nbad += bad
                    cells.append(f"{o:.3f} / {p:.3f}" + (" !" if bad else ""))
                stats.setdefault(name, [0, 0])
                stats[name][0] += 5 - nbad
                stats[name][1] += 5
                print(f"| {name} | {interp} | {k} | " + " | ".join(cells) + " |")
    print("\n| problem | entries within 0.0015 |\n|---|---|")
    for name, (ok, tot) in stats.items():
        print(f"| {name} | {ok}/{tot} |")

    # random field unknown: ensemble over seeds for the bilinear k = 1 row
    print("\nRandom diffusion, bilinear, k = 1, over 10 seeds of our field "
          "(paper: 0.900 0.891 0.900 0.895 0.866):\n")
    vals = np.array([factors(setup("Random diffusion", m, 1, A=pb.random_diffusion(m, seed=sd)))
                     for sd in range(10)])
    for lab, col in zip(("M", "M'", "S", "N", "rho(K)"), vals.T):
        print(f"- {lab}: min {col.min():.3f}, median {np.median(col):.3f}, max {col.max():.3f}")


if __name__ == "__main__":
    main()
