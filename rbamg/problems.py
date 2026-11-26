"""Model problems on an m x m interior grid of the unit square, in C/F ordering.

Coarse points are the grid points with even 0-based indices in both directions
((m+1)/2 per direction, m odd), so n_c = ((m+1)/2)^2; for m = 15 that is
n_c = 64, n_f = 161, matching the supplement.

Only Poisson is stated in the paper.  CD and the shifted Laplacian were
identified by matching Table SM1 (see results_sm1.md); the random-diffusion
field is a stand-in.
"""

import numpy as np


def lap1d(m):
    return 2 * np.eye(m) - np.eye(m, k=1) - np.eye(m, k=-1)


def poisson(m):
    """5-point Laplacian, Dirichlet, scaled by h^2."""
    I = np.eye(m)
    return np.kron(I, lap1d(m)) + np.kron(lap1d(m), I)


def shifted_laplacian(m, n_neg=3):
    """h^2-scaled Poisson minus a shift midway between its eigenvalues n_neg and
    n_neg + 1 (three negative modes).  Identified from Table SM1: this midpoint
    reproduces even the divergent P1 entries (M' = 2.380, N = 2.542 at k = 4),
    whereas s = 1/16, equal to three digits, gives 2.371 and 2.532."""
    A = poisson(m)
    lam = np.sort(np.linalg.eigvalsh(A))
    return A - 0.5 * (lam[n_neg - 1] + lam[n_neg]) * np.eye(m * m)


def convection_diffusion(m, eps, theta=np.pi / 6):
    """-eps Lap u + beta . grad u, beta = (cos theta, sin theta), first-order upwind,
    scaled by h^2.  theta = pi/6 and upwinding were identified from Table SM1; it
    also gives the mesh Peclet number max|beta_i| h / (2 eps) = 3.4 at m = 127
    quoted in SM2.1."""
    h = 1.0 / (m + 1)
    I = np.eye(m)
    D = np.eye(m) - np.eye(m, k=-1)            # (u_i - u_{i-1}), upwind for beta_i > 0
    Dx = np.kron(I, D) * np.cos(theta) * h      # x is the inner (fast) index
    Dy = np.kron(D, I) * np.sin(theta) * h
    return eps * poisson(m) + Dx + Dy


def random_diffusion(m, sigma=3.0, seed=0):
    """-div(kappa grad u) with kappa = 10^(sigma u), u ~ U(0,1) per node, and
    edge coefficients the arithmetic means of the two node values.

    STAND-IN: the paper's field is not specified, and none of the nine models
    tried (exp(3g) or 10^(3u), per node/cell/edge, harmonic/arithmetic, several
    seeds and MATLAB-compatible streams) reproduces Table SM1.  This one is the
    closest in behaviour (no scheme diverges, M' < N < M as in the paper)."""
    rng = np.random.default_rng(seed)
    kap = 10.0 ** (sigma * rng.random((m + 2, m + 2)))          # incl. boundary ring
    idx = lambda i, j: i * m + j
    A = np.zeros((m * m, m * m))
    for i in range(m):
        for j in range(m):
            p = idx(i, j)
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = kap[i + 1, j + 1], kap[i + 1 + di, j + 1 + dj]
                c = 0.5 * (a + b)
                A[p, p] += c
                ii, jj = i + di, j + dj
                if 0 <= ii < m and 0 <= jj < m:
                    A[p, idx(ii, jj)] -= c
    return A


def bilinear_cf(m):
    """Bilinear interpolation and the C/F permutation.

    Returns perm (natural index of each C/F position) and P = [I; W] in C/F order.
    """
    assert m % 2 == 1
    mc = (m + 1) // 2
    coarse = lambda i, j: i % 2 == 0 and j % 2 == 0
    nat = lambda i, j: i * m + j
    C = [nat(i, j) for i in range(m) for j in range(m) if coarse(i, j)]
    F = [nat(i, j) for i in range(m) for j in range(m) if not coarse(i, j)]
    cidx = {nat(2 * a, 2 * b): a * mc + b for a in range(mc) for b in range(mc)}

    Pnat = np.zeros((m * m, mc * mc))
    for i in range(m):
        for j in range(m):
            I_ = [i] if i % 2 == 0 else [i - 1, i + 1]      # coarse rows around i
            J_ = [j] if j % 2 == 0 else [j - 1, j + 1]
            w = 1.0 / (len(I_) * len(J_))
            for a in I_:
                for b in J_:
                    Pnat[nat(i, j), cidx[nat(a, b)]] = w
    perm = np.array(C + F)
    return perm, Pnat[perm]


def jacobi(A, omega=2.0 / 3.0):
    return omega * np.diag(1.0 / np.diag(A))


PROBLEMS = {
    "Poisson": lambda m: poisson(m),
    "Random diffusion*": lambda m: random_diffusion(m),
    "CD eps=1e-2": lambda m: convection_diffusion(m, 1e-2),
    "CD eps=1e-3": lambda m: convection_diffusion(m, 1e-3),
    "Shifted Laplacian": lambda m: shifted_laplacian(m),
}
