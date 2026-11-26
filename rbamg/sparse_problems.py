"""Sparse (scipy) convection-diffusion problems on (-1,1)^2, m x m interior grid,
Dirichlet boundary, natural ordering index = i*m + j (x = j fast, y = i).

    -eps Lap u + beta(x, y) . grad u,   h = 2/(m+1),  scaled by h^2.
Differencing of the advection term: 'upwind' (first order) or 'centered'.
"""

import numpy as np
import scipy.sparse as sp


def grid(m):
    h = 2.0 / (m + 1)
    t = -1 + h * np.arange(1, m + 1)
    X, Y = np.meshgrid(t, t)          # X[i, j] = t[j], Y[i, j] = t[i]
    return h, X.ravel(), Y.ravel()


def cd(m, beta, eps, scheme="centered"):
    h, x, y = grid(m)
    bx, by = beta(x, y)
    bx = np.broadcast_to(bx, x.shape).astype(float)
    by = np.broadcast_to(by, x.shape).astype(float)
    n = m * m
    idx = np.arange(n).reshape(m, m)
    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(r), cols.append(c), vals.append(v)

    diag = 4 * eps * np.ones(n)
    # neighbours: (di, dj, velocity component, sign of offset)
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        I, J = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
        ok = (0 <= I + di) & (I + di < m) & (0 <= J + dj) & (J + dj < m)
        p = idx[I, J].ravel()
        okr = ok.ravel()
        q = np.where(okr, idx[np.clip(I + di, 0, m - 1), np.clip(J + dj, 0, m - 1)].ravel(), -1)
        if scheme == "centered":
            coef = -eps + s * 0.5 * h * b          # h^2 * (b u_x) with (u_{+}-u_{-})/(2h)
        else:
            # upwind: b>0 uses the minus neighbour, b<0 the plus neighbour
            up = (s < 0) & (b > 0) | (s > 0) & (b < 0)
            coef = -eps - np.where(up, h * np.abs(b), 0.0)
            diag = diag + np.where(up, h * np.abs(b), 0.0) * 0 + 0
        add(p[okr], q[okr], coef[okr])
    if scheme == "upwind":
        diag = diag + h * (np.abs(bx) + np.abs(by))
    rows.append(np.arange(n)), cols.append(np.arange(n)), vals.append(diag)
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    A.sum_duplicates()
    A.eliminate_zeros()
    return A


def eps_for_peclet(m, beta, pe, norm="inf", factor2=True):
    """eps such that the maximal cell Peclet number max |beta| h / (2 eps) equals pe
    (or max |beta| h / eps when factor2 is False)."""
    h, x, y = grid(m)
    bx, by = beta(x, y)
    bx, by = np.broadcast_to(bx, x.shape), np.broadcast_to(by, x.shape)
    b = np.maximum(np.abs(bx), np.abs(by)) if norm == "inf" else np.hypot(bx, by)
    return b.max() * h / ((2.0 if factor2 else 1.0) * pe)


STAGNATION = lambda x, y: (x, -y)
SHEAR = lambda x, y: (-y, 0.3 + 0 * x)
RECIRC = lambda x, y: (2 * y * (1 - x ** 2), -2 * x * (1 - y ** 2))
