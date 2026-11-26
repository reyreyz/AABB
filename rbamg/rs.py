"""Ruge-Stueben coarsening and direct interpolation, ported from PyAMG's
amg_core (rs_cf_splitting, rs_direct_interpolation_pass1/2) so that the
tie-breaking, and hence the C/F splitting, follows the PyAMG code path.

Matrices are given as dense numpy arrays (n <= ~4000 here) and converted to
per-row sorted index lists, which is the CSR order PyAMG iterates in.
"""

import numpy as np

F_NODE, C_NODE, U_NODE = 0, 1, 2


def rows_of(A):
    return [np.nonzero(A[i])[0] for i in range(A.shape[0])]


def classical_strength(A, theta=0.25, norm="abs"):
    """S_ij = 1 if |a_ij| >= theta max_{k != i} |a_ik| (norm 'abs') or
    -a_ij >= theta max_{k != i} (-a_ik) (norm 'min'); diagonal excluded."""
    n = A.shape[0]
    S = np.zeros_like(A, dtype=bool)
    for i in range(n):
        r = A[i].copy()
        r[i] = 0.0
        v = np.abs(r) if norm == "abs" else -r
        mx = v.max()
        if mx <= 0:
            continue
        S[i] = (v >= theta * mx) & (r != 0)
    return S


def rs_cf_splitting(S):
    """PyAMG amg_core::rs_cf_splitting (first pass only). S boolean, no diagonal."""
    n = S.shape[0]
    Sr = rows_of(S)           # S_i : points i strongly depends on
    Tr = rows_of(S.T)         # S^T_i : points strongly depending on i
    lam = np.array([len(t) for t in Tr])
    interval_ptr = np.zeros(n + 1, dtype=int)
    interval_count = np.zeros(n + 1, dtype=int)
    index_to_node = np.zeros(n, dtype=int)
    node_to_index = np.zeros(n, dtype=int)
    for i in range(n):
        interval_count[lam[i]] += 1
    cumsum = 0
    for i in range(n):
        interval_ptr[i] = cumsum
        cumsum += interval_count[i]
        interval_count[i] = 0
    for i in range(n):
        li = lam[i]
        index = interval_ptr[li] + interval_count[li]
        index_to_node[index] = i
        node_to_index[i] = index
        interval_count[li] += 1
    split = np.full(n, U_NODE)
    for i in range(n):
        if lam[i] == 0 or (lam[i] == 1 and Tr[i][0] == i):
            split[i] = F_NODE
    for top in range(n - 1, -1, -1):
        i = index_to_node[top]
        li = lam[i]
        interval_count[li] -= 1
        if split[i] == F_NODE:
            continue
        split[i] = C_NODE
        for j in Tr[i]:
            if split[j] == U_NODE:
                split[j] = F_NODE
                for k in Sr[j]:
                    if split[k] == U_NODE:
                        if lam[k] >= n - 1:
                            continue
                        lk = lam[k]
                        old, new = node_to_index[k], interval_ptr[lk] + interval_count[lk] - 1
                        node_to_index[index_to_node[old]] = new
                        node_to_index[index_to_node[new]] = old
                        index_to_node[old], index_to_node[new] = index_to_node[new], index_to_node[old]
                        interval_count[lk] -= 1
                        interval_count[lk + 1] += 1
                        interval_ptr[lk + 1] = new
                        lam[k] += 1
        for j in Sr[i]:
            if split[j] == U_NODE:
                if lam[j] == 0:
                    continue
                lj = lam[j]
                old, new = node_to_index[j], interval_ptr[lj]
                node_to_index[index_to_node[old]] = new
                node_to_index[index_to_node[new]] = old
                index_to_node[old], index_to_node[new] = index_to_node[new], index_to_node[old]
                interval_count[lj] -= 1
                interval_count[lj - 1] += 1
                interval_ptr[lj] += 1
                interval_ptr[lj - 1] = interval_ptr[lj] - interval_count[lj - 1]
                lam[j] -= 1
    return split


def direct_interpolation(A, S, split):
    """PyAMG rs_direct_interpolation (natural ordering). Returns P (n x nc) and
    the list of C points in increasing order (their coarse numbering)."""
    n = A.shape[0]
    C = np.nonzero(split == C_NODE)[0]
    cnum = -np.ones(n, dtype=int)
    cnum[C] = np.arange(len(C))
    P = np.zeros((n, len(C)))
    for i in range(n):
        if split[i] == C_NODE:
            P[i, cnum[i]] = 1.0
            continue
        sj = [j for j in np.nonzero(S[i])[0] if split[j] == C_NODE and j != i]
        sp = sum(A[i, j] for j in sj if A[i, j] > 0)
        sn = sum(A[i, j] for j in sj if A[i, j] < 0)
        r = A[i].copy()
        diag = r[i]
        r[i] = 0
        ap, an = r[r > 0].sum(), r[r < 0].sum()
        alpha = an / sn if sn != 0 else 0.0
        beta = ap / sp if sp != 0 else 0.0
        if sp == 0:
            diag += ap
            beta = 0.0
        for j in sj:
            P[i, cnum[j]] = (-alpha if A[i, j] < 0 else -beta) * A[i, j] / diag
    return P, C


def ruge_stuben(A, theta=0.25, norm="abs"):
    """Returns perm (C first, then F, each in natural order), P = [I; W] in C/F order."""
    S = classical_strength(A, theta, norm)
    split = rs_cf_splitting(S)
    P, C = direct_interpolation(A, S, split)
    F = np.nonzero(split != C_NODE)[0]
    perm = np.concatenate([C, F])
    return perm, P[perm], S


# ---------------------------------------------------------------------------
# PyAMG gallery stencils
# ---------------------------------------------------------------------------
def diffusion_stencil_2d(epsilon=1.0, theta=0.0, type="FE"):
    eps, C, S = float(epsilon), np.cos(theta), np.sin(theta)
    CS, CC, SS = C * S, C * C, S * S
    if type == "FE":
        a = (-1 * eps - 1) * CC + (-1 * eps - 1) * SS + (3 * eps - 3) * CS
        b = (2 * eps - 4) * CC + (-4 * eps + 2) * SS
        c = (-1 * eps - 1) * CC + (-1 * eps - 1) * SS + (-3 * eps + 3) * CS
        d = (-4 * eps + 2) * CC + (2 * eps - 4) * SS
        e = (8 * eps + 8) * CC + (8 * eps + 8) * SS
        return np.array([[a, b, c], [d, e, d], [c, b, a]]) / 6.0
    a = 0.5 * (-1 * eps + 1) * CS
    b = -(eps * SS + CC)
    c = -a
    d = -(eps * CC + SS)
    e = 2.0 * (eps + 1)
    return np.array([[a, b, c], [d, e, d], [c, b, a]])


def stencil_grid(st, m):
    """Dirichlet m x m grid; st[di+1, dj+1] couples (i, j) to (i+di, j+dj), index i*m + j."""
    A = np.zeros((m * m, m * m))
    for i in range(m):
        for j in range(m):
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    ii, jj = i + di, j + dj
                    if 0 <= ii < m and 0 <= jj < m and st[di + 1, dj + 1] != 0:
                        A[i * m + j, ii * m + jj] += st[di + 1, dj + 1]
    return A
