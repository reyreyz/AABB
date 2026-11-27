"""Robust non-Markovian multilevel cycle, O(n) variant (simplest self-contained code).

Problem (section 7.3): -eps Lap u + beta . grad u on (-1,1)^2, centered differences,
beta = (2y(1-x^2), -2x(1-y^2)), cell Peclet |beta| h / eps = 20.  Damped Jacobi and the
Markovian cycle diverge on this problem.

Per level (C/F order, P0 = [I; 0], Q0 = [0; I], Jacobi omega = 1, T = I - D^{-1} A):
    F-points greedy:     sum_{j in F, j != i} |a_ij| <= delta |a_ii|   for i in F
    CR propagator:       K = Q0^T T Q0 = I - D_ff^{-1} A_ff
    coupling:            H = Q0^T T P0 = -D_ff^{-1} A_fc
    memory sum:          S = sum_{l<k} K^l H            (entries < tauP * row max dropped)
    effective prolong.:  P' = P0 + Q0 S = [I; S]
    restriction:         R  = P'^T                       (Galerkin with P'; Thm 4.5 allows any R)
    coarse operator:     A_c = R A P' = (A_cc + A_cf S) + S^T (A_fc + A_ff S)
    next level:          A_c filtered (row-max rule, lumped), tau = 1e-3 (levels < 8), 3e-4 below
    coarsest:            sparse LU once n_l < 16 sqrt(n_0)
Cycle (Algorithm 4.2 with the shift relation of Algorithm 7.3):
    x^0..x^k Jacobi;  d^l = (x^k - x^l)_C;  y = 0, y <- K y + H d^l (l < k);  h = Q0 y
    r   = b - A (x^k + h)
    r_c = R r + A_c d^0                                  (A_c exact, applied matrix-free)
    e   = cycle(r_c);  on levels 0-1 one refinement: e += cycle(r_c - A_c e)
    z   = e - d^0;     return x^k + h + P' z
"""

import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


# --------------------------------------------------------------------------- problem
def recirc(m, pe=20.0):
    h = 2.0 / (m + 1)
    t = -1 + h * np.arange(1, m + 1)
    X, Y = np.meshgrid(t, t)                                  # x = column index (fast)
    bx, by = (2 * Y * (1 - X ** 2)).ravel(), (-2 * X * (1 - Y ** 2)).ravel()
    eps = np.hypot(bx, by).max() * h / pe                     # cell Peclet |beta| h / eps = pe
    n = m * m
    I, J = np.divmod(np.arange(n), m)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [4 * eps * np.ones(n)]
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        rows.append(np.nonzero(ok)[0])
        cols.append(((I + di) * m + J + dj)[ok])
        vals.append((-eps + s * 0.5 * h * b)[ok])               # h^2 (-eps u_xx + b u_x), centered
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n))


# --------------------------------------------------------------------------- helpers
def greedy_f_points(A, delta):
    """F-set with ||I - D_ff^{-1} A_ff||_inf <= delta (greedy, natural order)."""
    Ac = abs(A).tocsc()
    d = np.abs(A.diagonal())
    isF, load = np.zeros(A.shape[0], bool), np.zeros(A.shape[0])
    for i in range(A.shape[0]):
        if load[i] > delta * d[i]:
            continue
        ks, vs = Ac.indices[Ac.indptr[i]:Ac.indptr[i + 1]], Ac.data[Ac.indptr[i]:Ac.indptr[i + 1]]
        off = ks != i
        nb = off & isF[ks]
        if np.any(load[ks[nb]] + vs[nb] > delta * d[ks[nb]]):
            continue
        isF[i] = True
        load[ks[off]] += vs[off]
    return isF


def drop_rel(X, tau):
    """drop entries below tau * (row max |x|)."""
    X = sp.coo_matrix(X)
    mx = np.zeros(X.shape[0])
    np.maximum.at(mx, X.row, np.abs(X.data))
    keep = np.abs(X.data) >= tau * mx[X.row]
    return sp.csr_matrix((X.data[keep], (X.row[keep], X.col[keep])), X.shape)


def filter_lumped(A, tau):
    """drop off-diagonals below tau * (row max off-diagonal), add them to the diagonal."""
    C = sp.coo_matrix(A)
    off = C.row != C.col
    mx = np.zeros(A.shape[0])
    np.maximum.at(mx, C.row[off], np.abs(C.data[off]))
    small = off & (np.abs(C.data) < tau * mx[C.row])
    lump = np.bincount(C.row[small], weights=C.data[small], minlength=A.shape[0])
    B = sp.csr_matrix((C.data[~small], (C.row[~small], C.col[~small])), A.shape)
    return (B + sp.diags(lump)).tocsr()


# --------------------------------------------------------------------------- hierarchy
class Level:
    pass


class RobustN:
    def __init__(self, A, k0=4, k1=8, delta=0.3, tauP=1e-5, tau_up=1e-3, tau_deep=3e-4,
                 deep_from=8, refine=1, coarsest_factor=16):
        self.nnz0, self.refine = A.nnz, refine
        n_coarsest = coarsest_factor * np.sqrt(A.shape[0])
        self.levels = []
        A = sp.csr_matrix(A)
        while True:
            L = Level()
            L.A, L.D, lev = A, A.diagonal(), len(self.levels)
            L.k = k0 if lev == 0 else k1
            self.levels.append(L)
            if A.shape[0] < n_coarsest:                                    # exact coarse solve
                L.lu = sla.splu(sp.csc_matrix(A))
                L.lu_nnz = L.lu.L.nnz + L.lu.U.nnz
                break
            isF = greedy_f_points(A, delta)
            C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
            Acc, Acf = A[C][:, C], A[C][:, F]
            Afc, Aff = A[F][:, C], A[F][:, F]
            Dfi = sp.diags(1.0 / Aff.diagonal())
            K = (sp.eye(len(F)) - Dfi @ Aff).tocsr()                       # K = I - D_ff^-1 A_ff
            H = (-Dfi @ Afc).tocsr()                                       # H = -D_ff^-1 A_fc
            S, V = H.copy(), H.copy()
            for _ in range(L.k - 1):                                       # S = sum_{l<k} K^l H
                V = (K @ V).tocsr()
                S = (S + V).tocsr()
            if lev >= 1:
                S = drop_rel(S, tauP)                                      # bounded P'
            APc = (Acc + Acf @ S).tocsr()                                  # (A P')_C
            APf = (Afc + Aff @ S).tocsr()                                  # (A P')_F
            Asig = (APc + S.T @ APf).tocsr()                               # A_c = P'^T A P'
            L.C, L.F, L.K, L.H, L.S, L.ST = C, F, K, H, S, S.T.tocsr()
            L.APc, L.APf = APc, APf
            tau = tau_deep if lev + 1 >= deep_from else tau_up
            A = filter_lumped(Asig, tau)                                   # next level operator

    def RAP(self, L, v):
        """A_c v = P'^T A P' v, exact and matrix-free."""
        return L.APc @ v + L.ST @ (L.APf @ v)

    def cycle(self, lev, b, x0, work):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            work[0] += L.lu_nnz / self.nnz0
            return L.lu.solve(b)
        A, C, F = L.A, L.C, L.F
        xs = [x0]
        for _ in range(L.k):                                               # Jacobi, omega = 1
            xs.append(xs[-1] + (b - A @ xs[-1]) / L.D)
        xk = xs[-1]
        y = np.zeros(len(F))
        for l in range(L.k):                                               # y <- K y + H d^l
            y = L.K @ y + L.H @ (xk[C] - xs[l][C])
        d0 = xk[C] - x0[C]
        xh = xk.copy()
        xh[F] += y                                                         # x^k + h, h = Q0 y
        r = b - A @ xh
        rc = r[C] + L.ST @ r[F] + self.RAP(L, d0)                          # R r + A_c d^0
        e = self.cycle(lev + 1, rc, np.zeros(len(C)), work)
        if lev <= self.refine and not hasattr(self.levels[lev + 1], "lu"):
            e = e + self.cycle(lev + 1, rc - self.RAP(L, e), np.zeros(len(C)), work)
            work[0] += (L.APc.nnz + L.APf.nnz + L.S.nnz) / self.nnz0
        z = e - d0
        xh[C] += z                                                         # + P' z = [z; S z]
        xh[F] += L.S @ z
        work[0] += ((L.k + 1) * A.nnz + L.k * (L.K.nnz + L.H.nnz)
                    + 3 * L.S.nnz + 2 * (L.APc.nnz + L.APf.nnz)) / self.nnz0
        return xh

    def solve(self, b, tol=1e-8, maxcyc=30):
        A = self.levels[0].A
        x, nb, work, hist = np.zeros(A.shape[0]), np.linalg.norm(b), [0.0], [1.0]
        for _ in range(maxcyc):
            x = self.cycle(0, b, x, work)
            hist.append(np.linalg.norm(b - A @ x) / nb)
            if hist[-1] < tol:
                break
        return x, np.array(hist), work[0]

    def complexity(self):
        op = sum(L.A.nnz for L in self.levels) / self.nnz0
        st = sum(L.A.nnz + (L.lu_nnz if hasattr(L, "lu") else L.K.nnz + L.H.nnz + 2 * L.S.nnz
                            + L.APc.nnz + L.APf.nnz) for L in self.levels) / self.nnz0
        return op, st


if __name__ == "__main__":
    import sys
    for m in [int(a) for a in sys.argv[1:]] or [127, 255]:
        A = recirc(m)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        t = time.time()
        ml = RobustN(A)
        ts = time.time() - t
        x, h, w = ml.solve(b)
        op, st = ml.complexity()
        cyc = len(h) - 1
        fac = (h[-1] / h[1]) ** (1 / max(cyc - 1, 1))
        print(f"m={m} n={A.shape[0]} levels={len(ml.levels)} factor={fac:.2e} cycles={cyc} "
              f"WU/cycle={w/cyc:.0f} WU={w:.0f} opcx={op:.2f} stored={st:.0f} setup={ts:.1f}s")
