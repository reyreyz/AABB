"""Multilevel cycles of section 7 (Algorithms 7.1, 7.2, 7.3), scipy.sparse.

Every coarse solve is a recursive call with zero initial guess; the coarsest
system (< 150 unknowns) is factored.  Work is counted in units of nnz(A) of the
finest operator (WU): one product with a sparse matrix X costs nnz(X)/nnz(A0).
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from pyamg.strength import classical_strength_of_connection
from pyamg.classical.split import RS
from pyamg.classical.interpolate import one_point_interpolation

COARSEST = 150


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def filter_diag(A, tau):
    """Drop |a_ij| < tau |a_ii| (j != i) and add them to a_ii (Alg. 7.1)."""
    A = sp.csr_matrix(A, copy=True)
    d = A.diagonal()
    C = A.tocoo()
    off = C.row != C.col
    small = off & (np.abs(C.data) < tau * np.abs(d[C.row]))
    lump = np.bincount(C.row[small], weights=C.data[small], minlength=A.shape[0])
    keep = ~small
    B = sp.csr_matrix((C.data[keep], (C.row[keep], C.col[keep])), shape=A.shape)
    return (B + sp.diags(lump.astype(float))).tocsr()


def filter_rowmax(A, tau):
    """Drop off-diagonal entries below tau * (largest off-diagonal magnitude in the
    row) and add them to the diagonal (Algs. 7.2, 7.3)."""
    if tau <= 0:
        return sp.csr_matrix(A)
    C = sp.coo_matrix(A)
    off = C.row != C.col
    mx = np.zeros(A.shape[0])
    np.maximum.at(mx, C.row[off], np.abs(C.data[off]))
    small = off & (np.abs(C.data) < tau * mx[C.row])
    lump = np.bincount(C.row[small], weights=C.data[small], minlength=A.shape[0])
    keep = ~small
    B = sp.csr_matrix((C.data[keep], (C.row[keep], C.col[keep])), shape=A.shape)
    return (B + sp.diags(lump.astype(float))).tocsr()


def greedy_f_points(A, delta=0.3):
    """F-points greedily (natural order) such that sum_{j in F\\{i}} |a_ij| <= delta |a_ii|
    for every F-point i, i.e. ||I - D_ff^{-1} A_ff||_inf <= delta."""
    A = sp.csr_matrix(A)
    n = A.shape[0]
    R = abs(A).tocsr()
    Cc = abs(A).tocsc()
    rp, ri, rv = R.indptr, R.indices, R.data
    cp, ci, cv = Cc.indptr, Cc.indices, Cc.data
    d = np.abs(A.diagonal())
    lim = delta * d
    isF = np.zeros(n, dtype=bool)
    load = np.zeros(n)                 # sum_{j in F, j != i} |a_ij|
    for i in range(n):
        if load[i] > lim[i]:
            continue
        s, e = cp[i], cp[i + 1]        # column i: entries a_ki
        ks, vs = ci[s:e], cv[s:e]
        nb = (ks != i) & isF[ks]
        if np.any(load[ks[nb]] + vs[nb] > lim[ks[nb]]):
            continue
        isF[i] = True
        off = ks != i
        load[ks[off]] += vs[off]
    return isF


def cf_blocks(A, isF):
    C = np.nonzero(~isF)[0]
    F = np.nonzero(isF)[0]
    A = sp.csr_matrix(A)
    return C, F, A[C][:, C], A[C][:, F], A[F][:, C], A[F][:, F]


# ---------------------------------------------------------------------------
# Algorithm 7.1: Markovian V-cycle with relaxation-generated prolongation
# ---------------------------------------------------------------------------
class Alg71:
    def __init__(self, A, theta=0.25, sweeps_P=3, sweeps_relax=3, tau=1e-4, nnz0=None, norm='abs', second_pass=False, free_first=True):
        self.A = sp.csr_matrix(A)
        self.n = A.shape[0]
        self.nnz0 = nnz0 or self.A.nnz
        self.sweeps = sweeps_relax
        self.free_first = free_first   # count the first sweep from a zero guess as free
        if self.n < COARSEST:
            self.lu = sla.splu(sp.csc_matrix(filter_diag(self.A, tau) if nnz0 else self.A))
            return
        At = filter_diag(self.A, tau)
        if nnz0 is not None:          # coarse levels relax with the filtered operator
            self.A = At
        S = classical_strength_of_connection(At, theta=theta, norm=norm)
        split = RS(S, second_pass=second_pass)
        P1 = one_point_interpolation(At, S, split)
        self.C = np.nonzero(split == 1)[0]
        self.F = np.nonzero(split == 0)[0]
        Atc = At.tocsr()
        Aff, Afc = Atc[self.F][:, self.F], Atc[self.F][:, self.C]
        Acf, Acc = Atc[self.C][:, self.F], Atc[self.C][:, self.C]
        dff = Aff.diagonal()
        W = sp.csr_matrix(P1[self.F])
        for _ in range(sweeps_P):                       # W <- W - D_ff^{-1}(A_ff W + A_fc)
            W = (W - sp.diags(1 / dff) @ (Aff @ W + Afc)).tocsr()
        self.W = W
        Asig = (Acc + Acf @ W).tocsr()                  # injection: A_sigma = A_cc + A_cf W
        # relaxation on this level uses A on the finest level, A~ below
        self.AF = self.A[self.F].tocsr()                # F rows for F-relaxation
        self.Aff_l = self.A[self.F][:, self.F].tocsr()
        self.Afc_l = self.A[self.F][:, self.C].tocsr()
        self.AC = self.A[self.C].tocsr()                # C rows for the residual
        self.dff = self.A.diagonal()[self.F]
        self.coarse = Alg71(Asig, theta, sweeps_P, sweeps_relax, tau, self.nnz0, norm, second_pass, free_first)
        self.coarse_uses_filtered = True

    def cycle(self, b, x, work):
        if self.n < COARSEST:
            return self.lu.solve(b)
        F, C = self.F, self.C
        zero = not np.any(x)
        # F-relaxation, omega = 1.  x_c is fixed during the sweeps, so
        # g = b_f - A_fc x_c is formed once and each sweep costs nnz(A_ff).
        g = b[F] - self.Afc_l @ x[C]
        if not zero:
            work[0] += self.Afc_l.nnz / self.nnz0
        for it in range(self.sweeps):
            x[F] += (g - self.Aff_l @ x[F]) / self.dff
            if not (self.free_first and zero and it == 0):
                work[0] += self.Aff_l.nnz / self.nnz0
        rc = b[C] - self.AC @ x                         # C-part of the residual only
        work[0] += self.AC.nnz / self.nnz0
        e = self.coarse.cycle(rc, np.zeros(len(C)), work)
        x[C] += e
        x[F] += self.W @ e
        work[0] += self.W.nnz / self.nnz0
        return x

    def levels(self):
        return [self] + (self.coarse.levels() if self.n >= COARSEST else [])

    def op_complexity(self):
        return sum(l.A.nnz for l in self.levels()) / self.nnz0


# ---------------------------------------------------------------------------
# Algorithm 7.2: semi-Markovian V-cycle
# ---------------------------------------------------------------------------
class Alg72:
    def __init__(self, A, k=4, omega=0.7, tau=1e-3, delta=0.3, nnz0=None, markov=False):
        self.markov = markov          # True: Petrov-Galerkin cycle with the same operators
        self.A = sp.csr_matrix(A)
        self.n = A.shape[0]
        self.nnz0 = nnz0 or self.A.nnz
        self.k, self.omega = k, omega
        if self.n < COARSEST:
            self.lu = sla.splu(sp.csc_matrix(self.A))
            return
        isF = greedy_f_points(self.A, delta)
        self.C, self.F, Acc, Acf, Afc, Aff = cf_blocks(self.A, isF)
        dff = Aff.diagonal()
        self.D = self.A.diagonal()
        Dfi = sp.diags(1 / dff)
        self.K = (sp.eye(len(self.F)) - omega * Dfi @ Aff).tocsr()      # K = I - w D_ff^{-1} A_ff
        self.H = (-omega * Dfi @ Afc).tocsr()                            # H = -w D_ff^{-1} A_fc
        G = (sp.eye(len(self.F)) - Dfi @ Aff).tocsr()
        Z = (sp.eye(len(self.F)) + G + G @ G) @ Dfi
        self.Rf = (-Acf @ Z).tocsr()                                     # Rhat = [I  -A_cf(I+G+G^2)D_ff^{-1}]
        self.Asig = (Acc + self.Rf @ Afc).tocsr()                        # Rhat A P0 (exact)
        self.coarse = Alg72(filter_rowmax(self.Asig, tau), k, omega, tau, delta, self.nnz0, markov)

    def cycle(self, b, x0, work):
        if self.n < COARSEST:
            return self.lu.solve(b)
        k, w, A, C, F = self.k, self.omega, self.A, self.C, self.F
        xs = [x0]
        for _ in range(k):                                               # Jacobi on all points
            xs.append(xs[-1] + w * (b - A @ xs[-1]) / self.D)
            work[0] += A.nnz / self.nnz0
        xk = xs[-1]
        if self.markov:
            r = b - A @ xk
            e = self.coarse.cycle(r[C] + self.Rf @ r[F], np.zeros(len(C)), work)
            out = xk.copy()
            out[C] += e
            return out
        d = [xk[C] - xs[l][C] for l in range(k)]
        y = np.zeros(len(F))
        for l in range(k):
            y = self.K @ y + self.H @ d[l]
            work[0] += (self.K.nnz + self.H.nnz) / self.nnz0
        r = b - A @ xk
        rc = r[C] + self.Rf @ r[F] + self.Asig @ d[0]                    # Rhat r + A_sigma d^(0)
        work[0] += (A.nnz + self.Rf.nnz + self.Asig.nnz) / self.nnz0
        e0 = self.coarse.cycle(rc, np.zeros(len(C)), work)
        z = e0 - d[0]
        out = xk.copy()
        out[C] += z
        v = self.H @ z
        s = v.copy()
        for _ in range(k - 1):                                           # sum_{l<k} K^l H z
            v = self.K @ v
            s += v
        work[0] += (self.H.nnz + (k - 1) * self.K.nnz) / self.nnz0
        out[F] += y + s
        return out

    def levels(self):
        return [self] + (self.coarse.levels() if self.n >= COARSEST else [])

    def op_complexity(self):
        return sum(l.A.nnz for l in self.levels()) / self.nnz0


# ---------------------------------------------------------------------------
# Algorithm 7.3: non-Markovian V-cycle
# ---------------------------------------------------------------------------
class Alg73:
    def __init__(self, A, level=0, k0=8, k1=12, delta=0.3, nnz0=None, markov=False):
        self.markov = markov
        self.A = sp.csr_matrix(A)
        self.n = A.shape[0]
        self.nnz0 = nnz0 or self.A.nnz
        self.level = level
        self.k = k0 if level == 0 else k1
        if self.n < COARSEST:
            self.lu = sla.splu(sp.csc_matrix(self.A))
            return
        isF = greedy_f_points(self.A, delta)
        self.C, self.F, Acc, Acf, Afc, Aff = cf_blocks(self.A, isF)
        self.D = self.A.diagonal()
        Dfi = sp.diags(1 / Aff.diagonal())
        self.K = (sp.eye(len(self.F)) - Dfi @ Aff).tocsr()               # omega = 1
        self.H = (-Dfi @ Afc).tocsr()
        W = self.H.copy()
        V = self.H.copy()
        for _ in range(self.k - 1):                                      # W = sum_{l<k} K^l H, exact
            V = (self.K @ V).tocsr()
            W = W + V
        self.W = W.tocsr()
        self.Acc, self.Acf = Acc, Acf
        Asig = (Acc + Acf @ self.W).tocsr()                              # (A [I; W])_c assembled
        tau = 1e-3 if level < 8 else 3e-5
        self.coarse = Alg73(filter_rowmax(Asig, tau), level + 1, k0, k1, delta, self.nnz0, markov)

    def calA(self, v, work):
        """matrix-free (A [I; W] v)_c"""
        work[0] += (self.W.nnz + self.Acc.nnz + self.Acf.nnz) / self.nnz0
        return self.Acc @ v + self.Acf @ (self.W @ v)

    def cycle(self, b, x0, work):
        if self.n < COARSEST:
            return self.lu.solve(b)
        k, A, C, F = self.k, self.A, self.C, self.F
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + (b - A @ xs[-1]) / self.D)
            work[0] += A.nnz / self.nnz0
        xk = xs[-1]
        if self.markov:
            r = b - A @ xk
            e = self.coarse.cycle(r[C], np.zeros(len(C)), work)
            out = xk.copy()
            out[C] += e
            out[F] += self.W @ e
            return out
        d = [xk[C] - xs[l][C] for l in range(k)]
        y = np.zeros(len(F))
        for l in range(k):
            y = self.K @ y + self.H @ d[l]
            work[0] += (self.K.nnz + self.H.nnz) / self.nnz0
        h = np.zeros(self.n)
        h[F] = y
        r = b - A @ (xk + h)
        work[0] += A.nnz / self.nnz0
        rc = r[C] + self.calA(d[0], work)
        e0 = self.coarse.cycle(rc, np.zeros(len(C)), work)
        if self.level <= 3:                                              # one refinement step
            e0 = e0 + self.coarse.cycle(rc - self.calA(e0, work), np.zeros(len(C)), work)
        z = e0 - d[0]
        out = xk + h
        out[C] += z
        out[F] += self.W @ z
        work[0] += self.W.nnz / self.nnz0
        return out

    def levels(self):
        return [self] + (self.coarse.levels() if self.n >= COARSEST else [])

    def op_complexity(self):
        return sum(l.A.nnz for l in self.levels()) / self.nnz0


def solve(ml, A, b, tol=1e-10, maxcyc=60):
    """Stationary iteration from x0 = 0; returns residual history and work."""
    x = np.zeros(A.shape[0])
    nb = np.linalg.norm(b)
    hist, works, work = [1.0], [0.0], [0.0]
    for _ in range(maxcyc):
        x = ml.cycle(b, x, work)
        r = np.linalg.norm(b - A @ x) / nb
        hist.append(r)
        works.append(work[0])
        if r < tol or not np.isfinite(r) or r > 1e30:
            break
    return np.array(hist), np.array(works)
