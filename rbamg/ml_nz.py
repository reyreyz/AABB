"""Multilevel non-Markovian V-cycle with an optimized complement Q = Q0 + P Z (N_Z).

Hierarchy: geometric bilinear interpolation on m = 2^L - 1 grids (coarse points
at even indices), C/F ordering per level, Galerkin restriction Rhat = P^T,
coarse operator A_sigma = P^T A P' (Petrov-Galerkin, as in (4.10)), filtered
(row-max rule, lumped) before it becomes the next level.  Relaxation: damped
Jacobi omega = 2/3, k steps per level.  The coarsest level (< 150) is factored.

Per level, with T = I - omega D^{-1} A in C/F blocks:
    K = T_ff - W T_cf,  H = T_fc + T_ff W - W T_cc - W T_cf W,  K_Z = K + H Z,  P' = P + Q sum_{l<k} K_Z^l H
Cycle (Algorithm 4.2 with complement Q, Qsh = Q0sh, Psh = P0sh - Z Q0sh):
    x^0..x^k by relaxation;  d^l = Psh (x^k - x^l)
    y = 0; for l < k: y <- K_Z y + H d^l;   h = Q y = Q0 y + P Z y
    A_sigma eps = Rhat (b - A x^k) - Rhat A h     (recursive cycle, zero guess)
    return x^k + P' eps + h
Z (pattern: coarse points interpolating to the fine point):
    'Q0'  Z = 0;   'LS' min ||K_Z||_F columnwise (the paper's N_Z);
    'ALS' local alternating LS for ||K_Z^k||_F;
    'Kk'  min ||K_Z^k||_F^2 over the pattern, L-BFGS with sparse gradients.
Markovian reference 'M': x^k + P eps with A_sigma = P^T A P.
"""

import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import minimize
import problems as pb
from ml_cycles import filter_rowmax

COARSEST = 150


def bilinear_1level(m, offset=0):
    """sparse bilinear P (natural ordering, n x nc) on an m x m grid with coarse
    points at indices offset, offset+2, ...; a missing neighbour is the Dirichlet
    boundary (weight dropped).  offset = 0 is the layout of the two-grid
    experiments, offset = 1 (m = 2^L - 1) keeps every level of that form."""
    cidx = list(range(offset, m, 2))
    mc = len(cidx)
    pos = {c: t for t, c in enumerate(cidx)}

    def nb(i):
        if (i - offset) % 2 == 0:
            return [(i, 1.0)]
        return [(c, 0.5) for c in (i - 1, i + 1) if 0 <= c < m]
    rows, cols, vals = [], [], []
    for i in range(m):
        for j in range(m):
            for a, wa in nb(i):
                for b, wb in nb(j):
                    rows.append(i * m + j)
                    cols.append(pos[a] * mc + pos[b])
                    vals.append(wa * wb)
    isC = np.zeros(m * m, dtype=bool)
    isC[[i * m + j for i in cidx for j in cidx]] = True
    return sp.csr_matrix((vals, (rows, cols)), shape=(m * m, mc * mc)), mc, isC


def cd_sparse(m, eps, theta=np.pi / 6):
    h = 1.0 / (m + 1)
    I = sp.eye(m)
    L1 = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(m, m))
    D = sp.diags([-1, 1], [-1, 0], shape=(m, m))
    return (eps * (sp.kron(I, L1) + sp.kron(L1, I)) + h * np.cos(theta) * sp.kron(I, D)
            + h * np.sin(theta) * sp.kron(D, I)).tocsr()


def mask_pattern(W):
    """pattern of Z (nc x nf) = transpose of the pattern of W (nf x nc)."""
    Wt = sp.csr_matrix(W.T)
    Wt.data[:] = 1.0
    return Wt


def z_ls(K, H, patt):
    """column j: min ||K e_j + H z||, z on the pattern."""
    Hc = sp.csc_matrix(H)
    Kc = sp.csc_matrix(K)
    pc = sp.csc_matrix(patt)
    rows, cols, vals = [], [], []
    for j in range(K.shape[1]):
        Cj = pc.indices[pc.indptr[j]:pc.indptr[j + 1]]
        Hj = Hc[:, Cj]
        r = np.unique(np.concatenate([Hj.indices, Kc.indices[Kc.indptr[j]:Kc.indptr[j + 1]]]))
        Hd = Hj[r].toarray()
        kd = Kc[r, j].toarray().ravel()
        z, *_ = np.linalg.lstsq(Hd, -kd, rcond=None)
        rows += list(Cj)
        cols += [j] * len(Cj)
        vals += list(z)
    return sp.csr_matrix((vals, (rows, cols)), shape=patt.shape)


def z_als(K, H, patt, k, sweeps=4, Z0=None):
    Z = z_ls(K, H, patt) if Z0 is None else Z0
    for _ in range(sweeps):
        KZ = (K + H @ Z).tocsr()
        L = sp.eye(K.shape[0], format="csr")
        for _ in range(k - 1):
            L = (L @ KZ).tocsr()
        Z = z_ls((L @ K).tocsr(), (L @ H).tocsr(), patt)
    return Z


def z_kk(K, H, patt, k, Z0, maxiter=100):
    """min ||(K + H Z)^k||_F^2 over Z on the pattern (sparse gradient)."""
    pc = sp.coo_matrix(patt)
    r_, c_ = pc.row, pc.col

    def build(z):
        return sp.csr_matrix((z, (r_, c_)), shape=patt.shape)

    def fg(z):
        KZ = (K + H @ build(z)).tocsr()
        pw = [sp.eye(K.shape[0], format="csr")]
        for _ in range(k):
            pw.append((pw[-1] @ KZ).tocsr())
        Kk = pw[k]
        f = float((Kk.multiply(Kk)).sum())
        G = sp.csr_matrix(K.shape)
        for i in range(k):            # d/dKZ ||KZ^k||^2 = 2 sum_i (KZ^T)^i KZ^k (KZ^T)^{k-1-i}
            G = G + pw[i].T @ Kk @ pw[k - 1 - i].T
        GZ = (2 * H.T @ G).tocsr()
        g = np.asarray(GZ[r_, c_]).ravel()
        return f, g
    z0 = np.asarray(sp.csr_matrix(Z0)[r_, c_]).ravel()
    res = minimize(fg, z0, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter))
    return build(res.x)


class Level:
    pass


class NZCycle:
    def __init__(self, A, m, k=2, omega=2 / 3, zmode="LS", tauA=1e-3, markov=False, nnz0=None, offset=1, max_levels=99,
                 omega_rule="fixed", z_levels=99, pattern="interp", gamma=1):
        self.levels = []
        self.k, self.omega, self.zmode, self.markov = k, omega, zmode, markov
        self.gamma = gamma            # recursive cycles per coarse solve (1: V, 2: W)
        self.nnz0 = nnz0 or A.nnz
        A = sp.csr_matrix(A)
        self.setup_time = 0.0
        t = time.time()
        while True:
            L = Level()
            L.A = A
            self.levels.append(L)
            if A.shape[0] < COARSEST or m < 3 or len(self.levels) >= max_levels:
                L.lu = sla.splu(sp.csc_matrix(A))
                break
            Pn, mc, isC = bilinear_1level(m, offset)
            n = A.shape[0]
            C, F = np.nonzero(isC)[0], np.nonzero(~isC)[0]
            # coarse numbering of C points = row-major on the coarse grid = order of C
            W = Pn[F].tocsr()
            D = A.diagonal()
            om = omega
            if omega_rule == "rho" or (omega_rule == "rho_coarse" and len(self.levels) > 1):
                DA = sp.diags(1 / D) @ A                    # omega_l = omega / rho(D^{-1} A_l)
                lam = abs(sla.eigs(DA, k=1, which="LM", return_eigenvectors=False, tol=1e-4, maxiter=20000)[0])
                om = omega * 1.5 / lam if omega_rule == "rho" else omega * 1.5 / max(lam, 1.5)
            L.omega = om
            T = (sp.eye(n) - om * sp.diags(1 / D) @ A).tocsr()
            Tcc, Tcf, Tfc, Tff = T[C][:, C], T[C][:, F], T[F][:, C], T[F][:, F]
            K = (Tff - W @ Tcf).tocsr()
            H = (Tfc + Tff @ W - W @ Tcc - W @ (Tcf @ W)).tocsr()   # Q0sh T P, P = [I; W]
            if pattern == "interp":
                patt = mask_pattern(W)
            else:                                        # coarse points within distance 3 in A
                G = abs(A).tocsr()
                G.data[:] = 1.0
                G3 = (G @ G @ G).tocsr()
                patt = sp.csr_matrix(G3[F][:, C].T)
                patt.data[:] = 1.0
            if markov or zmode == "Q0" or len(self.levels) > z_levels:
                Z = sp.csr_matrix(patt.shape)
            elif zmode == "LS":
                Z = z_ls(K, H, patt)
            elif zmode == "ALS":
                Z = z_als(K, H, patt, k)
            elif zmode == "Kk":
                Z = z_kk(K, H, patt, k, z_ls(K, H, patt))
            KZ = (K + H @ Z).tocsr()
            # P' = P + Q sum_{l<k} K_Z^l H, Q = Q0 + P Z  (C/F blocks)
            S = H.copy()
            V = H.copy()
            for _ in range(k - 1):
                V = (KZ @ V).tocsr()
                S = S + V
            ZS = (Z @ S).tocsr()                      # C-block of Q S
            WZS = (W @ ZS).tocsr()
            Ppc = sp.eye(len(C)) + ZS                 # P'_C
            Ppf = (W + S + WZS).tocsr()               # P'_F
            if markov:
                Ppc, Ppf = sp.eye(len(C), format="csr"), W
            perm = np.concatenate([C, F])
            inv = np.empty_like(perm)
            inv[perm] = np.arange(n)
            Pp = sp.vstack([Ppc, Ppf]).tocsr()[inv]   # natural ordering
            Pb = Pn                                   # bilinear P, natural
            Asig = (Pb.T @ A @ Pp).tocsr()            # Rhat A P', Rhat = P^T
            L.C, L.F, L.W, L.Z, L.KZ, L.H, L.Pp, L.Pb, L.D = C, F, W, Z, KZ, H, Pp, Pb, D
            L.T = T
            A = filter_rowmax(Asig, tauA)
            A.eliminate_zeros()
            m = mc
        self.setup_time = time.time() - t

    def op_complexity(self):
        return sum(L.A.nnz for L in self.levels) / self.nnz0

    def cycle(self, lev, b, x0, work):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            return L.lu.solve(b)
        A, k, w, nz = L.A, self.k, L.omega, self.nnz0
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + w * (b - A @ xs[-1]) / L.D)
            work[0] += A.nnz / nz
        xk = xs[-1]
        r = b - A @ xk
        work[0] += A.nnz / nz
        if self.markov:
            e = self.cycle(lev + 1, L.Pb.T @ r, np.zeros(len(L.C)), work)
            work[0] += 2 * L.Pb.nnz / nz
            return xk + L.Pb @ e
        C, F = L.C, L.F
        y = np.zeros(len(F))
        for l in range(k):
            dl = xk - xs[l]
            dF = dl[F] - L.W @ dl[C]                      # Q0sh (x^k - x^l)
            d = dl[C] - L.Z @ dF                          # Psh = P0sh - Z Q0sh
            y = L.KZ @ y + L.H @ d
            work[0] += (L.W.nnz + L.Z.nnz + L.KZ.nnz + L.H.nnz) / nz
        h = np.zeros(A.shape[0])
        h[F] = y
        zy = L.Z @ y
        h += L.Pb @ zy                                    # h = Q0 y + P Z y
        rhs = L.Pb.T @ (r - A @ h)
        work[0] += (L.Z.nnz + 2 * L.Pb.nnz + A.nnz) / nz
        Ac = self.levels[lev + 1].A
        e = self.cycle(lev + 1, rhs, np.zeros(len(C)), work)
        for _ in range(self.gamma - 1):
            if hasattr(self.levels[lev + 1], "lu"):
                break
            e = self.cycle(lev + 1, rhs, e, work)
        work[0] += L.Pp.nnz / nz
        return xk + L.Pp @ e + h

    def solve(self, b, tol=1e-10, maxcyc=100):
        A = self.levels[0].A
        x = np.zeros(A.shape[0])
        nb = np.linalg.norm(b)
        hist, works, work = [1.0], [0.0], [0.0]
        for _ in range(maxcyc):
            x = self.cycle(0, b, x, work)
            res = np.linalg.norm(b - A @ x) / nb
            hist.append(res)
            works.append(work[0])
            if res < tol or not np.isfinite(res) or res > 1e8:
                break
        return np.array(hist), np.array(works)


def two_level_check(m=15, k=3, eps=1e-2):
    """With an exact coarse solve, the cycle must reproduce E_N of nz_opt / rbamg."""
    import rbamg as rb
    from reproduce_sm1 import setup
    A = cd_sparse(m, eps)
    cyc = NZCycle(A, m, k=k, zmode="LS", tauA=0.0, offset=0)
    L0 = cyc.levels[0]
    Ac = cyc.levels[1].A
    cyc.levels[1].lu = sla.splu(sp.csc_matrix(Ac))   # exact coarse solve
    cyc.levels = cyc.levels[:2]
    n = A.shape[0]
    E = np.zeros((n, n))
    for i in range(n):
        x0 = -np.eye(n)[i]
        E[:, i] = -cyc.cycle(0, np.zeros(n), x0, [0.0])
    s = setup("CD eps=1e-2" if eps == 1e-2 else "CD eps=1e-3", m, k)
    perm, _ = pb.bilinear_cf(m)
    from reproduce_sm4 import ls_feedback, nonmarkov_with
    Zd = ls_feedback(s)
    KZ = s.K + s.H @ Zd
    Q = s.Q0 + s.P @ Zd
    Pp = s.P + Q @ sum(np.linalg.matrix_power(KZ, l) for l in range(k)) @ s.H
    RA = s.Rhat @ s.A
    Eref = (np.eye(n) - Pp @ np.linalg.solve(RA @ Pp, RA)) @ Q @ np.linalg.matrix_power(KZ, k) @ s.Q0sh
    Ecf = E[np.ix_(perm, perm)]
    return np.linalg.norm(Ecf - Eref, 2) / np.linalg.norm(Eref, 2)


if __name__ == "__main__":
    print("two-level check (multilevel code vs dense E_N):", two_level_check())
