"""N vs its Markovian twin in the multilevel, sparsified setting.

Same hierarchy for both (greedy split, P~ = [I; S~], R = P~^T, A_c = R A P~, filtered next level,
refinement, sparse LU).  Smoother: damped Jacobi T = I - w D^{-1} A.
  N    : k global Jacobi sweeps, exact history (K, H), coarse unknown = step-0 error (Alg 7.3)
  twin : k F-Jacobi sweeps  (T_Q = I - w Q0 Q0^T D^{-1} A), coarse correction with P~
Prediction (any coarse cycle B_c, any filtering):
  E_N e0 - E_twin e0 = -(I - P~ B_c R A) Q (S - S~) d0,   d0 = P^sharp (I - T^k) e0,
i.e. twin's truncation defect acts on c_0 = P^sharp e0, N's on c_k = P^sharp T^k e0.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from robust_n import recirc, greedy_f_points, drop_rel, filter_lumped, Level


def poisson(m, eps_y=1.0):
    T1 = sp.diags([-1, 2, -1], [-1, 0, 1], (m, m))
    return (sp.kron(sp.eye(m), T1) + eps_y * sp.kron(T1, sp.eye(m))).tocsr()


def truncate_lin(S, Bc, s=3):
    """keep the s largest |S_ij| per row, then the minimal change that matches S on the columns of Bc
    (test vectors on C, e.g. 1, x, y):  w = s_k + B_k^T (B_k B_k^T)^+ (B s_i - B_k s_k)."""
    S = sp.csr_matrix(S); rows, cols, vals = [], [], []
    for i in range(S.shape[0]):
        idx, v = S.indices[S.indptr[i]:S.indptr[i + 1]], S.data[S.indptr[i]:S.indptr[i + 1]]
        if len(idx) == 0: continue
        keep = np.argsort(-np.abs(v))[:s]; ik, vk = idx[keep], v[keep]
        Bk = Bc[ik].T                                        # p x s
        target = Bc[idx].T @ v                               # B s_i
        w = vk + Bk.T @ np.linalg.lstsq(Bk @ Bk.T, target - Bk @ vk, rcond=None)[0] if len(ik) > 1 else \
            np.array([target[0]])
        rows += [i] * len(ik); cols += list(ik); vals += list(w)
    return sp.csr_matrix((vals, (rows, cols)), S.shape)


def truncate(S, tau, mode):
    St = drop_rel(S, tau)
    if mode == "lump":                                   # keep row sums: S~ 1 = S 1
        lost = np.asarray(S.sum(1)).ravel() - np.asarray(St.sum(1)).ravel()
        St = St.tolil()
        Sa = abs(St).tocsr()
        for i in np.nonzero(np.abs(lost) > 0)[0]:
            row = Sa.indices[Sa.indptr[i]:Sa.indptr[i + 1]]
            if len(row):
                j = row[np.argmax(Sa.data[Sa.indptr[i]:Sa.indptr[i + 1]])]
                St[i, j] += lost[i]
        St = St.tocsr()
    return St


class ML:
    def __init__(self, A, omega=1.0, k0=4, k1=4, delta=0.3, tauP=0.0, mode="drop", trunc_from=0,
                 tau_filter=1e-3, tau_deep=None, shift=0.0, restrict="ppt", omega_c=None, xy=None, nkeep=3, lin_levels=99, tauP_deep=None, wc_levels=99, tv=None, split0=None, theta=None, omega_coarse=None, deep_from=8, refine=1, coarsest=None, max_levels=30):
        self.nnz0, self.refine, self.omega, self.shift, self.work = A.nnz, refine, omega, shift, 0.0
        self.restrict, self.omega_c = restrict, omega_c
        nco = coarsest if coarsest else 16 * np.sqrt(A.shape[0])
        self.levels = []
        A = sp.csr_matrix(A)
        while True:
            L = Level()
            L.A, L.D, lev = A, A.diagonal(), len(self.levels)
            L.k = k0 if lev == 0 else k1
            self.levels.append(L)
            if A.shape[0] < nco or lev + 1 == max_levels:
                L.lu = sla.splu(sp.csc_matrix(A))
                break
            isF = split0 if (split0 is not None and lev == 0) else greedy_f_points(A, delta)
            C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
            Acc, Acf, Afc, Aff = A[C][:, C], A[C][:, F], A[F][:, C], A[F][:, F]
            om = omega if (omega_coarse is None or lev == 0) else omega_coarse
            Dfi = sp.diags(om / Aff.diagonal())
            K = (sp.eye(len(F)) - Dfi @ Aff).tocsr()                 # K = Q^T T Q
            H = (-Dfi @ Afc).tocsr()                                 # H = Q^T T P
            S, V = H.copy(), H.copy()
            for _ in range(L.k - 1):
                V = (K @ V).tocsr(); S = (S + V).tocsr()
            L.Sex = S                                                # exact sum_{l<k} K^l H
            if mode == "lin" and lev >= lin_levels:
                if tauP_deep: S = drop_rel(S, tauP_deep)
            elif mode == "lin" and lev >= trunc_from:
                Bc = np.column_stack([np.ones(len(C)), xy[C]]) if tv is None else np.column_stack([np.ones(len(C)), tv[C]])
                if theta is None:
                    S = truncate_lin(S, Bc, nkeep)
                else:                                                # blend of 1-entry and 2-entry row-sum-exact S~
                    B1 = np.ones((len(C), 1))
                    S = ((1 - theta) * truncate_lin(S, B1, 1) + theta * truncate_lin(S, B1, 2)).tocsr()
            elif tauP > 0 and lev >= trunc_from:
                S = truncate(S, tauP, mode)
            if xy is not None: xy = xy[C]
            tv = None
            L.C, L.F, L.K, L.H, L.S, L.ST = C, F, K, H, S, S.T.tocsr()
            L.w = np.full(A.shape[0], om); L.w[C] = om if (omega_c is None or lev >= wc_levels) else omega_c
            L.om = om
            L.AFC, L.AFF = Afc.tocsr(), Aff.tocsr()
            L.APc, L.APf = (Acc + Acf @ S).tocsr(), (Afc + Aff @ S).tocsr()
            tf = tau_deep if (tau_deep and lev + 1 >= deep_from) else tau_filter
            Acn = (L.APc + L.ST @ L.APf).tocsr() if restrict == "ppt" else L.APc   # R A P~, R = P~^T or [I 0]
            A = filter_lumped(Acn, tf) if tf else Acn

    def RAP(self, L, v):
        return L.APc @ v + L.ST @ (L.APf @ v) if self.restrict == "ppt" else L.APc @ v

    def Rr(self, L, r):
        return r[L.C] + L.ST @ r[L.F] if self.restrict == "ppt" else r[L.C].copy()

    def coarse(self, lev, rc, which):
        e = self.cycle(lev + 1, rc, np.zeros(len(rc)), which)
        if lev <= self.refine and not hasattr(self.levels[lev + 1], "lu"):
            e = e + self.cycle(lev + 1, rc - self.RAP(self.levels[lev], e), np.zeros(len(rc)), which)
        return e

    def cycle(self, lev, b, x0, which):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            self.work += (L.lu.L.nnz + L.lu.U.nnz) / self.nnz0
            return L.lu.solve(b)
        A, C, F, w = L.A, L.C, L.F, getattr(L, "om", self.omega)
        ap = L.APc.nnz + L.APf.nnz
        self.work += {"twin": L.AFC.nnz + L.k * L.AFF.nnz + A.nnz + 2 * L.S.nnz,
                      "markovJ": (L.k + 1) * A.nnz + 2 * L.S.nnz,
                      "markovJF": (L.k + 2) * A.nnz + 2 * L.S.nnz,
                      "hat": (2 * L.k + 2) * A.nnz + L.k * (L.K.nnz + L.H.nnz) + 3 * L.S.nnz,
                      "N": (L.k + 1) * A.nnz + L.k * (L.K.nnz + L.H.nnz) + 3 * L.S.nnz + ap}[which] / self.nnz0
        if lev <= self.refine and lev + 1 < len(self.levels) and not hasattr(self.levels[lev + 1], "lu"):
            self.work += ap / self.nnz0
        if which == "twin":                                          # k F-Jacobi sweeps, then P~ correction
            x = x0.copy()
            g = b[F] - L.AFC @ x[C]
            for _ in range(L.k):
                x[F] += w * (g - L.AFF @ x[F]) / L.D[F]
            r = b - A @ x
            e = self.coarse(lev, self.Rr(L, r), which)
            x[C] += e; x[F] += L.S @ e
            return x
        if which == "hat":                                           # effective Markovian: x += M^ r, then P~ correction
            r = b - A @ x0
            uT = np.zeros(len(x0))                                   # sum_j T_Q^j Q Q^# M r  (k F-sweeps from 0)
            for _ in range(L.k):
                uT[F] += w * (r[F] - L.AFF @ uT[F]) / L.D[F]
            u = np.zeros(len(x0))                                    # sum_j T^j M r          (k global sweeps from 0)
            for _ in range(L.k):
                u = u + L.w * (r - A @ u) / L.D
            c = u[C]; s_ = np.zeros(len(F))
            for _ in range(L.k):                                     # S c = sum_{l<k} K^l H c
                s_ = L.K @ s_ + L.H @ c
            x = x0 + uT
            x[F] += s_ - L.S @ c                                     # + Q (S - S~) P^# sum_j T^j M r
            r2 = b - A @ x
            e = self.coarse(lev, self.Rr(L, r2), which)
            x[C] += e; x[F] += L.S @ e
            return x
        if which == "markovJ":                                       # k global Jacobi, then P~ correction
            x = x0.copy()
            for _ in range(L.k):
                x = x + L.w * (b - A @ x) / L.D
            r = b - A @ x
            e = self.coarse(lev, self.Rr(L, r), which)
            x[C] += e; x[F] += L.S @ e
            return x
        xs = [x0]                                                    # N (Alg 7.3)
        for _ in range(L.k):
            xs.append(xs[-1] + L.w * (b - A @ xs[-1]) / L.D)
        xk = xs[-1]
        y = np.zeros(len(F))
        for l in range(L.k):
            y = L.K @ y + L.H @ (xk[C] - xs[l][C])
        sh = self.shift[lev] if isinstance(self.shift, (list, tuple)) else self.shift
        j = int(round(sh * L.k))                             # coarse unknown = step-j error
        d0 = xk[C] - xs[j][C]                                        # d^j = c_j - c_k
        xh = xk.copy(); xh[F] += y
        r = b - A @ xh
        e = self.coarse(lev, self.Rr(L, r) + self.RAP(L, d0), which)
        z = e - d0
        xh[C] += z; xh[F] += L.S @ z
        return xh

    def solve(self, b, which, tol=1e-8, maxcyc=40):
        A = self.levels[0].A
        x, nb, h = np.zeros(A.shape[0]), np.linalg.norm(b), [1.0]
        self.work = 0.0
        for _ in range(maxcyc):
            x = self.cycle(0, b, x, which)
            h.append(np.linalg.norm(b - A @ x) / nb)
            if h[-1] < tol or not np.isfinite(h[-1]) or h[-1] > 1e8:
                break
        self.wpc = self.work / (len(h) - 1)
        return np.array(h)

    def opcx(self):
        return sum(L.A.nnz for L in self.levels) / self.nnz0, \
            sum(L.S.nnz for L in self.levels if hasattr(L, "S")) / self.nnz0


def factor(h):
    n = len(h) - 1
    if not np.isfinite(h[-1]) or h[-1] > 1:
        return float("inf")
    return (h[-1] / h[1]) ** (1 / max(n - 1, 1)) if n > 1 else h[-1]
