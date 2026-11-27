"""Level 0 with a non-injective complement Q = Q0 + P0 Z (Q^# = Q0^T, P^# = P0^T - Z Q0^T), coarse levels from ML.

C/F coordinates, P0 = [I; 0], Q0 = [0; I]:
    P^# x = x_C - Z x_F,   Q y = [Z y; y],   K = K0 + H Z,   K0 = I - w D_ff^-1 A_ff,   H = -w D_ff^-1 A_fc
    S = sum_{l<k} K^l H,   S~ = truncation of S,   P~ = P0 + Q S~ = [I + Z S~; S~],   R = P~^T,   A_c = R A P~
Cycles on level 0 (coarse solve = recursive cycle of ML on A_c, same method name):
    N    : k Jacobi sweeps (w_F on F, w_C on C), history y <- K y + H P^#(x^k - x^l), h = Q y,
           coarse rhs R r(x^k + h) + A_c d0, d0 = P^#(x^k - x^0), correction P~ (e - d0)
    twin : k steps x += Q (M r)_F  (= I - Q Q^# M A), then P~ correction
    markovJ : k Jacobi sweeps, then P~ correction
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from ml_twins import ML, truncate_lin, truncate
from robust_n import drop_rel


class MLZ:
    def __init__(self, A, isF, Z_fn, omega=1.0, omega_c=1.0, k=2, trunc=("lin", 1), xy=None, coarse_kw=None,
                 refine0=0):
        A = sp.csr_matrix(A)
        self.A, self.nnz0, self.k, self.refine0 = A, A.nnz, k, refine0
        C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
        self.C, self.F = C, F
        d = A.diagonal(); self.D = d
        self.w = np.full(A.shape[0], omega); self.w[C] = omega_c; self.omega = omega
        Acc, Acf, Afc, Aff = A[C][:, C], A[C][:, F], A[F][:, C], A[F][:, F]
        Dfi = sp.diags(omega / Aff.diagonal())
        K0 = (sp.eye(len(F)) - Dfi @ Aff).tocsr(); H = (-Dfi @ Afc).tocsr()
        Z = Z_fn(A, C, F, K0, H)
        self.Z = sp.csr_matrix(Z) if Z is not None else sp.csr_matrix((len(C), len(F)))
        K = (K0 + H @ self.Z).tocsr()
        S, V = H.copy(), H.copy()
        for _ in range(k - 1):
            V = (K @ V).tocsr(); S = (S + V).tocsr()
        self.K0, self.K, self.H, self.Sex = K0, K, H, S
        mode, par = trunc
        if mode == "lin":
            Bc = np.ones((len(C), 1)) if xy is None or par == 1 else np.column_stack([np.ones(len(C)), xy[C]])
            St = truncate_lin(S, Bc, par)
        elif mode == "blend":
            B1 = np.ones((len(C), 1))
            St = ((1 - par) * truncate_lin(S, B1, 1) + par * truncate_lin(S, B1, 2)).tocsr()
        elif mode == "drop":
            St = drop_rel(S, par)
        else:
            St = S
        self.S = sp.csr_matrix(St); self.ST = self.S.T.tocsr()
        self.IZS = (sp.eye(len(C)) + self.Z @ self.S).tocsr(); self.IZST = self.IZS.T.tocsr()
        self.AFC, self.AFF = Afc.tocsr(), Aff.tocsr()
        self.APc = (Acc @ self.IZS + Acf @ self.S).tocsr()
        self.APf = (Afc @ self.IZS + Aff @ self.S).tocsr()
        Ac = (self.IZST @ self.APc + self.ST @ self.APf).tocsr()
        self.sub = ML(Ac, **(coarse_kw or {}))
        self.levels = [None] + self.sub.levels

    # transfer helpers
    def Psh(self, x):
        return x[self.C] - self.Z @ x[self.F]

    def Rr(self, r):
        return self.IZST @ r[self.C] + self.ST @ r[self.F]

    def RAP(self, v):
        return self.IZST @ (self.APc @ v) + self.ST @ (self.APf @ v)

    def prolong_add(self, x, e):
        x[self.C] += self.IZS @ e; x[self.F] += self.S @ e

    def coarse(self, rc, which):
        sw = "N" if which == "N" else which
        e = self.sub.cycle(0, rc, np.zeros(len(rc)), sw)
        for _ in range(self.refine0):
            e = e + self.sub.cycle(0, rc - self.RAP(e), np.zeros(len(rc)), sw)
        return e

    def cycle(self, b, x0, which):
        A, C, F, k = self.A, self.C, self.F, self.k
        if which == "twin":
            x = x0.copy()
            for _ in range(k):
                g = self.omega * (b[F] - A[F] @ x) / self.D[F]
                x[F] += g; x[C] += self.Z @ g
            r = b - A @ x
            self.prolong_add(x, self.coarse(self.Rr(r), which))
            self.work += (k * (A.nnz * len(F) / A.shape[0] + 2 * self.Z.nnz) + A.nnz + 4 * self.S.nnz) / self.nnz0
            return x
        if which == "markovJ":
            x = x0.copy()
            for _ in range(k):
                x = x + self.w * (b - A @ x) / self.D
            r = b - A @ x
            self.prolong_add(x, self.coarse(self.Rr(r), which))
            self.work += ((k + 1) * A.nnz + 4 * self.S.nnz) / self.nnz0
            return x
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + self.w * (b - A @ xs[-1]) / self.D)
        xk = xs[-1]; dk = self.Psh(xk)
        y = np.zeros(len(F))
        for l in range(k):
            y = self.K @ y + self.H @ (dk - self.Psh(xs[l]))
        d0 = dk - self.Psh(x0)
        xh = xk.copy(); xh[C] += self.Z @ y; xh[F] += y
        r = b - A @ xh
        e = self.coarse(self.Rr(r) + self.RAP(d0), which)
        self.prolong_add(xh, e - d0)
        self.work += ((k + 1) * A.nnz + k * (self.K.nnz + self.H.nnz + 2 * self.Z.nnz) + 6 * self.S.nnz
                      + 2 * (self.APc.nnz + self.APf.nnz)) / self.nnz0
        return xh

    def solve(self, b, which, tol=1e-8, maxcyc=60):
        x, nb, h = np.zeros(self.A.shape[0]), np.linalg.norm(b), [1.0]
        self.work = 0.0; self.sub.work = 0.0
        for _ in range(maxcyc):
            x = self.cycle(b, x, which)
            h.append(np.linalg.norm(b - self.A @ x) / nb)
            if h[-1] < tol or not np.isfinite(h[-1]) or h[-1] > 1e8:
                break
        self.wpc = (self.work + self.sub.work * self.sub.nnz0 / self.nnz0) / (len(h) - 1)
        return np.array(h)

    def rho(self, which, maxiter=1500):
        n = self.A.shape[0]; z = np.zeros(n)
        self.work = 0.0; self.sub.work = 0.0
        op = sla.LinearOperator((n, n), matvec=lambda v: self.cycle(z, v, which))
        try:
            return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=maxiter, tol=1e-3, ncv=40)[0])
        except Exception:
            return float("nan")
