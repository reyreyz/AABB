"""Multilevel Markovian Petrov-Galerkin V-cycle with a good interpolation, same smoother and restriction as Alg. 7.2:
greedy split (delta 0.3), Rhat = [I, -A_cf (I+G+G^2) D_ff^{-1}], P = [I; W] with
  W1 : W = -D_ff^{-1} A_fc            (one-step Jacobi / direct interpolation)
  Pp : W = sum_{l<k} K^l H, K = I - w D_ff^{-1} A_ff, H = -w D_ff^{-1} A_fc, drop 1e-2   (S's interpolation, static)
coarse operator Rhat A P (filtered 1e-3, row-max), k Jacobi sweeps, coarse correction x += P A_c^{-1} Rhat r."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
import ml_cycles as mc
from robust_n import greedy_f_points, drop_rel
COARSEST = 150
class MLM:
    def __init__(self, A, k=2, omega=0.8, wmode="W1", tau=1e-3, nnz0=None):
        self.A = sp.csr_matrix(A); self.n = A.shape[0]; self.nnz0 = nnz0 or self.A.nnz; self.k, self.w = k, omega
        if self.n < COARSEST:
            self.lu = sla.splu(sp.csc_matrix(self.A)); self.lunnz = self.lu.L.nnz + self.lu.U.nnz; return
        isF = greedy_f_points(self.A, 0.3); C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]; self.C, self.F = C, F
        A = self.A; Acc, Acf, Afc, Aff = A[C][:, C], A[C][:, F], A[F][:, C], A[F][:, F]; self.D = A.diagonal()
        Dfi = sp.diags(1 / Aff.diagonal()); G = (sp.eye(len(F)) - Dfi @ Aff).tocsr()
        self.Rf = (-Acf @ ((sp.eye(len(F)) + G + G @ G) @ Dfi)).tocsr()
        if wmode == "W1":
            W = (-Dfi @ Afc).tocsr()
        else:
            K = (sp.eye(len(F)) - omega * Dfi @ Aff).tocsr(); H = (-omega * Dfi @ Afc).tocsr(); W, V = H.copy(), H.copy()
            for _ in range(k - 1): V = (K @ V).tocsr(); W = (W + V).tocsr()
            W = drop_rel(W, 1e-2)
        self.W = W
        APc, APf = (Acc + Acf @ W).tocsr(), (Afc + Aff @ W).tocsr()
        Ac = (APc + self.Rf @ APf).tocsr()                                # Rhat A P
        self.coarse = MLM(mc.filter_rowmax(Ac, tau), k, omega, wmode, tau, self.nnz0)
    def cycle(self, b, x0, work):
        if self.n < COARSEST:
            work[0] += self.lunnz / self.nnz0; return self.lu.solve(b)
        x = x0.copy()
        for _ in range(self.k):
            x = x + self.w * (b - self.A @ x) / self.D
        r = b - self.A @ x
        e = self.coarse.cycle(r[self.C] + self.Rf @ r[self.F], np.zeros(len(self.C)), work)
        x[self.C] += e; x[self.F] += self.W @ e
        work[0] += ((self.k + 1) * self.A.nnz + self.Rf.nnz + self.W.nnz) / self.nnz0
        return x
    def levels(self):
        return [self] + (self.coarse.levels() if self.n >= COARSEST else [])
    def op_complexity(self):
        return sum(l.A.nnz for l in self.levels()) / self.nnz0
