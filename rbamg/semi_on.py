"""Parametrized Algorithm 7.2 (semi-Markovian, injection + memory in the interpolation) for O(n) tuning.
Rhat = [I, -A_cf (I + G + ... + G^{L-1}) D_ff^{-1}],  coarse operator Rhat A P0 filtered (row-max, tau / tau_deep),
k0 sweeps on level 0, k1 on coarser levels, Jacobi w, greedy split delta, LU below 'coarsest' unknowns.
markov=True: the Markovian version (injection, same operators)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from ml_cycles import greedy_f_points, cf_blocks, filter_rowmax
from ml_z73 import cap_rows
from robust_n import drop_rel

class S72:
    def __init__(self, A, k0=4, k1=4, omega=0.7, L=3, tau=1e-3, tau_deep=None, deep_from=4, delta=0.3,
                 coarsest=150, markov=False, nnz0=None, lev=0, maxlev=40, delta_coarse=None, refine=-1, omega_coarse=None, tauR=0.0, cap=0, split_theta=0.0, nref=1):
        self.A = sp.csr_matrix(A); self.n = A.shape[0]; self.nnz0 = nnz0 or self.A.nnz; self.lev = lev
        self.k = k0 if lev == 0 else k1; self.omega = omega if (lev == 0 or omega_coarse is None) else omega_coarse; self.markov = markov; self.refine = refine; self.nref = nref
        if self.n < coarsest or lev + 1 >= maxlev:
            self.lu = sla.splu(sp.csc_matrix(self.A)); self.lunnz = self.lu.L.nnz + self.lu.U.nnz; return
        As = filter_rowmax(self.A, split_theta) if (split_theta > 0 and lev > 0) else self.A
        isF = greedy_f_points(As, delta if (lev == 0 or delta_coarse is None) else delta_coarse); omega = self.omega
        self.C, self.F, Acc, Acf, Afc, Aff = cf_blocks(self.A, isF)
        dff = Aff.diagonal(); self.D = self.A.diagonal(); Dfi = sp.diags(1 / dff)
        self.K = (sp.eye(len(self.F)) - omega * Dfi @ Aff).tocsr(); self.H = (-omega * Dfi @ Afc).tocsr()
        G = (sp.eye(len(self.F)) - Dfi @ Aff).tocsr(); Z = sp.eye(len(self.F), format="csr"); Gp = Z.copy()
        for _ in range(L - 1): Gp = (Gp @ G).tocsr(); Z = (Z + Gp).tocsr()
        self.Rf = (-Acf @ (Z @ Dfi)).tocsr()
        if tauR > 0: self.Rf = drop_rel(self.Rf, tauR)
        self.Asig = (Acc + self.Rf @ Afc).tocsr()
        t = tau_deep if (tau_deep is not None and lev + 1 >= deep_from) else tau
        Acn = filter_rowmax(self.Asig, t)
        if cap > 0: Acn = cap_rows(Acn, cap, lump=True)
        self.coarse = S72(Acn, k0, k1, omega, L, tau, tau_deep, deep_from, delta, coarsest,
                          markov, self.nnz0, lev + 1, maxlev, delta_coarse, refine, omega_coarse, tauR, cap, split_theta, nref)

    def cycle(self, b, x0, work):
        if hasattr(self, "lu"):
            work[0] += self.lunnz / self.nnz0; return self.lu.solve(b)
        k, w, A, C, F = self.k, self.omega, self.A, self.C, self.F
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + w * (b - A @ xs[-1]) / self.D)
        work[0] += k * A.nnz / self.nnz0
        xk = xs[-1]
        if self.markov:
            r = b - A @ xk; rc = r[C] + self.Rf @ r[F]; e = self.coarse.cycle(rc, np.zeros(len(C)), work)
            if self.lev <= self.refine and not hasattr(self.coarse, "lu"):
                e = e + self.coarse.cycle(rc - self.Asig @ e, np.zeros(len(C)), work)
            out = xk.copy(); out[C] += e; work[0] += (A.nnz + self.Rf.nnz) / self.nnz0; return out
        d = [xk[C] - xs[l][C] for l in range(k)]
        y = np.zeros(len(F))
        for l in range(k):
            y = self.K @ y + self.H @ d[l]
        r = b - A @ xk
        rc = r[C] + self.Rf @ r[F] + self.Asig @ d[0]
        e0 = self.coarse.cycle(rc, np.zeros(len(C)), work)
        if self.lev <= self.refine and not hasattr(self.coarse, "lu"):
            for _ in range(self.nref):
                e0 = e0 + self.coarse.cycle(rc - self.Asig @ e0, np.zeros(len(C)), work); work[0] += self.Asig.nnz / self.nnz0
        z = e0 - d[0]; out = xk.copy(); out[C] += z
        v = self.H @ z; s = v.copy()
        for _ in range(k - 1):
            v = self.K @ v; s += v
        out[F] += y + s
        work[0] += (k * (self.K.nnz + self.H.nnz) + A.nnz + self.Rf.nnz + self.Asig.nnz + self.H.nnz + (k - 1) * self.K.nnz) / self.nnz0
        return out

    def levels(self):
        return [self] + ([] if hasattr(self, "lu") else self.coarse.levels())
    def op_complexity(self):
        return sum(l.A.nnz for l in self.levels()) / self.nnz0
