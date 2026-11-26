"""O(N) two-level non-Markovian solver with a stencil decomposition (W, Z).

Given the sparse fine operator A, a C/F layout with stencil slots, and the fitted
stencil theta (O(1) numbers), the setup lays out W, Z (bounded nnz per row/col),
forms K, H_W, P' = P + Q sum_{l<k} K^l H_W and A_sigma = (A P')_C with sparse
products of bounded fill, and builds an O(N) coarse solver (PyAMG AIR on
A_sigma, a fixed number of cycles, or a sparse LU for reference).
Cycle (Alg. 4.2, injection):  x^0..x^k Jacobi (omega = 1);  d^l = Psh (x^k - x^l);
y <- K y + H_W d^l;  h = Q y;  A_sigma eps = (b - A x^k - A h)_C;  x^k + P' eps + h.
"""

import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


class TwoLevelWZ:
    def __init__(self, A, C, F, e, nslot, theta, k, mode="WZ", coarse="lu", coarse_cycles=2,
                 coarse_tol=None):
        t0 = time.time()
        A = sp.csr_matrix(A)
        self.A, self.k, self.mode = A, k, mode
        nc, nf = len(C), len(F)
        if mode == "WZ":
            vals = theta[np.concatenate([e[:, 2], nslot + e[:, 2]])]
            nW = len(e)
            W = sp.csr_matrix((vals[:nW], (e[:, 0], e[:, 1])), shape=(nf, nc))
            Z = sp.csr_matrix((vals[nW:], (e[:, 1], e[:, 0])), shape=(nc, nf))
        else:
            W, Z = sp.csr_matrix((nf, nc)), sp.csr_matrix((nc, nf))
        n = A.shape[0]
        D = A.diagonal()
        T = (sp.eye(n) - sp.diags(1 / D) @ A).tocsr()
        Tcc, Tcf, Tfc, Tff = T[C][:, C], T[C][:, F], T[F][:, C], T[F][:, F]
        HW = (Tfc + Tff @ W - W @ Tcc - W @ (Tcf @ W)).tocsr()
        K = (Tff - W @ Tcf + HW @ Z).tocsr()
        S, V = HW.copy(), HW.copy()
        for _ in range(k - 1):
            V = (K @ V).tocsr()
            S = (S + V).tocsr()
        ZS = (Z @ S).tocsr()
        self.PpC = (sp.eye(nc) + ZS).tocsr()
        self.PpF = (W + S + W @ ZS).tocsr()
        if mode == "M":
            self.PpC, self.PpF = sp.eye(nc, format="csr"), W
        Acc, Acf = A[C][:, C].tocsr(), A[C][:, F].tocsr()
        self.Asig = (Acc @ self.PpC + Acf @ self.PpF).tocsr()
        self.C, self.F, self.W, self.Z, self.K, self.HW, self.D = C, F, W, Z, K, HW, D
        self.coarse, self.coarse_cycles, self.coarse_tol = coarse, coarse_cycles, coarse_tol
        if coarse == "lu":
            self.lu = sla.splu(sp.csc_matrix(self.Asig))
        else:
            import pyamg
            self.amg = pyamg.air_solver(self.Asig, max_coarse=150)
        self.setup_time = time.time() - t0

    def coarse_solve(self, r, work):
        if self.coarse == "lu":
            return self.lu.solve(r)
        res = []
        if self.coarse_tol:
            x = self.amg.solve(r, tol=self.coarse_tol, maxiter=50, accel="gmres", residuals=res)
            work[0] += (len(res) - 1) * self.amg.cycle_complexity() * self.Asig.nnz / self.A.nnz
        else:
            x = self.amg.solve(r, maxiter=self.coarse_cycles, tol=1e-30, residuals=res)
            work[0] += self.coarse_cycles * self.amg.cycle_complexity() * self.Asig.nnz / self.A.nnz
        return x

    def cycle(self, b, x0, work):
        A, k, C, F = self.A, self.k, self.C, self.F
        nz = A.nnz
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + (b - A @ xs[-1]) / self.D)
            work[0] += 1.0
        xk = xs[-1]
        if self.mode == "M":
            r = b - A @ xk
            e = self.coarse_solve(r[C], work)
            out = xk.copy()
            out[C] += e
            return out
        y = np.zeros(len(F))
        for l in range(k):
            dl = xk - xs[l]
            qd = dl[F] - self.W @ dl[C]
            d = dl[C] - self.Z @ qd
            y = self.K @ y + self.HW @ d
            work[0] += (2 * self.W.nnz + self.Z.nnz + self.K.nnz + self.HW.nnz) / nz
        zy = self.Z @ y
        h = np.zeros(A.shape[0])
        h[C] = zy
        h[F] = y + self.W @ zy
        r = b - A @ (xk + h)
        work[0] += 1.0 + (self.Z.nnz + self.W.nnz) / nz
        e = self.coarse_solve(r[C], work)
        out = xk + h
        out[C] += self.PpC @ e
        out[F] += self.PpF @ e
        work[0] += (self.PpC.nnz + self.PpF.nnz) / nz
        return out

    def solve(self, b, tol=1e-10, maxcyc=40):
        x = np.zeros(self.A.shape[0])
        nb = np.linalg.norm(b)
        hist, works, work = [1.0], [0.0], [0.0]
        for _ in range(maxcyc):
            x = self.cycle(b, x, work)
            res = np.linalg.norm(b - self.A @ x) / nb
            hist.append(res)
            works.append(work[0])
            if res < tol or not np.isfinite(res) or res > 1e12:
                break
        return np.array(hist), np.array(works)

    def complexity(self):
        """stored nonzeros (W, Z, K, H_W, P', A_sigma) relative to nnz(A)"""
        tot = sum(X.nnz for X in (self.W, self.Z, self.K, self.HW, self.PpC, self.PpF, self.Asig))
        return tot / self.A.nnz
