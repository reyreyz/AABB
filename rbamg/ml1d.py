"""Multilevel non-Markovian cycle with stencil decompositions (W, Z), 1D, coarsening by 3.

Per level (sparse, C/F order inside the level):
  split   C = indices 2, 5, 8, ...   (n_l = 3 n_{l+1} - 1)
  T       Jacobi, omega = 1
  P = [I; W], Q = Q0 + P Z, Qsh = [-W I], Psh = [I 0] - Z Qsh
  K  = T_ff - W T_cf + H_W Z,   H_W = T_fc + T_ff W - W T_cc - W T_cf W
  S  = sum_{l<k} K^l H_W,  P' = [I + Z S; W + S + W Z S],  A_sigma = (A P')_C (injection)
  W, Z from a periodic stencil (w neighbours per side, period p blocks) fitted by
  min ||K^k||_F on a central window of the level operator (O(1) work per level),
  then laid out over the whole level (O(N)).
Cycle (Alg. 4.2, injection):  x^0..x^k Jacobi;  d^l = Psh (x^k - x^l);
  y <- K y + H_W d^l;  h = Q y;  A_sigma eps = (b - A x^k - A h)_C (recursive);
  return x^k + P' eps + h.   Coarsest level (< 60 unknowns) factored.
Optional bandwidth cap of A_sigma (lumped) keeps the stencils bounded.
"""

import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import jax
import jax.numpy as jnp
from scipy.optimize import minimize
from joint_wz import Joint
from stencil_wz import layout_general
from cd1d import cd1d

jax.config.update("jax_enable_x64", True)
COARSEST = 60


def fit_stencil(A_win, k, w, p, starts=12, seed=0, maxiter=600, lam=0.0):
    n = A_win.shape[0]
    C, F, e, nslot = layout_general(n, w, p)
    perm = np.concatenate([C, F])
    Ad = A_win[np.ix_(perm, perm)]
    pat = [e[e[:, 0] == j, 1] for j in range(len(F))]
    J = Joint(Ad, len(C), k, pat, pat)
    slots = np.concatenate([e[:, 2], nslot + e[:, 2]])
    Phi = np.zeros((len(slots), 2 * nslot))
    Phi[np.arange(len(slots)), slots] = 1.0
    Phi = jnp.array(Phi)
    f = lambda th: jnp.sum(jnp.linalg.matrix_power(J.KH(Phi @ th)[2], k) ** 2) + lam * len(F) * jnp.sum(th ** 2)
    vg = jax.jit(jax.value_and_grad(f))

    def fun(th):
        v, g = vg(jnp.array(th))
        return float(v), np.array(g)
    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        r = minimize(fun, rng.standard_normal(2 * nslot) * (0.5 if s else 0.0), jac=True,
                     method="L-BFGS-B", options=dict(maxiter=maxiter))
        if best is None or r.fun < best.fun:
            best = r
    Kk = jnp.linalg.matrix_power(J.KH(Phi @ jnp.array(best.x))[2], k)
    return best.x, float(np.sqrt(np.sum(np.array(Kk) ** 2) / len(F)))


def band_cap(A, bw):
    """keep |i - j| <= bw, lump the rest into the diagonal."""
    A = sp.coo_matrix(A)
    far = np.abs(A.row - A.col) > bw
    lump = np.bincount(A.row[far], weights=A.data[far], minlength=A.shape[0])
    B = sp.csr_matrix((A.data[~far], (A.row[~far], A.col[~far])), shape=A.shape)
    return (B + sp.diags(lump)).tocsr()


class Level:
    pass


class ML1D:
    def __init__(self, A, k=4, w=1, p=2, mode="WZ", window=95, bw=None, theta_cache=None,
                 max_levels=99, k_coarse=None, w_coarse=None, starts=12, p_coarse=None, lam=0.0, gamma=1):
        self.k, self.mode, self.gamma = k, mode, gamma
        self.nnz0 = A.nnz
        self.levels = []
        A = sp.csr_matrix(A)
        t0 = time.time()
        while True:
            L = Level()
            L.A = A
            self.levels.append(L)
            n = A.shape[0]
            lev = len(self.levels) - 1
            L.k = k if (lev == 0 or k_coarse is None) else k_coarse
            wl = w if (lev == 0 or w_coarse is None) else w_coarse
            pl = p if (lev == 0 or p_coarse is None) else p_coarse
            if n < COARSEST or len(self.levels) >= max_levels:
                L.lu = sla.splu(sp.csc_matrix(A))
                break
            C, F, e, nslot = layout_general(n, wl, pl)
            nc, nf = len(C), len(F)
            if mode in ("WZ", "W", "Z"):
                mid = (n // 2) - (n // 2) % 3 + 0       # window aligned to the 3-periodic split
                lo = max(0, mid - window // 2)
                lo -= lo % (3 * pl)
                hi = lo + window
                if hi > n:
                    lo, hi = 0, n
                th, res = fit_stencil(A[lo:hi][:, lo:hi].toarray(), L.k, wl, pl, starts=starts, lam=lam)
                L.fit_res = res
                vals = th[np.concatenate([e[:, 2], nslot + e[:, 2]])]
                nW = len(e)
                W = sp.csr_matrix((vals[:nW], (e[:, 0], e[:, 1])), shape=(nf, nc))
                Z = sp.csr_matrix((vals[nW:], (e[:, 1], e[:, 0])), shape=(nc, nf))
                if mode == "W":
                    Z = sp.csr_matrix((nc, nf))
                if mode == "Z":
                    W = sp.csr_matrix((nf, nc))
            else:                                           # Q0, P0 (mode 'Q0') or Markovian ('M')
                W = sp.csr_matrix((nf, nc))
                Z = sp.csr_matrix((nc, nf))
            D = A.diagonal()
            T = (sp.eye(n) - sp.diags(1 / D) @ A).tocsr()
            Tcc, Tcf, Tfc, Tff = T[C][:, C], T[C][:, F], T[F][:, C], T[F][:, F]
            HW = (Tfc + Tff @ W - W @ Tcc - W @ (Tcf @ W)).tocsr()
            K = (Tff - W @ Tcf + HW @ Z).tocsr()
            S, V = HW.copy(), HW.copy()
            for _ in range(L.k - 1):
                V = (K @ V).tocsr()
                S = (S + V).tocsr()
            ZS = (Z @ S).tocsr()
            PpC = (sp.eye(nc) + ZS).tocsr()
            PpF = (W + S + W @ ZS).tocsr()
            if mode == "M":
                PpC, PpF = sp.eye(nc, format="csr"), W
            Acc, Acf = A[C][:, C].tocsr(), A[C][:, F].tocsr()
            Asig = (Acc @ PpC + Acf @ PpF).tocsr()
            L.C, L.F, L.W, L.Z, L.K, L.HW, L.PpC, L.PpF, L.D = C, F, W, Z, K, HW, PpC, PpF, D
            L.AC = A[C].tocsr()
            A = band_cap(Asig, bw) if bw else Asig
            A.eliminate_zeros()
        self.setup_time = time.time() - t0

    def op_complexity(self):
        return sum(L.A.nnz for L in self.levels) / self.nnz0

    def cycle(self, lev, b, x0, work):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            return L.lu.solve(b)
        A, k, nz, C, F = L.A, L.k, self.nnz0, L.C, L.F
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + (b - A @ xs[-1]) / L.D)
            work[0] += A.nnz / nz
        xk = xs[-1]
        if self.mode == "M":
            r = b - A @ xk
            e = self.cycle(lev + 1, r[C], np.zeros(len(C)), work)
            out = xk.copy()
            out[C] += e
            return out
        y = np.zeros(len(F))
        for l in range(k):
            dl = xk - xs[l]
            qd = dl[F] - L.W @ dl[C]                        # Qsh dl
            d = dl[C] - L.Z @ qd                            # Psh dl
            y = L.K @ y + L.HW @ d
            work[0] += (2 * L.W.nnz + L.Z.nnz + L.K.nnz + L.HW.nnz) / nz
        zy = L.Z @ y
        h = np.zeros(A.shape[0])
        h[C] = zy
        h[F] = y + L.W @ zy
        r = b - A @ (xk + h)
        work[0] += (A.nnz + L.Z.nnz + L.W.nnz) / nz
        e = self.cycle(lev + 1, r[C], np.zeros(len(C)), work)
        for _ in range(self.gamma - 1):
            if hasattr(self.levels[lev + 1], "lu"):
                break
            e = self.cycle(lev + 1, r[C], e, work)
        out = xk + h
        out[C] += L.PpC @ e
        out[F] += L.PpF @ e
        work[0] += (L.PpC.nnz + L.PpF.nnz) / nz
        return out

    def solve(self, b, tol=1e-10, maxcyc=60):
        A = self.levels[0].A
        x = np.zeros(A.shape[0])
        nb = np.linalg.norm(b)
        hist, works, work = [1.0], [0.0], [0.0]
        for _ in range(maxcyc):
            x = self.cycle(0, b, x, work)
            res = np.linalg.norm(b - A @ x) / nb
            hist.append(res)
            works.append(work[0])
            if res < tol or not np.isfinite(res) or res > 1e12:
                break
        return np.array(hist), np.array(works)


def sizes(L, coarsest=26):
    n = [coarsest]
    for _ in range(L):
        n.append(3 * n[-1] - 1)
    return n[::-1]
