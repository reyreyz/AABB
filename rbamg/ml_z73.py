"""O(N) multilevel non-Markovian cycle with a sparse, bounded complement Q = Q0 + P0 Z.

Base: Algorithm 7.3 (greedy F-points with ||I - D_ff^{-1}A_ff||_inf <= delta, Jacobi
omega = 1, P0 = [I; 0], injection).  With the complement Q = Q0 + P0 Z
(Qsh = Q0sh = [0 I], Psh = [I -Z]):
    K = I - D_ff^{-1} A_ff,  H = -D_ff^{-1} A_fc,  K_Z = K + H Z
    S = sum_{l<k} K_Z^l H,   P' = P0 + Q S = [I + Z S; S]
    A_sigma = (A P')_C (injection), filtered (row-max, lumped) -> next level
Cycle (Alg. 4.2 with the shift relation of Alg. 7.3):
    x^0..x^k Jacobi;  d^l = Psh (x^k - x^l) = dC - Z dF
    y = 0; y <- K_Z y + H d^l (l < k);   h = Q y = [Z y; y]
    r_c = (b - A x^k - A h)_C + A_sigma_exact d^0      (A_sigma_exact matrix-free)
    e0  = recursive cycle(r_c, 0);  z = e0 - d^0;  return x^k + h + P' z
Z: at most q entries per column, chosen among the coarse points within graph
distance 2 of the fine point by the CR overlap |H[:,c]^T K[:,j]|, values from
min ||K_Z^k||_F over that pattern (sparse gradients, fixed iteration cap).
Everything is sparse with bounded fill, so setup and cycle are O(N) per level.
"""

import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import ml_cycles as mc
import ml_nz as mz

COARSEST = 150


def bounded_pattern(A, C, F, K, H, q):
    """for each fine j: <= q coarse points within distance 2, ranked by |H[:,c]^T K[:,j]|."""
    Af = sp.csr_matrix(A)
    G = (abs(Af) > 0).astype(float).tocsr()
    G2 = (G @ G).tocsr()
    cpos = -np.ones(A.shape[0], dtype=int)
    cpos[C] = np.arange(len(C))
    score = (abs(H).T @ abs(K)).tocsc()            # nc x nf: |H|^T |K| (overlap)
    rows, cols = [], []
    for jf, p in enumerate(F):
        nb = G2.indices[G2.indptr[p]:G2.indptr[p + 1]]
        cand = cpos[nb]
        cand = cand[cand >= 0]
        if len(cand) == 0:
            continue
        sc = np.asarray(score[cand, jf].todense()).ravel()
        keep = cand[np.argsort(-sc)[:q]]
        keep = keep[sc[np.argsort(-sc)[:q]] > 0]
        rows += list(keep)
        cols += [jf] * len(keep)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(C), len(F)))


def cap_rows(X, p, lump=False):
    """keep the p largest |entries| (off-diagonal when lump) per row; lump the rest
    into the diagonal when lump=True (row sums preserved)."""
    X = sp.csr_matrix(X)
    rows, cols, vals = [], [], []
    extra = np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        s0, e0 = X.indptr[i], X.indptr[i + 1]
        c, v = X.indices[s0:e0], X.data[s0:e0]
        if lump:
            dmask = c == i
            dval = v[dmask].sum()
            c, v = c[~dmask], v[~dmask]
        if len(v) > p:
            order = np.argsort(-np.abs(v))
            keep, drop = order[:p], order[p:]
            if lump:
                extra[i] += v[drop].sum()
            c, v = c[keep], v[keep]
        rows += [i] * len(c)
        cols += list(c)
        vals += list(v)
        if lump:
            rows.append(i), cols.append(i), vals.append(dval)
    Y = sp.csr_matrix((vals, (rows, cols)), shape=X.shape)
    if lump:
        Y = Y + sp.diags(extra)
    return Y.tocsr()


class Level:
    pass


class Alg73Z:
    def __init__(self, A, k=4, delta=1.5, q=4, zmode="Kk", tauA=1e-3, markov=False,
                 kk_iter=100, max_levels=99, capP=None, capA=None,
                 delta_coarse=None, z_levels=99, k_coarse=None, coarse="recursive", coarse_cycles=2,
                 theta=1.0, thetas=(0.0, 0.25, 0.5, 0.75, 1.0), pc_mask=None,
                 zpattern="dist2", z_from=0, refine_levels=-1):
        self.k, self.markov = k, markov
        self.nnz0 = A.nnz
        self.levels = []
        A = sp.csr_matrix(A)
        t = time.time()
        self.coarse, self.coarse_cycles = coarse, coarse_cycles
        self.refine_levels = refine_levels      # Alg. 7.3: one refinement step on levels <= this
        while True:
            L = Level()
            L.A = A
            self.levels.append(L)
            n = A.shape[0]
            lev = len(self.levels) - 1
            L.k = k if (lev == 0 or k_coarse is None) else k_coarse
            if n < COARSEST or len(self.levels) >= max_levels:
                L.lu = sla.splu(sp.csc_matrix(A))
                break
            if lev == 1 and coarse == "air":                # O(N) coarse solver: PyAMG AIR on A_sigma
                import pyamg
                L.amg = pyamg.air_solver(A, max_coarse=150)
                break
            dl = delta if (lev == 0 or delta_coarse is None) else delta_coarse
            isF = mc.greedy_f_points(A, dl)
            C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
            Acc, Acf = A[C][:, C].tocsr(), A[C][:, F].tocsr()
            Afc, Aff = A[F][:, C].tocsr(), A[F][:, F].tocsr()
            Dfi = sp.diags(1 / Aff.diagonal())
            K = (sp.eye(len(F)) - Dfi @ Aff).tocsr()
            H = (-Dfi @ Afc).tocsr()
            kl = L.k
            if zmode == "Q0" or markov or lev >= z_levels or lev < z_from:
                Z = sp.csr_matrix((len(C), len(F)))
            else:
                if zpattern == "acf":
                    # Z(i, f) only where A_cf(i, f) != 0: pattern(Z S) within pattern(A_cf S),
                    # i.e. no fill beyond N's coarse operator except through A_cc
                    patt = abs(Acf).tocsr().astype(float)
                    patt.data[:] = 1.0
                else:
                    patt = bounded_pattern(A, C, F, K, H, q)
                Z0 = mz.z_ls(K, H, patt)
                Z = mz.z_kk(K, H, patt, kl, Z0, maxiter=kk_iter) if zmode == "Kk" else Z0
                Z.eliminate_zeros()
                if theta == "auto":
                    # safeguard: theta minimizing the CR rate of the level pair
                    # max(rho(K_{theta Z}), rho(K_c)), K_c = I - D_ff^{-1} A_ff of the
                    # next level after its own greedy split; theta = 0 is plain N
                    best = None
                    for th in thetas:
                        Zt = (th * Z).tocsr()
                        KZt = (K + H @ Zt).tocsr()
                        St, Vt = H.copy(), H.copy()
                        for _ in range(kl - 1):
                            Vt = (KZt @ Vt).tocsr()
                            St = (St + Vt).tocsr()
                        Ac = mc.filter_rowmax((Acc @ (sp.eye(len(C)) + Zt @ St) + Acf @ St).tocsr(), tauA)
                        r_f = mz.cr_rate(KZt)
                        if Ac.shape[0] >= COARSEST:
                            isFc = mc.greedy_f_points(Ac, delta_coarse or delta)
                            Fc = np.nonzero(isFc)[0]
                            Affc = Ac[Fc][:, Fc]
                            Kc = (sp.eye(len(Fc)) - sp.diags(1 / Affc.diagonal()) @ Affc).tocsr()
                            r_c = mz.cr_rate(Kc)
                        else:
                            r_c = 0.0
                        if best is None or max(r_f, r_c) < best[0] - 1e-12:
                            best = (max(r_f, r_c), th)
                    L.theta = best[1]
                    Z = (best[1] * Z).tocsr()
            KZ = (K + H @ Z).tocsr()
            S, V = H.copy(), H.copy()
            for _ in range(kl - 1):
                V = (KZ @ V).tocsr()
                S = (S + V).tocsr()
            if capP:
                S = cap_rows(S, capP)                       # bounded P' rows
            ZS = (Z @ S).tocsr()
            if pc_mask == "diag":                           # no new fill: A_sigma has N's pattern
                ZS = sp.diags(ZS.diagonal()).tocsr()
            elif pc_mask == "acc":                          # pattern of A_cc (one ring of fill)
                Mk = Acc.copy().tocsr()
                Mk.data[:] = 1.0
                ZS = ZS.multiply(Mk).tocsr()
            PpC = (sp.eye(len(C)) + ZS).tocsr()
            PpF = S
            if markov:
                PpC, PpF = sp.eye(len(C), format="csr"), sp.csr_matrix((len(F), len(C)))
            Asig = (Acc @ PpC + Acf @ PpF).tocsr()           # (A P')_C, exact
            L.C, L.F, L.Z, L.KZ, L.H, L.PpC, L.PpF = C, F, Z, KZ, H, PpC, PpF
            L.Acc, L.Acf, L.Asig, L.D = Acc, Acf, Asig, A.diagonal()
            A = mc.filter_rowmax(Asig, tauA)
            if capA:
                A = cap_rows(A, capA, lump=True)            # bounded coarse stencils
            A.eliminate_zeros()
        self.setup_time = time.time() - t

    def op_complexity(self):
        return sum(L.A.nnz for L in self.levels) / self.nnz0

    def stored_complexity(self):
        """all stored sparse objects relative to nnz(A0): A, Z, K_Z, H, P'"""
        tot = 0
        for L in self.levels:
            tot += L.A.nnz
            if hasattr(L, "Z"):
                tot += L.Z.nnz + L.KZ.nnz + L.H.nnz + L.PpC.nnz + L.PpF.nnz + L.Asig.nnz
        return tot / self.nnz0

    def cycle(self, lev, b, x0, work):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            return L.lu.solve(b)
        if hasattr(L, "amg"):
            res = []
            x = L.amg.solve(b, x0=np.zeros_like(b), maxiter=self.coarse_cycles, tol=1e-30, residuals=res)
            work[0] += self.coarse_cycles * L.amg.cycle_complexity() * L.A.nnz / self.nnz0
            return x
        A, k, nz, C, F = L.A, L.k, self.nnz0, L.C, L.F
        xs = [x0]
        for _ in range(k):
            xs.append(xs[-1] + (b - A @ xs[-1]) / L.D)
            work[0] += A.nnz / nz
        xk = xs[-1]
        if self.markov:
            r = b - A @ xk
            work[0] += A.nnz / nz
            e = self.cycle(lev + 1, r[C], np.zeros(len(C)), work)
            out = xk.copy()
            out[C] += e
            return out
        y = np.zeros(len(F))
        d = []
        for l in range(k):
            dl = xk - xs[l]
            dd = dl[C] - L.Z @ dl[F]                       # Psh = [I  -Z]
            d.append(dd)
            y = L.KZ @ y + L.H @ dd
            work[0] += (L.Z.nnz + L.KZ.nnz + L.H.nnz) / nz
        h = np.zeros(A.shape[0])
        h[F] = y
        h[C] = L.Z @ y
        r = b - A @ (xk + h)
        rc = r[C] + L.Acc @ (L.PpC @ d[0]) + L.Acf @ (L.PpF @ d[0])
        work[0] += (L.Z.nnz + A.nnz + L.Acc.nnz + L.Acf.nnz + L.PpC.nnz + L.PpF.nnz) / nz
        e0 = self.cycle(lev + 1, rc, np.zeros(len(C)), work)
        if lev <= self.refine_levels and not hasattr(self.levels[lev + 1], "lu"):
            # refinement: e0 += cycle(r_c - A_sigma e0), A_sigma exact (matrix-free)
            Ae0 = L.Acc @ (L.PpC @ e0) + L.Acf @ (L.PpF @ e0)
            work[0] += (L.Acc.nnz + L.Acf.nnz + L.PpC.nnz + L.PpF.nnz) / nz
            e0 = e0 + self.cycle(lev + 1, rc - Ae0, np.zeros(len(C)), work)
        z = e0 - d[0]
        out = xk + h
        out[C] += L.PpC @ z
        out[F] += L.PpF @ z
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
            if res < tol or not np.isfinite(res) or res > 1e10:
                break
        return np.array(hist), np.array(works)
