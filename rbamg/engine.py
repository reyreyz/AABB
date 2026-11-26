"""One V-cycle engine for section 6.3 / Tables 3, SM13-SM15, built from PyAMG
components, with one work count for every family.

Configuration (keys of `cfg`):
  theta        classical strength threshold ('min' measure, as PyAMG's AIR)
  interp       'zero' (W = 0), '1pt' (one-point), 'dir' (direct)
  dP           number of F-Jacobi sweeps W <- W - D_ff^{-1}(A_ff W + A_fc)
  tauP         truncation of W (drop |w_ij| < tauP * row max, rescale row sums)
  restrict     'inj', 'air1', 'air2' (PyAMG local_air, degree 1/2), 'PT'
  pre, post    relaxation schedules, strings over {'F', 'C'} (Jacobi, omega = 1)
  tauA         lumped filter of the coarse operators
  filter_rule  'diag' (|a_ij| < tau |a_ii|) or 'rowmax'

Work (WU, units of nnz(A0)), counted on every level:
  F-sweep: nnz(A_FF), plus nnz(A_FC) once per block of consecutive F-sweeps
           (x_C is fixed during the block); the first sweep from a zero guess is free
  C-sweep: nnz(A_C:)
  residual: nnz(A) (only the C rows nnz(A_C:) for injection)
  restriction nnz(R) (0 for injection); interpolation nnz(W)
The coarsest system (< 150 unknowns) is factored; the hierarchy stops there.
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from pyamg.strength import classical_strength_of_connection
from pyamg.classical.split import RS
from pyamg.classical.interpolate import one_point_interpolation, direct_interpolation, local_air
from ml_cycles import filter_diag, filter_rowmax

COARSEST = 150


def truncate_rows(W, tau):
    if tau <= 0:
        return W
    W = sp.csr_matrix(W, copy=True)
    rs0 = np.asarray(W.sum(axis=1)).ravel()
    C = W.tocoo()
    mx = np.zeros(W.shape[0])
    np.maximum.at(mx, C.row, np.abs(C.data))
    keep = np.abs(C.data) >= tau * mx[C.row]
    Wt = sp.csr_matrix((C.data[keep], (C.row[keep], C.col[keep])), shape=W.shape)
    rs1 = np.asarray(Wt.sum(axis=1)).ravel()
    scale = np.where(np.abs(rs1) > 1e-14, rs0 / np.where(rs1 == 0, 1, rs1), 1.0)
    return (sp.diags(scale) @ Wt).tocsr()


class Level:
    pass


class VCycle:
    def __init__(self, A, cfg):
        self.cfg = cfg
        self.nnz0 = A.nnz
        self.levels = []
        A = sp.csr_matrix(A)
        while True:
            L = Level()
            L.A = A
            self.levels.append(L)
            if A.shape[0] < COARSEST or len(self.levels) > 30:
                L.lu = sla.splu(sp.csc_matrix(A))
                break
            S = classical_strength_of_connection(A, theta=cfg["theta"], norm="min")
            split = RS(S, second_pass=True)
            C = np.nonzero(split == 1)[0]
            F = np.nonzero(split == 0)[0]
            if len(C) == 0 or len(F) == 0:
                L.lu = sla.splu(sp.csc_matrix(A))
                break
            Aff, Afc = A[F][:, F].tocsr(), A[F][:, C].tocsr()
            if cfg["interp"] == "zero":
                W = sp.csr_matrix((len(F), len(C)))
            else:
                Pf = (one_point_interpolation(A, S, split) if cfg["interp"] == "1pt"
                      else direct_interpolation(A, S, split))
                W = sp.csr_matrix(Pf)[F]
            dinv = sp.diags(1.0 / Aff.diagonal())
            for _ in range(cfg.get("dP", 0)):
                W = (W - dinv @ (Aff @ W + Afc)).tocsr()
            W = truncate_rows(W, cfg.get("tauP", 0.0))
            W.eliminate_zeros()
            L.C, L.F, L.W = C, F, W
            L.AF, L.AC = A[F].tocsr(), A[C].tocsr()
            L.Aff, L.Afc = Aff, Afc
            L.dF = Aff.diagonal()
            L.dC = A[C][:, C].diagonal()
            P = sp.lil_matrix((A.shape[0], len(C)))
            Pfull = sp.vstack([sp.eye(len(C)), W]).tocsr()
            perm = np.concatenate([C, F])
            inv = np.empty_like(perm)
            inv[perm] = np.arange(len(perm))
            P = Pfull[inv]                                  # natural ordering
            r = cfg["restrict"]
            if r == "inj":
                L.R = None
                Ac = (A[C] @ P).tocsr()
            elif r.startswith("air"):
                R = local_air(A, split, theta=0.05, degree=int(r[3:]))
                L.R = sp.csr_matrix(R)
                Ac = (L.R @ A @ P).tocsr()
            else:
                L.R = P.T.tocsr()
                Ac = (L.R @ A @ P).tocsr()
            tau = cfg.get("tauA", 0.0)
            if tau > 0:
                Ac = filter_diag(Ac, tau) if cfg.get("filter_rule", "diag") == "diag" else filter_rowmax(Ac, tau)
            Ac.eliminate_zeros()
            A = Ac

    def op_complexity(self):
        return sum(L.A.nnz for L in self.levels) / self.nnz0

    def _relax(self, L, b, x, sched, work, zero):
        i = 0
        while i < len(sched):
            if sched[i] == "F":
                j = i
                while j < len(sched) and sched[j] == "F":
                    j += 1
                g = b[L.F] - L.Afc @ x[L.C]
                if not (zero and not np.any(x[L.C])):
                    work[0] += L.Afc.nnz / self.nnz0
                for s in range(j - i):
                    if zero and s == 0 and i == 0:
                        x[L.F] = g / L.dF               # free: from a zero guess
                    else:
                        x[L.F] += (g - L.Aff @ x[L.F]) / L.dF
                        work[0] += L.Aff.nnz / self.nnz0
                i = j
            else:
                x[L.C] += (b[L.C] - L.AC @ x) / L.dC
                work[0] += L.AC.nnz / self.nnz0
                i += 1
        return x

    def cycle(self, lev, b, x, work, zero=False):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            return L.lu.solve(b)
        x = self._relax(L, b, x, self.cfg["pre"], work, zero)
        if L.R is None:
            rc = b[L.C] - L.AC @ x
            work[0] += L.AC.nnz / self.nnz0
        else:
            r = b - L.A @ x
            rc = L.R @ r
            work[0] += (L.A.nnz + L.R.nnz) / self.nnz0
        e = self.cycle(lev + 1, rc, np.zeros(len(L.C)), work, zero=True)
        x[L.C] += e
        x[L.F] += L.W @ e
        work[0] += L.W.nnz / self.nnz0
        x = self._relax(L, b, x, self.cfg["post"], work, False)
        return x

    def solve(self, b, tol=1e-10, maxcyc=150):
        A = self.levels[0].A
        x = np.zeros(A.shape[0])
        nb = np.linalg.norm(b)
        hist, works, work = [1.0], [0.0], [0.0]
        for _ in range(maxcyc):
            x = self.cycle(0, b, x, work, zero=not np.any(x))
            r = np.linalg.norm(b - A @ x) / nb
            hist.append(r)
            works.append(work[0])
            if r < tol or not np.isfinite(r) or r > 1e8:
                break
        return np.array(hist), np.array(works)


def describe(cfg):
    s = {"zero": "W=0", "1pt": "1pt", "dir": "dir."}[cfg["interp"]]
    if cfg.get("dP"):
        s += f"+{cfg['dP']}F"
    if cfg.get("tauP"):
        s += f" tauP={cfg['tauP']}"
    s += " " + {"inj": "inj.", "air1": "lAIR1", "air2": "lAIR2", "PT": "P^T"}[cfg["restrict"]]
    s += f" {cfg['pre'] or '-'}|{cfg['post'] or '-'} theta={cfg['theta']} tauA={cfg.get('tauA', 0)}"
    return s
