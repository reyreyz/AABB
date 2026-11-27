"""Idea 3: constant-flow centered CD (no vortex), cell Pe fixed; two-level N / twin / MJ, Claim 11b level-0 setup."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from ml_twins import ML
def const_centered(m, pe, ang=0.4):
    h = 2.0 / (m + 1); n = m * m; bx, by = np.cos(ang), np.sin(ang)
    eps = h / pe; I, J = np.divmod(np.arange(n), m)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [4 * eps * np.ones(n)]
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        rows.append(np.nonzero(ok)[0]); cols.append(((I + di) * m + J + dj)[ok]); vals.append((-eps + s * 0.5 * h * b) * np.ones(ok.sum()))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n))
m = int(sys.argv[1]); pe = float(sys.argv[2])
A = const_centered(m, pe); n = A.shape[0]; zz = np.zeros(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-3, ncv=40)[0])
    except Exception: return float("nan")
for wc in (0.4, 0.6, 0.8):
    for nk in (1, 2):
        ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, mode="lin", xy=xy, nkeep=nk, max_levels=2, tau_filter=0)
        print(f"| const Pe {pe:g} | {m} | {wc} | {nk} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
