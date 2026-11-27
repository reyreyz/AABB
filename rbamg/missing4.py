"""recirc centered Pe 10, Jacobi w_F=1 / w_C, delta=0.3, k=1, P~ = 3 entries/row exact on 1, x, y; R = P~^T."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda x: ml.cycle(0, z, x, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-3, ncv=40)[0])
    except Exception: return float("nan")
for nk in (2,):
    for wc in (0.8,):
        ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, mode="lin", xy=xy, nkeep=nk, max_levels=2, tau_filter=0)
        print(f"| {m} | {nk} | {wc} | {ml.opcx()[1]:.2f} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
