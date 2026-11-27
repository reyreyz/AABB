import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, z, v, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-3, ncv=40)[0])
    except Exception: return float("nan")
for wc in (0.3, 0.4, 0.5, 0.6):
    ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, mode="lin", xy=xy, nkeep=2, max_levels=2, tau_filter=0)
    print(f"| {m} | {wc} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
