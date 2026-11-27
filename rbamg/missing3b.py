import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda x: ml.cycle(0, z, x, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-3, ncv=40)[0])
    except Exception: return float("nan")
for wc in (0.7, 0.8, 0.9):
    ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, tauP=1.0, mode="lump", max_levels=2, tau_filter=0)
    print(f"| {m} | {wc} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
