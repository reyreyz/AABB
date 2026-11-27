"""recirc centered Pe 10, Jacobi w_F = 1 with C-weight w_C, delta=0.3, k=1, single-entry lumped P~, R = P~^T."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda x: ml.cycle(0, z, x, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-4, ncv=60)[0])
    except Exception: return float("nan")
print(f"\nm={m}\n\n| w_C | k | tau | rho(E_N) | rho(E_twin) | rho(E_MJ) |\n|---|---|---|---|---|---|")
for k in (1, 2):
    for wc in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        ml = ML(A, omega=1.0, omega_c=wc, k0=k, delta=0.3, tauP=1.0, mode="lump", max_levels=2, tau_filter=0)
        print(f"| {wc} | {k} | 1 lump | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
