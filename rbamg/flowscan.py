"""Two-level N / twin / MJ (Claim 11b level-0 setup) for several flows, cell Pe, w_C."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
flow, m = sys.argv[1], int(sys.argv[2]); pes = [float(p) for p in sys.argv[3].split(",")]
for pe in pes:
    A, xy = cd(m, flow, pe); n = A.shape[0]; zz = np.zeros(n)
    def rho(ml, wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    for wc in (0.4, 0.6):
        ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, mode="lin", xy=xy, nkeep=2, max_levels=2, tau_filter=0)
        rn = rho(ml, "N")
        print(f"| {flow} | {pe:g} | {m} | {wc} | {rn:.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
