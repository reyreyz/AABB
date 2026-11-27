import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); pes = [float(p) for p in sys.argv[2].split(",")]; wcs = [float(w) for w in sys.argv[3].split(",")]
for pe in pes:
    A, xy = cd(m, "const", pe); n = A.shape[0]; zz = np.zeros(n)
    def rho(ml, wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    for wc in wcs:
        ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, mode="lin", xy=xy, nkeep=2, max_levels=2, tau_filter=0)
        print(f"| const | {pe:g} | {m} | {wc} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
