import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); cases = [c.split(":") for c in sys.argv[2].split(",")]   # flow:pe:k:wc:nk
for flow, pe, k, wc, nk in cases:
    pe, k, wc, nk = float(pe), int(k), float(wc), int(nk)
    A, xy = cd(m, flow, pe); n = A.shape[0]; zz = np.zeros(n)
    I, J = np.divmod(np.arange(n), m); isF = (I + J) % 2 == 1
    def rho(ml, wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    ml = ML(A, omega=1.0, omega_c=wc, k0=k, split0=isF, mode="lin", xy=xy, nkeep=nk, max_levels=2, tau_filter=0)
    print(f"| {flow} | {pe:g} | {m} | {k} | {wc} | {nk} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
