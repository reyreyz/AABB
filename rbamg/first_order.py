"""N vs Markovian 'k sweeps + 1 F-sweep' (first-order equivalent in w_C): two-level rho on the red-black candidates."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1])
for pe, th, k, wc in ((6, None, 2, 0.2), (6, None, 2, 0.05), (6, None, 2, 0.5), (12, 0.5, 8, 0.05), (12, 0.5, 8, 0.02), (12, 0.5, 4, 0.1)):
    A, xy = cd(m, "const", pe); n = A.shape[0]; zz = np.zeros(n); I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
    kw = dict(nkeep=1) if th is None else dict(theta=th)
    ml = ML(A, omega=1.0, omega_c=wc, k0=k, split0=rb, mode="lin", xy=xy, max_levels=2, tau_filter=0, **kw)
    def rho(wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    print(f"| {pe} | {th} | {k} | {wc} | {m} | {rho('N'):.3g} | {rho('markovJF'):.3g} | {rho('markovJ'):.3g} | {rho('twin'):.3g} |"); sys.stdout.flush()
