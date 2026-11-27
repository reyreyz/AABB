import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); cases = [c.split(":") for c in sys.argv[2].split(",")]   # pe:theta:k:wc
for pe, th, k, wc in cases:
    pe, th, k, wc = float(pe), float(th), int(k), float(wc)
    A, xy = cd(m, "const", pe); n = A.shape[0]; zz = np.zeros(n)
    I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
    def rho(ml, wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    ml = ML(A, omega=1.0, omega_c=wc, k0=k, split0=rb, mode="lin", xy=xy, theta=th, max_levels=2, tau_filter=0)
    print(f"| {pe:g} | {th} | {k} | {wc} | {m} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
