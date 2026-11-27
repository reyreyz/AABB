import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from hitgen import problem, H1
m = int(sys.argv[1]); A, xy = problem(m); n = A.shape[0]; zz = np.zeros(n); P = H1
ml = ML(A, omega=P["wf"], omega_c=P["wc"], k0=P["k"], delta=P["delta"], mode="lin", xy=xy, theta=P["theta"], max_levels=2, tau_filter=0)
def rho(wh):
    op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=2000, tol=1e-3, ncv=40)[0])
    except Exception: return float("nan")
print(f"| {m} | {rho('N'):.3g} | {rho('twin'):.3g} | {rho('markovJ'):.3g} | {rho('markovJF'):.3g} |"); sys.stdout.flush()
