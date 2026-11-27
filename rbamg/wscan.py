import sys, itertools, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, z, v, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1000, tol=1e-3, ncv=40)[0])
    except Exception: return float("nan")
print(f"\nm={m}\n\n| delta | w_F | w_C | N | twin |\n|---|---|---|---|---|")
for d, wf, wc in itertools.product((0.2, 0.3, 0.45), (0.9, 1.0), (0.6, 0.7, 0.8, 0.9)):
    ml = ML(A, omega=wf, omega_c=wc, k0=1, delta=d, mode="lin", xy=xy, nkeep=2, max_levels=2, tau_filter=0)
    rn = rho(ml, "N")
    rt = rho(ml, "twin") if rn < 0.6 else float("nan")
    print(f"| {d} | {wf} | {wc} | {rn:.3g} | {rt:.3g} |"); sys.stdout.flush()
