"""Two-level: which test vectors for S~ (finest level) reduce N's factor while the twin still diverges."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); x, y = X.ravel(), Y.ravel()
psi = (1 - x**2) * (1 - y**2)
def rho(ml, wh, vec=False):
    op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, z, v, wh))
    try:
        vals, vecs = sla.eigs(op, k=1, which="LM", maxiter=1000, tol=1e-3, ncv=40)
        return (abs(vals[0]), vecs[:, 0]) if vec else abs(vals[0])
    except Exception: return (float("nan"), None) if vec else float("nan")
base = dict(omega=1.0, omega_c=0.8, k0=1, delta=0.3, mode="lin", max_levels=2, tau_filter=0)
ml0 = ML(A, xy=np.column_stack([x, y]), nkeep=2, **base)
r0, v0 = rho(ml0, "N", vec=True)
sets = {"1,x,y": [x, y], "1,x,y,psi": [x, y, psi], "1,psi": [psi], "1,x,y,xy,x2-y2": [x, y, x * y, x**2 - y**2],
        "1,x,y + slow N mode": [x, y, v0.real, v0.imag]}
print(f"\nm={m}\n\n| test vectors | entries/row | P~ nnz/nnz(A) | N | twin | Markov-Jacobi |\n|---|---|---|---|---|---|")
for lab, vs in sets.items():
    for nk in (2, 3):
        ml = ML(A, xy=np.column_stack([x, y]), tv=np.column_stack(vs), nkeep=nk, **base)
        print(f"| {lab} | {nk} | {ml.opcx()[1]:.2f} | {rho(ml,'N'):.3g} | {rho(ml,'twin'):.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
