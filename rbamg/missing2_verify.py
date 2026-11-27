"""Shifted Laplacian -Lap - kappa^2, kappa^2 = 60 fixed: Jacobi w=1, greedy delta=0.3, k=4, lumped truncation."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from ml_twins import ML
m = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else "two-level"
h = 1.0 / (m + 1)
L2 = (sp.kron(sp.eye(m), sp.diags([-1, 2, -1], [-1, 0, 1], (m, m))) + sp.kron(sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)), sp.eye(m))).tocsr()
A = (L2 - 60 * h * h * sp.eye(m * m)).tocsr(); n = A.shape[0]; z = np.zeros(n)
kw = dict(max_levels=2, tau_filter=0) if mode == "two-level" else dict(tau_deep=3e-4)
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda x: ml.cycle(0, z, x, wh))
    try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=800, tol=1e-4, ncv=40)[0])
    except Exception: return float("nan")
out = []
for tau in (0.6, 1.0):
    ml = ML(A, omega=1.0, k0=4, k1=4, delta=0.3, tauP=tau, mode="lump", restrict="inj", **kw)
    out.append(f"tau={tau:g}: N {rho(ml,'N'):.3g}, twin {rho(ml,'twin'):.3g}, MJ {rho(ml,'markovJ'):.3g} (P~ {ml.opcx()[1]:.2f})")
print(f"| {m} | {mode} | " + " ; ".join(out) + " |"); sys.stdout.flush()
best = []
for R in ("inj", "ppt"):
    for w in (0.2, 0.4, 0.6, 0.8, 1.0):
        for k in (1, 2, 4, 8):
            for tau in (0.6, 1.0):
                mm = ML(A, omega=w, k0=k, k1=k, delta=0.3, tauP=tau, mode="lump", restrict=R, **kw)
                for wh in ("markovJ", "twin"):
                    best.append((rho(mm, wh), wh, R, w, k, tau))
best = sorted([b for b in best if np.isfinite(b[0])])
print("   best Markovian on the same truncation family (rho, method, R, w, k, tau):", [(round(b[0], 3),) + b[1:] for b in best[:3]]); sys.stdout.flush()
