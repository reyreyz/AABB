"""Idea 1: slow mode of two-level N (Claim 11b setup) - smoothness and localization."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); wc = float(sys.argv[2]) if len(sys.argv) > 2 else 0.6
A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); x, y = X.ravel(), Y.ravel()
ml = ML(A, omega=1.0, omega_c=wc, k0=1, delta=0.3, mode="lin", xy=np.column_stack([x, y]), nkeep=2, max_levels=2, tau_filter=0)
op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, z, v, "N"))
vals, vecs = sla.eigs(op, k=3, which="LM", maxiter=2000, tol=1e-4, ncv=60)
L = ml.levels[0]
for j in np.argsort(-abs(vals)):
    v = vecs[:, j]; p = abs(v) ** 2 / np.sum(abs(v) ** 2)
    dist_b = np.minimum.reduce([1 - abs(x), 1 - abs(y)])
    r = np.hypot(x, y)
    sm = np.linalg.norm(A @ v) / (abs(A).sum(0).max() * np.linalg.norm(v))
    ic = np.argmax(p); 
    print(f"m={m} wC={wc}: |lam|={abs(vals[j]):.3f} arg={np.angle(vals[j]):.2f} | smooth ‖Av‖/‖A‖‖v‖={sm:.2e} | C-frac={p[L.C].sum():.2f} | "
          f"mass near boundary (<5h)={p[dist_b < 5*h].sum():.2f}, r<0.2={p[r < 0.2].sum():.2f}, 0.2<r<0.7={p[(r>=0.2)&(r<0.7)].sum():.2f} | "
          f"peak at ({x[ic]:.2f},{y[ic]:.2f}) | mass within 10h of peak={p[np.hypot(x-x[ic], y-y[ic]) < 10*h].sum():.2f}")
    np.save(f"/tmp/claude-0/-home-user-AABB/72d75a69-cd12-59c4-ba81-ce7a8de908cb/scratchpad/slow_{m}_{j}.npy", v)
