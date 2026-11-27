import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
def rho(ml, wh):
    op = sla.LinearOperator((n, n), matvec=lambda x: ml.cycle(0, z, x, wh))
    vals, vecs = sla.eigs(op, k=1, which="LM", maxiter=600, tol=1e-4, ncv=40)
    return abs(vals[0]), vecs[:, 0]
print(f"\nm={m}\n\n| omega | delta | tau | k | rho(E_N) | eigvec: ‖v_C‖/‖v‖ | ‖A v‖/(‖A‖_1 ‖v‖) | rho(E_twin) |\n|---|---|---|---|---|---|---|---|")
for w in (1.0, 0.8):
    for tau in (0.0, 0.1, 0.3):
        ml = ML(A, omega=w, k0=1, delta=0.3, tauP=tau, mode="lump", max_levels=2, tau_filter=0)
        try:
            r, v = rho(ml, "N"); rt, _ = rho(ml, "twin")
            L = ml.levels[0]; vC = np.linalg.norm(v[L.C]) / np.linalg.norm(v)
            sm = np.linalg.norm(A @ v) / (abs(A).sum(0).max() * np.linalg.norm(v))
            print(f"| {w} | 0.3 | {tau:g} | 1 | {r:.3g} | {vC:.2f} | {sm:.2e} | {rt:.3g} |")
        except Exception as e:
            print(f"| {w} | 0.3 | {tau:g} | 1 | fail {e} |")
        sys.stdout.flush()
