import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]; z = np.zeros(n)
print(f"\nm={m}\n\n| delta | trunc | k | P~ nnz/nnz(A) | rho(E_N) | rho(E_twin) | rho(E_MJ) |\n|---|---|---|---|---|---|---|")
for delta in (0.3, 0.6):
    for tau in (1.0, 0.6, 0.3):
        for k in (1, 2):
            ml = ML(A, omega=1.0, k0=k, delta=delta, tauP=tau, mode="lump", max_levels=2, tau_filter=0)
            r = []
            for wh in ("N", "twin", "markovJ"):
                op = sla.LinearOperator((n, n), matvec=lambda x, wh=wh: ml.cycle(0, z, x, wh))
                try: r.append(abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=400, tol=1e-4, ncv=40)[0]))
                except Exception as e: r.append(np.nan)
            print(f"| {delta} | {tau:g} lump | {k} | {ml.opcx()[1]:.2f} | {r[0]:.3g} | {r[1]:.3g} | {r[2]:.3g} |"); sys.stdout.flush()
