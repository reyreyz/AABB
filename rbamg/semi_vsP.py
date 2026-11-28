"""M with a good interpolation vs S, two-level dense.  Same smoother T (damped Jacobi), same greedy split, same
restriction Rhat = [I, -A_cf (I+G+G^2) D_ff^{-1}].  Interpolations P = [I; W]:
  P0 (injection), W1 = -D_ff^{-1} A_fc (one-step Jacobi), P' = P0 + Q0 sum_{l<k} K^l H (S's interpolation, static),
  ideal W = -A_ff^{-1} A_fc.
  E^M(P) = (I - P (Rhat A P)^{-1} Rhat A) T^k,   E^S(P) = (Pi_P T Pi_P)^k."""
import sys, numpy as np, scipy.sparse as sp
import sparse_problems as spp
from robust_n import greedy_f_points
from flows import cd as cdflow
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 31
def poisson(m, ey=1.0):
    T1 = sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)); return (sp.kron(sp.eye(m), T1) + ey * sp.kron(T1, sp.eye(m))).tocsr()
cases = [("shear Pe 20", lambda: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 20, "2", False), "centered"), 0.7, 4),
         ("const Pe 6", lambda: cdflow(m, "const", 6)[0], 0.8, 2), ("shear Pe 6", lambda: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 6, "2", False), "centered"), 0.8, 2),
         ("shear Pe 2", lambda: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 2, "2", False), "centered"), 0.8, 2),
         ("aniso 0.01", lambda: poisson(m, 0.01), 0.8, 4), ("Poisson", lambda: poisson(m), 0.8, 4)]
print(f"m={m}\n| problem | w, k | rho(T) | M(P0) | M(W1) | M(P') | M(ideal) | S(P0) | S(W1) |\n|---|---|---|---|---|---|---|---|---|")
for name, mk, w, k in cases:
    A = sp.csr_matrix(mk()); isF = greedy_f_points(A, 0.3); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad)
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nc, nf = len(C), len(F)
    T = np.eye(n) - w * Ad / d[:, None]
    Aff, Acf, Afc = Ad[np.ix_(F, F)], Ad[np.ix_(C, F)], Ad[np.ix_(F, C)]
    Dfi = np.diag(1 / np.diag(Aff)); G = np.eye(nf) - Dfi @ Aff
    Rh = np.zeros((nc, n)); Rh[np.arange(nc), C] = 1; Rh[:, F] = -Acf @ (np.eye(nf) + G + G @ G) @ Dfi
    def P_of(W):
        P = np.zeros((n, nc)); P[C, np.arange(nc)] = 1; P[F] = W; return P
    K0, H0 = T[np.ix_(F, F)], T[np.ix_(F, C)]                    # local complement Q0 with P0 (W = 0)
    Sk = sum(np.linalg.matrix_power(K0, l) for l in range(k)) @ H0
    Ws = {"P0": np.zeros((nf, nc)), "W1": -Dfi @ Afc, "P'": Sk, "ideal": -np.linalg.solve(Aff, Afc)}
    Tk = np.linalg.matrix_power(T, k); out = []
    Pis = {}
    for key, W in Ws.items():
        P = P_of(W); Ac = Rh @ Ad @ P
        Pi = np.eye(n) - P @ np.linalg.solve(Ac, Rh @ Ad); Pis[key] = Pi
        out.append(rho(Pi @ Tk))
    for key in ("P0", "W1"):
        out.append(rho(np.linalg.matrix_power(Pis[key] @ T @ Pis[key], k)))
    f = lambda r: f"{r:.3g}" + (" div" if r > 1 else "")
    print(f"| {name} | {w}, {k} | {rho(T):.2f} | " + " | ".join(f(r) for r in out) + " |"); sys.stdout.flush()
