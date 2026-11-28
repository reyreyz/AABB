"""Is the history shift in the interpolation necessary?  Two-level, dense, exact coarse solve, same smoother T
(damped Jacobi w), same greedy split, same Rhat = [I, -A_cf (I+G+G^2) D_ff^{-1}].
  M(P0)   = Pi_R T^k                                   (injection)
  M'      = (I - P'(Rhat A P')^{-1} Rhat A) T^k,  P' = P0 + Q0 sum_{l<k} K^l H   (static memory interpolation)
  M(ideal)= same with W = -A_ff^{-1} A_fc
  S       = (Pi_R T Pi_R)^k                            (Alg. 4.1: shift relation + Horner replay)
  S-noshift: eps_l = eps_k for all l  ->  e <- e_k - [sum_{j=0}^{k} (Pi_R T)^j P0] (Rhat A P0)^{-1} Rhat A e_k
  rho(K_can) = rho(Pi_R T Pi_R) (CR rate of the canonical complement; 0 would mean trivial memory)."""
import sys, itertools, numpy as np, scipy.sparse as sp
import sparse_problems as spp
from robust_n import greedy_f_points
from flows import cd as cdflow
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 31
PROBS_ALL = {"shear Pe 20": lambda: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 20, "2", False), "centered"),
         "shear Pe 10": lambda: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 10, "2", False), "centered"),
         "const Pe 20": lambda: cdflow(m, "const", 20)[0], "recirc Pe 20": lambda: cdflow(m, "recirc", 20)[0]}
probs = {k_: v for k_, v in PROBS_ALL.items() if k_ in sys.argv[2].split(",")} if len(sys.argv) > 2 else PROBS_ALL
WS = [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else (0.3, 0.5, 0.7, 0.9)
KS = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else (2, 4)
print(f"m={m}\n| problem | w | k | rho(T) | rho(K_can) | M(P0) | M' | M(ideal) | S | S-noshift |\n|---|---|---|---|---|---|---|---|---|---|")
for (name, mk), w, k in itertools.product(probs.items(), WS, KS):
    A = sp.csr_matrix(mk()); isF = greedy_f_points(A, 0.3); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad)
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nc, nf = len(C), len(F)
    M = np.diag(w / d); T = np.eye(n) - M @ Ad
    Aff, Acf, Afc = Ad[np.ix_(F, F)], Ad[np.ix_(C, F)], Ad[np.ix_(F, C)]
    Dfi = np.diag(1 / np.diag(Aff)); G = np.eye(nf) - Dfi @ Aff
    Rh = np.zeros((nc, n)); Rh[np.arange(nc), C] = 1; Rh[:, F] = -Acf @ (np.eye(nf) + G + G @ G) @ Dfi
    def Pof(W):
        P = np.zeros((n, nc)); P[C, np.arange(nc)] = 1; P[F] = W; return P
    def proj(P): return np.eye(n) - P @ np.linalg.solve(Rh @ Ad @ P, Rh @ Ad)
    K0, H0 = T[np.ix_(F, F)], T[np.ix_(F, C)]; Sk = sum(np.linalg.matrix_power(K0, l) for l in range(k)) @ H0
    P0 = Pof(np.zeros((nf, nc))); Pi = proj(P0); Tk = np.linalg.matrix_power(T, k)
    rM0 = rho(Pi @ Tk); rMp = rho(proj(Pof(Sk)) @ Tk); rMi = rho(proj(Pof(-np.linalg.solve(Aff, Afc))) @ Tk)
    PiT = Pi @ T; rS = rho(np.linalg.matrix_power(Pi @ T @ Pi, k))
    Rep = sum(np.linalg.matrix_power(PiT, j) for j in range(k + 1)) @ P0
    rN = rho((np.eye(n) - Rep @ np.linalg.solve(Rh @ Ad @ P0, Rh @ Ad)) @ Tk)
    f = lambda r: f"{r:.3g}"
    print(f"| {name} | {w} | {k} | {rho(T):.2f} | {rho(Pi @ T @ Pi):.3f} | {f(rM0)} | {f(rMp)} | {f(rMi)} | {f(rS)} | {f(rN)} |"); sys.stdout.flush()
