"""M vs S (semi-Markovian: memory in the interpolation only), two-level, dense, natural problems.
Same smoother (damped Jacobi, uniform w), same greedy split, same restriction
Rhat = [I, -A_cf (I + G + ... + G^{L-1}) D_ff^{-1}],  G = I - D_ff^{-1} A_ff,  same coarse operator Rhat A P0.
  E^M = Pi_R T^k,   E^S = (Pi_R T Pi_R)^k,   Pi_R = I - P0 (Rhat A P0)^{-1} Rhat A."""
import sys, itertools, numpy as np, scipy.sparse as sp
import sparse_problems as spp
from robust_n import greedy_f_points
from flows import cd as cdflow
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1])
def poisson(m, ey=1.0):
    T1 = sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)); return (sp.kron(sp.eye(m), T1) + ey * sp.kron(T1, sp.eye(m))).tocsr()
probs = [("Poisson", lambda: poisson(m)), ("aniso 0.01", lambda: poisson(m, 0.01))]
for pe in (2, 6, 20):
    probs.append((f"shear Pe {pe}", lambda pe=pe: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, pe, "2", False), "centered")))
    probs.append((f"recirc Pe {pe}", lambda pe=pe: cdflow(m, "recirc", pe)[0]))
    probs.append((f"const Pe {pe}", lambda pe=pe: cdflow(m, "const", pe)[0]))
print(f"m={m}\n| problem | w | L | k | rho(T) | rho(K) | M | S | S work/digit vs M |\n|---|---|---|---|---|---|---|---|---|")
for (pn, mk), w, L, k in itertools.product(probs, (0.6, 0.8, 1.0), (1, 3), (1, 2, 4)):
    A = sp.csr_matrix(mk()); isF = greedy_f_points(A, 0.3); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad)
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nc = len(C)
    T = np.eye(n) - w * Ad / d[:, None]
    Aff, Acf = Ad[np.ix_(F, F)], Ad[np.ix_(C, F)]; Dfi = np.diag(1 / np.diag(Aff)); G = np.eye(len(F)) - Dfi @ Aff
    Z = sum(np.linalg.matrix_power(G, j) for j in range(L)) @ Dfi
    Rh = np.zeros((nc, n)); Rh[np.arange(nc), C] = 1; Rh[:, F] = -Acf @ Z
    P0 = np.zeros((n, nc)); P0[C, np.arange(nc)] = 1
    Ac = Rh @ Ad @ P0
    if np.linalg.cond(Ac) > 1e12: continue
    Pi = np.eye(n) - P0 @ np.linalg.solve(Ac, Rh @ Ad)
    rM = rho(Pi @ np.linalg.matrix_power(T, k)); rS = rho(np.linalg.matrix_power(Pi @ T @ Pi, k))
    K = (Pi @ T @ Pi)                       # CR operator in the canonical complement (restricted to range Pi)
    # work per digit: M = k sweeps + residual + coarse(~3); S = M + history (~k * 1.5) + replay (~k)
    wM = (k + 1 + 3) / -np.log10(rM) if rM < 1 else np.inf; wS = (k + 1 + 3 + 2.5 * k) / -np.log10(rS) if rS < 1 else np.inf
    tag = "RESCUE" if (rM >= 1 and rS < 0.7) else ("COST x%.1f" % (wM / wS) if (rS < 1 and wM / wS > 2) else "")
    if tag:
        print(f"| {pn} | {w} | {L} | {k} | {rho(T):.2f} | {rho(Pi @ T @ Pi):.2f} | {rM:.3g} | {rS:.3g} | {tag} |"); sys.stdout.flush()
