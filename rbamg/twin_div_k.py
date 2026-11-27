"""Pe 10 centered recirc, damped Jacobi, single-entry lumped P~ (tau = 1): rho of N, twin, Markov-Jacobi vs k."""
import numpy as np, scipy.sparse as sp, sys
from robust_n import recirc, greedy_f_points
from ml_twins import truncate
m, pe = int(sys.argv[1]), float(sys.argv[2])
A = recirc(m, pe); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad)
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
for w in (0.05, 0.1, 0.15, 0.2):
    print(f"omega={w}: rho(T) = {rho(np.eye(n) - w * Ad / d[:, None]):.5f}")
print("\n| omega | delta | k | rho(K) | rho(E_N) | rho(E_twin) | rho(E_MJ) | rho(E_twin) k->inf |\n|---|---|---|---|---|---|---|---|")
for w in (0.1,):
    T = np.eye(n) - w * Ad / d[:, None]
    for delta in (0.3, 0.6):
        isF = greedy_f_points(A, delta); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]
        K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]; nf, nc = len(F), len(C)
        Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1
        Q = np.zeros((n, nf)); Q[F, np.arange(nf)] = 1
        def Pfrom(S):
            P = np.zeros((n, nc)); P[C, np.arange(nc)] = 1; P[F] = S; return P
        Sinf = np.linalg.solve(np.eye(nf) - K, H)
        Pt = Pfrom(truncate(sp.csr_matrix(Sinf), 1.0, "lump").toarray())
        lim = rho((np.eye(n) - Pt @ np.linalg.solve(Pt.T @ Ad @ Pt, Pt.T @ Ad)) @ Pfrom(Sinf) @ Psh)
        S, V, Tk, Kk = H.copy(), H.copy(), T.copy(), K.copy()
        for k in range(1, 257):
            if k > 1: V = K @ V; S = S + V; Tk = T @ Tk; Kk = K @ Kk
            if k not in (4, 16, 64, 256): continue
            St = truncate(sp.csr_matrix(S), 1.0, "lump").toarray()
            Pt = Pfrom(St); I_Pi = np.eye(n) - Pt @ np.linalg.solve(Pt.T @ Ad @ Pt, Pt.T @ Ad)
            Pex = Pfrom(S)
            # N with truncated correction: (I - Pi~)[P' P^sharp T^k ... ] = (I-Pi~)(P' Psh T^k + Q K^k Q^T) with exact history
            EN = I_Pi @ (Pex @ Psh @ Tk + Q @ Kk @ Q.T) + 0
            # careful: N's pre-correction error is P' c_k + Q K^k q0 where c_k = Psh T^k e0 and q0 = Q^T e0
            ET = I_Pi @ (Pex @ Psh + Q @ Kk @ Q.T)
            EM = I_Pi @ Tk
            print(f"| {w} | {delta} | {k} | {rho(K):.3f} | {rho(EN):.3g} | {rho(ET):.3g} | {rho(EM):.3g} | {lim:.3g} |"); sys.stdout.flush()
