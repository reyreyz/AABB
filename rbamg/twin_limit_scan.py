"""Dense scan: rho(T) < 1 and rho(E_twin, k->inf) > 1 ?  E_twin_inf = (I - Pi~)(P~ + Q (S_inf - S~)) P^sharp."""
import numpy as np, scipy.sparse as sp, sys
from robust_n import recirc, greedy_f_points
from ml_twins import truncate
m = int(sys.argv[1]) if len(sys.argv) > 1 else 31
print("| Pe | omega | rho(T) | delta | trunc | rho(K) | rho(E_twin,inf) |\n|---|---|---|---|---|---|---|")
for pe in (3, 5, 10, 20):
    A = recirc(m, pe); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad)
    for w in (0.1, 0.3, 0.5, 0.8):
        T = np.eye(n) - w * Ad / d[:, None]
        rT = np.abs(np.linalg.eigvals(T)).max()
        if rT >= 1: continue
        for delta in (0.3, 0.6):
            isF = greedy_f_points(A, delta); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]
            K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
            rK = np.abs(np.linalg.eigvals(K)).max()
            Sinf = np.linalg.solve(np.eye(len(F)) - K, H)          # sum_{l>=0} K^l H
            for tau, mode in ((0.1, "drop"), (0.3, "drop"), (0.3, "lump"), (1.0, "lump")):
                St = truncate(sp.csr_matrix(Sinf), tau, mode).toarray()
                Pt = np.zeros((n, len(C))); Pt[C, np.arange(len(C))] = 1; Pt[F] = St
                Pi = Pt @ np.linalg.solve(Pt.T @ Ad @ Pt, Pt.T @ Ad)
                Pinf = np.zeros((n, len(C))); Pinf[C, np.arange(len(C))] = 1; Pinf[F] = Sinf
                Psh = np.zeros((len(C), n)); Psh[np.arange(len(C)), C] = 1
                E = (np.eye(n) - Pi) @ Pinf @ Psh
                print(f"| {pe} | {w} | {rT:.3f} | {delta} | {tau:g} {mode} | {rK:.3f} | {np.abs(np.linalg.eigvals(E)).max():.3g} |"); sys.stdout.flush()
