"""E_N - E_MJF_k = (I - Pi~) Q K^k Q^T (I - T^k)  (two-level, exact coarse solve, P0/Q0, any S~). Dense check."""
import numpy as np, scipy.sparse as sp, itertools
from robust_n import greedy_f_points
from ml_twins import truncate_lin
from flows import cd
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
for pe, delta, k, wc, nk in itertools.product((6, 20), (0.3, 0.6), (2, 4), (0.3, 1.0), (1, 2)):
    m = 15; A, xy = cd(m, "recirc", pe); isF = greedy_f_points(A, delta); A = A.toarray(); n = len(A); d = np.diag(A)
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    W = np.where(isF, 1.0, wc); T = np.eye(n) - (W / d)[:, None] * A
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]; S = H.copy(); V = H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    St = truncate_lin(sp.csr_matrix(S), np.column_stack([np.ones(nc), xy[C]]), nk).toarray()
    Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
    I_Pi = np.eye(n) - Pt @ np.linalg.solve(Pt.T @ A @ Pt, Pt.T @ A)
    Q = np.zeros((n, nf)); Q[F, np.arange(nf)] = 1; Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1
    Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
    EN = I_Pi @ Q @ ((S - St) @ Psh @ Tk + Kk @ Q.T)
    TQ = np.eye(n) - Q @ Q.T @ ((np.where(isF, 1.0, 0.0) / d)[:, None] * A)      # F-sweep, w_F = 1
    EMJF = I_Pi @ np.linalg.matrix_power(TQ, k) @ Tk
    pred = I_Pi @ Q @ Kk @ Q.T @ (np.eye(n) - Tk)
    err = np.linalg.norm(EN - EMJF - pred) / max(np.linalg.norm(EN), 1)
    print(f"| {pe} | {delta} | {k} | {wc} | {nk} | rho(T)={rho(T):.2f} rho(K)={rho(K):.2f} | N {rho(EN):.3g} | MJF_k {rho(EMJF):.3g} | identity err {err:.1e} |")
