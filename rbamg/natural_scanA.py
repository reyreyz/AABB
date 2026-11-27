"""Class A only (N converges, twin diverges, Markov-Jacobi converges): damped Jacobi inside its stability range."""
import sys, itertools, numpy as np, scipy.sparse as sp
from robust_n import greedy_f_points, drop_rel
from ml_twins import truncate, truncate_lin
from flows import cd
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1])
print(f"m={m}\n| flow | Pe | w | split | k | trunc | rho(T) | N | twin | MJ |\n|---|---|---|---|---|---|---|---|---|---|")
for (flow, pe), w, split, k, tr in itertools.product([("const", 6), ("const", 8), ("const", 10), ("recirc", 10), ("recirc", 20), ("poiseuille", 10), ("poiseuille", 20)],
                                                     (0.05, 0.1, 0.2, 0.3), ("rb", "greedy"), (2, 4, 8, 16), ("onepoint", "lump0.5")):
    A, xy = cd(m, flow, pe); n = A.shape[0]; I, J = np.divmod(np.arange(n), m)
    isF = ((I + J) % 2 == 1) if split == "rb" else greedy_f_points(A, 0.3)
    Ad = A.toarray(); d = np.diag(Ad); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    T = np.eye(n) - w * Ad / d[:, None]; rT = rho(T)
    if rT >= 1: continue
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]; S = H.copy(); V = H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    Ss = sp.csr_matrix(S)
    St = (truncate_lin(Ss, np.ones((nc, 1)), 1) if tr == "onepoint" else truncate(Ss, 0.5, "lump")).toarray()
    Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St; Ac = Pt.T @ Ad @ Pt
    if np.linalg.cond(Ac) > 1e12: continue
    I_Pi = np.eye(n) - Pt @ np.linalg.solve(Ac, Pt.T @ Ad); I_PiQ = I_Pi[:, F]
    Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
    Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
    rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh))
    if rN >= 0.8: continue
    rTw = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
    if rTw <= 1: continue
    print(f"| {flow} | {pe} | {w} | {split} | {k} | {tr} | {rT:.4f} | {rN:.3g} | {rTw:.3g} | {rho(I_Pi @ Tk):.3g} |"); sys.stdout.flush()
