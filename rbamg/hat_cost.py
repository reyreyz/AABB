"""Cost of the effective Markovian smoother M^ vs N: (1) applied matrix-free (work per cycle), (2) as an explicit matrix."""
import numpy as np, scipy.sparse as sp
from ml_twins import ML
from flows import cd
from robust_n import recirc
# (1) work per cycle, same hierarchy, multilevel
for lab, A, kw in (("Claim 13 (const Pe 6, red-black, k=2)", None, None), ("Claim 11 (recirc Pe 10, k=1)", None, None)):
    pass
m = 127
A, xy = cd(m, "const", 6.0); n = A.shape[0]; I, J = np.divmod(np.arange(n), m); b = np.random.default_rng(0).standard_normal(n)
cfgs = [("const Pe 6, red-black, k=2 (Claim 13)", A, dict(omega=1.0, omega_c=0.2, wc_levels=1, k0=2, split0=(I + J) % 2 == 1, mode="lin", xy=xy, nkeep=1,
         lin_levels=1, k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0))]
A2 = recirc(m, 10.0); h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t)
cfgs.append(("recirc Pe 10, greedy, k=1 (Claim 11)", A2, dict(omega=1.0, omega_c=0.6, wc_levels=1, k0=1, k1=8, delta=0.3, mode="lin",
             xy=np.column_stack([X.ravel(), Y.ravel()]), nkeep=2, lin_levels=1, tauP_deep=1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=0)))
cfgs.append(("recirc Pe 20, greedy, k=4/8, untruncated level 0 (robust N)", recirc(m, 20.0), dict(omega=1.0, k0=4, k1=8, delta=0.3, tauP=1e-5,
             trunc_from=1, mode="drop", tau_filter=1e-3, tau_deep=3e-4, refine=1)))
print("| hierarchy | N: WU/cycle | M^ cycle: WU/cycle | iterates equal |\n|---|---|---|---|")
for lab, AA, kw in cfgs:
    ml = ML(AA, **kw); bb = np.random.default_rng(0).standard_normal(AA.shape[0])
    hN = ml.solve(bb, "N", maxcyc=3); wN = ml.wpc; xN = np.zeros(AA.shape[0])
    hH = ml.solve(bb, "hat", maxcyc=3); wH = ml.wpc
    x1 = ml.cycle(0, bb, np.zeros(AA.shape[0]), "N"); x2 = ml.cycle(0, bb, np.zeros(AA.shape[0]), "hat")
    print(f"| {lab} | {wN:.0f} | {wH:.0f} | {np.linalg.norm(x1 - x2)/np.linalg.norm(x1):.0e} |")
# (2) explicit matrix: nnz of M^ vs M, two-level, small grid
print("\n| case | k | nnz(M)/n | nnz(twin smoother)/n | nnz(M^)/n | nnz(A)/n |\n|---|---|---|---|---|---|")
for k in (1, 2, 4):
    m = 31; A, xy = cd(m, "const", 6.0); n = A.shape[0]; I, J = np.divmod(np.arange(n), m); isF = (I + J) % 2 == 1
    Ad = A.toarray(); d = np.diag(Ad); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    W = np.where(isF, 1.0, 0.2); M = np.diag(W / d); T = np.eye(n) - M @ Ad
    Qm = np.zeros((n, nf)); Qm[F, np.arange(nf)] = 1; Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1
    TQ = np.eye(n) - Qm @ Qm.T @ M @ Ad
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]; S = H.copy(); V = H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    from ml_twins import truncate_lin
    St = truncate_lin(sp.csr_matrix(S), np.ones((nc, 1)), 1).toarray()
    Mtw = sum(np.linalg.matrix_power(TQ, j) for j in range(k)) @ Qm @ Qm.T @ M
    Mhat = Mtw + Qm @ (S - St) @ Psh @ sum(np.linalg.matrix_power(T, j) for j in range(k)) @ M
    nz = lambda X: np.count_nonzero(np.abs(X) > 1e-14 * np.abs(X).max()) / n
    print(f"| const Pe 6, red-black, m=31 | {k} | {nz(M):.1f} | {nz(Mtw):.1f} | {nz(Mhat):.1f} | {nz(Ad):.1f} |")
