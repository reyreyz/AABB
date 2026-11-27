"""Idea 19: potential well A = L + c(x), c = -gamma on a disk of radius r (cells): few localized bound states
amplified by every damped Jacobi.  Dense two-level: N (S~ exact on constants + bound-state C-parts), twin,
best damped Markov-Jacobi over (w, w_C ratio, k)."""
import numpy as np, scipy.sparse as sp, itertools, sys
from ml_twins import truncate_lin
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 23
T1 = sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)); L2 = (sp.kron(sp.eye(m), T1) + sp.kron(T1, sp.eye(m))).toarray()
n = m * m; I, J = np.divmod(np.arange(n), m); ctr = (m - 1) / 2
isF = (I + J) % 2 == 1; F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
print(f"m={m}\n| gamma | r | #neg eig | w_C | k | fit | N | twin | best damped MJ (w, ratio, k) |\n|---|---|---|---|---|---|---|---|---|")
for gamma, r in itertools.product((2.0, 4.0, 6.0), (1.5, 2.5)):
    well = np.hypot(I - ctr, J - ctr) <= r
    A = L2 - gamma * np.diag(well.astype(float)); d = np.diag(A).copy()
    if np.any(d <= 0): continue
    ev = np.linalg.eigvalsh(A); nneg = int((ev < 0).sum())
    for wc, k in itertools.product((0.3, 0.6, 1.0), (1, 2, 4)):
        W = np.where(isF, 1.0, wc); T = np.eye(n) - (W / d)[:, None] * A
        K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
        S = H.copy(); V = H.copy()
        for _ in range(k - 1): V = K @ V; S = S + V
        Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
        lam, vec = np.linalg.eig(T); amp = np.argsort(-np.abs(lam))[:max(nneg, 1) * 2]
        for fit in ("const", "const+bound"):
            Bc = np.ones((nc, 1)) if fit == "const" else np.column_stack([np.ones(nc), vec[C][:, amp].real, vec[C][:, amp].imag])
            St = truncate_lin(sp.csr_matrix(S), Bc, 2 if fit == "const" else 3).toarray()
            Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
            I_Pi = np.eye(n) - Pt @ np.linalg.solve(Pt.T @ A @ Pt, Pt.T @ A); I_PiQ = I_Pi[:, F]
            rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh)); rT = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
            best = (np.inf, None)
            if rN < 0.8:
                for w, rr, kk in itertools.product((0.1, 0.2, 0.4, 0.7, 1.0), (0.0, 0.1, 0.3, 1.0), (1, 2, 4, 8, 16)):
                    Wm = np.where(isF, w, rr * w); Tm = np.eye(n) - (Wm / d)[:, None] * A
                    rm = rho(I_Pi @ np.linalg.matrix_power(Tm, kk))
                    if rm < best[0]: best = (rm, (w, rr, kk))
            print(f"| {gamma} | {r} | {nneg} | {wc} | {k} | {fit} | {rN:.3g} | {rT:.3g} | {best[0]:.3g} {best[1]} |"); sys.stdout.flush()
