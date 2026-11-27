"""Idea 21: constant-flow centered CD (red-black) + localized well A - gamma*diag(A)*chi_disk.  Dense two-level:
N with S~ fitted to constants + bound-state C-parts, twin, best damped MJ (incl. w_F=1 with small w_C)."""
import numpy as np, scipy.sparse as sp, itertools, sys
from ml_twins import truncate_lin
from flows import cd
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 23
n = m * m; I, J = np.divmod(np.arange(n), m); ctr = (m - 1) / 2
isF = (I + J) % 2 == 1; F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
rows = []
print(f"m={m}\n| Pe | gamma | r | #Re<0 | w_C | k | nk | p | N | twin | best damped MJ (wF, wC, k) |\n|---|---|---|---|---|---|---|---|---|---|---|")
for pe, gamma, r in itertools.product((6.0, 10.0), (0.7, 0.85), (1.5, 2.5)):
    A0, _ = cd(m, "const", pe); A0 = A0.toarray(); well = np.hypot(I - ctr, J - ctr) <= r
    A = A0 - gamma * np.diag(np.diag(A0) * well); d = np.diag(A)
    ev = np.linalg.eigvals(A); nneg = int((ev.real < 0).sum())
    for wc, k in itertools.product((0.05, 0.2, 0.5), (1, 2, 4, 8)):
        W = np.where(isF, 1.0, wc); T = np.eye(n) - (W / d)[:, None] * A
        K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
        S = H.copy(); V = H.copy()
        for _ in range(k - 1): V = K @ V; S = S + V
        Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
        lam, vec = np.linalg.eig(T); amp = np.argsort(-np.abs(lam))
        for nk, p in itertools.product((1, 2, 3), (0, 2, 6)):
            Bc = np.ones((nc, 1)) if p == 0 else np.column_stack([np.ones(nc), vec[C][:, amp[:p]].real, vec[C][:, amp[:p]].imag])
            St = truncate_lin(sp.csr_matrix(S), Bc, nk).toarray()
            Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
            Ac = Pt.T @ A @ Pt
            if np.linalg.cond(Ac) > 1e12: continue
            I_Pi = np.eye(n) - Pt @ np.linalg.solve(Ac, Pt.T @ A); I_PiQ = I_Pi[:, F]
            rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh))
            if rN > 0.75: continue
            rT = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
            if rT < 1: continue
            best = (np.inf, None)
            for wf, wcc, kk in itertools.product((0.1, 0.3, 0.6, 1.0), (0.0, 0.01, 0.03, 0.1, 0.3, 1.0), (1, 2, 4, 8, 16)):
                Wm = np.where(isF, wf, wcc * wf if wf < 1 else wcc); Tm = np.eye(n) - (Wm / d)[:, None] * A
                rm = rho(I_Pi @ np.linalg.matrix_power(Tm, kk))
                if rm < best[0]: best = (rm, (wf, wcc, kk))
            print(f"| {pe:g} | {gamma} | {r} | {nneg} | {wc} | {k} | {nk} | {p} | {rN:.3g} | {rT:.3g} | {best[0]:.3g} {best[1]} |"); sys.stdout.flush()
