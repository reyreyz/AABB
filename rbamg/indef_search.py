"""Idea 17: indefinite operators, red-black, S~ fitted to the C-parts of the amplified modes of T.
Dense two-level: rho(E_N), rho(E_twin), and the best damped Markov-Jacobi over (w, w_C ratio, k <= 16)."""
import numpy as np, scipy.sparse as sp, itertools, sys
from ml_twins import truncate_lin
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 15
T1 = sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)); L2 = (sp.kron(sp.eye(m), T1) + sp.kron(T1, sp.eye(m))).toarray()
n = m * m; I, J = np.divmod(np.arange(n), m); isF = (I + J) % 2 == 1; F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]
nf, nc = len(F), len(C)
h = 1.0 / (m + 1); conv = np.zeros((n, n))
for i in range(n):                                            # small centered convection in x to break symmetry
    if J[i] > 0: conv[i, i - 1] = -0.5
    if J[i] < m - 1: conv[i, i + 1] = 0.5
def two_level(A, T, St, S, k):
    Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
    I_Pi = np.eye(n) - Pt @ np.linalg.solve(Pt.T @ A @ Pt, Pt.T @ A)
    return Pt, I_Pi
Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
allrows = []
print(f"m={m}\n| shift s | conv | w_C | k | test vecs p | entries | N | twin | best damped MJ (w, ratio, k) |\n|---|---|---|---|---|---|---|---|---|")
for s, cv in itertools.product((0.5, 1.0, 1.5, 2.5), (0.0, 0.5)):
    A = L2 - s * np.eye(n) + cv * conv; d = np.diag(A)
    for wc, k in itertools.product((0.3, 0.6, 1.0), (1, 2, 3)):
        W = np.where(isF, 1.0, wc); T = np.eye(n) - (W / d)[:, None] * A
        K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
        S = H.copy(); V = H.copy()
        for _ in range(k - 1): V = K @ V; S = S + V
        ev, evec = np.linalg.eig(T); amp = np.argsort(-np.abs(ev))
        Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
        for p, nk in itertools.product((2, 4, 8), (1, 2, 3)):
            V_ = evec[:, amp[:p]]; tv = np.column_stack([V_.real, V_.imag])[C]
            # fit S~ rows to reproduce S on tv (not to reproduce the modes' F-parts)
            Bc = np.column_stack([np.ones(nc), tv])
            # truncate_lin matches S on Bc columns using the nk largest entries of S per row
            St = truncate_lin(sp.csr_matrix(S), Bc, nk).toarray()
            Pt, I_Pi = two_level(A, T, St, S, k)
            I_PiQ = I_Pi[:, F]
            rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh))
            rT = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
            allrows.append((rN, rT, s, cv, wc, k, p, nk))
            if rN > 0.7 or rT < 1: continue
            best = (np.inf, None)
            for w, r, kk in itertools.product((0.1, 0.2, 0.4, 0.7, 1.0), (0.0, 0.3, 1.0), (1, 2, 4, 8, 16)):
                Wm = np.where(isF, w, r * w); Tm = np.eye(n) - (Wm / d)[:, None] * A
                rm = rho(I_Pi @ np.linalg.matrix_power(Tm, kk))
                if rm < best[0]: best = (rm, (w, r, kk))
            print(f"| {s} | {cv} | {wc} | {k} | {p} | {nk} | {rN:.3g} | {rT:.3g} | {best[0]:.3g} {best[1]} |"); sys.stdout.flush()

allrows.sort()
print("\nbest N overall (N, twin, s, conv, wC, k, p, nk):")
for r in allrows[:10]: print("  ", tuple(round(x, 3) if isinstance(x, float) else x for x in r))
print("best twin-divergent by N:")
for r in [r for r in allrows if r[1] > 1][:10]: print("  ", tuple(round(x, 3) if isinstance(x, float) else x for x in r))
