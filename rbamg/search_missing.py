"""Search: twin diverges AND Markov-Jacobi diverges AND N converges (two-grid, dense, truncated P~)."""
import numpy as np, scipy.sparse as sp, itertools, sys
from robust_n import recirc, greedy_f_points
from ml_twins import truncate
from twin_diverges import upwind

def const_centered(m, pe):
    h = 2.0 / (m + 1); n = m * m; bx, by = np.cos(0.4), np.sin(0.4)
    eps = h / pe; I, J = np.divmod(np.arange(n), m)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [4 * eps * np.ones(n)]
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        rows.append(np.nonzero(ok)[0]); cols.append(((I + di) * m + J + dj)[ok]); vals.append((-eps + s * 0.5 * h * b) * np.ones(ok.sum()))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n))

rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 15
probs = {f"recirc centered Pe {p}": recirc(m, p) for p in (5, 10, 20, 50)}
probs.update({f"const centered Pe {p}": const_centered(m, p) for p in (5, 20)})
probs.update({f"upwind recirc Pe {p}": upwind(m, p) for p in (20, 100)})
smoothers = [("Jacobi", 0.5), ("Jacobi", 0.8), ("Jacobi", 1.0), ("GS", 1.0)]
truncs = [(0.3, "lump"), (0.6, "lump"), (1.0, "lump"), (0.3, "drop"), (1.0, "drop")]
hits, near = [], []
for (pn, A), (sm, w), delta in itertools.product(probs.items(), smoothers, (0.3, 0.6)):
    Ad = A.toarray(); n = len(Ad)
    Minv = w * np.diag(1 / np.diag(Ad)) if sm == "Jacobi" else np.linalg.inv(np.tril(Ad))
    T = np.eye(n) - Minv @ Ad
    isF = greedy_f_points(A, delta); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
    rT, rK = rho(T), rho(K)
    if rK >= 1: continue
    Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1
    Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
    S, V, Tk, Kk = H.copy(), H.copy(), T.copy(), K.copy()
    for k in range(1, 5):
        if k > 1: V = K @ V; S = S + V; Tk = T @ Tk; Kk = K @ Kk
        for tau, mode in truncs:
            St = truncate(sp.csr_matrix(S), tau, mode).toarray()
            Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
            I_Pi = np.eye(n) - Pt @ np.linalg.solve(Pt.T @ Ad @ Pt, Pt.T @ Ad)
            I_PiQ = I_Pi[:, F]
            rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh))
            rTw = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
            rM = rho(I_Pi @ Tk)
            row = (pn, f"{sm} w={w}", delta, k, f"{tau:g} {mode}", rT, rK, rN, rTw, rM)
            if rN < 1 and rTw > 1 and rM > 1: hits.append(row)
            elif rN < 1 and (rTw > 1 or rM > 1): near.append(row)
    sys.stdout.flush()
fmt = lambda r: f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} | {r[5]:.3g} | {r[6]:.3g} | {r[7]:.3g} | {r[8]:.3g} | {r[9]:.3g} |"
hdr = "| problem | smoother | delta | k | trunc | rho(T) | rho(K) | N | twin | Markov-Jacobi |\n|---|---|---|---|---|---|---|---|---|---|"
print(f"m={m}\n\n## Hits: N < 1, twin > 1, Markov-Jacobi > 1 ({len(hits)})\n\n{hdr}")
for r in sorted(hits, key=lambda r: r[7]): print(fmt(r))
print(f"\n## Near misses: N < 1 and one of the two > 1 ({len(near)}), best 30 by N\n\n{hdr}")
for r in sorted(near, key=lambda r: r[7])[:30]: print(fmt(r))
