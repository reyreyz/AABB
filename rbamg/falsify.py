"""Randomized falsification search (runs until killed): look for rho(E_N) < 0.7 while every simple Markovian cycle on
the same T, Q, P~ has rho >= 0.95: twin (k F-sweeps), MJ (k T-sweeps), MJF_m (k T-sweeps + m F-sweeps, m <= 2k).
Two-level, exact coarse solve, dense, small grids, random flows/splits/weights/truncations."""
import sys, time, numpy as np, scipy.sparse as sp
from robust_n import greedy_f_points
from ml_twins import truncate_lin
from robust_n import drop_rel
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
out = open(sys.argv[2] if len(sys.argv) > 2 else "results_falsify.md", "a")
def problem():
    m = int(rng.integers(9, 15)); h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t)
    a = rng.normal(size=6)
    bx = a[0] + a[1] * Y + a[2] * np.sin(np.pi * (X + Y) * rng.uniform(0.5, 2)); by = a[3] + a[4] * X + a[5] * np.cos(np.pi * X * rng.uniform(0.5, 2))
    bx, by = bx.ravel(), by.ravel(); pe = rng.uniform(2, 25); eps = np.hypot(bx, by).max() * h / pe
    n = m * m; I, J = np.divmod(np.arange(n), m); A = np.zeros((n, n)); A[np.arange(n), np.arange(n)] = 4 * eps
    shift = rng.choice([0.0, 0.0, rng.uniform(0, 3)]) * eps
    A -= shift * np.eye(n)
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        A[np.nonzero(ok)[0], ((I + di) * m + J + dj)[ok]] = (-eps + s * 0.5 * h * b)[ok]
    return m, A, np.column_stack([X.ravel(), Y.ravel()]), pe, shift / eps
count, hits, t0 = 0, 0, time.time()
CC = 3.0
def wpd(r, sweeps):
    return (sweeps + CC) / (-np.log10(r)) if r < 1 else np.inf
while True:
    m, A, xy, pe, sh = problem(); n = len(A); d = np.diag(A)
    if np.any(d <= 0): continue
    I, J = np.divmod(np.arange(n), m)
    sk = rng.integers(0, 3)
    isF = ((I + J) % 2 == 1) if sk == 0 else greedy_f_points(sp.csr_matrix(A), rng.uniform(0.1, 0.9)) if sk == 1 else (rng.random(n) < 0.6)
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    if nc < 4 or nf < 4: continue
    wf, wc, k = rng.uniform(0.5, 1.0), rng.uniform(0, 1.0), int(rng.integers(1, 7))
    W = np.where(isF, wf, wc); T = np.eye(n) - (W / d)[:, None] * A
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
    if rho(K) >= 1: continue
    S = H.copy(); V = H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    tk = rng.integers(0, 3)
    if tk == 0: St = truncate_lin(sp.csr_matrix(S), np.column_stack([np.ones(nc), xy[C]]) if rng.random() < 0.5 else np.ones((nc, 1)), int(rng.integers(1, 4))).toarray()
    elif tk == 1: St = drop_rel(sp.csr_matrix(S), rng.uniform(0.05, 0.9)).toarray()
    else:
        th = rng.uniform(0, 1); B1 = np.ones((nc, 1)); St = ((1 - th) * truncate_lin(sp.csr_matrix(S), B1, 1) + th * truncate_lin(sp.csr_matrix(S), B1, 2)).toarray()
    Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
    Ac = Pt.T @ A @ Pt
    if np.linalg.cond(Ac) > 1e10: continue
    I_Pi = np.eye(n) - Pt @ np.linalg.solve(Ac, Pt.T @ A); I_PiQ = I_Pi[:, F]
    Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
    Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
    rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh)); count += 1
    if rN < 0.7:
        rT = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
        if rT >= 0.95:
            rM = rho(I_Pi @ Tk)
            TQ = np.eye(n) - Qsh.T @ Qsh @ ((np.where(isF, wf, 0.0) / d)[:, None] * A)
            best = min(rho(I_Pi @ np.linalg.matrix_power(TQ, mm) @ Tk) for mm in range(1, 2 * k + 1))
            if rM >= 0.95 and best >= 0.95:
                # extended Markovian family on the same P~ (two-level, exact coarse solve)
                ext = []
                FW = lambda wF_: np.eye(n) - Qsh.T @ Qsh @ ((np.where(isF, wF_, 0.0) / d)[:, None] * A)
                for fwt in (wf, 1.0):
                    TQw = FW(fwt); TQp = [np.eye(n)]
                    for _ in range(12): TQp.append(TQw @ TQp[-1])
                    Tp = np.eye(n)
                    for kk in range(0, 2 * k + 3):
                        if kk > 0: Tp = T @ Tp
                        for mm in range(0, 13):
                            if kk + mm == 0: continue
                            r_ = rho(I_Pi @ TQp[mm] @ Tp)
                            ext.append((wpd(r_, kk + mm * nf / n), r_, f"MJF T^{kk} F^{mm} (Fw {fwt:.2f})"))
                for w_, r_, kk in [(w_, r_, kk) for w_ in (0.1, 0.2, 0.4, 0.7, 1.0) for r_ in (0.0, 0.1, 0.3, 1.0) for kk in (1, 2, 4, 8, 16)]:
                    Tm = np.eye(n) - (np.where(isF, w_, r_ * w_) / d)[:, None] * A; Tmk = np.linalg.matrix_power(Tm, kk)
                    for mm in (0, 1, 3):
                        rr = rho(I_Pi @ np.linalg.matrix_power(FW(1.0), mm) @ Tmk)
                        ext.append((wpd(rr, kk + mm * nf / n), rr, f"dampedMJ w {w_} wC {r_*w_:.2f} k {kk} F^{mm}"))
                for w_ in (0.3, 0.5, 0.7, 1.0):
                    TQd = FW(w_)
                    for kk in (1, 2, 4, 8, 16):
                        rr = rho(I_Pi @ np.linalg.matrix_power(TQd, kk)); ext.append((wpd(rr, kk * nf / n), rr, f"twin wF {w_} k {kk}"))
                bestext = min(ext); wN = wpd(rN, k + 1 + 0.5 * k)
                if bestext[0] < 2 * wN:
                    out.write(f"near (best competitor {bestext[2]}: rho {bestext[1]:.3g}, work/digit {bestext[0]:.1f} vs N {wN:.1f}): N={rN:.3g} m={m} Pe={pe:.1f} split={sk} wF={wf:.2f} wC={wc:.2f} k={k}\n"); out.flush()
                    continue
                hits += 1
                out.write(f"STRONG HIT (N work/digit {wN:.1f}; best Markovian {bestext[2]}: rho {bestext[1]:.3g}, work/digit {bestext[0]:.1f})  count={count}: m={m} Pe={pe:.2f} shift={sh:.2f} split={sk} wF={wf:.2f} wC={wc:.2f} k={k} trunc={tk} "
                          f"rho(T)={rho(T):.2f} rho(K)={rho(K):.2f} N={rN:.3g} twin={rT:.3g} MJ={rM:.3g} bestMJF={best:.3g}\n"); out.flush()
    if count % 500 == 0:
        out.write(f"progress: {count} samples, {hits} hits, {time.time()-t0:.0f}s\n"); out.flush()
