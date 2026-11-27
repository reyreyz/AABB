"""Randomized falsification search with a persistence test (runs until killed).
Strong hit at size m (two-level, dense): rho(N) < 0.7; twin, MJ at N's weights, MJF_m (m <= 2k) all >= 0.95; and N's
work per digit <= 1/2 of every extended Markovian variant (MJF T^a F^b, damped MJ + F-sweeps, damped twins).
Persistent hit: the same flow/parameters rebuilt at m2 = 2m+1 is again a strong hit."""
import sys, time, numpy as np, scipy.sparse as sp
from robust_n import greedy_f_points, drop_rel
from ml_twins import truncate_lin
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
seed = int(sys.argv[1]); rng = np.random.default_rng(seed); out = open(sys.argv[2], "a"); CC = 3.0
wpd = lambda r, sw: (sw + CC) / (-np.log10(r)) if r < 1 else np.inf

def draw():
    return dict(m=int(rng.integers(9, 15)), a=rng.normal(size=6), f1=rng.uniform(0.5, 2), f2=rng.uniform(0.5, 2),
                pe=rng.uniform(2, 25), shift=float(rng.choice([0.0, 0.0, rng.uniform(0, 3)])), sk=int(rng.integers(0, 2)),
                delta=rng.uniform(0.1, 0.9), wf=rng.uniform(0.5, 1.0), wc=rng.uniform(0, 1.0), k=int(rng.integers(1, 7)),
                tk=int(rng.integers(0, 3)), nk=int(rng.integers(1, 4)), lin=bool(rng.random() < 0.5), tau=rng.uniform(0.05, 0.9),
                theta=rng.uniform(0, 1))

def build(P, m):
    h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); a = P["a"]
    bx = (a[0] + a[1] * Y + a[2] * np.sin(np.pi * (X + Y) * P["f1"])).ravel(); by = (a[3] + a[4] * X + a[5] * np.cos(np.pi * X * P["f2"])).ravel()
    eps = np.hypot(bx, by).max() * h / P["pe"]; n = m * m; I, J = np.divmod(np.arange(n), m)
    A = np.zeros((n, n)); A[np.arange(n), np.arange(n)] = 4 * eps - P["shift"] * eps
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        A[np.nonzero(ok)[0], ((I + di) * m + J + dj)[ok]] = (-eps + s * 0.5 * h * b)[ok]
    isF = ((I + J) % 2 == 1) if P["sk"] == 0 else greedy_f_points(sp.csr_matrix(A), P["delta"])
    return A, isF, np.column_stack([X.ravel(), Y.ravel()])

def evaluate(P, m, full=True):
    A, isF, xy = build(P, m); n = len(A); d = np.diag(A)
    if np.any(d <= 0): return None
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    if nc < 4 or nf < 4: return None
    wf, wc, k = P["wf"], P["wc"], P["k"]
    T = np.eye(n) - (np.where(isF, wf, wc) / d)[:, None] * A; K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
    if rho(K) >= 1: return None
    S = H.copy(); V = H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    B = np.column_stack([np.ones(nc), xy[C]]) if P["lin"] else np.ones((nc, 1)); B1 = np.ones((nc, 1))
    St = (truncate_lin(sp.csr_matrix(S), B, P["nk"]) if P["tk"] == 0 else drop_rel(sp.csr_matrix(S), P["tau"]) if P["tk"] == 1 else
          ((1 - P["theta"]) * truncate_lin(sp.csr_matrix(S), B1, 1) + P["theta"] * truncate_lin(sp.csr_matrix(S), B1, 2))).toarray()
    Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St; Ac = Pt.T @ A @ Pt
    if np.linalg.cond(Ac) > 1e10: return None
    I_Pi = np.eye(n) - Pt @ np.linalg.solve(Ac, Pt.T @ A); I_PiQ = I_Pi[:, F]
    Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
    Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
    rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh)); res = dict(N=rN, rhoT=rho(T), rhoK=rho(K))
    if rN >= 0.7: return res
    res["twin"] = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh)); res["MJ"] = rho(I_Pi @ Tk)
    if res["twin"] < 0.95 or res["MJ"] < 0.95 or not full: return res
    FW = lambda w_: np.eye(n) - Qsh.T @ Qsh @ ((np.where(isF, w_, 0.0) / d)[:, None] * A)
    ext = []
    for fwt in (wf, 1.0):
        TQw = FW(fwt); TQp = [np.eye(n)]
        for _ in range(12): TQp.append(TQw @ TQp[-1])
        Tp = np.eye(n)
        for kk in range(0, 2 * k + 3):
            if kk > 0: Tp = T @ Tp
            for mm in range(0, 13):
                if kk + mm: r_ = rho(I_Pi @ TQp[mm] @ Tp); ext.append((wpd(r_, kk + mm * nf / n), r_, f"MJF T^{kk} F^{mm} Fw {fwt:.2f}"))
    for w_ in (0.1, 0.2, 0.4, 0.7, 1.0):
        for r_ in (0.0, 0.1, 0.3, 1.0):
            Tm = np.eye(n) - (np.where(isF, w_, r_ * w_) / d)[:, None] * A
            for kk in (1, 2, 4, 8, 16):
                Tmk = np.linalg.matrix_power(Tm, kk)
                for mm in (0, 1, 3):
                    rr = rho(I_Pi @ np.linalg.matrix_power(FW(1.0), mm) @ Tmk); ext.append((wpd(rr, kk + mm * nf / n), rr, f"dMJ w {w_} wC {r_*w_:.2f} k {kk} F^{mm}"))
    for w_ in (0.3, 0.5, 0.7, 1.0):
        for kk in (1, 2, 4, 8, 16):
            rr = rho(I_Pi @ np.linalg.matrix_power(FW(w_), kk)); ext.append((wpd(rr, kk * nf / n), rr, f"twin wF {w_} k {kk}"))
    res["best"] = min(ext); res["wN"] = wpd(rN, k + 1 + 0.5 * k); res["strong"] = res["best"][0] >= 2 * res["wN"]
    return res

count, hits, t0 = 0, 0, time.time()
while True:
    P = draw(); r = evaluate(P, P["m"]); count += 1
    if r and r.get("strong"):
        m2 = 2 * P["m"] + 1; r2 = evaluate(P, m2)
        tag = "PERSISTENT HIT" if (r2 and r2.get("strong")) else "transient hit"
        out.write(f"{tag} seed={seed} count={count}: m={P['m']} -> {m2}: N {r['N']:.3g} -> {r2['N'] if r2 else float('nan'):.3g}; "
                  f"twin {r['twin']:.3g}, MJ {r['MJ']:.3g}; best Markovian {r['best'][2]} rho {r['best'][1]:.3g} wpd {r['best'][0]:.1f} vs N {r['wN']:.1f}; "
                  f"P={ {k_: (round(float(v), 4) if np.isscalar(v) else [round(float(x), 4) for x in v]) for k_, v in P.items()} }\n"); out.flush()
        if tag == "PERSISTENT HIT": np.save(f"/tmp/claude-0/-home-user-AABB/72d75a69-cd12-59c4-ba81-ce7a8de908cb/scratchpad/persistent_{seed}_{count}.npy", P, allow_pickle=True)
    if count % 500 == 0: out.write(f"progress: {count} samples, {time.time()-t0:.0f}s\n"); out.flush()
