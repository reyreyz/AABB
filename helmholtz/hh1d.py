"""1D Helmholtz: pseudomodes of damped Jacobi and sparsity-constrained interpolation.

Pure NumPy (no scipy). Everything is dense and exact: small grids only.

A = tridiag(-1, 2, -1) - (kh)^2 I        (scaled by h^2, Dirichlet, n interior points)
S = I - (omega/d) A,  d = 2 - (kh)^2     (damped Jacobi; constant diagonal => S = poly(A), symmetric)
coarse points: fine indices 1, 3, ..., n-2 (0-based), n_c = (n-1)/2
column j of P is supported on S_j = {c_j - r, ..., c_j + r}, width w = 2r+1
two-grid: E = S^nu (I - P (P^T A P)^{-1} P^T A) S^nu, exact coarse solve, rho = max|eig(E)|
"""
import numpy as np

OMEGA = 2.0 / 3.0


# ---------------------------------------------------------------- operators
def helmholtz(n, kh):
    A = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1) - kh**2 * np.eye(n)
    d = 2 - kh**2
    S = np.eye(n) - (OMEGA / d) * A
    return A, S


def supports(n, w):
    r = (w - 1) // 2
    cs = np.arange(1, n - 1, 2)
    return cs, [np.arange(max(c - r, 0), min(c + r, n - 1) + 1) for c in cs]


def assemble(n, sups, cols):
    P = np.zeros((n, len(sups)))
    for j, (s, v) in enumerate(zip(sups, cols)):
        P[s, j] = v
    return P


def unpack(x, sups):
    out, k = [], 0
    for s in sups:
        out.append(x[k:k + len(s)]); k += len(s)
    return out


def twogrid(A, S, P, nu=1):
    Ac = P.T @ A @ P
    K = np.eye(len(A)) - P @ np.linalg.solve(Ac, P.T @ A)
    Snu = np.linalg.matrix_power(S, nu)
    return Snu @ K @ Snu


def rho(E):
    return np.max(np.abs(np.linalg.eigvals(E)))


def tg_rho(A, S, P, nu=1):
    try:
        Ac = P.T @ A @ P
        if np.linalg.cond(Ac) > 1e13:
            return np.inf
        return rho(twogrid(A, S, P, nu))
    except np.linalg.LinAlgError:
        return np.inf


# ---------------------------------------------------------------- interpolations
def hat(n, w_hat=3):
    cs, sups = supports(n, w_hat)
    r = (w_hat - 1) // 2
    cols = [1 - np.abs(s - c) / (r + 1) for c, s in zip(cs, sups)]
    return assemble(n, sups, cols)


def local_pseudomode(S, s_idx, zgrid):
    """min over z (real) and unit p with supp p in s_idx of ||(zI-S)p||."""
    n = len(S)
    best = (np.inf, None, None)
    for z in zgrid:
        B = (z * np.eye(n) - S)[:, s_idx]
        _, sv, Vt = np.linalg.svd(B, full_matrices=False)
        if sv[-1] < best[0]:
            best = (sv[-1], z, Vt[-1])
    # golden refine around best z
    z0 = best[1]; dz = zgrid[1] - zgrid[0]
    a, b = z0 - dz, z0 + dz
    f = lambda z: np.linalg.svd((z * np.eye(n) - S)[:, s_idx], compute_uv=False)[-1]
    g = (np.sqrt(5) - 1) / 2
    for _ in range(60):
        c, e = b - g * (b - a), a + g * (b - a)
        if f(c) < f(e): b = e
        else: a = c
    z = (a + b) / 2
    _, sv, Vt = np.linalg.svd((z * np.eye(n) - S)[:, s_idx], full_matrices=False)
    v = Vt[-1]
    v = v * np.sign(v[np.argmax(np.abs(v))])
    return sv[-1], z, v


def constrained_min(n, sups, cs, M, B):
    """Exact solution of  min sum_j p_j^T M p_j  s.t. supp p_j in S_j,  P B_c = B
    (B_c = B at coarse points). Null-space method; returns P and the constraint
    residual (0 if feasible, else the least-squares-feasible optimum)."""
    nv = sum(len(s) for s in sups)
    r = B.shape[1]
    Bc = B[cs]
    C = np.zeros((n * r, nv)); f = np.zeros(n * r)
    H = np.zeros((nv, nv))
    k = 0
    for j, s in enumerate(sups):
        H[k:k + len(s), k:k + len(s)] = M[np.ix_(s, s)]
        for t in range(r):
            C[t * n + s, k + np.arange(len(s))] = Bc[j, t]
        k += len(s)
    f = B.T.reshape(-1)
    U, sv, Vt = np.linalg.svd(C)
    tol = sv[0] * 1e-10
    rk = int(np.sum(sv > tol))
    x0 = Vt[:rk].T @ ((U[:, :rk].T @ f) / sv[:rk])
    res = np.linalg.norm(C @ x0 - f) / np.linalg.norm(f)
    N = Vt[rk:].T
    if N.shape[1]:
        y = np.linalg.solve(N.T @ H @ N, -N.T @ H @ x0)
        x = x0 + N @ y
    else:
        x = x0
    return assemble(n, sups, unpack(x, sups)), res


def spectral_opt(S, nc, nu=1):
    s = np.sort(np.abs(np.linalg.eigvalsh(S)))[::-1]
    return s[nc] ** (2 * nu)


# ---------------------------------------------------------------- true optimum
def _bfgs(f, x, iters=300, h=1e-6):
    def grad(x):
        g = np.zeros_like(x); fx = f(x)
        for i in range(len(x)):
            e = np.zeros_like(x); e[i] = h
            g[i] = (f(x + e) - f(x - e)) / (2 * h)
        return g
    fx, g = f(x), grad(x)
    Hinv = np.eye(len(x))
    for _ in range(iters):
        p = -Hinv @ g
        if g @ p >= 0:
            Hinv = np.eye(len(x)); p = -g
        t = 1.0
        while t > 1e-12:
            xn = x + t * p; fn = f(xn)
            if fn <= fx + 1e-4 * t * (g @ p): break
            t *= 0.5
        if t <= 1e-12: break
        gn = grad(xn); s_, y_ = xn - x, gn - g
        if y_ @ s_ > 1e-14:
            rho_ = 1 / (y_ @ s_); I = np.eye(len(x))
            Hinv = (I - rho_ * np.outer(s_, y_)) @ Hinv @ (I - rho_ * np.outer(y_, s_)) + rho_ * np.outer(s_, s_)
        if abs(fx - fn) < 1e-13 * max(1, abs(fx)): x, fx = xn, fn; break
        x, fx, g = xn, fn, gn
    return x, fx


def true_opt(A, S, sups, starts, nu=1, ms=(4, 16, 64)):
    """Minimise the two-grid spectral radius over all P with the given sparsity.
    Smooth surrogate ||E^m||_F^(1/m) with growing m, then report rho.
    Multi-start; returns best (rho, P)."""
    n = len(A)
    best = (np.inf, None)
    for x in starts:
        for m in ms:
            def f(x, m=m):
                P = assemble(n, sups, unpack(x, sups))
                Ac = P.T @ A @ P
                if np.linalg.cond(Ac) > 1e12: return 1e3
                E = twogrid(A, S, P, nu)
                return np.log(np.linalg.norm(np.linalg.matrix_power(E, m)) + 1e-300) / m
            x, _ = _bfgs(f, x.copy())
        P = assemble(n, sups, unpack(x, sups))
        r = tg_rho(A, S, P, nu)
        if r < best[0]: best = (r, P)
    return best


def pack(P, sups):
    return np.concatenate([P[s, j] for j, s in enumerate(sups)])


def invariant_scan(A, S, n, w, nu=1, grid=None):
    """Global grid search over translation-invariant symmetric stencils
    [.., b, a, 1, a, b, ..] (truncated at boundaries)."""
    cs, sups = supports(n, w)
    r = (w - 1) // 2
    grid = np.linspace(-2, 2, 81) if grid is None else grid
    import itertools
    best = (np.inf, None)
    for par in itertools.product(grid, repeat=r):
        st = np.concatenate([par[::-1], [1.0], par])
        cols = [st[s - c + r] for c, s in zip(cs, sups)]
        P = assemble(n, sups, cols)
        v = tg_rho(A, S, P, nu)
        if v < best[0]: best = (v, np.array(par))
    return best


def constrained_min_masked(n, sups, cs, M, B, mask):
    """As constrained_min, but constraint (row i, vector t) is imposed only where mask[i, t]."""
    nv = sum(len(s) for s in sups)
    r = B.shape[1]
    Bc = B[cs]
    rows, f = [], []
    H = np.zeros((nv, nv))
    k = 0; offs = []
    for j, s in enumerate(sups):
        H[k:k + len(s), k:k + len(s)] = M[np.ix_(s, s)]; offs.append(k); k += len(s)
    for t in range(r):
        for i in range(n):
            if not mask[i, t]: continue
            row = np.zeros(nv)
            for j, s in enumerate(sups):
                pos = np.where(s == i)[0]
                if len(pos): row[offs[j] + pos[0]] = Bc[j, t]
            rows.append(row); f.append(B[i, t])
    C = np.array(rows); f = np.array(f)
    U, sv, Vt = np.linalg.svd(C)
    rk = int(np.sum(sv > sv[0] * 1e-10))
    x0 = Vt[:rk].T @ ((U[:, :rk].T @ f) / sv[:rk])
    res = np.linalg.norm(C @ x0 - f) / np.linalg.norm(f)
    N = Vt[rk:].T
    x = x0 + N @ np.linalg.solve(N.T @ H @ N, -N.T @ H @ x0) if N.shape[1] else x0
    return assemble(n, sups, unpack(x, sups)), res


def waves_bc(n, kh, r):
    """Local Helmholtz near-kernel with Dirichlet walls: sin(th0 (i+1)) vanishes at the left ghost
    point, sin(th0 (n-i)) at the right one; each is imposed except within r+1 rows of the wall it violates."""
    th0 = np.arccos(1 - kh**2 / 2); i = np.arange(n)
    if th0 == 0:   # Poisson limit of the two one-sided sines: linear functions
        B = np.c_[i + 1.0, n - i + 0.0]
    else:
        B = np.c_[np.sin(th0 * (i + 1)), np.sin(th0 * (n - i))]
    mask = np.ones((n, 2), bool)
    mask[n - (r + 1):, 0] = False
    mask[:r + 1, 1] = False
    return B, mask


def true_opt_norm(A, S, sups, starts, nu=1, ps=(8, 32, 128)):
    """Minimise ||E||_2 over P with the given sparsity (Schatten-p surrogate, growing p)."""
    n = len(A)
    best = (np.inf, None)
    for x in starts:
        for p in ps:
            def f(x, p=p):
                P = assemble(n, sups, unpack(x, sups))
                if np.linalg.cond(P.T @ A @ P) > 1e12: return 1e3
                sv = np.linalg.svd(twogrid(A, S, P, nu), compute_uv=False)
                return np.log(sv[0]) + np.log(np.sum((sv / sv[0]) ** p)) / p
            x, _ = _bfgs(f, x.copy())
        P = assemble(n, sups, unpack(x, sups))
        try:
            v = np.linalg.norm(twogrid(A, S, P, nu), 2)
        except np.linalg.LinAlgError:
            v = np.inf
        if v < best[0]: best = (v, P)
    return best
