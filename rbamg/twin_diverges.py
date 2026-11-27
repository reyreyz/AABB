"""Twin diverges, N converges: convergent smoother + heavily truncated P~ (row sums kept).
k -> inf: E_twin -> (I - Pi~) Q (S_inf - S~) P^sharp (fixed);  E_N -> that times P^sharp T^k -> 0."""
import sys, numpy as np
from ml_twins import ML, poisson
from robust_n import recirc

def rate(h):
    h = np.asarray(h)
    if len(h) < 4: return h[-1] ** (1 / (len(h) - 1))
    j = max(1, len(h) - 6)
    return (h[-1] / h[j]) ** (1 / (len(h) - 1 - j))

def run(ml, b, wh, maxcyc=40):
    A = ml.levels[0].A; x = np.zeros(A.shape[0]); nb = np.linalg.norm(b); h = [1.0]
    for _ in range(maxcyc):
        x = ml.cycle(0, b, x, wh); h.append(np.linalg.norm(b - A @ x) / nb)
        if h[-1] < 1e-10 or h[-1] > 1e6 or not np.isfinite(h[-1]): break
    return rate(h), len(h) - 1, h[-1]

import scipy.sparse as sp
def upwind(m, pe, flow="recirc"):
    h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t)
    if flow == "recirc": bx, by = (2 * Y * (1 - X ** 2)).ravel(), (-2 * X * (1 - Y ** 2)).ravel()
    else: bx, by = np.cos(0.4) * np.ones(m * m), np.sin(0.4) * np.ones(m * m)
    eps = np.hypot(bx, by).max() * h / pe; n = m * m; I, J = np.divmod(np.arange(n), m)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [4 * eps + h * (np.abs(bx) + np.abs(by))]
    for di, dj, b in ((0, -1, np.maximum(bx, 0)), (0, 1, np.maximum(-bx, 0)), (-1, 0, np.maximum(by, 0)), (1, 0, np.maximum(-by, 0))):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        rows.append(np.nonzero(ok)[0]); cols.append(((I + di) * m + J + dj)[ok]); vals.append((-eps - h * b)[ok])
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n))
def _main():
    mode, m = sys.argv[1], int(sys.argv[2])
    sel = sys.argv[3] if len(sys.argv) > 3 else "base"
    probs = [("recirc Pe 2", recirc(m, 2.0)), ("Poisson", poisson(m))] if sel in ("base", "drop") else \
            [(f"upwind {fl} Pe {pe}", upwind(m, pe, fl)) for fl in ("recirc", "const") for pe in (20, 100)]
    W = 0.8 if sel in ("base", "drop") else 1.0
    if sel == "drop": probs.append(("upwind recirc Pe 20", upwind(m, 20)))
    TM = "drop" if sel == "drop" else "lump"
    TAUS = (0.3, 0.6, 1.0) if sel == "drop" else (0.5, 1.0)
    print(f"\n{mode}, m={m}\n\n| problem | trunc | k | P~ nnz/nnz(A) | N rate | twin rate | Markov-Jacobi rate |\n|---|---|---|---|---|---|---|")
    for pn, A in probs:
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for tau in TAUS:
            for k in ((2, 4, 8) if sel == 'base' else (2, 4, 8, 16) if sel == 'upwind' else (2, 4, 8, 16, 32)):
                kw = dict(max_levels=2, tau_filter=0) if mode == "two-level" else dict(k1=8, tau_deep=3e-4)
                ml = ML(A, omega=W, k0=k, tauP=tau, mode=TM, **kw)
                pc = ml.opcx()[1]
                r = [run(ml, b, wh) for wh in ("N", "twin", "markovJ")]
                f = lambda t: f"{t[0]:.3g}" + (" (div.)" if t[0] > 1 else "")
                print(f"| {pn} | {tau:g} {TM} | {k} | {pc:.2f} | {f(r[0])} | {f(r[1])} | {f(r[2])} |"); sys.stdout.flush()

if __name__ == "__main__":
    _main()
