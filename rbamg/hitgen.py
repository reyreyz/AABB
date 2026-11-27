"""Problem and level-0 setup of falsification hit 1 (seed 12, sample 635), at any m."""
import numpy as np, scipy.sparse as sp
H1 = dict(a=[-0.98489481, -0.00789497, -1.06081079, -0.78232303, -0.36845338, 0.96149543], f1=0.9550906842488585,
          f2=1.3439989079159567, pe=20.525384441816787, delta=0.3946463251653395, wf=0.51, wc=0.51, k=1, theta=0.5213219631577246)
def problem(m, P=H1, pe=None):
    h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); a = P["a"]
    bx = (a[0] + a[1] * Y + a[2] * np.sin(np.pi * (X + Y) * P["f1"])).ravel()
    by = (a[3] + a[4] * X + a[5] * np.cos(np.pi * X * P["f2"])).ravel()
    eps = np.hypot(bx, by).max() * h / (pe or P["pe"]); n = m * m; I, J = np.divmod(np.arange(n), m)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [4 * eps * np.ones(n)]
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        rows.append(np.nonzero(ok)[0]); cols.append(((I + di) * m + J + dj)[ok]); vals.append((-eps + s * 0.5 * h * b)[ok])
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n))
    return A, np.column_stack([X.ravel(), Y.ravel()])
