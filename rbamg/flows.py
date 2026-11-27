"""Centered CD on (-1,1)^2 with a general velocity field; cell Peclet defined with max |beta| (as recirc)."""
import numpy as np, scipy.sparse as sp
FLOWS = {
    "recirc": lambda X, Y: (2 * Y * (1 - X**2), -2 * X * (1 - Y**2)),
    "poiseuille": lambda X, Y: (1 - Y**2, 0 * X),
    "hyperbolic": lambda X, Y: (X, -Y),
    "const": lambda X, Y: (np.cos(0.4) + 0 * X, np.sin(0.4) + 0 * Y),
    "half-vortex": lambda X, Y: (2 * (Y + 1) * (1 - X**2) / 2, -2 * X * (1 - ((Y + 1) / 2 * 2 - 1)**2) * 0 - X * (1 - Y**2)),
    "inflow-bend": lambda X, Y: ((1 + Y) * (1 - X**2) * 0 + (1 - Y**2), -X * (1 - Y**2) * 0 + 0.5 * (1 - X**2) * np.sign(X + 1e-12) * 0 + 0.3 * (1 + X)),
}
def cd(m, flow, pe):
    h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t)
    bx, by = FLOWS[flow](X, Y); bx, by = np.asarray(bx).ravel(), np.asarray(by).ravel()
    eps = np.hypot(bx, by).max() * h / pe; n = m * m; I, J = np.divmod(np.arange(n), m)
    rows, cols, vals = [np.arange(n)], [np.arange(n)], [4 * eps * np.ones(n)]
    for di, dj, b, s in ((0, -1, bx, -1), (0, 1, bx, 1), (-1, 0, by, -1), (1, 0, by, 1)):
        ok = (I + di >= 0) & (I + di < m) & (J + dj >= 0) & (J + dj < m)
        rows.append(np.nonzero(ok)[0]); cols.append(((I + di) * m + J + dj)[ok]); vals.append((-eps + s * 0.5 * h * b)[ok])
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), (n, n)), np.column_stack([X.ravel(), Y.ravel()])
