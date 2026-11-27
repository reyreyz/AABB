"""Two-grid convergence: N (= (I-Pi) Q0 K^k Q0^T) vs Markovian M (= (I-Pi) T^k), same T, same P', R = P'^T.
Sweep problems x splits x smoothers x k; count who converges (rho < 1)."""
import numpy as np, scipy.sparse as sp, itertools, sys
from robust_n import recirc, greedy_f_points

def cd(m, pe, scheme, flow):
    h = 1.0 / (m + 1); n = m * m
    t = h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t)
    if flow == "const": bx, by = np.cos(0.4) * np.ones(n), np.sin(0.4) * np.ones(n)
    else: bx, by = (2 * (Y - .5) * (1 - (2 * X - 1) ** 2)).ravel(), (-2 * (X - .5) * (1 - (2 * Y - 1) ** 2)).ravel()
    eps = np.hypot(bx, by).max() * h / pe
    I, J = np.divmod(np.arange(n), m); A = sp.lil_matrix((n, n))
    for i in range(n):
        A[i, i] = 4 * eps
        for di, dj, b in ((0, -1, -bx[i]), (0, 1, bx[i]), (-1, 0, -by[i]), (1, 0, by[i])):
            ii, jj = I[i] + di, J[i] + dj
            if 0 <= ii < m and 0 <= jj < m:
                c = -eps + 0.5 * h * b if scheme == "centered" else -eps + h * max(b, 0) * (1 if b > 0 else 0)
                if scheme == "upwind":
                    c = -eps - h * max(-b, 0)
                    A[i, i] += h * abs(b) / 2 * 0 + 0
                A[i, ii * m + jj] = c
        if scheme == "upwind": A[i, i] = -A[i].sum() + A[i, i] + (4 * eps - (-A[i].sum() + A[i, i])) * 0
    A = A.tocsr()
    if scheme == "upwind":       # make rows of interior points sum like eps-Laplacian + upwind convection
        A = A + sp.diags(-np.asarray(A.sum(1)).ravel() + np.asarray(A.sum(1)).ravel() * 0)
        A = A.tolil()
        for i in range(n):
            A[i, i] = 4 * eps + h * (abs(bx[i]) + abs(by[i]))
        A = A.tocsr()
    return A

def splits(A, m):
    I, J = np.divmod(np.arange(m * m), m)
    out = {f"greedy d={d}": greedy_f_points(A, d) for d in (0.3, 0.6, 0.9)}
    out["red-black (C = black)"] = (I + J) % 2 == 0
    out["full coarsening (C = even,even)"] = ~((I % 2 == 1) & (J % 2 == 1))
    out["aggressive (C = every 3rd, 3rd)"] = ~((I % 3 == 1) & (J % 3 == 1))
    return out

def rates(A, isF, smoother, w, k):
    A = A.toarray(); n = len(A); D = np.diag(np.diag(A))
    Minv = (w * np.linalg.inv(D)) if smoother == "Jacobi" else np.linalg.inv(np.tril(A))
    T = np.eye(n) - Minv @ A
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
    S, V = H.copy(), H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    Pp = np.zeros((n, len(C))); Pp[C, np.arange(len(C))] = 1; Pp[F] = S
    Ac = Pp.T @ A @ Pp
    Pi = Pp @ np.linalg.solve(Ac, Pp.T @ A)
    QK = np.zeros((n, n)); QK[np.ix_(F, F)] = np.linalg.matrix_power(K, k)
    EN, EM = (np.eye(n) - Pi) @ QK, (np.eye(n) - Pi) @ np.linalg.matrix_power(T, k)
    r = lambda E: np.abs(np.linalg.eigvals(E)).max()
    return r(EN), r(EM), np.abs(np.linalg.eigvals(T)).max(), np.abs(np.linalg.eigvals(K)).max() if len(F) else 0

m = 15
probs = {"Poisson": (sp.kron(sp.eye(m), sp.diags([-1, 2, -1], [-1, 0, 1], (m, m))) + sp.kron(sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)), sp.eye(m))).tocsr()}
probs["aniso 0.01"] = (sp.kron(sp.eye(m), sp.diags([-1, 2, -1], [-1, 0, 1], (m, m))) + 0.01 * sp.kron(sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)), sp.eye(m))).tocsr()
for pe in (2, 20, 100):
    probs[f"recirc centered Pe {pe}"] = recirc(m, pe)
    probs[f"const-flow upwind Pe {pe}"] = cd(m, pe, "upwind", "const")
rows, tally = [], {}
for (pn, A), (sm, w) in itertools.product(probs.items(), (("Jacobi", 2 / 3), ("Jacobi", 1.0), ("GS", 1.0))):
    for sn, isF in splits(A, m).items():
        for k in (1, 2, 4):
            rN, rM, rT, rK = rates(A, isF, sm, w, k)
            cls = ("both" if rN < 1 and rM < 1 else "only N" if rN < 1 else "only M" if rM < 1 else "neither")
            tally[cls] = tally.get(cls, 0) + 1
            rows.append((pn, f"{sm} w={w:.2g}", sn, k, rT, rK, rN, rM, cls))
print(f"m={m}, {len(rows)} cases\n\n| outcome | cases |\n|---|---|")
for c in ("both", "only N", "only M", "neither"): print(f"| {c} | {tally.get(c, 0)} |")
both = [r for r in rows if r[-1] == "both"]
print(f"\nWhen both converge: rho(E_N) <= rho(E_M) in {sum(r[6] <= r[7] * (1 + 1e-9) for r in both)} of {len(both)}; "
      f"median rho(E_N)/rho(E_M) = {np.median([r[6] / max(r[7], 1e-300) for r in both]):.3g}")
print("\n## Cases where M converges and N does not, or M is faster by more than 1.5x\n")
print("| problem | smoother | split | k | rho(T) | rho(K) | rho(E_N) | rho(E_M) |\n|---|---|---|---|---|---|---|---|")
for r in rows:
    if r[-1] == "only M" or (r[-1] == "both" and r[6] > 1.5 * r[7] and r[6] > 1e-3):
        print(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]:.3g} | {r[5]:.3g} | {r[6]:.3g} | {r[7]:.3g} |")
print("\n## Full table\n\n| problem | smoother | split | k | rho(T) | rho(K) | rho(E_N) | rho(E_M) | outcome |\n|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]:.3g} | {r[5]:.3g} | {r[6]:.3g} | {r[7]:.3g} | {r[8]} |")
