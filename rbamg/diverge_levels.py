import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML, factor
from robust_n import recirc
m = int(sys.argv[1])
A = recirc(m); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
ml = ML(A, omega=1.0, k1=8, tau_deep=3e-4)
nl = len(ml.levels) - 1
print(f"m={m}, {nl} smoothing levels + LU\n\n| level | n | k | rho(T) | growth of C-error over k sweeps (random e) |\n|---|---|---|---|---|")
rng = np.random.default_rng(1)
for l, L in enumerate(ml.levels[:-1]):
    Tm = sla.LinearOperator(L.A.shape, matvec=lambda v, L=L: v - (L.A @ v) / L.D)
    rho = abs(sla.eigs(Tm, k=1, which="LM", return_eigenvectors=False, maxiter=3000, tol=1e-3)[0])
    gs = []
    for _ in range(5):
        e = rng.standard_normal(L.A.shape[0])
        for _ in range(5): e = e - (L.A @ e) / L.D          # warm: dominant modes
        e /= np.linalg.norm(e); ek = e.copy()
        for _ in range(L.k): ek = ek - (L.A @ ek) / L.D
        gs.append(np.linalg.norm(ek[L.C]) / np.linalg.norm(e[L.C]))
    print(f"| {l} | {L.A.shape[0]} | {L.k} | {rho:.2f} | {np.median(gs):.3g} |")
print("\n| shift j=k applied on | factor | cycles |\n|---|---|---|")
for lab, sh in (("no level (Alg 7.3, j=0 everywhere)", [0.0] * nl),
                ("level 0 only", [1.0] + [0.0] * (nl - 1)),
                ("levels >= 1 only", [0.0] + [1.0] * (nl - 1)),
                ("all levels", [1.0] * nl)):
    ml.shift = sh
    h = ml.solve(b, "N")
    print(f"| {lab} | {factor(h):.3g} | {len(h)-1 if h[-1] < 1e-8 else 'no conv.'} |"); sys.stdout.flush()
