"""Rebuild every construction (fast, exact) + load the optimised P's; print a markdown table of
two-grid rho and ||E||_2 (n = 31, damped Jacobi omega = 2/3, nu = 1 + 1, exact coarse solve)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hh1d import *

n = 31
KH = (0.0, 0.3, 0.6, 0.9)
W = (3, 5, 7)

def metr(A, S, P):
    r = tg_rho(A, S, P)
    if not np.isfinite(r): return "sing."
    return f"{r:.3f} / {np.linalg.norm(twogrid(A, S, P), 2):.2f}"

def build(kh, w, A, S):
    lam, V = np.linalg.eigh(A); order = np.argsort(np.abs(lam))
    cs, sups = supports(n, w); r = (w - 1) // 2
    out = {}
    if w == 3: out["linear (hat)"] = hat(n)
    zs = np.linspace(0.5, np.linalg.eigvalsh(S).max() + 0.05, 120)
    out["local pseudomode opt (each column)"] = assemble(n, sups, [local_pseudomode(S, sj, zs)[2] for sj in sups])
    out["constr. min A^2, B = 1 near-kernel eigvec"] = constrained_min(n, sups, cs, A @ A, V[:, order[:1]])[0]
    out["constr. min A^2, B = 2 near-kernel eigvecs"] = constrained_min(n, sups, cs, A @ A, V[:, order[:2]])[0]
    B, mask = waves_bc(n, kh, r)
    out["constr. min A^2, B = one-sided sines (BC-aware)"] = constrained_min_masked(n, sups, cs, A @ A, B, mask)[0]
    f = f"Popt_n{n}_kh{kh}_w{w}.npy"
    if os.path.exists(f): out["**min rho, all P with this support**"] = np.load(f)
    f = f"Pnorm_n{n}_kh{kh}_w={w}.npy"
    if os.path.exists(f): out["**min ||E||2, all P with this support**"] = np.load(f)
    return out

names = None; table = {}
for kh in KH:
    A, S = helmholtz(n, kh)
    for w in W:
        for name, P in build(kh, w, A, S).items():
            table.setdefault((name, w), {})[kh] = metr(A, S, P)
print(f"| interpolation | width | " + " | ".join(f"kh={kh}" for kh in KH) + " |")
print("|---|---|" + "---|" * len(KH))
for w in W:
    for (name, ww), d in table.items():
        if ww != w: continue
        print(f"| {name} | {w} | " + " | ".join(d.get(kh, "–") for kh in KH) + " |")
for kh in KH:
    A, S = helmholtz(n, kh); nc = (n - 1) // 2
    lam, V = np.linalg.eigh(S); idx = np.argsort(-np.abs(lam))[:nc]
    print(f"eigenvector coarse space (dense), kh={kh}: {metr(A, S, V[:, idx])}")
