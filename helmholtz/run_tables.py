"""Driver: hat pseudomode stats, local optima, constrained optima, true two-grid optima."""
import sys, time, json
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hh1d import *

n = int(sys.argv[1]); kh = float(sys.argv[2]); widths = [int(w) for w in sys.argv[3].split(',')]
rng = np.random.default_rng(0)
A, S = helmholtz(n, kh)
lam, V = np.linalg.eigh(A)
s = 1 - OMEGA / (2 - kh**2) * lam
cs3, s3 = supports(n, 3)
nc = len(cs3)
out = {"n": n, "kh": kh, "k": kh * (n + 1), "spectral": spectral_opt(S, nc)}

def colstats(P):
    eps, zs = [], []
    for j in range(P.shape[1]):
        p = P[:, j] / np.linalg.norm(P[:, j]); z = p @ S @ p
        eps.append(np.linalg.norm(S @ p - z * p)); zs.append(z)
    return float(np.median(zs)), float(np.median(eps))

# hat as pseudomode (centre column)
Ph = hat(n)
p = Ph[:, nc // 2] / np.linalg.norm(Ph[:, nc // 2]); c = V.T @ p
out["hat_zp"] = float(p @ S @ p); out["hat_eps"] = float(np.linalg.norm(S @ p - (p @ S @ p) * p))
out["hat_w_amplified"] = float(np.sum(c[s > 1] ** 2))
out["hat_w_slow"] = float(np.sum(c[np.abs(s) > 0.5] ** 2))
inb = np.argsort(np.abs(lam))[0]
out["nk_wavenumber_index"] = int(inb)

rows = []
def add(name, w, P, extra=""):
    r = tg_rho(A, S, P)
    nE = float(np.linalg.norm(twogrid(A, S, P), 2)) if np.isfinite(r) else np.inf
    zm, em = colstats(P)
    rows.append(dict(name=name, w=w, rho=r, normE=nE, condP=float(np.linalg.cond(P)), z_med=zm, eps_med=em, extra=extra))
    print(f"{name:30s} w={w} rho={r:9.4f} |E|2={nE:9.3f} cond(P)={np.linalg.cond(P):9.2f} z_med={zm:.3f} eps_med={em:.4f} {extra}", flush=True)

add("linear (hat)", 3, Ph)
B_nk1 = V[:, [inb]]
B_lonk = V[:, [0, inb]]
B_nk2 = V[:, np.argsort(np.abs(lam))[:2]]
zslow = np.linspace(0.5, s.max() + 0.05, 120)
for w in widths:
    cs, sups = supports(n, w)
    # individually optimal local pseudomodes
    for tag, zg in (("local opt, z free>=0.5", zslow), ("local opt, z=1", None)):
        cols = []
        for sj in sups:
            if zg is None:
                _, _, Vt = np.linalg.svd((np.eye(n) - S)[:, sj], full_matrices=False); cols.append(Vt[-1])
            else:
                cols.append(local_pseudomode(S, sj, zg)[2])
        add(tag, w, assemble(n, sups, cols))
    # constrained, exact
    starts = [pack(assemble(n, sups, [np.interp(sj, [c - 2, c, c + 2], [0, 1, 0]) for c, sj in zip(cs, sups)]), sups)]  # linear hat, padded
    for bname, B in (("nk1", B_nk1), ("lo+nk", B_lonk), ("nk2", B_nk2)):
        P, res = constrained_min(n, sups, cs, A @ A, B)
        add(f"constr M=A^2 B={bname}", w, P, f"constraint resid={res:.1e}")
        starts.append(pack(P, sups))
    if kh == 0:
        P, res = constrained_min(n, sups, cs, A, B_nk1)
        add("constr M=A B=nk1", w, P, f"constraint resid={res:.1e}")
        starts.append(pack(P, sups))
    # translation-invariant global scan
    t = time.time()
    grid = {1: np.linspace(-3, 3, 1201), 2: np.linspace(-2, 2, 161), 3: np.linspace(-1.5, 1.5, 41)}[(w - 1) // 2]
    rinv, par = invariant_scan(A, S, n, w, grid=grid)
    st = np.concatenate([par[::-1], [1.0], par]); r_ = (w - 1) // 2
    Pinv = assemble(n, sups, [st[sj - c + r_] for c, sj in zip(cs, sups)])
    add("invariant stencil, global scan", w, Pinv, f"stencil={np.round(st, 4).tolist()} ({time.time()-t:.0f}s)")
    starts.append(pack(Pinv, sups))
    for _ in range(3):
        starts.append(starts[0] + 0.3 * rng.standard_normal(len(starts[0])))
    t = time.time()
    ropt, Popt = true_opt(A, S, sups, starts)
    add("TRUE OPT (all P, this support)", w, Popt, f"({len(starts)} starts, {time.time()-t:.0f}s)")
    np.save(f"Popt_n{n}_kh{kh}_w{w}.npy", Popt)
print(f"eigenvector coarse space (dense P; optimal only for SPD) rho = {out['spectral']:.4f}")
out["rows"] = rows
json.dump(out, open(f"res_n{n}_kh{kh}.json", "w"), indent=1)
