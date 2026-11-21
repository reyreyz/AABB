"""Pseudospectra of the two-grid operator E for four interpolations (n = 31, damped Jacobi 2/3):
linear hat, rho-optimal (w=3 search result), ||E||_2-optimal (w=3), eigenvector coarse space (dense).
Norms: Euclidean, and |A|-norm (|A| = (A^2)^(1/2); equals the energy norm for Poisson).
Outputs: psE.json (grids for plotting), psE.txt (pseudospectral radii, transient norms, Kreiss bounds)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hh1d import *

n = 31
EPS = [1e-1, 1e-2, 1e-3, 1e-4]
out = {}
lines = []
def say(s): print(s, flush=True); lines.append(s)

def smin_grid(E, X, Y):
    I = np.eye(len(E)); G = np.empty((len(Y), len(X)))
    for a, y in enumerate(Y):
        for b, x in enumerate(X):
            G[a, b] = np.linalg.svd((x + 1j * y) * I - E, compute_uv=False)[-1]
    return G

def ps_radius(E, eps_list, nphi=181, nr=500):
    """alpha_eps = max{|z| : sigma_min(zI-E) <= eps} on a polar grid (E real: upper half suffices).
    Lambda_eps lies in |z| <= ||E|| + eps, so the grid covers it."""
    I = np.eye(len(E)); R = np.linalg.norm(E, 2) + max(eps_list)
    phis = np.linspace(0, np.pi, nphi); rs = np.linspace(0, R, nr)
    S = np.array([[np.linalg.svd(r * np.exp(1j * p) * I - E, compute_uv=False)[-1] for r in rs] for p in phis])
    ev = np.linalg.eigvals(E)
    f = lambda z: np.linalg.svd(z * I - E, compute_uv=False)[-1]
    res = []
    for e in eps_list:
        m = S <= e
        best = float(rs[np.where(m.any(0))[0].max()]) if m.any() else 0.0
        # march outward from every eigenvalue along its own ray and along +-90 deg rotations of it
        for lam in ev:
            if abs(lam) < best - 2 * R / nr: continue
            for d in (lam / abs(lam) if abs(lam) > 1e-14 else 1.0,):
                t0, step = 0.0, max(e, 1e-6) * 1e-2
                while f(lam + (t0 + step) * d) <= e: t0 += step; step *= 1.5
                lo, hi = t0, t0 + step
                for _ in range(50):
                    mid = (lo + hi) / 2
                    if f(lam + mid * d) <= e: lo = mid
                    else: hi = mid
                best = max(best, abs(lam + lo * d))
        res.append(best)
    return res

for kh in (0.0, 0.3, 0.6, 0.9):
    A, S = helmholtz(n, kh)
    lam, V = np.linalg.eigh(A)
    absA_h = V @ np.diag(np.abs(lam) ** 0.5) @ V.T
    absA_mh = V @ np.diag(np.abs(lam) ** -0.5) @ V.T
    sl, sv = np.linalg.eigh(S); idx = np.argsort(-np.abs(sl))[:(n - 1) // 2]
    Ps = {"linear": hat(n),
          "rho-opt (w=3)": np.load(f"Popt_n{n}_kh{kh}_w3.npy"),
          "norm-opt (w=3)": np.load(f"Pnorm_n{n}_kh{kh}_w=3.npy"),
          "eigvec (dense)": sv[:, idx]}
    say(f"\n=== kh = {kh} ===")
    say(f"{'P':16s} {'norm':5s} {'rho':>6s} {'||E||':>7s} " + " ".join(f"a(1e{int(np.log10(e))})" for e in EPS)
        + "   ||E^m|| m=1,2,5,10,20          cyc(1e-8) rho-pred  dep(E)")
    out[str(kh)] = {}
    for name, P in Ps.items():
        E2 = twogrid(A, S, P)
        for nm, E in (("2", E2), ("|A|", absA_h @ E2 @ absA_mh)):
            ev = np.linalg.eigvals(E); r = np.max(np.abs(ev))
            a = ps_radius(E, EPS)
            pw = [np.linalg.norm(np.linalg.matrix_power(E, m), 2) for m in range(0, 41)]
            hit = [m for m in range(1, 41) if pw[m] <= 1e-8]
            cyc = str(hit[0]) if hit else (">40" if r < 1 else "div.")
            pred = f"{np.log(1e-8) / np.log(r):.1f}" if r < 1 else "div."
            # Henrici departure from normality, relative
            dep = np.sqrt(max(np.linalg.norm(E, 'fro')**2 - np.sum(np.abs(ev)**2), 0)) / np.linalg.norm(E, 'fro')
            say(f"{name:16s} {nm:5s} {r:6.3f} {pw[1]:7.3f} " + " ".join(f"{x:7.3f}" for x in a)
                + "   " + " ".join(f"{pw[m]:8.2g}" for m in (1, 2, 5, 10, 20)) + f"   {cyc:>6s}   {pred:>6s}   {dep:.2f}")
            if nm == "|A|" or True:
                R = max(1.3, 1.15 * (a[0] if np.isfinite(a[0]) else r))
                X = np.linspace(-R, R, 161); Y = np.linspace(0, R, 81)
                G = smin_grid(E, X, Y)
                out[str(kh)].setdefault(name, {})[nm] = dict(R=R, X=X.tolist(), Y=Y.tolist(), logs=np.log10(G).tolist(),
                                                           ev=[[float(z.real), float(z.imag)] for z in ev], rho=float(r),
                                                           normE=float(pw[1]), alpha=a, powers=[float(p) for p in pw])
json.dump(out, open("psE.json", "w"))
open("psE.txt", "w").write("\n".join(lines) + "\n")
