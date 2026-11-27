"""Verify the missing example with sparse iterative runs: recirc centered Pe 10, Jacobi w=1, delta=0.3,
k=1, single-entry lumped P~ (tau=1).  Competitors on the same P~: twin (k F-sweeps), Markov-Jacobi with
any omega and k, and their best case."""
import sys, numpy as np
from ml_twins import ML
from robust_n import recirc
from twin_diverges import run
m, ml_mode = int(sys.argv[1]), sys.argv[2]
A = recirc(m, 10.0); b = np.random.default_rng(0).standard_normal(A.shape[0])
kw = dict(max_levels=2, tau_filter=0) if ml_mode == "two-level" else dict(k1=1, tau_deep=3e-4)
fmt = lambda t: f"{t[0]:.3g}" + (" div." if t[2] > 1 else f" / {t[1]}" if t[2] < 1e-10 else " no conv.")
ml = ML(A, omega=1.0, k0=1, delta=0.3, tauP=1.0, mode="lump", **kw)
print(f"| {m} | {ml_mode} | {len(ml.levels)} | {ml.opcx()[1]:.2f} | N k=1: {fmt(run(ml, b, 'N', 60))} | "
      f"twin k=1: {fmt(run(ml, b, 'twin', 60))} | MJ k=1 w=1: {fmt(run(ml, b, 'markovJ', 60))} |"); sys.stdout.flush()
best = []
for w in (0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0):
    for k in (1, 2, 4, 8):
        mm = ML(A, omega=w, k0=k, k1=k if ml_mode != "two-level" else 1, delta=0.3, tauP=1.0, mode="lump", **{x: y for x, y in kw.items() if x != "k1"})
        for wh in ("markovJ", "twin"):
            best.append((run(mm, b, wh, 40)[0], wh, w, k))
best.sort()
print("   best Markovian competitors on the same truncation (rate, method, omega, k):", [(round(r, 3), wh, w, k) for r, wh, w, k in best[:4]])
