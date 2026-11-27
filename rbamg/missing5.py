"""Missing example, checks: iterative runs, Markovian competitor sweep on the same P~ family, multilevel."""
import sys, numpy as np, itertools
from ml_twins import ML
from robust_n import recirc
from twin_diverges import run
m, what = int(sys.argv[1]), sys.argv[2]
A = recirc(m, 10.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
fmt = lambda t: f"{t[0]:.3g}" + (" div." if t[2] > 1 else f" / {t[1]} cyc" if t[2] < 1e-10 else " no conv.")
if what == "runs":
    for lab, kw in (("two-level", dict(max_levels=2, tau_filter=0)), ("multilevel k1=1", dict(k1=1, tau_deep=3e-4)),
                    ("multilevel k1=2", dict(k1=2, tau_deep=3e-4))):
        ml = ML(A, omega=1.0, omega_c=0.8, k0=1, delta=0.3, mode="lin", xy=xy, nkeep=2, **kw)
        print(f"| {m} | {lab} ({len(ml.levels)} levels) | N {fmt(run(ml, b, 'N', 80))} | twin {fmt(run(ml, b, 'twin', 80))} | MJ {fmt(run(ml, b, 'markovJ', 80))} |"); sys.stdout.flush()
else:
    res = []
    for w, k in itertools.product((0.2, 0.4, 0.6, 0.8, 1.0), (1, 2, 4, 8)):
        for wc_rel in (1.0, 0.5, 0.0):
            ml = ML(A, omega=w, omega_c=w * wc_rel, k0=k, delta=0.3, mode="lin", xy=xy, nkeep=2, max_levels=2, tau_filter=0)
            res.append((run(ml, b, "markovJ", 30)[0], "MJ", w, w * wc_rel, k))
            if wc_rel == 1.0: res.append((run(ml, b, "twin", 30)[0], "twin", w, "-", k))
            if wc_rel == 1.0 and k <= 2: res.append((run(ml, b, "N", 30)[0], "N", w, w, k))
    res.sort(key=lambda r: r[0])
    print(f"m={m}: best 8 of {len(res)} (rate, method, omega_F, omega_C, k):")
    for r in res[:8]: print("  ", (round(r[0], 3),) + r[1:])
    print("   best Markovian:", next((round(r[0], 3),) + r[1:] for r in res if r[1] != "N"))
