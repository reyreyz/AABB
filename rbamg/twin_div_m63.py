import sys, numpy as np
from ml_twins import ML
from robust_n import recirc
from twin_diverges import run
m = int(sys.argv[1])
A = recirc(m, 10.0); b = np.random.default_rng(0).standard_normal(A.shape[0])
print(f"\nm={m}, sparse two-level runs (exact coarse LU), Pe 10, Jacobi w=0.1, tau=1 lump\n")
print("| delta | k | N: rate / cycles | twin | Markov-Jacobi |\n|---|---|---|---|---|")
for delta in (0.3, 0.6):
    for k in (16, 64):
        ml = ML(A, omega=0.1, k0=k, delta=delta, tauP=1.0, mode="lump", max_levels=2, tau_filter=0)
        out = []
        for wh in ("N", "twin", "markovJ"):
            r, c, last = run(ml, b, wh, maxcyc=60)
            out.append(f"{r:.3g} / {c}" + (" (div.)" if last > 1 else "" if last < 1e-10 else " (no conv.)"))
        print(f"| {delta} | {k} | " + " | ".join(out) + " |"); sys.stdout.flush()
