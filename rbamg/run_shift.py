import sys, numpy as np
from ml_twins import ML, factor
from robust_n import recirc
m, pe, w = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
taus = [float(t) for t in sys.argv[4].split(",")]
A = recirc(m, pe); b = np.random.default_rng(0).standard_normal(A.shape[0])
for tau in taus:
    out = []
    for wh, sh in (("twin", 0), ("N", 0.0), ("N", 0.5), ("N", 1.0)):
        ml = ML(A, omega=w, tauP=tau, mode="lump", k1=8, tau_deep=3e-4, shift=sh)
        h = ml.solve(b, wh); conv = h[-1] < 1e-8
        f = factor(h)
        out.append(f"{'div.' if f == float('inf') else f'{f:.3g}'} / {len(h)-1 if conv else '-'} / {ml.wpc*(len(h)-1) if conv else float('inf'):.0f}")
    print(f"| {pe:g} | {m} | {'none' if tau == 0 else f'{tau:g} lump'} | " + " | ".join(out) + " |"); sys.stdout.flush()
