import sys, numpy as np
from ml_twins import ML, factor
from robust_n import recirc
m, pe = int(sys.argv[1]), float(sys.argv[2])
A = recirc(m, pe); b = np.random.default_rng(0).standard_normal(A.shape[0])
for tau in (0.03, 0.1, 0.2, 0.3):
    ml = ML(A, omega=0.8, tauP=tau, mode="lump", k1=8, tau_deep=3e-4)
    oc, pc = ml.opcx(); out = []
    for wh in ("N", "twin", "markovJ"):
        h = ml.solve(b, wh); conv = h[-1] < 1e-8
        out.append(f"{factor(h):.3g} / {len(h)-1 if conv else '-'} / {ml.wpc*(len(h)-1) if conv else float('inf'):.0f}")
    print(f"| {pe:g} | {m} | {tau:g} lump | {pc:.2f} | {oc:.2f} | " + " | ".join(out) + " |"); sys.stdout.flush()
