import sys, numpy as np
from ml_twins import ML, factor
from robust_n import recirc
m = int(sys.argv[1]); pe = float(sys.argv[2])
A = recirc(m, pe); b = np.random.default_rng(0).standard_normal(A.shape[0])
for tau, mode in ((0, "drop"), (0.3, "lump"), (0.5, "lump")):
    ml = ML(A, omega=0.8, tauP=tau, mode=mode)
    oc, pc = ml.opcx()
    out = []
    for wh in ("N", "twin", "markovJ"):
        h = ml.solve(b, wh); conv = h[-1] < 1e-8
        out.append(f"{factor(h):.3g} / {len(h)-1 if conv else '-'} / {ml.wpc*(len(h)-1) if conv else float('inf'):.0f}")
    print(f"| {pe:g} | {m} | {'none' if tau == 0 else f'{tau:g} lump'} | {pc:.2f} | {oc:.2f} | " + " | ".join(out) + " |")
    sys.stdout.flush()
