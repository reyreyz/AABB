"""Multilevel natural examples for any flow.
A-type: damped Jacobi w (inside its stability range) and k sweeps on all levels, split on level 0, one-point or
        lump-0.5 S~ on level 0, S truncated at 1e-2 below, filter 1e-2.
B-type: Jacobi w and k sweeps on all levels, greedy split, S truncated at tau on all levels, filter 1e-3."""
import sys, time, numpy as np
from ml_twins import ML
from flows import cd
m, flow, pe, case, w, k = int(sys.argv[1]), sys.argv[2], float(sys.argv[3]), sys.argv[4], float(sys.argv[5]), int(sys.argv[6])
opt = sys.argv[7] if len(sys.argv) > 7 else ""
methods = sys.argv[8].split(",") if len(sys.argv) > 8 else ["N", "twin", "markovJ"]
A, xy = cd(m, flow, pe); n = A.shape[0]; I, J = np.divmod(np.arange(n), m)
if case == "A":
    split, tr = opt.split("/") if opt else ("greedy", "onepoint")
    kw = dict(omega=w, k0=k, k1=k, delta=0.3, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0, lin_levels=1)
    if split == "rb": kw["split0"] = (I + J) % 2 == 1
    kw.update(dict(mode="lin", xy=xy, nkeep=1) if tr == "onepoint" else dict(mode="lump", tauP=0.5, trunc_from=0))
    if tr != "onepoint": kw.pop("tauP_deep"); kw.pop("lin_levels")
else:
    tau = float(opt) if opt else 0.1
    kw = dict(omega=w, k0=k, k1=k, delta=0.3, mode="drop", tauP=tau, tau_filter=1e-3, tau_deep=1e-3, refine=0)
for extra in (sys.argv[9].split(",") if len(sys.argv) > 9 else []):
    key, val = extra.split("="); kw[key] = float(val) if "." in val or "e" in val else int(val)
b = np.random.default_rng(0).standard_normal(n); t0 = time.time(); ml = ML(A, **kw); ts = time.time() - t0
for wh in methods:
    hist = ml.solve(b, wh, tol=1e-8, maxcyc=80); nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {flow} Pe {pe:g} {case} (w {w}, k {k}, {opt} {sys.argv[9] if len(sys.argv) > 9 else ""}) | {m} | {n} | {len(ml.levels)} | {wh} | {f:.3g}{'' if conv else ' (no conv.)'} | "
          f"{nc if conv else '-'} | {ml.wpc:.0f} | {ml.wpc*nc if conv else float('inf'):.0f} | {ml.opcx()[0]:.2f} | {ts:.0f} |"); sys.stdout.flush()
