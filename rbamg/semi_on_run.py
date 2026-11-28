import sys, time, numpy as np, scipy.sparse as sp
import sparse_problems as spp
from flows import cd as cdflow
from semi_on import S72
m, prob = int(sys.argv[1]), sys.argv[2]; cfgs = sys.argv[3].split(";"); markov = len(sys.argv) > 4 and sys.argv[4] == "M"
name, pe = prob.split(":"); pe = float(pe)
A = sp.csr_matrix(spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, pe, "2", False), "centered") if name == "shear" else cdflow(m, "const", pe)[0])
b = np.random.default_rng(0).standard_normal(A.shape[0]); nb = np.linalg.norm(b)
for cfg in cfgs:
    kw = {}
    for item in cfg.split(","):
        key, val = item.split("="); kw[key] = (float(val) if ("." in val or "e" in val) else int(val))
    if kw.get("coarsest", 0) < 0: kw["coarsest"] = int(-kw["coarsest"] * np.sqrt(A.shape[0]))
    t0 = time.time(); ml = S72(A, markov=markov, **kw); ts = time.time() - t0
    x, work, hist = np.zeros(A.shape[0]), [0.0], [1.0]
    for c in range(60):
        x = ml.cycle(b, x, work); hist.append(np.linalg.norm(b - A @ x) / nb)
        if hist[-1] < 1e-8 or not np.isfinite(hist[-1]) or hist[-1] > 1e8: break
    nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {prob} | {m} | {'M ' if markov else 'S '}{cfg} | {len(ml.levels())} | {f:.3g}{'' if conv else ' (no conv.)'} | {nc if conv else '-'} | "
          f"{work[0]/nc:.0f} | {work[0] if conv else float('inf'):.0f} | {ml.op_complexity():.2f} | {ts:.0f} |"); sys.stdout.flush()
