"""Multilevel M vs S: Algorithm 7.2 (semi-Markovian: memory in the interpolation only) vs its Markovian version
(same smoother, split, Rhat = [I, -A_cf (I+G+G^2) D_ff^{-1}], coarse operator Rhat A P0, P0 = injection)."""
import sys, time, numpy as np, scipy.sparse as sp
import ml_cycles as mc, sparse_problems as spp
from mlm_goodP import MLM
from flows import cd as cdflow
m, prob, w, k = int(sys.argv[1]), sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
def poisson(m, ey=1.0):
    T1 = sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)); return (sp.kron(sp.eye(m), T1) + ey * sp.kron(T1, sp.eye(m))).tocsr()
name, pe = (prob.split(":") + ["0"])[:2]; pe = float(pe)
A = {"poisson": lambda: poisson(m), "aniso": lambda: poisson(m, 0.01),
     "shear": lambda: spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, pe, "2", False), "centered"),
     "recirc": lambda: cdflow(m, "recirc", pe)[0], "const": lambda: cdflow(m, "const", pe)[0]}[name]()
A = sp.csr_matrix(A); b = np.random.default_rng(0).standard_normal(A.shape[0])
which = sys.argv[5].split(",") if len(sys.argv) > 5 else ["S", "M0"]
builders = {"S": ("S (Alg. 7.2, injection + memory)", lambda: mc.Alg72(A, k=k, omega=w)),
            "M0": ("M, injection P0", lambda: mc.Alg72(A, k=k, omega=w, markov=True)),
            "MW1": ("M, P = [I; -D_ff^-1 A_fc]", lambda: MLM(A, k=k, omega=w, wmode="W1")),
            "MPp": ("M, P = P' (S's interpolation, static)", lambda: MLM(A, k=k, omega=w, wmode="Pp"))}
for key in which:
    lab = builders[key][0]
    t0 = time.time(); ml = builders[key][1](); ts = time.time() - t0
    x, work, hist = np.zeros(A.shape[0]), [0.0], [1.0]; nb = np.linalg.norm(b)
    for c in range(80):
        x = ml.cycle(b, x, work); hist.append(np.linalg.norm(b - A @ x) / nb)
        if hist[-1] < 1e-8 or not np.isfinite(hist[-1]) or hist[-1] > 1e8: break
    nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {prob} | w {w}, k {k} | {m} | {lab} | {f:.3g}{'' if conv else (' (diverges)' if hist[-1] > 1 else ' (no conv. in 80)')} | "
          f"{nc if conv else '-'} | {work[0]/nc:.0f} | {work[0] if conv else float('inf'):.0f} | {ml.op_complexity():.2f} |"); sys.stdout.flush()
