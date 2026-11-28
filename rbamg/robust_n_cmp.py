"""Example 3b: robust N on the recirculating flow (Pe 20) vs Markovian cycles on the SAME hierarchy (robust_n.py):
  twin : k_l F-sweeps, then x += P' A_c^{-1} P'^T r      (N's Markovian twin with Q0)
  MJ   : k_l Jacobi sweeps (the smoother N uses), then the same correction
  MJF  : k_l Jacobi sweeps + k_l F-sweeps, then the same correction."""
import sys, time, numpy as np
from robust_n import recirc, RobustN

class Cmp(RobustN):
    def cyc(self, lev, b, x0, work, kind):
        L = self.levels[lev]
        if hasattr(L, "lu"):
            work[0] += L.lu_nnz / self.nnz0; return L.lu.solve(b)
        x = x0.copy()
        if kind in ("MJ", "MJF", "MJF1"):
            for _ in range(L.k):
                x = x + (b - L.A @ x) / L.D
            work[0] += L.k * L.A.nnz / self.nnz0
        if kind in ("twin", "MJF"):
            x = self.f_sweeps(L, b, x, L.k, work)
        if kind == "MJF1":
            x = self.f_sweeps(L, b, x, 1, work)
        r = b - L.A @ x; rc = r[L.C] + L.ST @ r[L.F]
        e = self.cyc(lev + 1, rc, np.zeros(len(L.C)), work, kind)
        if lev <= self.refine and not hasattr(self.levels[lev + 1], "lu"):
            e = e + self.cyc(lev + 1, rc - self.RAP(L, e), np.zeros(len(L.C)), work, kind)
            work[0] += (L.APc.nnz + L.APf.nnz + L.S.nnz) / self.nnz0
        x[L.C] += e; x[L.F] += L.S @ e
        work[0] += (L.A.nnz + 2 * L.S.nnz) / self.nnz0
        return x

m = int(sys.argv[1]); A = recirc(m); b = np.random.default_rng(0).standard_normal(A.shape[0]); nb = np.linalg.norm(b)
t0 = time.time(); ml = Cmp(A); ts = time.time() - t0
op, st = ml.complexity()
for kind in (sys.argv[2].split(",") if len(sys.argv) > 2 else ("N", "twin", "MJF", "MJF1", "MJ")):
    x, work, hist = np.zeros(A.shape[0]), [0.0], [1.0]
    for c in range(40):
        x = ml.cycle(0, b, x, work) if kind == "N" else ml.cyc(0, b, x, work, kind)
        hist.append(np.linalg.norm(b - A @ x) / nb)
        if hist[-1] < 1e-8 or not np.isfinite(hist[-1]) or hist[-1] > 1e8: break
    nc = len(hist) - 1; conv = hist[-1] < 1e-8
    f = (hist[-1] / hist[1]) ** (1 / max(nc - 1, 1)) if conv else hist[-1] ** (1 / nc)
    print(f"| {m} | {A.shape[0]} | {len(ml.levels)} | {kind} | {f:.3g}{'' if conv else ' (diverges)' if hist[-1] > 1 else ' (no conv.)'} | "
          f"{nc if conv else '-'} | {work[0]/nc:.0f} | {work[0] if conv else float('inf'):.0f} | {op:.2f} | {ts:.0f} |", flush=True)
