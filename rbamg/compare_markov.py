"""Honest Markovian competitors for the robust N cycle (robust_n.py).

Same hierarchy for N and the Markovian cycles: greedy split, P' = [I; S],
R = P'^T, A_c = P'^T A P', truncation, filters, refinement on levels 0-1,
sparse LU below 16 sqrt(n).
  N        : non-Markovian cycle, k undamped Jacobi steps on all points (rho(T) ~ 4.5)
  M_J      : Markovian, k Jacobi steps then coarse correction with P' (same operators)
  M_F(a,b) : Markovian, a F-Jacobi sweeps, coarse correction with P', b F-sweeps
             (the Corollary 4.7 member: with F-relaxation the history vanishes)
External: PyAMG AIR and Ruge-Stueben (defaults), stand-alone and as GMRES preconditioners.
Work: WU = nnz(A) units; for PyAMG its cycle_complexity per cycle; GMRES orthogonalization
not counted for anyone.
"""

import sys
import time
import numpy as np
import scipy.sparse.linalg as sla
import pyamg
from robust_n import recirc, RobustN


def stat(hist):
    h = np.asarray(hist)
    n = len(h) - 1
    return (h[-1] / h[1]) ** (1 / max(n - 1, 1)) if n > 1 else h[-1], n


def markov_jacobi(ml, b, tol=1e-8, maxcyc=30):
    """M_J: k Jacobi steps on all points, then x += P' A_c^{-1} R r (recursive, same hierarchy)."""
    def cyc(lev, rhs, x0, work):
        L = ml.levels[lev]
        if hasattr(L, "lu"):
            work[0] += L.lu_nnz / ml.nnz0
            return L.lu.solve(rhs)
        x = x0.copy()
        for _ in range(L.k):
            x = x + (rhs - L.A @ x) / L.D
        r = rhs - L.A @ x
        rc = r[L.C] + L.ST @ r[L.F]
        e = cyc(lev + 1, rc, np.zeros(len(L.C)), work)
        if lev <= ml.refine and not hasattr(ml.levels[lev + 1], "lu"):
            e = e + cyc(lev + 1, rc - ml.RAP(L, e), np.zeros(len(L.C)), work)
        x[L.C] += e
        x[L.F] += L.S @ e
        work[0] += ((L.k + 1) * L.A.nnz + 2 * L.S.nnz) / ml.nnz0
        return x
    A = ml.levels[0].A
    x, nb, work, hist = np.zeros(A.shape[0]), np.linalg.norm(b), [0.0], [1.0]
    for _ in range(maxcyc):
        x = cyc(0, b, x, work)
        hist.append(np.linalg.norm(b - A @ x) / nb)
        if hist[-1] < tol or not np.isfinite(hist[-1]) or hist[-1] > 1e10:
            break
    return hist, work[0]


def gmres_with(precond, A, b, tol=1e-8, maxiter=200):
    """right-preconditioned GMRES (restart 50); returns iterations and convergence flag."""
    n = A.shape[0]
    count = [0]

    def apply(v):
        count[0] += 1
        return precond(v)
    M = sla.LinearOperator((n, n), matvec=apply)
    res = []
    x, info = sla.gmres(A, b, M=M, rtol=tol, restart=50, maxiter=maxiter,
                        callback=lambda r: res.append(r), callback_type="pr_norm")
    ok = np.linalg.norm(b - A @ x) / np.linalg.norm(b) < 10 * tol
    return count[0], ok


def row(m, lab, fac, cyc, wuc, wu, conv, extra=""):
    print(f"| {m} | {lab} | {fac:.3g} | {cyc}{'' if conv else ' (no conv.)'} | {wuc:.0f} | "
          f"{wu if conv else float('inf'):.0f} | {extra} |")
    sys.stdout.flush()


def main(ms=(127, 255, 511)):
    print("# Markovian competitors vs N (section 7.3 problem, Pe 20)\n")
    print("| m | method | factor | cycles / its to 1e-8 | WU per cycle | WU total | note |")
    print("|---|---|---|---|---|---|---|")
    for m in ms:
        A = recirc(m)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        ml = RobustN(A)
        _, h, w = ml.solve(b)
        f, c = stat(h)
        row(m, "N (ours)", f, c, w / c, w, h[-1] < 1e-8)
        h, w = markov_jacobi(ml, b)
        f, c = stat(h)
        row(m, "M_J (Jacobi, same operators)", f, c, w / c, w, h[-1] < 1e-8)
        best = None
        for pre, post in ((1, 0), (2, 0), (4, 0), (8, 0), (2, 2), (4, 4)):
            _, h, w = ml.solve(b, markov=(pre, post))
            f, c = stat(h)
            conv = h[-1] < 1e-8
            row(m, f"M_F({pre},{post}) (F-relaxation, same operators)", f, c, w / c, w, conv)
            if conv and (best is None or w < best[0]):
                best = (w, pre, post, w / c)
        for kind in ("AIR", "RS"):
            t = time.time()
            amg = pyamg.air_solver(A) if kind == "AIR" else pyamg.ruge_stuben_solver(A)
            res = []
            amg.solve(b, tol=1e-8, maxiter=30, residuals=res)
            res = np.array(res) / res[0]
            f, c = stat(res)
            cc = amg.cycle_complexity()
            row(m, f"PyAMG {kind} stand-alone", f, c, cc, cc * c, res[-1] < 1e-8)
            its, ok = gmres_with(lambda v: amg.solve(v, maxiter=1, tol=1e-30), A, b)
            row(m, f"PyAMG {kind} + GMRES", float("nan"), its, cc, cc * its, ok, "orthog. not counted")
        wN = [0.0]
        its, ok = gmres_with(lambda v: ml.cycle(0, v, np.zeros_like(v), wN), A, b)
        row(m, "N + GMRES", float("nan"), its, wN[0] / max(its, 1), wN[0], ok, "orthog. not counted")
        if best:
            wM = [0.0]
            its, ok = gmres_with(lambda v: ml.cycle_markov(0, v, np.zeros_like(v), wM, best[1], best[2]), A, b)
            row(m, f"M_F({best[1]},{best[2]}) + GMRES", float("nan"), its, wM[0] / max(its, 1), wM[0], ok,
                "orthog. not counted")


if __name__ == "__main__":
    main(tuple(int(a) for a in sys.argv[1:]) or (127, 255, 511))
