"""Diverging smoother: coarse solve on the step-0 vs step-k C-error.
E_j e0 = P~ (I - B_c A_c) c_j + (I - P~ B_c R A) Q K^k q0   (no truncation)
=> one-cycle growth ~ delta(c_j) * ||c_j|| / ||e0||, with ||c_k|| = g ||c_0||."""
import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from robust_n import recirc
for m in (127, 255):
    A = recirc(m); n = A.shape[0]
    rng = np.random.default_rng(0); b = rng.standard_normal(n); xs = sla.spsolve(A.tocsc(), b)
    ml = ML(A, omega=1.0, k1=8, tau_deep=3e-4)
    L = ml.levels[0]
    x = np.zeros(n)
    for _ in range(2): x = ml.cycle(0, b, x, "N")          # representative error
    e0 = xs - x
    ek = e0.copy()
    for _ in range(L.k): ek = ek - (A @ ek) / L.D            # T^k e0
    c0, ck = e0[L.C], ek[L.C]
    def defect(c):                                           # ||(I - B_c A_c) c|| / ||c||
        ml.shift = 0.0
        return np.linalg.norm(c - ml.coarse(0, ml.RAP(L, c), "N")) / np.linalg.norm(c)
    out = []
    for sh in (0.0, 1.0):
        ml.shift = sh
        out.append(np.linalg.norm(xs - ml.cycle(0, b, x, "N")) / np.linalg.norm(e0))
    ml.shift = 0.0
    eMJ = np.linalg.norm(xs - ml.cycle(0, b, x, "markovJ")) / np.linalg.norm(e0)
    eT = np.linalg.norm(xs - ml.cycle(0, b, x, "twin")) / np.linalg.norm(e0)
    g = np.linalg.norm(ck) / np.linalg.norm(c0)
    d0, dk = defect(c0), defect(ck)
    print(f"| {m} | {g:.3g} | {d0:.2e} | {dk:.2e} | {d0*np.linalg.norm(c0)/np.linalg.norm(e0):.2e} | {out[0]:.2e} | "
          f"{dk*np.linalg.norm(ck)/np.linalg.norm(e0):.2e} | {out[1]:.2e} | {eT:.2e} | {eMJ:.2e} |"); sys.stdout.flush()
