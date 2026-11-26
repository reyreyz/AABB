"""Multilevel tests of the N_Z V-cycle (ml_nz.py) on the CD problems.

CD: -eps Lap u + beta . grad u, beta = (cos pi/6, sin pi/6), upwind, unit square.
For each problem, size and k: cycles to relative residual 1e-10 (x0 = 0, random
b), asymptotic factor, work per cycle and total (WU = nnz(A0)), operator
complexity, setup time.  Baselines: Markovian cycle with the same hierarchy (M),
and PyAMG's AIR and Ruge-Stueben solvers (PyAMG's cycle_complexity as work).
"""

import sys
import time
import numpy as np
import pyamg
import ml_nz as mz


def rate(h):
    h = np.asarray(h)
    a = max(1, len(h) // 2)
    return (h[-1] / h[a]) ** (1.0 / (len(h) - 1 - a)) if len(h) - 1 > a else h[-1]


def run_mode(A, m, k, mode, tauA, b, max_levels=99):
    t = time.time()
    cyc = mz.NZCycle(A, m, k=k, zmode=("Q0" if mode == "M" else mode), tauA=tauA, markov=(mode == "M"), max_levels=max_levels)
    ts = time.time() - t
    h, w = cyc.solve(b)
    conv = h[-1] < 1e-10
    return dict(cyc=len(h) - 1, rho=rate(h), wuc=w[-1] / (len(h) - 1), wu=w[-1] if conv else np.inf,
                oc=cyc.op_complexity(), ts=ts, conv=conv, levels=len(cyc.levels))


def run_pyamg(A, b, kind):
    t = time.time()
    ml = pyamg.air_solver(A, max_coarse=150) if kind == "AIR" else pyamg.ruge_stuben_solver(A, max_coarse=150)
    ts = time.time() - t
    res = []
    ml.solve(b, tol=1e-10, maxiter=100, residuals=res, accel=None)
    res = np.array(res) / res[0]
    conv = res[-1] < 1e-10
    cc = ml.cycle_complexity()
    return dict(cyc=len(res) - 1, rho=rate(res), wuc=cc, wu=cc * (len(res) - 1) if conv else np.inf,
                oc=ml.operator_complexity(), ts=ts, conv=conv, levels=len(ml.levels))


VARIANTS = [  # label, NZCycle keyword arguments
    ("M (Markovian)", dict(zmode="Q0", markov=True)),
    ("N (Q0)", dict(zmode="Q0")),
    ("N_Z LS, all levels", dict(zmode="LS")),
    ("N_Z LS, finest only", dict(zmode="LS", z_levels=1)),
    ("N_Z ALS, finest only", dict(zmode="ALS", z_levels=1)),
    ("N_Z Kk, finest only", dict(zmode="Kk", z_levels=1)),
]


def main(ms=(63, 127, 255), ks=(2, 4), eps_list=(1e-2, 1e-3), tauA=1e-3):
    print("# Multilevel N_Z V-cycles on CD (bilinear hierarchy, Galerkin Rhat = P^T, Jacobi 2/3, "
          f"filter tauA = {tauA})\n")
    for eps in eps_list:
        print(f"\n## CD eps = {eps:g}\n")
        print("| m | k | cycle | levels | op.cx | cycles | factor | WU/cycle | WU total | setup [s] |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for m in ms:
            A = mz.cd_sparse(m, eps)
            b = np.random.default_rng(0).standard_normal(A.shape[0])
            for k in ks:
                for lab, kw in VARIANTS:
                    t = time.time()
                    cyc = mz.NZCycle(A, m, k=k, tauA=tauA, **kw)
                    ts = time.time() - t
                    h, w = cyc.solve(b, maxcyc=80)
                    conv = h[-1] < 1e-10
                    print(f"| {m} | {k} | {lab} | {len(cyc.levels)} | {cyc.op_complexity():.2f} | "
                          f"{len(h)-1}{'' if conv else ' (no conv.)'} | {rate(h):.2e} | "
                          f"{w[-1]/(len(h)-1):.1f} | {w[-1] if conv else float('inf'):.1f} | {ts:.1f} |")
                    sys.stdout.flush()
            for kind in ("AIR", "RS"):
                r = run_pyamg(A, b, kind)
                print(f"| {m} | – | PyAMG {kind} | {r['levels']} | {r['oc']:.2f} | {r['cyc']}{'' if r['conv'] else ' (no conv.)'} | "
                      f"{r['rho']:.2e} | {r['wuc']:.1f} | {r['wu']:.1f} | {r['ts']:.1f} |")
                sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (63, 127, 255)
    main(ms)
