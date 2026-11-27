"""2D multilevel: start from a convergent N (Algorithm 7.3 setting) and optimize to N_Z.

Problem (section 7.3): centered recirculating flow on (-1,1)^2, cell Peclet 20
(Jacobi diverges, rho(T) ~ 4.5; the Markovian cycle diverges).
Cycle: Alg73Z (ml_z73.py): greedy F-points ||I - D_ff^{-1}A_ff||_inf <= delta,
Jacobi omega = 1, P0, injection, k steps per level, Petrov-Galerkin coarse
operators (A P')_C filtered, recursive V-cycle.  N uses Z = 0 (complement Q0);
N_Z uses a bounded Z (<= q entries per column, coarse points within distance 2,
values from min ||K_Z^k||_F).  Reported: cycles and total work to 1e-8.
"""

import sys
import time
import numpy as np
import sparse_problems as spp
import ml_z73 as mzz
from run_ml_nz import rate


def problem(m):
    return spp.cd(m, spp.RECIRC, spp.eps_for_peclet(m, spp.RECIRC, 20, "2", False), "centered")


def run(A, b, **kw):
    t = time.time()
    c = mzz.Alg73Z(A, **kw)
    ts = time.time() - t
    h, w = c.solve(b, tol=1e-8, maxcyc=40)
    conv = h[-1] < 1e-8
    return dict(levels=len(c.levels), rho=rate(h), cyc=len(h) - 1, conv=conv,
                wu=w[-1] if conv else np.inf, wuc=w[-1] / max(len(h) - 1, 1),
                oc=c.op_complexity(), sc=c.stored_complexity(), ts=ts)


CONFIGS = [  # label, kwargs
    ("N k=4, dc=0.3, tauA=1e-3", dict(zmode="Q0", k=4, delta=0.3)),
    ("N k=4, dc=0.3, tauA=1e-2", dict(zmode="Q0", k=4, delta=0.3, tauA=1e-2)),
    ("N_Z k=4, dc=0.3, tauA=1e-2", dict(k=4, delta=0.3, q=4, theta="auto", tauA=1e-2)),
    ("N_Z k=4, dc=0.3, tauA=1e-2, capP=12", dict(k=4, delta=0.3, q=4, theta="auto", tauA=1e-2, capP=12)),
    ("N k=4, dc=0.6, tauA=1e-2", dict(zmode="Q0", k=4, delta=0.3, delta_coarse=0.6, tauA=1e-2)),
    ("N_Z k=4, dc=0.6, tauA=1e-2", dict(k=4, delta=0.3, delta_coarse=0.6, q=4, theta="auto", tauA=1e-2)),
    ("N_Z k=4, dc=0.6, tauA=3e-2", dict(k=4, delta=0.3, delta_coarse=0.6, q=4, theta="auto", tauA=3e-2)),
    ("N_Z k=4, dc=0.6, tauA=1e-2, capP=12", dict(k=4, delta=0.3, delta_coarse=0.6, q=4, theta="auto", tauA=1e-2, capP=12)),
]


def main(ms=(63,), tauA=1e-3):
    print("# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)\n")
    print("| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = problem(m)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for lab, kw in CONFIGS:
            kw = dict(kw)
            r = run(A, b, tauA=kw.pop('tauA', tauA), **kw)
            print(f"| {m} | {lab} | {r['levels']} | {r['rho']:.3g} | {r['cyc']}{'' if r['conv'] else ' (no conv.)'} | "
                  f"{r['wuc']:.0f} | {r['wu']:.0f} | {r['oc']:.2f} | {r['sc']:.1f} | {r['ts']:.1f} |")
            sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (63,)
    main(ms)
