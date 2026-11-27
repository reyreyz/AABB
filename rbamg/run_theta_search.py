"""Convergent N first, then N_Z by a safeguarded line search on theta.

Z (bounded, CR-optimized) on levels 1..L_Z is scaled by theta; theta is chosen
from {0, 1/8, 1/4, 1/2, 3/4, 1} by the measured contraction of the actual
multilevel cycle (4 test cycles on b = 0 from a random start).  theta = 0 is
plain N, so the selected N_Z is never worse than N on the test.
Problem: section 7.3 recirculating flow, Pe 20; cycle with refinement on levels 0-3.
"""

import sys
import time
import numpy as np
import ml_z73 as mzz
from run_base_nz import problem
from run_ml_nz import rate

THETAS = (0.0, 0.125, 0.25, 0.5, 0.75, 1.0)


def contraction(c, n, cycles=4, seed=3):
    x = np.random.default_rng(seed).standard_normal(n)
    b = np.zeros(n)
    nr = [np.linalg.norm(c.levels[0].A @ x)]
    for _ in range(cycles):
        x = c.cycle(0, b, x, [0.0])
        nr.append(np.linalg.norm(c.levels[0].A @ x))
        if not np.isfinite(nr[-1]) or nr[-1] > 1e30 * nr[0]:
            return np.inf
    return (nr[-1] / nr[1]) ** (1.0 / (cycles - 1))


def run(A, b, kw):
    t = time.time()
    c = mzz.Alg73Z(A, **kw)
    ts = time.time() - t
    h, w = c.solve(b, tol=1e-8, maxcyc=40)
    conv = h[-1] < 1e-8
    return c, dict(rho=rate(h), cyc=len(h) - 1, conv=conv, wu=w[-1] if conv else np.inf,
                   wuc=w[-1] / max(len(h) - 1, 1), oc=c.op_complexity(), ts=ts, levels=len(c.levels))


def main(ms=(63, 127, 255)):
    print("# Convergent N, then N_Z by a safeguarded theta search (recirculating flow, Pe 20, refinement on levels 0-3)\n")
    print("| m | cycle | theta | test contraction | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | setup [s] |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = problem(m)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for kk, kc in ((4, 6), (3, 4)):
            base = dict(k=kk, k_coarse=kc, delta=0.3, refine_levels=3)
            c0, r0 = run(A, b, dict(base, zmode="Q0"))
            print(f"| {m} | N k={kk}/{kc} | – | {contraction(c0, A.shape[0]):.3g} | {r0['rho']:.3g} | "
                  f"{r0['cyc']}{'' if r0['conv'] else ' (no conv.)'} | {r0['wuc']:.0f} | {r0['wu']:.0f} | {r0['oc']:.2f} | {r0['ts']:.1f} |")
            sys.stdout.flush()
            t = time.time()
            best = None
            for th in THETAS:
                kw = dict(base, q=4, z_from=1, z_levels=3, theta=th) if th > 0 else dict(base, zmode="Q0")
                c = mzz.Alg73Z(A, **kw)
                ct = contraction(c, A.shape[0])
                if best is None or ct < best[0]:
                    best = (ct, th, kw)
            tsearch = time.time() - t
            c1, r1 = run(A, b, best[2])
            print(f"| {m} | N_Z k={kk}/{kc}, Z on levels 1-2 | {best[1]:g} | {best[0]:.3g} | {r1['rho']:.3g} | "
                  f"{r1['cyc']}{'' if r1['conv'] else ' (no conv.)'} | {r1['wuc']:.0f} | {r1['wu']:.0f} | {r1['oc']:.2f} | "
                  f"{r1['ts']:.1f} (+{tsearch:.0f} search) |")
            sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (63, 127, 255)
    main(ms)
