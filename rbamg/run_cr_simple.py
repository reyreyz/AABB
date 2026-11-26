"""CR-centric, simple complements in the multilevel N_Z V-cycle (CD problems).

Z per level (pattern: the coarse points interpolating to each fine point):
  LS   column j: min ||K_Z e_j||          (one-step CR, the paper's N_Z)
  LSk  column j: min ||K^{k-1} K_Z e_j||  (k-step CR with the Q0 propagator frozen)
optionally damped, Z <- theta Z, theta in {0, 1/4, 1/2, 3/4, 1} chosen per level to
minimize the CR rate of the level pair max(rho(K_{theta Z}), rho(K_c(theta))),
K_c = local CR propagator of the next coarse operator.
"""

import sys
import time
import numpy as np
import ml_nz as mz
from run_ml_nz import rate, run_pyamg

VARIANTS = [
    ("N (Q0)", dict(zmode="Q0")),
    ("LS", dict(zmode="LS")),
    ("LS, theta auto", dict(zmode="LS", theta="auto")),
    ("LSk", dict(zmode="LSk")),
    ("LSk, theta auto", dict(zmode="LSk", theta="auto")),
    ("LSk, finest only", dict(zmode="LSk", z_levels=1)),
]


def main(ms=(63,), eps_list=(1e-3, 1e-2), ks=(2, 4), tauA=1e-3, with_pyamg=True):
    print("# CR-centric complements, multilevel V-cycles (bilinear hierarchy, Galerkin, Jacobi 2/3, "
          f"tauA = {tauA})\n")
    for eps in eps_list:
        print(f"\n## CD eps = {eps:g}\n")
        print("| m | k | Z | theta per level | factor | cycles | WU/cycle | WU total | op.cx | setup [s] |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for m in ms:
            A = mz.cd_sparse(m, eps)
            b = np.random.default_rng(0).standard_normal(A.shape[0])
            for k in ks:
                for lab, kw in VARIANTS:
                    t = time.time()
                    c = mz.NZCycle(A, m, k=k, tauA=tauA, **kw)
                    ts = time.time() - t
                    h, w = c.solve(b, maxcyc=80)
                    conv = h[-1] < 1e-10
                    th = ",".join(f"{getattr(L, 'theta', 1.0):g}" for L in c.levels[:-1]) if kw.get("zmode") != "Q0" else "–"
                    print(f"| {m} | {k} | {lab} | {th} | {rate(h):.3f} | {len(h)-1}{'' if conv else ' (no conv.)'} | "
                          f"{w[-1]/(len(h)-1):.1f} | {w[-1] if conv else float('inf'):.0f} | {c.op_complexity():.2f} | {ts:.1f} |")
                    sys.stdout.flush()
            if with_pyamg:
                for kind in ("AIR", "RS"):
                    r = run_pyamg(A, b, kind)
                    print(f"| {m} | – | PyAMG {kind} | – | {r['rho']:.3f} | {r['cyc']}{'' if r['conv'] else ' (no conv.)'} | "
                          f"{r['wuc']:.1f} | {r['wu']:.0f} | {r['oc']:.2f} | {r['ts']:.1f} |")
                    sys.stdout.flush()


if __name__ == "__main__":
    ms = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (63,)
    main(ms)
