"""Section 7: multilevel examples (Algorithms 7.1-7.3), Table 4 and Fig. 4, with
PyAMG comparisons.

Problems on (-1,1)^2, m x m interior points, x0 = 0, random b (seed 0):
  7.1  upwind beta . grad u = f, stagnation flow beta = (x, -y), m = 180
  7.2  centered -eps Lap u + beta . grad u, shear beta = (-y, 0.3)
  7.3  centered recirculating beta = (2y(1-x^2), -2x(1-y^2))
Maximal cell Peclet number |beta| h / eps = 20 (identified: gives rho(T) = 4.09
for damped Jacobi omega = 0.7 on the shear flow; paper "rho(T) ~ 4").
RS settings for Alg. 7.1 (identified: operator complexity 3.72 vs 3.7):
classical strength theta = 0.25, 'min' measure, RS with second pass.
"""

import sys
import time
import numpy as np
import scipy.sparse.linalg as sla
import pyamg
import sparse_problems as spp
import ml_cycles as mc


def rate(h):
    """asymptotic factor: geometric mean over the second half of the history"""
    h = np.asarray(h)
    a = len(h) // 2
    return (h[-1] / h[a]) ** (1.0 / (len(h) - 1 - a)) if len(h) - 1 > a else h[-1]


def pyamg_run(ml, A, b, tol=1e-10, maxiter=200, accel=None):
    res = []
    ml.solve(b, tol=tol, maxiter=maxiter, residuals=res, accel=accel)
    res = np.array(res) / res[0]
    it = len(res) - 1
    cc = ml.cycle_complexity()
    conv = res[-1] < tol
    return it, cc, cc * it, conv, res


def air_configs():
    """A small configuration search for PyAMG's AIR (F/C- and F-relaxation)."""
    cfgs = []
    for theta in (0.1, 0.25, 0.3, 0.5):
        for deg in (1, 2):
            for f, c in ((1, 1), (2, 1), (3, 1), (1, 0), (2, 0), (3, 0)):
                cfgs.append(dict(theta=theta, degree=deg, f=f, c=c))
    return cfgs


def build_air(A, cfg):
    return pyamg.air_solver(
        A, strength=("classical", {"theta": cfg["theta"], "norm": "min"}),
        restrict=("air", {"theta": 0.05, "degree": cfg["degree"]}),
        postsmoother=("fc_jacobi", {"omega": 1.0, "iterations": 1, "withrho": False,
                                    "f_iterations": cfg["f"], "c_iterations": cfg["c"]}),
        max_coarse=150)


def sec71(m=180):
    print(f"\n## 7.1 Markovian V-cycle, stagnation flow, m = {m}\n")
    A = spp.cd(m, spp.STAGNATION, 0.0, "upwind")
    b = np.random.default_rng(0).standard_normal(A.shape[0])
    ml = mc.Alg71(A, theta=0.25, norm="min", second_pass=True)
    h, w = mc.solve(ml, A, b)
    print(f"Algorithm 7.1: {len(h)-1} cycles, {w[-1]:.1f} WU, op. complexity "
          f"{ml.op_complexity():.2f}  (paper: 4 cycles, 19.3 WU, 3.7)\n")
    best = {"F/C": None, "F": None}
    for cfg in air_configs():
        try:
            it, cc, wu, conv, _ = pyamg_run(build_air(A, cfg), A, b)
        except Exception:
            continue
        fam = "F/C" if cfg["c"] else "F"
        if conv and (best[fam] is None or wu < best[fam][0]):
            best[fam] = (wu, it, cfg)
    for fam, v in best.items():
        print(f"best PyAMG AIR, {fam}-relaxation ({len(air_configs())//2} configs): "
              f"{v[0]:.1f} WU, {v[1]} cycles, {v[2]}")
    print("(paper: AIR F/C best of 41 configs 149.6 WU / 12 cycles; F-relaxation, 64 configs: "
          "57.4 WU / 7 cycles; Alg. 7.1 cheaper by 7.8 and 3.0)")
    print(f"ratio vs best F/C: {best['F/C'][0]/w[-1]:.1f}, vs best F: {best['F'][0]/w[-1]:.1f}")
    print("\nWU of PyAMG solvers = ml.cycle_complexity() x cycles (PyAMG's own count).")
    return h, w


def table4(ms=(127, 255, 511, 1023)):
    print("\n## Table 4: Alg. 7.2 on the shear flow, Alg. 7.3 on the recirculating flow\n")
    print("| m | n | 7.2 rho | WU/cyc | WU/digit | Markov 7.2 | 7.3 cycles to 1e-8 | WU/cyc | WU/digit | Markov 7.3 (1 cycle) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    hists = {}
    for m in ms:
        n = m * m
        A2 = spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 20, "2", False), "centered")
        b = np.random.default_rng(0).standard_normal(n)
        ml2 = mc.Alg72(A2)
        h2, w2 = mc.solve(ml2, A2, b, maxcyc=40)
        r2 = rate(h2)
        wc2 = w2[-1] / (len(h2) - 1)
        hm, _ = mc.solve(mc.Alg72(A2, markov=True), A2, b, maxcyc=5)
        A3 = spp.cd(m, spp.RECIRC, spp.eps_for_peclet(m, spp.RECIRC, 20, "2", False), "centered")
        ml3 = mc.Alg73(A3)
        h3, w3 = mc.solve(ml3, A3, b, tol=1e-8, maxcyc=10)
        wc3 = w3[-1] / (len(h3) - 1)
        hm3, _ = mc.solve(mc.Alg73(A3, markov=True), A3, b, tol=1e-8, maxcyc=1)
        hists[m] = (h2, w2, h3, w3)
        print(f"| {m} | {n:.1e} | {r2:.3f} | {wc2:.0f} | {wc2/-np.log10(r2):.0f} | "
              f"{'diverges' if hm[-1] > 1 else f'{hm[-1]:.1e}'} | {len(h3)-1} | {wc3:.0f} | "
              f"{w3[-1]/-np.log10(h3[-1]):.0f} | {hm3[-1]:.1e} |")
        sys.stdout.flush()
    print("\nPaper: 7.2 rho 0.063/0.090/0.108/0.157, WU/cyc 214/225/225/231 (m = 127..1023); "
          "7.3 two cycles at every size, WU/cyc 3324/3952/4398/4694; Markovian cycles diverge "
          "(7.3: residual 1e41 after one cycle).")
    return hists


def pyamg_comparison(m=511):
    print(f"\n## Fig. 4 comparison with PyAMG stand-alone cycles, m = {m}\n")
    print("| problem | solver | cycles | WU (PyAMG count) | final rel. residual |\n|---|---|---|---|---|")
    for name, beta in (("shear", spp.SHEAR), ("recirc", spp.RECIRC)):
        A = spp.cd(m, beta, spp.eps_for_peclet(m, beta, 20, "2", False), "centered")
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        solvers = {
            "AIR F/C": lambda: pyamg.air_solver(A, max_coarse=150),
            "AIR F": lambda: build_air(A, dict(theta=0.3, degree=2, f=2, c=0)),
            "SA": lambda: pyamg.smoothed_aggregation_solver(A, max_coarse=150),
            "RS": lambda: pyamg.ruge_stuben_solver(A, max_coarse=150),
        }
        for sname, mk in solvers.items():
            try:
                it, cc, wu, conv, res = pyamg_run(mk(), A, b, maxiter=60)
                status = f"{res[-1]:.1e}" if np.isfinite(res[-1]) else "nan"
                print(f"| {name} | {sname} | {it} | {wu:.0f} | {status}{'' if conv else ' (no conv.)'} |")
            except Exception as e:
                print(f"| {name} | {sname} | – | – | failed: {type(e).__name__} |")
            sys.stdout.flush()
        if name == "recirc":
            ml = pyamg.air_solver(A, max_coarse=150)
            it, cc, wu, conv, _ = pyamg_run(ml, A, b, maxiter=200, accel="gmres")
            print(f"| recirc | AIR F/C + GMRES | {it} | {wu:.0f} (no orthog.) | {'conv' if conv else 'no conv.'} |")
    print("\nPaper (m = 511): shear: AIR F/C diverges, AIR F converges 3-5x cheaper than Alg. 7.2 "
          "(738 vs 2698 WU); recirc: AIR (F/C, F), SA, RS all diverge; AIR + GMRES 35 its, 570 WU.")


if __name__ == "__main__":
    # usage: reproduce_sec7.py [71] [t4 m1 m2 ...] [cmp m]
    t = time.time()
    args = sys.argv[1:] or ["71", "t4", "127", "255", "511", "1023", "cmp", "511"]
    print("# Section 7 reproduction")
    i = 0
    while i < len(args):
        if args[i] == "71":
            h71, w71 = sec71()
            np.savez("sec7_alg71.npz", h=h71, w=w71)
            i += 1
        elif args[i] == "t4":
            j = i + 1
            while j < len(args) and args[j].isdigit():
                j += 1
            hists = table4(tuple(int(a) for a in args[i + 1:j]))
            for m, hv in hists.items():
                np.savez(f"sec7_table4_m{m}.npz", **dict(zip(("h72", "w72", "h73", "w73"), hv)))
            i = j
        elif args[i] == "cmp":
            pyamg_comparison(int(args[i + 1]))
            i += 2
    print(f"\n(total time {time.time()-t:.0f} s)")
