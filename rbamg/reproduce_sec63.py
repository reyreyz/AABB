"""Section 6.3: V-cycles from PyAMG components, one work count (engine.py).

Table SM15: one fixed configuration per method and problem, total WU of a
stationary solve to relative residual 1e-10, m = 127, 255, 511, random b.
Table 3 / Fig. 3: best configuration per family over a (reduced) search.

Problems: PyAMG's upwind advection_2d (theta = pi/4), and the recirculating
flow with eps = 1e-4, 1e-3 (STAND-IN discretization: node upwind on the
unit square with mapped coordinates; the paper's discretization is not given).
"""

import sys
import itertools
import numpy as np
import pyamg
import engine
from reproduce_table1 import recirc as recirc_dense
import scipy.sparse as sp
import sparse_problems as spp


def advection(m):
    A, _ = pyamg.gallery.advection_2d((m, m))
    return sp.csr_matrix(A)


def recirc(m, eps):
    """sparse version of the stand-in recirculating flow (upwind, mapped coords)."""
    beta = spp.RECIRC
    A = spp.cd(m, beta, eps * 2.0, "upwind")   # step 2h on (-1,1)^2: eps*2 keeps eps/(h|beta|) of the dense stand-in
    return A


PROBLEMS = {"Advection": advection,
            "Recirc. 1e-4": lambda m: recirc(m, 1e-4),
            "Recirc. 1e-3": lambda m: recirc(m, 1e-3)}

SM15 = {  # problem: [(method, cfg, paper WU m=127, 255, 511)]
    "Advection": [
        ("lAIR", dict(theta=0.25, interp="1pt", restrict="air2", pre="", post="FF", tauA=1e-4), (29.7, 30.5, 31.2)),
        ("N (P')", dict(theta=0.25, interp="dir", dP=2, restrict="inj", pre="FFF", post="", tauA=1e-4), (17.4, 17.9, 17.9)),
    ],
    "Recirc. 1e-4": [
        ("lAIR", dict(theta=0.25, interp="1pt", restrict="air2", pre="", post="FC", tauA=1e-3), (77.5, 79.2, 97.2)),
        ("N (P')", dict(theta=0.5, interp="dir", dP=2, restrict="inj", pre="FF", post="", tauA=1e-4), (57.4, 72.9, 86.2)),
    ],
    "Recirc. 1e-3": [
        ("lAIR", dict(theta=0.5, interp="1pt", restrict="air2", pre="", post="FFC", tauA=1e-4), (118.4, 104.4, 124.6)),
        ("N (P') [dir.+2F instead of dir. Z2]", dict(theta=0.5, interp="dir", dP=2, restrict="inj", pre="", post="FF", tauA=1e-4), (82.5, 105.1, 118.3)),
    ],
}


def run(A, cfg, b):
    V = engine.VCycle(A, cfg)
    h, w = V.solve(b)
    ok = h[-1] < 1e-10
    return (w[-1] if ok else np.inf), len(h) - 1, V.op_complexity()


def table_sm15(ms=(127, 255, 511)):
    print("\n## Table SM15 (fixed configurations): total WU to 1e-10, ours / paper (cycles, oc)\n")
    print("| problem | method | config | " + " | ".join(f"m = {m}" for m in ms) + " |")
    print("|---|---|---|" + "---|" * len(ms))
    for pname, rows in SM15.items():
        mats = {m: PROBLEMS[pname](m) for m in ms}
        for meth, cfg, paper in rows:
            cells = []
            for m, p in zip(ms, paper):
                A = mats[m]
                b = np.random.default_rng(0).standard_normal(A.shape[0])
                wu, cyc, oc = run(A, cfg, b)
                cells.append(f"{wu:.1f} / {p} ({cyc}, {oc:.2f})")
            print(f"| {pname} | {meth} | {engine.describe(cfg)} | " + " | ".join(cells) + " |")
            sys.stdout.flush()


def family_space(fam):
    th = (0.25, 0.5)
    tA = (1e-4, 1e-3)
    if fam == "lAIR":
        for t, ta, deg, post in itertools.product(th, tA, ("air1", "air2"), ("F", "FF", "FC", "FFC", "FFF")):
            yield dict(theta=t, interp="1pt", restrict=deg, pre="", post=post, tauA=ta)
    elif fam == "lAIR+P'":
        for t, ta, dP, post in itertools.product(th, tA, (1, 2, 3), ("F", "FF", "FC", "FFC")):
            yield dict(theta=t, interp="1pt", dP=dP, restrict="air2", pre="", post=post, tauA=ta)
    elif fam == "RS":
        for t, ta, sch in itertools.product(th, tA, (("F", "F"), ("FC", "FC"), ("FF", "FF"))):
            yield dict(theta=t, interp="dir", restrict="PT", pre=sch[0], post=sch[1], tauA=ta)
    elif fam == "N (P')":
        for t, ta, ip, dP, sch in itertools.product(th, tA, ("1pt", "dir"), (1, 2, 3),
                                                     (("F", ""), ("FF", ""), ("FFF", ""), ("F", "F"), ("FF", "C"))):
            yield dict(theta=t, interp=ip, dP=dP, restrict="inj", pre=sch[0], post=sch[1], tauA=ta)


PAPER_T3 = {  # (problem, m): lAIR, lAIR+P', RS, N(P')
    ("Advection", 127): (28.9, 29.6, 166.1, 12.8), ("Advection", 255): (30.5, 34.0, 171.5, 17.9),
    ("Recirc. 1e-4", 127): (58.8, 57.8, 153.4, 51.1), ("Recirc. 1e-4", 255): (66.7, 70.2, 148.2, 64.7),
    ("Recirc. 1e-3", 127): (77.3, 81.4, 175.9, 61.8), ("Recirc. 1e-3", 255): (104.4, 125.0, 232.3, 95.4),
}


def table3(ms=(127, 255)):
    print("\n## Table 3 (reduced search): best total WU to 1e-10 per family, ours / paper\n")
    fams = ("lAIR", "lAIR+P'", "RS", "N (P')")
    print("| problem | m | " + " | ".join(fams) + " |\n|---|---|" + "---|" * len(fams))
    best_cfgs = {}
    for pname in PROBLEMS:
        for m in ms:
            A = PROBLEMS[pname](m)
            b = np.random.default_rng(0).standard_normal(A.shape[0])
            cells = []
            for fi, fam in enumerate(fams):
                best = (np.inf, None, 0)
                for cfg in family_space(fam):
                    try:
                        wu, cyc, oc = run(A, cfg, b)
                    except Exception:
                        continue
                    if wu < best[0]:
                        best = (wu, cfg, cyc)
                best_cfgs[(pname, m, fam)] = best
                cells.append(f"{best[0]:.1f} / {PAPER_T3[(pname, m)][fi]}")
            print(f"| {pname} | {m} | " + " | ".join(cells) + " |")
            sys.stdout.flush()
    print("\nBest configurations:\n")
    for (p, m, f), (wu, cfg, cyc) in best_cfgs.items():
        print(f"- {p}, m = {m}, {f}: {wu:.1f} WU, {cyc} cycles, {engine.describe(cfg) if cfg else '-'}")


if __name__ == "__main__":
    print("# Section 6.3 reproduction (PyAMG components, own work count)")
    print("Recirculating flow: stand-in discretization (not identified).")
    table_sm15()
    table3()
