"""Check of the prediction E_N - E_twin = -(I - P~ B_c R A) Q (S - S~) P^sharp (I - T^k)."""
import numpy as np
from ml_twins import ML, poisson
from robust_n import recirc
rng = np.random.default_rng(0)
print("| problem | setting | levels | ‖x_N - x_twin‖/‖x_twin‖ (one cycle) | predicted difference vs measured |")
print("|---|---|---|---|---|")
for pname, A, w in (("recirc Pe 20, w=1", recirc(63), 1.0), ("Poisson, w=0.8", poisson(63), 0.8)):
    for lab, kw in (("no truncation, filtered A_c, recursive, refinement", dict(tauP=0)),
                    ("no truncation, unfiltered", dict(tauP=0, tau_filter=0)),
                    ("truncation 1e-2 (drop), two-level, filtered A_c + LU", dict(tauP=1e-2, max_levels=2)),
                    ("truncation 1e-2 (drop), recursive", dict(tauP=1e-2))):
        ml = ML(A, omega=w, **kw)
        n = A.shape[0]
        b, x0 = np.zeros(n), rng.standard_normal(n)
        xN, xT = ml.cycle(0, b, x0, "N"), ml.cycle(0, b, x0, "twin")
        rel = np.linalg.norm(xN - xT) / np.linalg.norm(xT)
        pred = ""
        if len(ml.levels) == 2:                      # B_c = exact LU of the filtered coarse operator
            L = ml.levels[0]
            xk = x0.copy()
            for _ in range(L.k):
                xk = xk + w * (b - A @ xk) / L.D
            d0 = xk[L.C] - x0[L.C]                   # = c_0 - c_k  (errors are -x since x*=0)
            v = np.zeros(n); v[L.F] = (L.Sex - L.S) @ d0       # Q (S - S~) d0
            Av = A @ v
            v2 = v.copy(); ec = ml.levels[1].lu.solve(Av[L.C] + L.ST @ Av[L.F])
            v2[L.C] -= ec; v2[L.F] -= L.S @ ec                 # (I - P~ B_c R A) Q (S-S~) d0
            # errors e = x* - x = -x  =>  x_N - x_T = +(I - P~B_cRA) Q (S-S~) d0
            pred = f"{np.linalg.norm((xN - xT) - v2) / np.linalg.norm(xN - xT):.1e}"
        print(f"| {pname} | {lab} | {len(ml.levels)} | {rel:.1e} | {pred} |")
