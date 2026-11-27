import sys, numpy as np, pyamg
from flows import cd
from compare_markov import gmres_with, stat
for pe in (6.0, 12.0):
    for m in (127, 255):
        A, xy = cd(m, "const", pe); b = np.random.default_rng(0).standard_normal(A.shape[0])
        for kind in ("AIR", "RS"):
            try:
                amg = pyamg.air_solver(A) if kind == "AIR" else pyamg.ruge_stuben_solver(A)
                res = []; amg.solve(b, tol=1e-8, maxiter=60, residuals=res); res = np.array(res) / res[0]
                f, c = stat(res); cc = amg.cycle_complexity()
                its, ok = gmres_with(lambda v: amg.solve(v, maxiter=1, tol=1e-30), A, b, maxiter=300)
                print(f"| {pe:g} | {m} | {kind} | {f:.3g} | {c if res[-1] < 1e-8 else 'no conv.'} | {cc:.1f} | {cc*c if res[-1] < 1e-8 else float('inf'):.0f} | GMRES its {its}{'' if ok else ' (no conv.)'} | {cc*its if ok else float('inf'):.0f} |"); sys.stdout.flush()
            except Exception as e:
                print(f"| {pe:g} | {m} | {kind} | error {e} |")
