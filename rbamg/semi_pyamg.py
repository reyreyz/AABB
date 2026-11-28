import sys, numpy as np, scipy.sparse as sp, pyamg
import sparse_problems as spp
from flows import cd as cdflow
from compare_markov import gmres_with, stat
m = int(sys.argv[1])
for name, A in (("shear:20", spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 20, "2", False), "centered")),
                ("const:6", cdflow(m, "const", 6)[0]), ("shear:6", spp.cd(m, spp.SHEAR, spp.eps_for_peclet(m, spp.SHEAR, 6, "2", False), "centered"))):
    A = sp.csr_matrix(A); b = np.random.default_rng(0).standard_normal(A.shape[0])
    for kind in ("AIR", "RS"):
        try:
            amg = pyamg.air_solver(A) if kind == "AIR" else pyamg.ruge_stuben_solver(A)
            res = []; amg.solve(b, tol=1e-8, maxiter=80, residuals=res); res = np.array(res) / res[0]; f, c = stat(res); cc = amg.cycle_complexity()
            its, ok = gmres_with(lambda v: amg.solve(v, maxiter=1, tol=1e-30), A, b, maxiter=200)
            print(f"| {name} | {m} | PyAMG {kind} | {f:.3g} | {c if res[-1] < 1e-8 else 'no conv.'} | {cc:.1f} | {cc*c if res[-1] < 1e-8 else float('inf'):.0f} | GMRES: {its if ok else 'no conv.'} its, {cc*its if ok else float('inf'):.0f} WU |"); sys.stdout.flush()
        except Exception as e:
            print(f"| {name} | {m} | PyAMG {kind} | error {e} |")
