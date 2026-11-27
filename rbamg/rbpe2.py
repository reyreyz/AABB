"""Idea 15: S~ with 2 entries/row exact on constants only (no x, y), red-black or greedy, Pe 8-12."""
import sys, itertools, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); flow = sys.argv[2]; pes = [float(p) for p in sys.argv[3].split(",")]
for pe in pes:
    A, xy = cd(m, flow, pe); n = A.shape[0]; zz = np.zeros(n)
    I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
    def rho(ml, wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    for split, k, wc in itertools.product(("rb", "greedy"), (1, 2, 3, 4), (0.05, 0.1, 0.2, 0.4, 0.6)):
        ml = ML(A, omega=1.0, omega_c=wc, k0=k, split0=rb if split == "rb" else None, delta=0.3, mode="lin", xy=xy,
                tv=np.zeros((n, 0)), nkeep=2, max_levels=2, tau_filter=0)
        rn = rho(ml, "N")
        if rn < 0.8:
            rt = rho(ml, "twin")
            print(f"| {flow} | {pe:g} | {m} | {split} | {k} | {wc} | {rn:.3g} | {rt:.3g} | {rho(ml,'markovJ') if rt > 1 else float('nan'):.3g} |"); sys.stdout.flush()
