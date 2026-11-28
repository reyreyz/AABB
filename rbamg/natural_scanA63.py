"""Class A (N converges, twin diverges, MJ converges) for recirc and Poiseuille, two-level at m (sparse eigs)."""
import sys, itertools, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); flows = sys.argv[2].split(",")
print(f"m={m}\n| flow | Pe | w | split | k | trunc | N | twin | MJ |\n|---|---|---|---|---|---|---|---|---|")
for flow, pe, w, split, k, tr in itertools.product(flows, (4, 6, 8, 10), (0.05, 0.1, 0.2, 0.3, 0.4), ("rb", "greedy"), (8, 16, 32), ("onepoint", "lump0.5")):
    A, xy = cd(m, flow, pe); n = A.shape[0]; zz = np.zeros(n); I, J = np.divmod(np.arange(n), m)
    kw = dict(split0=((I + J) % 2 == 1)) if split == "rb" else dict(delta=0.3)
    kw.update(dict(mode="lin", xy=xy, nkeep=1) if tr == "onepoint" else dict(mode="lump", tauP=0.5))
    ml = ML(A, omega=w, k0=k, max_levels=2, tau_filter=0, **kw)
    def rho(wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=30)[0])
        except Exception: return float("nan")
    rN = rho("N")
    if not rN < 0.75: continue
    rT = rho("twin")
    if not rT > 1: continue
    rM = rho("markovJ")
    print(f"| {flow} | {pe} | {w} | {split} | {k} | {tr} | {rN:.3g} | {rT:.3g} | {rM:.3g} |"); sys.stdout.flush()
