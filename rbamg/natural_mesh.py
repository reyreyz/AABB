import sys, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); cases = [c.split(":") for c in sys.argv[2].split(",")]   # flow:pe:w:split:k:trunc
for flow, pe, w, split, k, tr in cases:
    pe, w, k = float(pe), float(w), int(k)
    A, xy = cd(m, flow, pe); n = A.shape[0]; zz = np.zeros(n); I, J = np.divmod(np.arange(n), m)
    kw = dict(split0=((I + J) % 2 == 1)) if split == "rb" else dict(delta=0.3)
    kw.update(dict(mode="lin", xy=xy, nkeep=1) if tr == "onepoint" else dict(mode="drop", tauP=0.3) if tr == "drop0.3" else dict(mode="lump", tauP=0.5) if tr == "lump0.5" else dict(mode="drop", tauP=0.0))
    ml = ML(A, omega=w, k0=k, max_levels=2, tau_filter=0, **kw)
    def rho(wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=2000, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    print(f"| {flow} | {pe:g} | {w} | {split} | {k} | {tr} | {m} | {rho('N'):.3g} | {rho('twin'):.3g} | {rho('markovJ'):.3g} |"); sys.stdout.flush()
