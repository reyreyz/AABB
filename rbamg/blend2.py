import sys, itertools, numpy as np, scipy.sparse.linalg as sla
from ml_twins import ML
from flows import cd
m = int(sys.argv[1]); flow = sys.argv[2]; pes = [float(p) for p in sys.argv[3].split(",")]; ths = [float(t) for t in sys.argv[4].split(",")]
for pe in pes:
    A, xy = cd(m, flow, pe); n = A.shape[0]; zz = np.zeros(n)
    I, J = np.divmod(np.arange(n), m); rb = (I + J) % 2 == 1
    def rho(ml, wh):
        op = sla.LinearOperator((n, n), matvec=lambda v: ml.cycle(0, zz, v, wh))
        try: return abs(sla.eigs(op, k=1, which="LM", return_eigenvectors=False, maxiter=1500, tol=1e-3, ncv=40)[0])
        except Exception: return float("nan")
    for th in ths:
        rt = rho(ML(A, omega=1.0, omega_c=0.05, k0=1, split0=rb, mode="lin", xy=xy, theta=th, max_levels=2, tau_filter=0), "twin")
        for k, wc in itertools.product((4, 6, 8), (0.02, 0.05, 0.1)):
            ml = ML(A, omega=1.0, omega_c=wc, k0=k, split0=rb, mode="lin", xy=xy, theta=th, max_levels=2, tau_filter=0)
            rn = rho(ml, "N")
            if rn < 0.6:
                print(f"| {flow} | {pe:g} | {m} | {th} | {k} | {wc} | {rn:.3g} | {rt:.3g} | {rho(ml,'markovJ'):.3g} |"); sys.stdout.flush()
