"""Idea 20, step 1: splits where pointwise F-relaxation diverges (rho(K0) > 1) but K0 + H Z is stable with sparse Z."""
import sys, numpy as np, scipy.sparse.linalg as sla
import ml_nz as mz, ml_z73 as mzz
from flows import cd
def rho(M):
    try: return abs(sla.eigs(M, k=1, which="LM", return_eigenvectors=False, maxiter=3000, tol=1e-4)[0])
    except Exception: return float("nan")
m = int(sys.argv[1])
for flow in ("const", "recirc"):
    for pe in (6, 10, 20):
        A, xy = cd(m, flow, pe); n = A.shape[0]; I, J = np.divmod(np.arange(n), m)
        for sname, isC in (("full (odd,odd)", (I % 2 == 1) & (J % 2 == 1)), ("semi-x (odd cols)", J % 2 == 1)):
            isF = ~isC; C, F = np.nonzero(isC)[0], np.nonzero(isF)[0]
            for w in (1.0, 0.7):
                Aff, Afc = A[F][:, F], A[F][:, C]
                import scipy.sparse as sp
                Dfi = sp.diags(w / Aff.diagonal()); K0 = (sp.eye(len(F)) - Dfi @ Aff).tocsr(); H = (-Dfi @ Afc).tocsr()
                r0 = rho(K0); out = [f"rho(K0)={r0:.3g}"]
                for q in (2, 4):
                    patt = mzz.bounded_pattern(A, C, F, K0, H, q)
                    Zl = mz.z_ls(K0, H, patt); rl = rho((K0 + H @ Zl).tocsr())
                    Zk = mz.z_kk(K0, H, patt, 2, Zl, maxiter=60); rk = rho((K0 + H @ Zk).tocsr())
                    out.append(f"q={q}: LS {rl:.3g}, KK {rk:.3g}")
                print(f"| {flow} | {pe} | {sname} | w={w} | " + " | ".join(out) + " |"); sys.stdout.flush()
