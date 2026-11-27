import sys, numpy as np
from ml_twins import ML
from robust_n import recirc
m = int(sys.argv[1]); A = recirc(m, 10.0); n = A.shape[0]
h = 2.0 / (m + 1); t = -1 + h * np.arange(1, m + 1); X, Y = np.meshgrid(t, t); xy = np.column_stack([X.ravel(), Y.ravel()])
kw = dict(omega=1.0, omega_c=0.8, wc_levels=1, k0=1, k1=8, delta=0.3, mode="lin", xy=xy, nkeep=2,
          lin_levels=1, tauP_deep=float(sys.argv[2]) if len(sys.argv) > 2 else 1e-5, tau_filter=1e-3, tau_deep=3e-4, refine=1)
ml = ML(A, **kw); n0 = A.nnz
print(f"m={m}, tauP_deep={kw['tauP_deep']}\n| lev | n | nnz(A)/nnz0 | nnz/row | nnz(K+H)/nnz0 | nnz(S)/nnz0 | k |")
for l, L in enumerate(ml.levels):
    if hasattr(L, "lu"): print(f"| {l} | {L.A.shape[0]} | LU {(L.lu.L.nnz + L.lu.U.nnz)/n0:.2f} | | | | |"); continue
    print(f"| {l} | {L.A.shape[0]} | {L.A.nnz/n0:.3f} | {L.A.nnz/L.A.shape[0]:.1f} | {(L.K.nnz+L.H.nnz)/n0:.3f} | {L.S.nnz/n0:.3f} | {L.k} |")
