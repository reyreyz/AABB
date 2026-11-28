import numpy as np, scipy.sparse as sp
import sparse_problems as spp
from robust_n import greedy_f_points
A = sp.csr_matrix(spp.cd(23, spp.SHEAR, spp.eps_for_peclet(23, spp.SHEAR, 20, "2", False), "centered"))
isF = greedy_f_points(A, 0.3); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nc, nf = len(C), len(F)
for w, k in ((0.7, 4), (0.3, 2), (0.05, 16)):
    T = np.eye(n) - w * Ad / d[:, None]
    Aff, Acf = Ad[np.ix_(F, F)], Ad[np.ix_(C, F)]; Dfi = np.diag(1 / np.diag(Aff)); G = np.eye(nf) - Dfi @ Aff
    Rh = np.zeros((nc, n)); Rh[np.arange(nc), C] = 1; Rh[:, F] = -Acf @ (np.eye(nf) + G + G @ G) @ Dfi
    P0 = np.zeros((n, nc)); P0[C, np.arange(nc)] = 1; Pi = np.eye(n) - P0 @ np.linalg.solve(Rh @ Ad @ P0, Rh @ Ad)
    K0, H0 = T[np.ix_(F, F)], T[np.ix_(F, C)]; Pp = P0.copy(); Pp[F] = sum(np.linalg.matrix_power(K0, l) for l in range(k)) @ H0
    Tk = np.linalg.matrix_power(T, k)
    EMp = (np.eye(n) - Pp @ np.linalg.solve(Rh @ Ad @ Pp, Rh @ Ad)) @ Tk
    Rep = sum(np.linalg.matrix_power(Pi @ T, j) for j in range(k + 1)) @ P0
    ENS = (np.eye(n) - Rep @ np.linalg.solve(Rh @ Ad @ P0, Rh @ Ad)) @ Tk
    print(f"w={w}, k={k}: ||E_noshift - E_M'|| / ||E_M'|| = {np.linalg.norm(ENS - EMp) / np.linalg.norm(EMp):.1e}")
