"""E^S = (Pi_R T Pi_R)^k = (I - Pi_R M A)^k Pi_R : S's Markovian twin = coarse correction, then k projected relaxation
steps x <- x + Pi_R M (b - A x), each needing one coarse solve (k+1 coarse solves per cycle vs 1 for S)."""
import numpy as np, scipy.sparse as sp
import sparse_problems as spp
from robust_n import greedy_f_points
from flows import cd as cdflow
for name, A in (("shear Pe 20", spp.cd(23, spp.SHEAR, spp.eps_for_peclet(23, spp.SHEAR, 20, "2", False), "centered")),
                ("const Pe 6", cdflow(23, "const", 6)[0])):
    A = sp.csr_matrix(A); isF = greedy_f_points(A, 0.3); Ad = A.toarray(); n = len(Ad); d = np.diag(Ad)
    F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nc = len(C)
    for w, k in ((0.7, 4), (0.8, 2)):
        M = np.diag(w / d); T = np.eye(n) - M @ Ad
        Aff, Acf = Ad[np.ix_(F, F)], Ad[np.ix_(C, F)]; Dfi = np.diag(1 / np.diag(Aff)); G = np.eye(len(F)) - Dfi @ Aff
        Rh = np.zeros((nc, n)); Rh[np.arange(nc), C] = 1; Rh[:, F] = -Acf @ (np.eye(len(F)) + G + G @ G) @ Dfi
        P0 = np.zeros((n, nc)); P0[C, np.arange(nc)] = 1
        Pi = np.eye(n) - P0 @ np.linalg.solve(Rh @ Ad @ P0, Rh @ Ad)
        ES = np.linalg.matrix_power(Pi @ T @ Pi, k); Etw = np.linalg.matrix_power(np.eye(n) - Pi @ M @ Ad, k) @ Pi
        print(f"{name}, w={w}, k={k}: ||E^S - (I - Pi M A)^k Pi|| / ||E^S|| = {np.linalg.norm(ES - Etw) / np.linalg.norm(ES):.1e}")
