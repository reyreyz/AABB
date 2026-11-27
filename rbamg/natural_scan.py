"""Natural examples, two-level dense (exact coarse solve): uniform damped Jacobi (same w on all points), standard
splits (red-black, greedy 0.3), standard truncations of S (one-point = single entry carrying the row sum; drop 0.3;
lump 0.5).  Classes: A = N converges, twin diverges; B = N converges, Markov-Jacobi diverges."""
import sys, itertools, numpy as np, scipy.sparse as sp
from robust_n import greedy_f_points, drop_rel
from ml_twins import truncate, truncate_lin
from flows import cd
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]); probs = [("const", 4), ("const", 6), ("const", 8), ("recirc", 4), ("recirc", 6), ("recirc", 10), ("recirc", 20), ("poiseuille", 10), ("poiseuille", 20)]
print(f"m={m}\n| flow | Pe | w | split | k | trunc | rho(T) | N | twin | MJ | class |\n|---|---|---|---|---|---|---|---|---|---|---|")
for (flow, pe), w, split, k, tr in itertools.product(probs, (0.6, 0.8, 1.0), ("rb", "greedy"), (1, 2, 3), ("onepoint", "drop0.3", "lump0.5", "exact")):
    A, xy = cd(m, flow, pe); n = A.shape[0]; I, J = np.divmod(np.arange(n), m)
    isF = ((I + J) % 2 == 1) if split == "rb" else greedy_f_points(A, 0.3)
    Ad = A.toarray(); d = np.diag(Ad); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    T = np.eye(n) - w * Ad / d[:, None]; K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]
    S = H.copy(); V = H.copy()
    for _ in range(k - 1): V = K @ V; S = S + V
    Ss = sp.csr_matrix(S)
    St = {"onepoint": lambda: truncate_lin(Ss, np.ones((nc, 1)), 1), "drop0.3": lambda: drop_rel(Ss, 0.3),
          "lump0.5": lambda: truncate(Ss, 0.5, "lump"), "exact": lambda: Ss}[tr]().toarray()
    Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St; Ac = Pt.T @ Ad @ Pt
    if np.linalg.cond(Ac) > 1e12: continue
    I_Pi = np.eye(n) - Pt @ np.linalg.solve(Ac, Pt.T @ Ad); I_PiQ = I_Pi[:, F]
    Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1; Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
    Tk = np.linalg.matrix_power(T, k); Kk = np.linalg.matrix_power(K, k)
    rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh)); rT = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh)); rM = rho(I_Pi @ Tk)
    cls = ("A" if rT > 1 else "") + ("B" if rM > 1 else "") if rN < 0.7 else ""
    if cls:
        print(f"| {flow} | {pe} | {w} | {split} | {k} | {tr} | {rho(T):.2f} | {rN:.3g} | {rT:.3g} | {rM:.3g} | {cls} |"); sys.stdout.flush()
