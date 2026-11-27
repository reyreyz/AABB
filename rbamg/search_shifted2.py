"""Search with shifted (indefinite) operators, heavy truncation, R = P~^T or injection."""
import numpy as np, scipy.sparse as sp, itertools, sys
from robust_n import recirc, greedy_f_points
from ml_twins import truncate
rho = lambda X: np.abs(np.linalg.eigvals(X)).max()
m = int(sys.argv[1]) if len(sys.argv) > 1 else 15
h = 1.0 / (m + 1)
L2 = (sp.kron(sp.eye(m), sp.diags([-1, 2, -1], [-1, 0, 1], (m, m))) + sp.kron(sp.diags([-1, 2, -1], [-1, 0, 1], (m, m)), sp.eye(m))).tocsr()
probs = {}
for kap2 in (30, 60, 120):
    probs[f"shifted Lap k2={kap2}"] = (L2 - kap2 * h * h * sp.eye(m * m)).tocsr()
    for pe in (1, 2, 5):
        Rc = recirc(m, pe); Rc = Rc / Rc.diagonal().max() * 4
        probs[f"recirc Pe {pe} shifted k2={kap2}"] = (Rc - kap2 * h * h * sp.eye(m * m)).tocsr()
hits, near = [], []
for (pn, A), w, delta, Rk in itertools.product(probs.items(), (0.6, 0.8, 1.0), (0.3, 0.6), ("PtT", "inj")):
    Ad = A.toarray(); n = len(Ad); T = np.eye(n) - w * Ad / np.diag(Ad)[:, None]
    isF = greedy_f_points(A, delta); F, C = np.nonzero(isF)[0], np.nonzero(~isF)[0]; nf, nc = len(F), len(C)
    K, H = T[np.ix_(F, F)], T[np.ix_(F, C)]; rT, rK = rho(T), rho(K)
    if rK >= 1: continue
    Psh = np.zeros((nc, n)); Psh[np.arange(nc), C] = 1
    Qsh = np.zeros((nf, n)); Qsh[np.arange(nf), F] = 1
    S, V, Tk, Kk = H.copy(), H.copy(), T.copy(), K.copy()
    for k in range(1, 5):
        if k > 1: V = K @ V; S = S + V; Tk = T @ Tk; Kk = K @ Kk
        for tau in (0.6, 1.0):
            St = truncate(sp.csr_matrix(S), tau, "lump").toarray()
            Pt = np.zeros((n, nc)); Pt[C, np.arange(nc)] = 1; Pt[F] = St
            R = Pt.T if Rk == "PtT" else Psh
            Ac = R @ Ad @ Pt
            if abs(np.linalg.det(Ac)) == 0 or np.linalg.cond(Ac) > 1e12: continue
            I_Pi = np.eye(n) - Pt @ np.linalg.solve(Ac, R @ Ad); I_PiQ = I_Pi[:, F]
            rN = rho(I_PiQ @ ((S - St) @ Psh @ Tk + Kk @ Qsh))
            rTw = rho(I_PiQ @ ((S - St) @ Psh + Kk @ Qsh))
            rM = rho(I_Pi @ Tk)
            row = (pn, f"Jacobi w={w}", delta, Rk, k, f"{tau:g} lump", rT, rK, rN, rTw, rM)
            (hits if (rN < 1 and rTw > 1 and rM > 1) else near if rN < 1 and (rTw > 1 or rM > 1) else []).append(row)
    sys.stdout.flush()
fmt = lambda r: "| " + " | ".join(str(x) if not isinstance(x, float) else f"{x:.3g}" for x in r) + " |"
hdr = "| problem | smoother | delta | R | k | trunc | rho(T) | rho(K) | N | twin | Markov-Jacobi |\n|" + "---|" * 11
print(f"m={m}\n\n## Hits ({len(hits)}), sorted by max(N, 1/twin, 1/MJ)\n\n{hdr}")
for r in sorted(hits, key=lambda r: max(r[8], 1 / r[9], 1 / r[10])): print(fmt(r))
