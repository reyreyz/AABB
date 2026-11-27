"""Solid examples for the claims about N and its Markovian twin.

Notation: decomposition B = [P Q], duals Psh, Qsh; relaxation T = I - M A;
K = Qsh T Q, H = Qsh T P, P' = P + Q sum_{l<k} K^l H, Pi' = P'(R A P')^{-1} R A.
N      : Algorithm 4.2 on iterates (k relaxation steps, history y <- K y + H d, h = Q y,
         coarse solve R A P' eps = R(b - A x^k) - R A h, return x^k + P' eps + h)
Twin   : k steps x <- x + Q Qsh M (b - A x)  (relaxation I - M~ A, M~ = Q Qsh M),
         then the Petrov-Galerkin correction x <- x + P'(R A P')^{-1} R (b - A x)

A  equivalence (propagators and iterates), several problems, complements and relaxations
B  cost with Jacobi: sparse two-level, operation counts and wall time, Q0 and Q0 + P0 Z
C  cost with Gauss-Seidel: K is dense, N needs a matrix-free history
D  black-box smoother (iterate generator, not restartable): only N converges, O(n) per cycle
"""

import sys
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

mpow = np.linalg.matrix_power


# =========================================================================== A
def dense_alg42(A, M, P, Q, Psh, Qsh, R, k, x0, b):
    """Algorithm 4.2 on iterates (dense, any M)."""
    n = A.shape[0]
    T = np.eye(n) - M @ A
    K, H = Qsh @ T @ Q, Qsh @ T @ P
    Pp = P + Q @ sum(mpow(K, l) for l in range(k)) @ H
    xs = [x0]
    for _ in range(k):
        xs.append(xs[-1] + M @ (b - A @ xs[-1]))
    xk = xs[-1]
    y = np.zeros(Q.shape[1])
    for l in range(k):
        y = K @ y + H @ (Psh @ (xk - xs[l]))
    h = Q @ y
    eps = np.linalg.solve(R @ A @ Pp, R @ (b - A @ xk) - R @ A @ h)
    return xk + Pp @ eps + h


def dense_twin(A, M, P, Q, Psh, Qsh, R, k, x0, b):
    n = A.shape[0]
    T = np.eye(n) - M @ A
    K, H = Qsh @ T @ Q, Qsh @ T @ P
    Pp = P + Q @ sum(mpow(K, l) for l in range(k)) @ H
    x = x0.copy()
    for _ in range(k):
        x = x + Q @ (Qsh @ (M @ (b - A @ x)))                 # x <- x + Q Qsh M r
    return x + Pp @ np.linalg.solve(R @ A @ Pp, R @ (b - A @ x))


def propagators(A, M, P, Q, Psh, Qsh, R, k):
    n = A.shape[0]
    T = np.eye(n) - M @ A
    K, H = Qsh @ T @ Q, Qsh @ T @ P
    Pp = P + Q @ sum(mpow(K, l) for l in range(k)) @ H
    CC = np.eye(n) - Pp @ np.linalg.solve(R @ A @ Pp, R @ A)
    EN = CC @ Q @ mpow(K, k) @ Qsh
    EM = CC @ mpow(np.eye(n) - Q @ Qsh @ M @ A, k)
    rho = lambda X: np.max(np.abs(np.linalg.eigvals(X)))
    scale = max(np.linalg.norm(EN, 2), np.linalg.norm(EM, 2), 1.0)       # E_N can be exactly 0 (K = 0)
    return np.linalg.norm(EN - EM, 2) / scale, rho(EN), rho(T)


def section_A():
    import problems as pb
    from reproduce_sm1 import setup
    from reproduce_sm4 import ls_feedback
    import ml_cycles as mc
    import sparse_problems as spp
    print("## A. Equivalence of N and its Markovian twin (propagators and iterates)\n")
    print("| problem | complement Q | relaxation M | rho(T) | rho(E_N) | ‖E_N - E_twin‖ / max(‖E_N‖, ‖E_twin‖, 1) | ‖x_N - x_twin‖/‖x_N‖ (one cycle) |")
    print("|---|---|---|---|---|---|---|")
    rng = np.random.default_rng(0)

    def cases():
        for name in ("CD eps=1e-2", "CD eps=1e-3", "Shifted Laplacian"):
            s = setup(name, 15, 3)
            yield name + " (bilinear P, R = P^T)", s.A, s.P, s.P.T, s.nc, lambda s=s: ls_feedback(s), s
        m = 15
        Asp = spp.cd(m, spp.RECIRC, spp.eps_for_peclet(m, spp.RECIRC, 20, "2", False), "centered")
        isF = mc.greedy_f_points(Asp, 0.3)
        C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
        perm = np.concatenate([C, F])
        A = Asp.toarray()[np.ix_(perm, perm)]
        nc = len(C)
        P0 = np.vstack([np.eye(nc), np.zeros((len(F), nc))])
        yield "recirc. Pe 20 (greedy split, P0, injection)", A, P0, P0.T.copy(), nc, None, None

    for lab, A, P, R, nc, lsz, s in cases():
        n = A.shape[0]
        nf = n - nc
        W = P[nc:]
        Q0 = np.vstack([np.zeros((nc, nf)), np.eye(nf)])
        Q0sh = np.hstack([-W, np.eye(nf)])
        D = np.diag(np.diag(A))
        relax = {"Jacobi 2/3": (2 / 3) * np.linalg.inv(D), "Jacobi 1": np.linalg.inv(D),
                 "Gauss-Seidel": np.linalg.inv(np.tril(A))}
        Zs = {"Q0 (injective)": np.zeros((nc, nf))}
        if lsz is not None:
            Zs["Q0 + P Z, LS feedback"] = lsz()
        Zr = np.zeros((nc, nf))
        for j in range(nf):                                  # random sparse Z, 2 entries per column
            Zr[rng.choice(nc, 2, replace=False), j] = 0.3 * rng.standard_normal(2)
        Zs["Q0 + P Z, random sparse"] = Zr
        for qlab, Z in Zs.items():
            Q = Q0 + P @ Z
            Qsh = Q0sh
            Psh = np.hstack([np.eye(nc), np.zeros((nc, nf))]) - Z @ Q0sh
            for mlab, M in relax.items():
                k = 3
                dE, rE, rT = propagators(A, M, P, Q, Psh, Qsh, R, k)
                b, x0 = rng.standard_normal(n), rng.standard_normal(n)
                xn = dense_alg42(A, M, P, Q, Psh, Qsh, R, k, x0, b)
                xt = dense_twin(A, M, P, Q, Psh, Qsh, R, k, x0, b)
                print(f"| {lab} | {qlab} | {mlab} | {rT:.3g} | {rE:.3g} | {dE:.1e} | "
                      f"{np.linalg.norm(xn - xt) / np.linalg.norm(xn):.1e} |")
                sys.stdout.flush()


# =========================================================================== B
class TwoLevelP0:
    """Sparse two-level N and twin for P0 = [I; 0], Q = Q0 + P0 Z = [Z; I], injection,
    Jacobi omega = 1, exact coarse LU. Qsh = [0 I], Psh = [I -Z]. Operation counts in nnz(A)."""

    def __init__(self, A, isF, Z, k):
        A = sp.csr_matrix(A)
        C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
        self.A, self.C, self.F, self.k, self.nnzA = A, C, F, k, A.nnz
        self.D = A.diagonal()
        Acc, Acf = A[C][:, C].tocsr(), A[C][:, F].tocsr()
        Afc, Aff = A[F][:, C].tocsr(), A[F][:, F].tocsr()
        Dfi = sp.diags(1 / Aff.diagonal())
        T_fc, T_ff = (-Dfi @ Afc).tocsr(), (sp.eye(len(F)) - Dfi @ Aff).tocsr()
        self.Z = sp.csr_matrix(Z)
        self.K = (T_ff + T_fc @ self.Z).tocsr()               # K = Qsh T Q = T_ff + T_fc Z
        self.H = T_fc                                        # H = Qsh T P0 = T_fc
        S, V = self.H.copy(), self.H.copy()
        for _ in range(k - 1):
            V = (self.K @ V).tocsr()
            S = (S + V).tocsr()
        self.PpC = (sp.eye(len(C)) + self.Z @ S).tocsr()     # P' = P0 + Q S = [I + Z S; S]
        self.PpF = S
        self.Asig = (Acc @ self.PpC + Acf @ self.PpF).tocsr()  # (A P')_C, injection
        self.lu = sla.splu(sp.csc_matrix(self.Asig))
        self.AF = A[F].tocsr()                               # F rows (twin's residual)
        self.AC = A[C].tocsr()

    def coarse(self, rc, ops):
        ops[0] += (self.lu.L.nnz + self.lu.U.nnz) / self.nnzA
        return self.lu.solve(rc)

    def prolong(self, x, e, ops):
        x[self.C] += self.PpC @ e
        x[self.F] += self.PpF @ e
        ops[0] += (self.PpC.nnz + self.PpF.nnz) / self.nnzA
        return x

    def cycle_N(self, x0, b, ops):
        A, C, F, k = self.A, self.C, self.F, self.k
        xs = [x0]
        for _ in range(k):                                   # k full Jacobi sweeps
            xs.append(xs[-1] + (b - A @ xs[-1]) / self.D)
        xk = xs[-1]
        ops[0] += k * (A.nnz + A.shape[0]) / self.nnzA
        y = np.zeros(len(F))
        for l in range(k):                                   # history: d = Psh dx, y <- K y + H d
            dx = xk - xs[l]
            d = dx[C] - self.Z @ dx[F]
            y = self.K @ y + self.H @ d
        ops[0] += k * (self.Z.nnz + self.K.nnz + self.H.nnz) / self.nnzA
        h = np.zeros(A.shape[0])
        h[F] = y
        h[C] = self.Z @ y
        r = b - A @ (xk + h)                                 # R(b - A x^k) - R A h, injection
        ops[0] += (self.Z.nnz + A.nnz) / self.nnzA
        e = self.coarse(r[C], ops)
        return self.prolong(xk + h, e, ops)

    def cycle_twin(self, x0, b, ops):
        A, C, F, k = self.A, self.C, self.F, self.k
        x = x0.copy()
        for _ in range(k):                                   # x <- x + Q Qsh M r, Qsh = [0 I]: F rows only
            g = (b[F] - self.AF @ x) / self.D[F]
            x[F] += g
            x[C] += self.Z @ g
        ops[0] += k * (self.AF.nnz + len(F) + self.Z.nnz) / self.nnzA
        rc = b[C] - self.AC @ x                              # injection: C rows of the residual
        ops[0] += self.AC.nnz / self.nnzA
        e = self.coarse(rc, ops)
        return self.prolong(x, e, ops)

    def run(self, which, b, cycles):
        x, ops = np.zeros(self.A.shape[0]), [0.0]
        f = self.cycle_N if which == "N" else self.cycle_twin
        t = time.perf_counter()
        for _ in range(cycles):
            x = f(x, b, ops)
        return x, ops[0] / cycles, (time.perf_counter() - t) / cycles


def section_B(ms=(63, 127, 255)):
    import ml_cycles as mc
    import ml_nz as mz
    import ml_z73 as mzz
    import sparse_problems as spp
    print("\n## B. Cost with Jacobi: sparse two-level N vs twin (recirculating flow, Pe 20, greedy split, P0, injection, exact coarse LU)\n")
    print("| m | complement | k | ‖x_N - x_twin‖/‖x_N‖ after 3 cycles | rel. residual | N: WU per cycle | twin: WU per cycle | N: ms per cycle | twin: ms per cycle |")
    print("|---|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = spp.cd(m, spp.RECIRC, spp.eps_for_peclet(m, spp.RECIRC, 20, "2", False), "centered")
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for delta, qlab in ((0.3, "Q0 (twin = F-relaxation)"), (0.6, "Q0 + P0 Z (bounded, CR-optimized)")):
            isF = mc.greedy_f_points(A, delta)
            C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
            k = 4
            if delta == 0.3:
                Z = sp.csr_matrix((len(C), len(F)))
            else:
                Aff, Afc = A[F][:, F], A[F][:, C]
                Dfi = sp.diags(1 / Aff.diagonal())
                K0, H0 = (sp.eye(len(F)) - Dfi @ Aff).tocsr(), (-Dfi @ Afc).tocsr()
                patt = mzz.bounded_pattern(A, C, F, K0, H0, 4)
                Z = mz.z_kk(K0, H0, patt, k, mz.z_ls(K0, H0, patt), maxiter=100)
            tl = TwoLevelP0(A, isF, Z, k)
            xn, wn, tn = tl.run("N", b, 3)
            xt, wt, tt = tl.run("twin", b, 3)
            res = np.linalg.norm(b - A @ xn) / np.linalg.norm(b)
            print(f"| {m} | {qlab} | {k} | {np.linalg.norm(xn - xt) / np.linalg.norm(xn):.1e} | {res:.1e} | "
                  f"{wn:.2f} | {wt:.2f} | {1e3*tn:.1f} | {1e3*tt:.1f} |")
            sys.stdout.flush()


# =========================================================================== C
def section_C(ms=(15, 31)):
    from reproduce_sm1 import setup
    from reproduce_sm4 import ls_feedback
    print("\n## C. Gauss-Seidel: K = Qsh T Q is dense, N needs a matrix-free history\n")
    print("| problem | m | nnz(K)/n_f^2 (Jacobi) | nnz(K)/n_f^2 (GS) | ‖x_N - x_twin‖/‖x_N‖ | N: GS sweeps per cycle (matrix-free K, H) | twin: GS sweeps per cycle |")
    print("|---|---|---|---|---|---|---|")
    rng = np.random.default_rng(1)
    for name in ("CD eps=1e-2", "CD eps=1e-3"):
        for m in ms:
            k = 3
            s = setup(name, m, k)
            n, nc, nf = s.n, s.nc, s.nf
            Z = ls_feedback(s)
            Q = s.Q0 + s.P @ Z
            Qsh, Psh = s.Q0sh, s.P0sh - Z @ s.Q0sh
            Mgs = np.linalg.inv(np.tril(s.A))
            Kj = Qsh @ (np.eye(n) - s.M @ s.A) @ Q
            Kg = Qsh @ (np.eye(n) - Mgs @ s.A) @ Q
            dens = lambda X: np.count_nonzero(np.abs(X) > 1e-14 * np.abs(X).max()) / X.size
            b, x0 = rng.standard_normal(n), rng.standard_normal(n)
            xn = dense_alg42(s.A, Mgs, s.P, Q, Psh, Qsh, s.P.T, k, x0, b)
            xt = dense_twin(s.A, Mgs, s.P, Q, Psh, Qsh, s.P.T, k, x0, b)
            # N: k relaxation sweeps + history y <- K y + H d with K y = Qsh(Q y - M A Q y),
            #    H d = Qsh(P d - M A P d): two matrix-free GS applications per history step
            sweeps_N = k + 2 * k
            sweeps_twin = k
            print(f"| {name} | {m} | {dens(Kj):.3f} | {dens(Kg):.3f} | {np.linalg.norm(xn - xt)/np.linalg.norm(xn):.1e} | "
                  f"{sweeps_N} | {sweeps_twin} |")
            sys.stdout.flush()


# =========================================================================== D
class BlackBoxJacobi:
    """Opaque smoother: given (b, x0) it runs k undamped Jacobi steps and returns the
    iterate history [x^0, ..., x^k].  It cannot be modified or restarted mid-run; every
    call is counted."""

    def __init__(self, A, k):
        self._A, self._D, self._k = A, A.diagonal(), k
        self.calls = 0

    def run(self, b, x0):
        self.calls += 1
        xs = [x0]
        for _ in range(self._k):
            xs.append(xs[-1] + (b - self._A @ xs[-1]) / self._D)
        return xs


def section_D(ms=(127, 255, 511)):
    from robust_n import recirc, RobustN
    print("\n## D. Black-box smoother (iterate generator, not restartable): only N converges\n")
    print("| m | n | method | uses only the black-box iterates? | factor | cycles to 1e-8 | WU per cycle | note |")
    print("|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = recirc(m)
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        ml = RobustN(A)
        boxes = [BlackBoxJacobi(L.A, L.k) for L in ml.levels]

        def cycle_N(lev, rhs, x0, work):
            L = ml.levels[lev]
            if hasattr(L, "lu"):
                work[0] += L.lu_nnz / ml.nnz0
                return L.lu.solve(rhs)
            C, F = L.C, L.F
            xs = boxes[lev].run(rhs, x0)                     # the ONLY use of the smoother
            xk = xs[-1]
            y = np.zeros(len(F))
            for l in range(L.k):
                y = L.K @ y + L.H @ (xk[C] - xs[l][C])
            d0 = xk[C] - x0[C]
            xh = xk.copy()
            xh[F] += y
            r = rhs - L.A @ xh
            rc = r[C] + L.ST @ r[F] + ml.RAP(L, d0)
            e = cycle_N(lev + 1, rc, np.zeros(len(C)), work)
            if lev <= ml.refine and not hasattr(ml.levels[lev + 1], "lu"):
                e = e + cycle_N(lev + 1, rc - ml.RAP(L, e), np.zeros(len(C)), work)
            z = e - d0
            xh[C] += z
            xh[F] += L.S @ z
            work[0] += ((L.k + 1) * L.A.nnz + L.k * (L.K.nnz + L.H.nnz) + 3 * L.S.nnz
                        + 2 * (L.APc.nnz + L.APf.nnz)) / ml.nnz0
            return xh

        def cycle_M(lev, rhs, x0, work):
            """non-intrusive Markovian: black-box x^k, then x^k + P' A_c^{-1} R r."""
            L = ml.levels[lev]
            if hasattr(L, "lu"):
                work[0] += L.lu_nnz / ml.nnz0
                return L.lu.solve(rhs)
            x = boxes[lev].run(rhs, x0)[-1]
            r = rhs - L.A @ x
            rc = r[L.C] + L.ST @ r[L.F]
            e = cycle_M(lev + 1, rc, np.zeros(len(L.C)), work)
            if lev <= ml.refine and not hasattr(ml.levels[lev + 1], "lu"):
                e = e + cycle_M(lev + 1, rc - ml.RAP(L, e), np.zeros(len(L.C)), work)
            x[L.C] += e
            x[L.F] += L.S @ e
            work[0] += ((L.k + 1) * L.A.nnz + 2 * L.S.nnz) / ml.nnz0
            return x

        for lab, cyc, note in (("N (Algorithm 4.2 on the iterates)", cycle_N, ""),
                               ("Markovian, non-intrusive (x^k + coarse correction)", cycle_M, "")):
            for bx in boxes:
                bx.calls = 0
            x, nb, work, hist = np.zeros(A.shape[0]), np.linalg.norm(b), [0.0], [1.0]
            for _ in range(30):
                x = cyc(0, b, x, work)
                hist.append(np.linalg.norm(b - A @ x) / nb)
                if hist[-1] < 1e-8 or not np.isfinite(hist[-1]) or hist[-1] > 1e10:
                    break
            h = np.array(hist)
            c = len(h) - 1
            fac = (h[-1] / h[1]) ** (1 / max(c - 1, 1)) if c > 1 else h[-1]
            conv = h[-1] < 1e-8
            print(f"| {m} | {A.shape[0]} | {lab} | yes | {fac:.3g} | {c}{'' if conv else ' (no conv.)'} | "
                  f"{work[0]/c:.0f} | black-box calls on level 0: {boxes[0].calls} |")
            sys.stdout.flush()
        # the twin needs x <- x + Q Qsh M r after EVERY step, i.e. k restarts per cycle
        _, h, w = ml.solve(b, markov=(4, 0))
        c = len(h) - 1
        fac = (h[-1] / h[1]) ** (1 / max(c - 1, 1))
        print(f"| {m} | {A.shape[0]} | twin M_F(4,0) (needs k = 4 single-step restarts per cycle) | no | {fac:.3g} | "
              f"{c}{'' if h[-1] < 1e-8 else ' (no conv.)'} | {w/c:.0f} | not available under the black-box constraint |")
        sys.stdout.flush()


# =========================================================================== E
def section_E():
    """Basis change Q -> Q S (sparse S, dense S^{-1}): Q^sharp dense, but the twin uses only
    Q Q^sharp = I - P P^sharp (basis independent); N's K_S = S^{-1} K S, H_S = S^{-1} H densify."""
    from reproduce_sm1 import setup
    print("\n## E. Dense Q^sharp from a basis change Q -> Q S: hurts N, not the twin\n")
    print("| problem | m | fill Q^sharp: Q0 / Q0 S | fill K: Q0 / Q0 S | fill H: Q0 / Q0 S | ‖QQ^sharp(Q0 S) - QQ^sharp(Q0)‖ | iterates N(Q0 S) vs N(Q0) | twin(Q0 S) vs twin(Q0) |")
    print("|---|---|---|---|---|---|---|---|")
    rng = np.random.default_rng(2)
    dens = lambda X: np.count_nonzero(np.abs(X) > 1e-12 * np.abs(X).max()) / X.size
    for name in ("CD eps=1e-2", "CD eps=1e-3"):
        for m in (15, 31):
            k = 3
            s = setup(name, m, k)
            n, nc, nf = s.n, s.nc, s.nf
            S = np.eye(nf) + 0.9 * np.eye(nf, k=-1)          # sparse bidiagonal, dense inverse
            Sinv = np.linalg.inv(S)
            Q0, Q0sh, P0sh = s.Q0, s.Q0sh, s.P0sh
            Q1, Q1sh = Q0 @ S, Sinv @ Q0sh                   # Prop. 3.6: Qsh -> S^{-1} Qsh, Psh unchanged
            T = np.eye(n) - s.M @ s.A
            K0, H0 = Q0sh @ T @ Q0, Q0sh @ T @ s.P
            K1, H1 = Q1sh @ T @ Q1, Q1sh @ T @ s.P
            b, x0 = rng.standard_normal(n), rng.standard_normal(n)
            R = s.P.T
            xn0 = dense_alg42(s.A, s.M, s.P, Q0, P0sh, Q0sh, R, k, x0, b)
            xn1 = dense_alg42(s.A, s.M, s.P, Q1, P0sh, Q1sh, R, k, x0, b)
            xt0 = dense_twin(s.A, s.M, s.P, Q0, P0sh, Q0sh, R, k, x0, b)
            xt1 = dense_twin(s.A, s.M, s.P, Q1, P0sh, Q1sh, R, k, x0, b)
            print(f"| {name} | {m} | {dens(Q0sh):.3f} / {dens(Q1sh):.3f} | {dens(K0):.3f} / {dens(K1):.3f} | "
                  f"{dens(H0):.3f} / {dens(H1):.3f} | {np.linalg.norm(Q1 @ Q1sh - Q0 @ Q0sh):.1e} | "
                  f"{np.linalg.norm(xn1 - xn0)/np.linalg.norm(xn0):.1e} | {np.linalg.norm(xt1 - xt0)/np.linalg.norm(xt0):.1e} |")
            sys.stdout.flush()


# =========================================================================== F
def section_F(ms=(127, 255)):
    """Latency model: critical-path sparse-matvec rounds (one halo exchange each) vs flops.
    N pipelined:  y = S Psh x^k - u,  u <- K u + H Psh x^l computed alongside the sweeps
                  (off the critical path).  Twin: the projected update must finish before the
                  next residual (on the critical path)."""
    import ml_cycles as mc
    import ml_nz as mz
    import ml_z73 as mzz
    import sparse_problems as spp
    print("\n## F. Latency model: N (pipelined history) vs twin, recirculating flow Pe 20, P0, injection\n")
    print("| m | complement | k | pipelined N = N (iterates) | N: rounds on critical path | twin: rounds | N: flops / nnz(A) | twin: flops / nnz(A) |")
    print("|---|---|---|---|---|---|---|---|")
    for m in ms:
        A = spp.cd(m, spp.RECIRC, spp.eps_for_peclet(m, spp.RECIRC, 20, "2", False), "centered")
        b = np.random.default_rng(0).standard_normal(A.shape[0])
        for delta, qlab in ((0.3, "Q0"), (0.6, "Q0 + P0 Z")):
            isF = mc.greedy_f_points(A, delta)
            C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
            for k in (2, 4, 8):
                if delta == 0.3:
                    Z = sp.csr_matrix((len(C), len(F)))
                else:
                    Aff, Afc = A[F][:, F], A[F][:, C]
                    Dfi = sp.diags(1 / Aff.diagonal())
                    K0, H0 = (sp.eye(len(F)) - Dfi @ Aff).tocsr(), (-Dfi @ Afc).tocsr()
                    patt = mzz.bounded_pattern(A, C, F, K0, H0, 4)
                    Z = mz.z_kk(K0, H0, patt, k, mz.z_ls(K0, H0, patt), maxiter=60)
                tl = TwoLevelP0(A, isF, Z, k)
                # pipelined N on one cycle
                x0 = np.zeros(A.shape[0])
                xs = [x0]
                u = np.zeros(len(F))
                Psh = lambda v: v[C] - tl.Z @ v[F]
                for l in range(k):
                    u = tl.K @ u + tl.H @ Psh(xs[-1])        # off the critical path
                    xs.append(xs[-1] + (b - A @ xs[-1]) / tl.D)
                xk = xs[-1]
                y = tl.PpF @ Psh(xk) - u                    # S = P'_F
                h = np.zeros(A.shape[0])
                h[F], h[C] = y, tl.Z @ y
                r = b - A @ (xk + h)
                xp = tl.prolong(xk + h, tl.coarse(r[C], [0.0]), [0.0])
                xn = tl.cycle_N(x0, b, [0.0])
                same = np.linalg.norm(xp - xn) / np.linalg.norm(xn)
                zr = 1 if tl.Z.nnz else 0
                # critical path (fine level, coarse solve excluded):
                #  N:    k sweeps; then Psh x^k (zr), S (1), h_C = Z y (zr), residual (1), prolongation (1)
                #  twin: per step residual F-rows (1) + Q g (zr); then residual C-rows (1), prolongation (1)
                rN = k + zr + 1 + zr + 1 + 1
                rT = k * (1 + zr) + 1 + 1
                fN = (k * A.nnz + k * (tl.K.nnz + tl.H.nnz + tl.Z.nnz) + tl.PpF.nnz + 2 * tl.Z.nnz + A.nnz
                      + tl.PpC.nnz + tl.PpF.nnz) / A.nnz
                fT = (k * (tl.AF.nnz + tl.Z.nnz) + tl.AC.nnz + tl.PpC.nnz + tl.PpF.nnz) / A.nnz
                print(f"| {m} | {qlab} | {k} | {same:.1e} | {rN} | {rT} | {fN:.2f} | {fT:.2f} |")
                sys.stdout.flush()


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "ABCD"
    print("# Claims: N and its Markovian twin\n")
    if "A" in which:
        section_A()
    if "B" in which:
        section_B()
    if "C" in which:
        section_C()
    if "D" in which:
        section_D()
    if "E" in which:
        section_E()
    if "F" in which:
        section_F()
