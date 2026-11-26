"""Table SM5: flow of transfer operators (random diffusion, m = 15).

Flow with stationary dual (Lemma 5.2), coordinates (P0, Q0), P_tau = P0 + Q0 X:
    X_{tau+1} = X_tau + S_k(K_tau) R(X_tau),   K_tau = T_ff - X_tau T_cf,
    R(X) = T_fc + T_ff X - X T_cc - X T_cf X,  S_k(K) = sum_{l<k} K^l,
    S_inf(K) = (I - K)^{-1}.
Fixed point X*: invariant subspace of the n_c eigenvalues of Ahat of smallest modulus.
Theorem 5.5: local rate gamma_k = max_ij |g_k(mu_i, lambda_j)|,
    g_k(mu, lam) = mu^k + (1 - lam) (1 - mu^k)/(1 - mu),
lambda_j = eig(T*^(0)) = eig(T_cc + T_cf X*), mu_i = eig(K*).
"""

import numpy as np
import problems as pb

mpow = np.linalg.matrix_power


def blocks(m, field_seed=0):
    A = pb.random_diffusion(m, seed=field_seed)
    perm, _ = pb.bilinear_cf(m)
    A = A[np.ix_(perm, perm)]
    nc = ((m + 1) // 2) ** 2
    T = np.eye(len(A)) - pb.jacobi(A) @ A
    return A, T, nc


def fixed_point(T, nc):
    Ahat = np.eye(len(T)) - T
    alpha, V = np.linalg.eig(Ahat)
    order = np.argsort(np.abs(alpha))
    Vs = V[:, order[:nc]]
    X = Vs[nc:] @ np.linalg.inv(Vs[:nc])        # X* = Q0^sharp V (P0^sharp V)^{-1}
    gap = abs(alpha[order[nc - 1]]) / abs(alpha[order[nc]])
    return X.real, gap, np.linalg.cond(V)


def step(X, T, nc, k):
    Tcc, Tcf, Tfc, Tff = T[:nc, :nc], T[:nc, nc:], T[nc:, :nc], T[nc:, nc:]
    K = Tff - X @ Tcf
    R = Tfc + Tff @ X - X @ Tcc - X @ Tcf @ X
    if k == np.inf:
        return X + np.linalg.solve(np.eye(len(K)) - K, R)
    S = sum(mpow(K, l) for l in range(k))
    return X + S @ R


def gamma(T, nc, X, k):
    Tcc, Tcf, Tff = T[:nc, :nc], T[:nc, nc:], T[nc:, nc:]
    lam = np.linalg.eigvals(Tcc + Tcf @ X)
    mu = np.linalg.eigvals(Tff - X @ Tcf)
    L, Mu = np.meshgrid(lam, mu)
    if k == np.inf:
        g = (1 - L) / (1 - Mu)
    else:
        g = Mu ** k + (1 - L) * (1 - Mu ** k) / (1 - Mu)
    return np.abs(g).max()


def observed_rate(T, nc, Xs, k, steps=1500, pert=2e-4, seed=1):
    rng = np.random.default_rng(seed)
    D = rng.standard_normal(Xs.shape)
    X = Xs + pert * D / np.linalg.norm(D, 2)
    errs = []
    for _ in range(steps):
        X = step(X, T, nc, k)
        errs.append(np.linalg.norm(X - Xs, 2))
    e = np.array(errs)
    # local rate: least-squares slope of log e over the stretch between the
    # initial transient (e below the start value) and the roundoff floor
    floor = e.min()
    ok = np.nonzero((e > 1e3 * floor) & (e < 0.5 * pert))[0]
    ok = ok[len(ok) // 4:]
    slope = np.polyfit(ok, np.log(e[ok]), 1)[0]
    return np.exp(slope), e


def global_run(T, nc, Xs, X0, k, steps=1500):
    X = X0.copy()
    for _ in range(steps):
        X = step(X, T, nc, k)
        if not np.all(np.isfinite(X)) or np.linalg.norm(X) > 1e8:
            return np.inf
    return np.linalg.norm(X - Xs, 2)


PAPER = {1: (0.9955, 0.9947, "diverges", 3.2e1), 2: (0.9927, 0.9922, "diverges", 4.6e-1),
         4: (0.9898, 0.9894, "diverges", 5.8e-3), np.inf: (0.9879, 0.9886, 1.9e-4, 3.1e-4)}


def main(m=15):
    A, T, nc = blocks(m)
    Xs, gap, condV = fixed_point(T, nc)
    _, Pb = pb.bilinear_cf(m)
    Wbil = Pb[nc:]
    print(f"# Table SM5 (m = {m}, stand-in random field): ours / paper\n")
    print(f"gap |alpha_64/alpha_65| = {gap:.4f} (paper 0.9879), cond(V) = {condV:.1e}\n")
    print("| k | gamma_k (Thm 5.5) | observed local | bilinear start | P0 start |")
    print("|---|---|---|---|---|")
    for k in (1, 2, 4, np.inf):
        g = gamma(T, nc, Xs, k)
        r, _ = observed_rate(T, nc, Xs, k)
        eb = global_run(T, nc, Xs, Wbil, k)
        e0 = global_run(T, nc, Xs, np.zeros_like(Xs), k)
        pg, pr, pb_, p0 = PAPER[k]
        fmt = lambda v: "diverges" if v == np.inf else f"{v:.1e}"
        pf = lambda v: v if isinstance(v, str) else f"{v:.1e}"
        print(f"| {k} | {g:.4f} / {pg} | {r:.4f} / {pr} | {fmt(eb)} / {pf(pb_)} | {fmt(e0)} / {pf(p0)} |")
    print("\nTheorem 5.5 check is the agreement of columns 2 and 3 (independent of the field).")


if __name__ == "__main__":
    main()
