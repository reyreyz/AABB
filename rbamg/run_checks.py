"""Checks of the relaxation-based AMG theory against an independent implementation.

  1. Table SM9: algorithms vs closed-form propagators (m = 15, k = 3, bilinear P),
     and the exact scheme after one cycle.
  2. Structural statements: (4.4) splitting, (4.8) rho(E^S) = rho(K_can)^k,
     k = 1 spectra of S and M, Proposition 4.4 (Hermitian case).
  3. Table SM1, Poisson rows (m = 31, bilinear P, damped Jacobi 2/3, Galerkin),
     compared with the published values.
"""

import sys
import numpy as np
import rbamg as rb
import problems as pb

np.set_printoptions(precision=3)


def build(name, m, k):
    A = pb.PROBLEMS[name](m)
    perm, P = pb.bilinear_cf(m)
    A = A[np.ix_(perm, perm)]                 # C/F ordering
    M = pb.jacobi(A)
    return rb.Setup(A, M, P, P.T, k)          # Galerkin: Rhat = P^T


def alg_matrix(alg, s):
    """Error propagator of an algorithm: b = 0, x0 = -I  =>  x_out = -E."""
    return -alg(s, -np.eye(s.n), np.zeros((s.n, s.n)))


def reldev(X, Y):
    return np.linalg.norm(X - Y, 2) / np.linalg.norm(Y, 2)


def table_sm9(m=15, k=3):
    print(f"\n## 1. Algorithms vs closed forms (m={m}, k={k}, bilinear P)\n")
    print("| problem | M | M' | S | N | exact: ||e_out||/||e0|| |")
    print("|---|---|---|---|---|---|")
    worst = 0.0
    for name in pb.PROBLEMS:
        s = build(name, m, k)
        devs = [
            reldev(alg_matrix(rb.alg_markov, s), rb.E_markov(s)),
            reldev(alg_matrix(lambda s, X, B: rb.alg_markov(s, X, B, s.Pp), s), rb.E_markov_pp(s)),
            reldev(alg_matrix(rb.alg_semi, s), rb.E_semi(s)),
            reldev(alg_matrix(rb.alg_nonmarkov, s), rb.E_nonmarkov(s)),
        ]
        ex = np.linalg.norm(alg_matrix(rb.alg_exact, s), 2)   # ||E_exact||_2, should be ~eps
        worst = max(worst, *devs, ex)
        print(f"| {name} | " + " | ".join(f"{d:.1e}" for d in devs) + f" | {ex:.1e} |")
    return worst


def structural(m=15):
    print(f"\n## 2. Structural statements (m={m})\n")
    print("| problem | k | (4.4) split | rho(E^S) vs rho(K_can)^k | k=1: rho(E^S)-rho(E^M) | rho(E^S) <= rho(E^M) |")
    print("|---|---|---|---|---|---|")
    worst = 0.0
    for name in pb.PROBLEMS:
        for k in (1, 2, 4):
            s = build(name, m, k)
            # canonical complement Q = Pi_R Q0, duals from [P Q]^{-1}
            Q = s.PiR @ s.Q0
            Binv = np.linalg.inv(np.hstack([s.P, Q]))
            Psh, Qsh = Binv[: s.nc], Binv[s.nc:]
            K, H = s.blocks(Q, Qsh)
            # (4.4): E^M = sum_{l=1}^k P^(l) Psh T^{k-l} + Q K^k Qsh,  P^(l) = Q K^{l-1} H
            EM = rb.E_markov(s)
            split = sum(Q @ rb.mpow(K, l - 1) @ H @ Psh @ rb.mpow(s.T, k - l)
                        for l in range(1, k + 1)) + Q @ rb.mpow(K, k) @ Qsh
            d44 = reldev(split, EM)
            rS, rM, rK = rb.rho(rb.E_semi(s)), rb.rho(EM), rb.rho(K) ** k
            d48 = abs(rS - rK) / rK
            worst = max(worst, d44, d48)
            k1 = f"{rS - rM:.1e}" if k == 1 else ""
            if k == 1:
                worst = max(worst, abs(rS - rM))
            p44 = ""
            if name == "Poisson" and k % 2 == 0:          # Hermitian, Galerkin: Prop. 4.4
                p44 = f"{rS:.4f} <= {rM:.4f}: {rS <= rM + 1e-12}"
                worst = max(worst, 0.0 if rS <= rM + 1e-12 else 1.0)
            print(f"| {name} | {k} | {d44:.1e} | {d48:.1e} | {k1} | {p44} |")
    return worst


PAPER_SM1_POISSON = {      # Table SM1, Poisson, bilinear P, m = 31: M, M', S, N, rho(K)^k
    1: (0.662, 0.661, 0.662, 0.662, 0.662),
    2: (0.439, 0.437, 0.438, 0.438, 0.438),
    4: (0.195, 0.192, 0.192, 0.192, 0.192),
}


def table_sm1_poisson(m=31):
    print(f"\n## 3. Table SM1, Poisson, bilinear P, m={m}\n")
    print("| k | | M | M' | S | N | rho(K)^k |")
    print("|---|---|---|---|---|---|---|")
    worst = 0.0
    for k in (1, 2, 4):
        s = build("Poisson", m, k)
        ours = (rb.rho(rb.E_markov(s)), rb.rho(rb.E_markov_pp(s)), rb.rho(rb.E_semi(s)),
                rb.rho(rb.E_nonmarkov(s)), rb.rho(s.K) ** k)
        paper = PAPER_SM1_POISSON[k]
        worst = max(worst, max(abs(a - b) for a, b in zip(ours, paper)))
        print(f"| {k} | ours | " + " | ".join(f"{v:.3f}" for v in ours) + " |")
        print(f"| {k} | paper | " + " | ".join(f"{v:.3f}" for v in paper) + " |")
    return worst


if __name__ == "__main__":
    print("Problems marked * use our own discretization (not specified in the paper).")
    w1 = table_sm9()
    w2 = structural()
    w3 = table_sm1_poisson()
    print(f"\nworst identity deviation (1, 2): {max(w1, w2):.1e}")
    print(f"worst |ours - paper| on Table SM1 Poisson: {w3:.4f}")
    ok = max(w1, w2) < 1e-10 and w3 < 1.5e-3
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
