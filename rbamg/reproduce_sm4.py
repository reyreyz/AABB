"""Table SM4 (m = 15, k = 4, bilinear P): effect of the complement.

Complements Q of P = [I; W] (C/F order), all written as Q = Q0 + P Z, so that
  Qsh = Q0sh,  Psh = P0sh - Z Q0sh,  K_Z = K + H Z,  H_Z = H     (Prop. 3.6, S = I).
  * Q0              : Z = 0
  * canonical       : col(Q) = ker(Rhat A), Z = -(Rhat A P)^{-1} Rhat A Q0
  * LS feedback     : column j of Z supported on the coarse points interpolating
                      to fine point j, minimizing ||K e_j + H Z e_j||_2
  * pole placement  : K + H Z with prescribed eigenvalues in |mu| <= 0.05 (KNV)
Non-Markovian propagator with complement Q (4.13), Galerkin Rhat = P^T:
  E^N = (I - P'(Rhat A P')^{-1} Rhat A) Q K_Z^k Qsh,  P' = P + Q sum_{l<k} K_Z^l H.
Semi-Markovian scheme with the restriction induced by Q (Lemma 4.1(iii)):
  rows of R span the left annihilator of Ahat Q;  E^S = Q K_Z^k Qsh, rho = rho(K_Z)^k.
"""

import numpy as np
import rbamg as rb
from reproduce_sm1 import setup

mpow = np.linalg.matrix_power

PAPER = {  # problem: dict of rows (rho(K), ||K^4||, rho(E), nnz(Z)/nf); PBH value
    "Poisson": dict(PBH=1.8e-4, M=0.186, Q0=(0.648, 1.8e-1, 0.176, 0.0),
                    can=(0.647, 1.8e-1, 0.176), LS=(0.552, 9.7e-2, 0.115, 2.6), place=None),
    "CD eps=1e-2": dict(PBH=1.4e-1, M=0.107, Q0=(0.563, 2.6e-1, 0.089, 0.0),
                        can=(0.486, 1.5e-1, 0.056), LS=(0.347, 3.9e-2, 0.017, 2.6),
                        place=(0.050, 2.1e-3, 5.2e-5, 64.0), placeS=6.2e-6),
    "CD eps=1e-3": dict(PBH=6.6e-1, M=0.133, Q0=(0.388, 3.1e-1, 0.025, 0.0),
                        can=(0.464, 1.6e-1, 0.047), LS=(0.445, 8.3e-2, 0.045, 2.6),
                        place=(0.050, 3.8e-3, 5.7e-5, 64.0), placeS=6.2e-6),
}


def pbh(s):
    """min over left eigenvectors w of T of ||w P|| / ||w||."""
    _, VL = np.linalg.eig(s.T.T)          # columns: right eigvecs of T^T = left eigvecs of T
    W = VL.T
    return np.min(np.linalg.norm(W @ s.P, axis=1) / np.linalg.norm(W, axis=1))


def nonmarkov_with(s, Z):
    """rho(K_Z), ||K_Z^k||_2, rho(E^N) for the complement Q0 + P Z."""
    k = s.k
    KZ = s.K + s.H @ Z
    Q = s.Q0 + s.P @ Z
    Pp = s.P + Q @ sum(mpow(KZ, l) for l in range(k)) @ s.H
    RA = s.Rhat @ s.A
    E = (np.eye(s.n) - Pp @ np.linalg.solve(RA @ Pp, RA)) @ Q @ mpow(KZ, k) @ s.Q0sh
    return rb.rho(KZ), np.linalg.norm(mpow(KZ, k), 2), rb.rho(E)


def semi_induced(s, Z):
    """Semi-Markovian scheme with the restriction induced by Q = Q0 + P Z."""
    Q = s.Q0 + s.P @ Z
    U, sv, Vt = np.linalg.svd((s.Ahat @ Q).T)          # left annihilator of Ahat Q
    Rrows = Vt[Q.shape[1]:].conj()                      # n_c x n, R Ahat Q = 0
    Psh = np.linalg.solve(Rrows @ s.Ahat @ s.P, Rrows @ s.Ahat)
    PiR = np.eye(s.n) - s.P @ Psh
    return rb.rho(mpow(PiR @ s.T @ PiR, s.k))


def ls_feedback(s):
    nc = s.nc
    W = s.P[nc:, :]
    Z = np.zeros((nc, s.nf))
    for j in range(s.nf):
        Cj = np.nonzero(W[j])[0]
        z, *_ = np.linalg.lstsq(s.H[:, Cj], -s.K[:, j], rcond=None)
        Z[Cj, j] = z
    return Z


def place_knv(A, B, poles, sweeps=3):
    """Kautsky-Nichols-Van Dooren method 0: F with eig(A + B F) = poles.

    Returns F (real if poles come in conjugate pairs with conjugate ordering).
    Raises LinAlgError if an eigenvector subspace is empty or X is singular.
    """
    # H is numerically rank deficient here (cond ~ 1e16): place with the reduced
    # input B V_r over the numerical range and map back, F = V_r F_r.
    _, sv, Vt = np.linalg.svd(B, full_matrices=False)
    r = int(np.sum(sv > 1e-10 * sv[0]))
    Vr = Vt[:r].T
    Fr, cond = _place_full_rank(A, B @ Vr, poles, sweeps)
    return Vr @ Fr, cond


def _place_full_rank(A, B, poles, sweeps):
    n, mB = B.shape
    Qb, Rb = np.linalg.qr(B, mode="complete")
    U0, U1, Zb = Qb[:, :mB], Qb[:, mB:], Rb[:mB]
    S = []
    for lam in poles:   # S_i = null(U1^T (A - lam I))
        M = U1.conj().T @ (A - lam * np.eye(n))
        _, sv, Vh = np.linalg.svd(M)
        S.append(Vh[n - mB:].conj().T)          # dim mB (generic)
    X = np.column_stack([Si[:, 0] for Si in S]).astype(complex)
    X /= np.linalg.norm(X, axis=0)
    for _ in range(sweeps):
        for i in range(n):
            Xo = np.delete(X, i, axis=1)
            q, _ = np.linalg.qr(Xo, mode="complete")
            y = q[:, -1]                              # orthogonal to the other columns
            x = S[i] @ (S[i].conj().T @ y)
            if np.linalg.norm(x) < 1e-14:
                raise np.linalg.LinAlgError("empty projection")
            X[:, i] = x / np.linalg.norm(x)
    # enforce conjugate pairs
    for i in range(n - 1):
        if np.imag(poles[i]) > 0 and np.isclose(poles[i + 1], np.conj(poles[i])):
            X[:, i + 1] = X[:, i].conj()
    cond = np.linalg.cond(X)
    if not np.isfinite(cond) or cond > 1e14:
        raise np.linalg.LinAlgError(f"cond(X) = {cond:.1e}")
    Acl = X @ np.diag(poles) @ np.linalg.inv(X)
    F = np.linalg.solve(Zb, U0.conj().T @ (Acl - A))
    return F.real, cond


def random_poles(n, r, rng):
    """n poles in the disc |mu| <= r, conjugate pairs, at least one on |mu| = r."""
    poles = []
    while len(poles) < n:
        rad, ang = r * np.sqrt(rng.random()), np.pi * rng.random()
        if len(poles) == 0:
            rad = r
        if len(poles) <= n - 2 and rng.random() < 0.5:
            mu = rad * np.exp(1j * ang)
            poles += [mu, np.conj(mu)]
        else:
            poles.append(rad * np.sign(np.cos(ang)) + 0j)
    return np.array(poles)


def main(m=15, k=4):
    print(f"# Table SM4 (m = {m}, k = {k}, bilinear P): ours / paper\n")
    rng = np.random.default_rng(0)
    for name, pap in PAPER.items():
        s = setup(name, m, k)
        print(f"\n## {name}: PBH {pbh(s):.1e} / {pap['PBH']:.1e}, "
              f"Markovian rho(E) {rb.rho(rb.E_markov(s)):.3f} / {pap['M']:.3f}\n")
        print("| complement | rho(K) | ‖K^4‖₂ | rho(E) | nnz(Z)/nf |\n|---|---|---|---|---|")
        nc, nf = s.nc, s.nf
        r = nonmarkov_with(s, np.zeros((nc, nf)))
        print("| Q0 | " + " | ".join(f"{a:.3g} / {b:.3g}" for a, b in zip(r, pap["Q0"])) + " | 0 |")
        Zc = -np.linalg.solve(s.Rhat @ s.A @ s.P, s.Rhat @ s.A @ s.Q0)
        Kc = s.K + s.H @ Zc
        rc = (rb.rho(Kc), np.linalg.norm(mpow(Kc, k), 2), rb.rho(rb.E_semi(s)))
        print("| ker(Rhat A), S scheme | " + " | ".join(f"{a:.3g} / {b:.3g}" for a, b in zip(rc, pap["can"])) + " | – |")
        Zl = ls_feedback(s)
        rl = nonmarkov_with(s, Zl) + (np.count_nonzero(Zl) / nf,)
        print("| Q0 + PZ, LS | " + " | ".join(f"{a:.3g} / {b:.3g}" for a, b in zip(rl, pap["LS"])) + " |")
        try:
            poles = random_poles(nf, 0.05, rng)
            F, cond = place_knv(s.K, s.H, poles)
            Zp = F
            achieved = np.max(np.abs(np.linalg.eigvals(s.K + s.H @ Zp)))
            rp = nonmarkov_with(s, Zp) + (np.count_nonzero(Zp) / nf,)
            rs = semi_induced(s, Zp)
            pp = pap["place"]
            print("| Q0 + PZ, placement (cond X %.1e) | " % cond +
                  " | ".join(f"{a:.3g} / {b:.3g}" for a, b in zip(rp, pp)) + " |")
            print(f"| same, semi-Markovian | {achieved:.3g} | | {rs:.3g} / {pap['placeS']:.3g} | |")
        except (np.linalg.LinAlgError, TypeError) as e:
            print(f"| Q0 + PZ, placement | failed ({e}); paper: "
                  f"{'failed' if pap['place'] is None else pap['place']} | | | |")


if __name__ == "__main__":
    main()
