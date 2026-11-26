"""Relaxation-based two-level AMG (Moussa & Kahl), dense reference implementation.

Notation follows the paper:
    A x = b,  relaxation  x <- x + M (b - A x),  Ahat = M A,  T = I - Ahat.
    Decomposition B = [P Q],  B^{-1} = [Psh; Qsh]  (Psh = P^sharp, Qsh = Q^sharp).
    K = Qsh T Q  (CR propagator),  H = Qsh T P  (coarse-to-fine coupling).
    Rhat = R M; every coarse operator R Ahat P' is written Rhat A P'.

Two kinds of objects are provided for every scheme:
  * closed-form error propagators  E  (Theorems 4.2, 4.3, 4.5, 4.11);
  * the algorithms themselves (Algorithms 4.1, 4.2, the exact scheme), run on
    iterates.  They are affine in (x0, b); with b = 0 and x0 = -I (all unit
    vectors at once) they return -E, which is how the two are compared.

All matrices are dense numpy arrays in C/F ordering (coarse points first).
"""

import numpy as np

solve = np.linalg.solve
mpow = np.linalg.matrix_power


def rho(E):
    """Spectral radius."""
    return np.max(np.abs(np.linalg.eigvals(E)))


# ---------------------------------------------------------------------------
# Decomposition objects
# ---------------------------------------------------------------------------
class Setup:
    """A, M, k, P = [I; W] in C/F order, Rhat, and the local complement Q0."""

    def __init__(self, A, M, P, Rhat, k):
        n, nc = P.shape
        self.A, self.M, self.P, self.Rhat, self.k = A, M, P, Rhat, k
        self.n, self.nc, self.nf = n, nc, n - nc
        I = np.eye(n)
        self.Ahat = M @ A                                  # Ahat = M A
        self.T = I - self.Ahat                             # T = I - Ahat   (1.5)

        # local complement Q0 = [0; I], Psh = [I 0], Qsh = [-W I]   (3.11), (5.1)
        W = P[nc:, :]
        self.Q0 = np.vstack([np.zeros((nc, self.nf)), np.eye(self.nf)])
        self.P0sh = np.hstack([np.eye(nc), np.zeros((nc, self.nf))])
        self.Q0sh = np.hstack([-W, np.eye(self.nf)])
        self.K, self.H = self.blocks(self.Q0, self.Q0sh)   # K = T_ff - W T_cf

        # effective prolongation P' = P + Q sum_{l<k} K^l H            (4.11)
        S = sum(mpow(self.K, l) for l in range(k))
        self.Pp = P + self.Q0 @ S @ self.H

        # canonical complement induced by Rhat: Pi_R = Q Q^sharp      (4.1)
        RA = Rhat @ A                                      # = R Ahat
        self.Psh_can = solve(RA @ P, RA)                   # P^sharp = (R Ahat P)^{-1} R Ahat
        self.PiR = I - P @ self.Psh_can                    # Pi_R = I - P (R Ahat P)^{-1} R Ahat

    def blocks(self, Q, Qsh):
        """K = Qsh T Q, H = Qsh T P   (3.5)."""
        return Qsh @ self.T @ Q, Qsh @ self.T @ self.P


# ---------------------------------------------------------------------------
# Closed-form propagators
# ---------------------------------------------------------------------------
def coarse_proj(s, Pc):
    """I - Pc (Rhat A Pc)^{-1} Rhat A."""
    RA = s.Rhat @ s.A
    return np.eye(s.n) - Pc @ solve(RA @ Pc, RA)


def E_markov(s):          # (4.3)  E^M = Pi_R T^k
    return coarse_proj(s, s.P) @ mpow(s.T, s.k)


def E_markov_pp(s):       # Markovian with P' in place of P (scheme M')
    return coarse_proj(s, s.Pp) @ mpow(s.T, s.k)


def E_semi(s):            # (4.8)  E^S = Q K^k Q^sharp = (Pi_R T Pi_R)^k, canonical complement
    return mpow(s.PiR @ s.T @ s.PiR, s.k)


def E_nonmarkov(s):       # (4.13) E^N = (I - P'(R Ahat P')^{-1} R Ahat) Q K^k Q^sharp
    return coarse_proj(s, s.Pp) @ s.Q0 @ mpow(s.K, s.k) @ s.Q0sh


# ---------------------------------------------------------------------------
# Algorithms (columns of X are independent iterates)
# ---------------------------------------------------------------------------
def relax(s, X, B):
    """x^(0..k+1) from (1.2); returns the list of iterates."""
    xs = [X]
    for _ in range(s.k + 1):
        xs.append(xs[-1] + s.M @ (B - s.A @ xs[-1]))       # x <- x + M (b - A x)
    return xs


def alg_markov(s, X, B, Pc=None):
    """Petrov-Galerkin two-grid: Rhat A Pc eps = R (x^(k+1) - x^(k)) = Rhat r^(k)."""
    Pc = s.P if Pc is None else Pc
    xs = relax(s, X, B)
    xk = xs[s.k]
    rhat = s.Rhat @ (B - s.A @ xk)                         # R M r^(k)
    eps = solve(s.Rhat @ s.A @ Pc, rhat)                   # (4.2)
    return xk + Pc @ eps


def alg_semi(s, X, B):
    """Algorithm 4.1 (semi-Markovian, canonical complement)."""
    k, P, Psh = s.k, s.P, s.Psh_can
    xs = relax(s, X, B)
    xk = xs[k]
    eps_k = solve(s.Rhat @ s.A @ P, s.Rhat @ (B - s.A @ xk))   # line 4, exact coarse solve
    eps = [eps_k + Psh @ (xk - xs[l]) for l in range(k)]        # line 5, shift relation
    eps.append(eps_k)
    QT = s.PiR @ s.T                                             # 𝒬 T,  𝒬 = I - P P^sharp
    z = P @ eps[0]                                               # line 6
    for l in range(1, k + 1):                                    # lines 7-9, Horner replay
        z = P @ eps[l] + QT @ z
    return xk + z


def alg_nonmarkov(s, X, B):
    """Algorithm 4.2 (non-Markovian, local complement Q0, injection P^sharp = [I 0])."""
    k = s.k
    xs = relax(s, X, B)
    xk = xs[k]
    d = [s.P0sh @ (xk - xs[l]) for l in range(k)]               # line 4
    y = np.zeros((s.nf, X.shape[1]))
    for l in range(k):                                           # lines 6-8
        y = s.K @ y + s.H @ d[l]
    h = s.Q0 @ y                                                 # line 9, history correction
    rhs = s.Rhat @ (B - s.A @ xk) - s.Rhat @ s.A @ h             # line 10
    eps = solve(s.Rhat @ s.A @ s.Pp, rhs)
    return xk + s.Pp @ eps + h                                   # line 11


def alg_exact(s, X, B):
    """Theorem 4.11 with the local complement Q0: exact in one cycle."""
    k, P, Q, Qsh, K, H = s.k, s.P, s.Q0, s.Q0sh, s.K, s.H
    xs = relax(s, X, B)
    xk = xs[k]
    I_Kk = np.eye(s.nf) - mpow(K, k)                             # I - K^k
    Wt = [solve(I_Kk, mpow(K, l - 1) @ H) for l in range(1, k + 1)]   # W~(l), l = 1..k
    d = [s.P0sh @ (xk - xs[l]) for l in range(k + 1)]            # d^(l), d^(k) = 0
    xi = solve(I_Kk, mpow(K, k) @ (s.Q0sh @ (xk - xs[0])))       # xi^(k)
    QAQ = Qsh @ s.Ahat @ Q                                       # = I - K
    Pt = P - Q @ solve(QAQ, Qsh @ s.Ahat @ P)                    # P~ (4.18)
    RAQ = s.Rhat @ s.A @ Q
    st = -sum(RAQ @ Wt[l - 1] @ d[k - l] for l in range(1, k + 1)) - RAQ @ xi
    e_s = solve(s.Rhat @ s.A @ Pt, s.Rhat @ (B - s.A @ xk) + st)        # (4.16)
    e_f = sum(Wt[l - 1] @ (e_s + d[k - l]) for l in range(1, k + 1)) + xi  # (4.15)
    return xk + P @ e_s + Q @ e_f
