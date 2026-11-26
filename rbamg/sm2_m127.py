"""Table SM2 at m = 127 (n = 16129), CD eps = 1e-3, bilinear P: matrix-free.

A x is applied with the 5-point upwind stencil on the grid, P and P' are kept
dense (n x n_c), coarse solves use the inverse of the n_c x n_c coarse matrix,
and rho(E) is estimated from the Ritz values of an Arnoldi run (dimension `kdim`;
the estimate is reported for two dimensions to show convergence).

The C/F ordering of problems.bilinear_cf is used, so every formula is the one
of rbamg.py; only the way the operators are applied changes.
"""

import sys
import numpy as np
import problems as pb


def make(m, eps=1e-3, theta=np.pi / 6):
    h = 1.0 / (m + 1)
    bx, by = np.cos(theta) * h, np.sin(theta) * h
    diag = 4 * eps + bx + by
    perm, P = pb.bilinear_cf(m)
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))

    def A_nat(U):                      # U: (m, m, c) natural grid, x = second index
        V = diag * U
        V[:, 1:] -= (eps + bx) * U[:, :-1]     # west (upwind)
        V[:, :-1] -= eps * U[:, 1:]            # east
        V[1:, :] -= (eps + by) * U[:-1, :]     # south (upwind)
        V[:-1, :] -= eps * U[1:, :]            # north
        return V

    def A(X):                          # X in C/F order, shape (n,) or (n, c)
        X2 = X.reshape(len(perm), -1)
        U = np.empty_like(X2)
        U[perm] = X2
        V = A_nat(U.reshape(m, m, -1)).reshape(len(perm), -1)[perm]
        return V.reshape(X.shape)

    return A, diag, P


def run(m=127, kdim=(200, 300), ks=(2, 4)):
    A, diag, P = make(m)
    n, nc = P.shape
    omega = 2.0 / 3.0
    T = lambda X: X - (omega / diag) * A(X)             # T = I - M A, M = omega D^{-1}
    W = P[nc:]
    Q0 = lambda Y: np.concatenate([np.zeros((nc,) + Y.shape[1:]), Y])
    Q0sh = lambda X: X[nc:] - W @ X[:nc]
    K = lambda Y: Q0sh(T(Q0(Y)))
    TP = T(P)
    H = Q0sh(TP)                                        # n_f x n_c
    AP = A(P)
    Ac = P.T @ AP
    Acinv = np.linalg.inv(Ac)
    PiR = lambda X: X - P @ (Acinv @ (P.T @ A(X)))      # Pi_R = I - P (P^T A P)^{-1} P^T A
    out = {}
    for k in ks:
        Y = H.copy()
        S = H.copy()
        for _ in range(k - 1):                          # sum_{l<k} K^l H
            Y = K(Y)
            S += Y
        Pp = P + Q0(S)
        APp = A(Pp)
        Acp_inv = np.linalg.inv(P.T @ APp)
        CPp = lambda X: X - Pp @ (Acp_inv @ (P.T @ A(X)))

        def Tk(x):
            for _ in range(k):
                x = T(x)
            return x

        def Kk(y):
            for _ in range(k):
                y = K(y)
            return y

        ops = {
            "M": lambda x: PiR(Tk(x)),
            "M'": lambda x: CPp(Tk(x)),
            "S": lambda x: _pow(lambda z: PiR(T(PiR(z))), k, x),
            "N": lambda x: CPp(Q0(Kk(Q0sh(x)))),
        }
        for name, op in ops.items():
            out[(k, name)] = [arnoldi_rho(op, n, d) for d in kdim]
            print(k, name, out[(k, name)])
            sys.stdout.flush()
    return out


def _pow(f, k, x):
    for _ in range(k):
        x = f(x)
    return x


def arnoldi_rho(op, n, kdim, seed=0):
    rng = np.random.default_rng(seed)
    V = np.zeros((kdim + 1, n))
    Hm = np.zeros((kdim + 1, kdim))
    v = rng.standard_normal(n)
    V[0] = v / np.linalg.norm(v)
    for j in range(kdim):
        w = op(V[j])
        for _ in range(2):                              # reorthogonalized Gram-Schmidt
            c = V[: j + 1] @ w
            w -= c @ V[: j + 1]
            Hm[: j + 1, j] += c
        Hm[j + 1, j] = np.linalg.norm(w)
        if Hm[j + 1, j] < 1e-14:
            kdim = j + 1
            break
        V[j + 1] = w / Hm[j + 1, j]
    return float(np.max(np.abs(np.linalg.eigvals(Hm[:kdim, :kdim]))))


PAPER = {(2, "M"): 0.275, (2, "M'"): 0.345, (2, "S"): 0.313, (2, "N"): 0.368,
         (4, "M"): 0.186, (4, "M'"): 0.087, (4, "S"): 0.096, (4, "N"): 0.146}

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "check":
        # validate the matrix-free path against the dense code at m = 31
        import rbamg as rb
        from reproduce_sm1 import setup
        res = run(31, kdim=(150, 250), ks=(2,))
        s = setup("CD eps=1e-3", 31, 2)
        print("dense:", rb.rho(rb.E_markov(s)), rb.rho(rb.E_markov_pp(s)),
              rb.rho(rb.E_semi(s)), rb.rho(rb.E_nonmarkov(s)))
        sys.exit()
    res = run()
    print("\n## Table SM2, CD eps=1e-3, m = 127 (Arnoldi dim 200 / 300): ours / paper\n")
    print("| k | scheme | ours (200) | ours (300) | paper |\n|---|---|---|---|---|")
    for (k, name), v in res.items():
        print(f"| {k} | {name} | {v[0]:.3f} | {v[1]:.3f} | {PAPER[(k, name)]:.3f} |")
