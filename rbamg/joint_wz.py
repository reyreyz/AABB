"""Joint local (W, Z): CR-centric optimization of the whole decomposition.

Decomposition (C/F order): P = [I; W], Q = Q0 + P Z, Qsh = [-W I], Psh = [I 0] - Z Qsh.
    K_{W,Z} = Qsh T Q = (T_ff - W T_cf) + H_W Z,
    H_W     = Qsh T P = T_fc + T_ff W - W T_cc - W T_cf W.
W (interpolation, output injection on rows) and Z (feedback on columns) both on
bounded local patterns.  Objective: ||K_{W,Z}^k||_F^2 (or the rho surrogate),
L-BFGS with JAX gradients.  Two-level error propagator (4.13) with restriction
'inj' (Rhat = [I 0]) or 'gal' (Rhat = P^T):
    P' = P + Q sum_{l<k} K^l H_W,   E_N = (I - P'(Rhat A P')^{-1} Rhat A) Q K^k Qsh.
"""

import sys
import numpy as np
import jax
import jax.numpy as jnp
from scipy.optimize import minimize
import rbamg as rb

jax.config.update("jax_enable_x64", True)


class Joint:
    def __init__(self, A, nc, k, patW, patZ, omega=1.0, restrict="inj"):
        self.A = jnp.array(A)
        n = A.shape[0]
        self.n, self.nc, self.nf, self.k = n, nc, n - nc, k
        T = np.eye(n) - omega * A / np.diag(A)[:, None]
        self.T = jnp.array(T)
        self.Tcc, self.Tcf = self.T[:nc, :nc], self.T[:nc, nc:]
        self.Tfc, self.Tff = self.T[nc:, :nc], self.T[nc:, nc:]
        self.wr = np.concatenate([np.full(len(p), j) for j, p in enumerate(patW)])   # rows of W (fine)
        self.wc = np.concatenate(patW)                                                # cols of W (coarse)
        self.zr = np.concatenate(patZ)                                                # rows of Z (coarse)
        self.zc = np.concatenate([np.full(len(p), j) for j, p in enumerate(patZ)])   # cols of Z (fine)
        self.nw, self.nz = len(self.wr), len(self.zr)
        self.restrict = restrict

    def unpack(self, x):
        W = jnp.zeros((self.nf, self.nc)).at[self.wr, self.wc].set(x[:self.nw])
        Z = jnp.zeros((self.nc, self.nf)).at[self.zr, self.zc].set(x[self.nw:])
        return W, Z

    def KH(self, x):
        W, Z = self.unpack(x)
        HW = self.Tfc + self.Tff @ W - W @ self.Tcc - W @ (self.Tcf @ W)
        K = self.Tff - W @ self.Tcf + HW @ Z
        return W, Z, K, HW

    def EN(self, x):
        W, Z, K, HW = self.KH(x)
        nc, nf = self.nc, self.nf
        P = jnp.vstack([jnp.eye(nc), W])
        Q = jnp.vstack([jnp.zeros((nc, nf)), jnp.eye(nf)]) + P @ Z
        Qsh = jnp.hstack([-W, jnp.eye(nf)])
        S, V = HW, HW
        for _ in range(self.k - 1):
            V = K @ V
            S = S + V
        Pp = P + Q @ S
        R = jnp.hstack([jnp.eye(nc), jnp.zeros((nc, nf))]) if self.restrict == "inj" else P.T
        RA = R @ self.A
        CC = jnp.eye(self.n) - Pp @ jnp.linalg.solve(RA @ Pp, RA)
        return CC @ Q @ jnp.linalg.matrix_power(K, self.k) @ Qsh

    def fit(self, x0, objective="Kk", maxiter=400, ms=(2, 4, 8, 16)):
        x = np.array(x0, dtype=float)
        stages = ms if objective == "Krho" else (None,)
        for m in stages:
            if objective == "Krho":
                f = lambda y, m=m: jnp.sum(jnp.linalg.matrix_power(self.KH(y)[2], m) ** 2) ** (1.0 / m)
            else:
                f = lambda y: jnp.sum(jnp.linalg.matrix_power(self.KH(y)[2], self.k) ** 2)
            vg = jax.jit(jax.value_and_grad(f))

            def fun(y):
                v, g = vg(jnp.array(y))
                return float(v), np.array(g)
            x = minimize(fun, x, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter)).x
        return x

    def report(self, x):
        W, Z, K, HW = self.KH(jnp.array(x))
        K = np.array(K)
        E = np.array(self.EN(jnp.array(x)))
        return dict(rhoK=rb.rho(K), nKk=np.linalg.norm(np.linalg.matrix_power(K, self.k), 2),
                    rhoE=rb.rho(E), nW=np.count_nonzero(np.abs(np.array(W)) > 1e-14) / self.nf,
                    nZ=np.count_nonzero(np.abs(np.array(Z)) > 1e-14) / self.nf)


def run_1d(n=191, pes=(5, 8, 20), ks=(2, 3, 4), width=2):
    from cd1d import cd1d, split3
    print(f"# 1D centered CD, coarsening by 3, joint local (W, Z), n = {n}, width {width}\n")
    print("| Pe | k | variant | nnz W/nf | nnz Z/nf | rho(K) | ‖K^k‖₂ | rho(E_N) inj | rho(E_N) gal |")
    print("|---|---|---|---|---|---|---|---|---|")
    for pe in pes:
        A = cd1d(n, pe).toarray()
        isF = split3(n)
        C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
        perm = np.concatenate([C, F])
        A = A[np.ix_(perm, perm)]
        nc = len(C)
        near = [np.argsort(np.abs(C - f))[:width] for f in F]
        for k in ks:
            for lab, pW, pZ in (("Q0, P0", [np.array([], int)] * len(F), [np.array([], int)] * len(F)),
                                ("Z only", [np.array([], int)] * len(F), near),
                                ("W only", near, [np.array([], int)] * len(F)),
                                ("W and Z", near, near)):
                out = []
                for rs in ("inj", "gal"):
                    J = Joint(A, nc, k, pW, pZ, restrict=rs)
                    x0 = np.zeros(J.nw + J.nz)
                    x = J.fit(x0) if len(x0) else x0
                    out.append(J.report(x))
                r = out[0]
                f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                print(f"| {pe} | {k} | {lab} | {r['nW']:.1f} | {r['nZ']:.1f} | {f(r['rhoK'])} | {r['nKk']:.2g} | "
                      f"{f(r['rhoE'])} | {f(out[1]['rhoE'])} |")
                sys.stdout.flush()


if __name__ == "__main__":
    run_1d()
