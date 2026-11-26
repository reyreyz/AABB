"""Optimizing the complement Q = Q0 + P Z of the non-Markovian scheme (N_Z) for CD.

Objects (C/F order, P = [I; W], Q0 = [0; I], Q0sh = [-W I]):
    K_Z = K + H Z                       CR propagator of the complement      (Prop. 3.6)
    P'  = P + Q sum_{l<k} K_Z^l H       effective prolongation              (4.11)
    E_N = (I - P'(Rhat A P')^{-1} Rhat A) Q K_Z^k Q0sh                      (4.13)
Z is sparse with a fixed pattern S (column j = fine point j):
    'interp' : the coarse points interpolating to j (the paper's N_Z, ~2.6/col)
    'dist3'  : coarse points within graph distance 3 of j in A (distance 2 adds
               nothing for this coarsening: no coarse point is two steps from a fine one)
    'dense'  : all coarse points
Objectives f(Z), minimized over the entries of Z in S (L-BFGS, JAX gradients),
started from the LS solution:
    LS      ||K_Z||_F^2                       (closed form, columnwise; = paper)
    Kk      ||K_Z^k||_F^2                     (CR noise without Q)
    noise   ||Q K_Z^k Q0sh||_F^2              (the CR noise of (4.13))
    EN_F    ||E_N||_F^2                       (whole two-grid error)
    EN_rho  ||E_N^m||_F^(2/m), m = 2, 4, 8, 16 (spectral-radius surrogate)
    K_rho   ||K_Z^m||_F^(2/m), m = 2, 4, 8, 16 (sparse "pole placement")
    ALS     local alternating LS for ||K_Z^k||: column j of the update minimizes
            ||K_{Z_t}^{k-1} (K e_j + H z)||_2 over z in S_j (cheap; multilevel-ready)
"""

import sys
import time
import numpy as np
import jax
import jax.numpy as jnp
from scipy.optimize import minimize
import rbamg as rb
from reproduce_sm1 import setup

jax.config.update("jax_enable_x64", True)


def pattern(s, kind):
    nc, nf = s.nc, s.nf
    W = s.P[nc:]
    if kind == "interp":
        return [np.nonzero(W[j])[0] for j in range(nf)]
    if kind == "dense":
        return [np.arange(nc) for _ in range(nf)]
    d = int(kind[4:]) if kind.startswith("dist") else 3
    G = (np.abs(s.A) > 0).astype(float)
    Gd = np.linalg.matrix_power(G, d) > 0
    return [np.nonzero(Gd[nc + j, :nc])[0] for j in range(nf)]


class Problem:
    def __init__(self, s, pat):
        self.s = s
        self.pat = pat
        self.rows = np.concatenate([p for p in pat])
        self.cols = np.concatenate([np.full(len(p), j) for j, p in enumerate(pat)])
        self.nvar = len(self.rows)
        self.K, self.H = jnp.array(s.K), jnp.array(s.H)
        self.P, self.Q0, self.Q0sh = jnp.array(s.P), jnp.array(s.Q0), jnp.array(s.Q0sh)
        self.RA = jnp.array(s.Rhat @ s.A)
        self.k = s.k

    def Z(self, z):
        return jnp.zeros((self.s.nc, self.s.nf)).at[self.rows, self.cols].set(z)

    def parts(self, z):
        Z = self.Z(z)
        KZ = self.K + self.H @ Z
        Q = self.Q0 + self.P @ Z
        return Z, KZ, Q

    def EN(self, z):
        Z, KZ, Q = self.parts(z)
        S = self.H
        acc = self.H
        for _ in range(self.k - 1):
            S = KZ @ S
            acc = acc + S
        Pp = self.P + Q @ acc
        CC = jnp.eye(self.s.n) - Pp @ jnp.linalg.solve(self.RA @ Pp, self.RA)
        return CC @ Q @ jnp.linalg.matrix_power(KZ, self.k) @ self.Q0sh

    def obj(self, name, m=2):
        k = self.k
        if name in ("LS", "K"):
            return lambda z: jnp.sum(self.parts(z)[1] ** 2)
        if name == "Kk":
            return lambda z: jnp.sum(jnp.linalg.matrix_power(self.parts(z)[1], k) ** 2)
        if name == "noise":
            def f(z):
                Z, KZ, Q = self.parts(z)
                return jnp.sum((Q @ jnp.linalg.matrix_power(KZ, k) @ self.Q0sh) ** 2)
            return f
        if name == "EN_F":
            return lambda z: jnp.sum(self.EN(z) ** 2)
        if name == "EN_rho":
            return lambda z: jnp.sum(jnp.linalg.matrix_power(self.EN(z), m) ** 2) ** (1.0 / m)
        if name == "K_rho":
            return lambda z: jnp.sum(jnp.linalg.matrix_power(self.parts(z)[1], m) ** 2) ** (1.0 / m)
        raise ValueError(name)

    def ls_start(self):
        s = self.s
        z = []
        for j, Cj in enumerate(self.pat):
            sol, *_ = np.linalg.lstsq(s.H[:, Cj], -s.K[:, j], rcond=None)
            z.append(sol)
        return np.concatenate(z)

    def als(self, z0, sweeps=6):
        s = self.s
        z = z0.copy()
        for _ in range(sweeps):
            Z = np.array(self.Z(z))
            KZ = s.K + s.H @ Z
            Lk = np.linalg.matrix_power(KZ, s.k - 1)
            LK, LH = Lk @ s.K, Lk @ s.H
            out = []
            for j, Cj in enumerate(self.pat):
                sol, *_ = np.linalg.lstsq(LH[:, Cj], -LK[:, j], rcond=None)
                out.append(sol)
            z = np.concatenate(out)
        return z

    def optimize(self, name, z0, maxiter=300):
        ms = (2, 4, 8, 16) if name.endswith("rho") else (None,)
        z = z0.copy()
        for m in ms:
            f = self.obj(name, m or 2)
            vg = jax.jit(jax.value_and_grad(f))

            def fun(x):
                v, g = vg(jnp.array(x))
                return float(v), np.array(g)
            res = minimize(fun, z, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter))
            z = res.x
        return z

    def report(self, z):
        s = self.s
        Z = np.array(self.Z(z))
        KZ = s.K + s.H @ Z
        E = np.array(self.EN(jnp.array(z)))
        return dict(rhoK=rb.rho(KZ), nKk=np.linalg.norm(np.linalg.matrix_power(KZ, s.k), 2),
                    rhoE=rb.rho(E), nE=np.linalg.norm(E, 2),
                    nnz=np.count_nonzero(np.abs(Z) > 1e-14) / s.nf)


def study(name, m, k, pats=("interp", "dist3"), objs=("LS", "ALS", "Kk", "noise", "EN_F", "EN_rho", "K_rho")):
    s = setup(name, m, k)
    print(f"\n### {name}, m = {m}, k = {k} (bilinear P, Galerkin, Jacobi 2/3)\n")
    print(f"References: Markovian rho(E_M) = {rb.rho(rb.E_markov(s)):.3e}, "
          f"N with Q0: rho = {rb.rho(rb.E_nonmarkov(s)):.3e}, "
          f"S (canonical): rho = {rb.rho(rb.E_semi(s)):.3e}\n")
    print("| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |")
    print("|---|---|---|---|---|---|---|---|")
    out = {}
    for pk in pats:
        pr = Problem(s, pattern(s, pk))
        z0 = pr.ls_start()
        for ob in objs:
            t = time.time()
            if ob == "LS":
                z = z0
            elif ob == "ALS":
                z = pr.als(z0)
            else:
                z = pr.optimize(ob, z0)
            r = pr.report(z)
            out[(pk, ob)] = (r, z)
            print(f"| {pk} | {ob} | {r['rhoK']:.3e} | {r['nKk']:.2e} | {r['rhoE']:.3e} | "
                  f"{r['nE']:.2e} | {r['nnz']:.1f} | {time.time()-t:.1f} |")
            sys.stdout.flush()
    return out


if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    pats = tuple(sys.argv[2].split(",")) if len(sys.argv) > 2 else ("interp", "dist3")
    print(f"# N_Z complement optimization (CD, m = {m})")
    for name in ("CD eps=1e-2", "CD eps=1e-3"):
        for k in (2, 4):
            study(name, m, k, pats=pats)
