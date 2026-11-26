"""Coarse-aware complement: choose Z so that the CR propagator is small AND the
next coarse operator stays relaxable.

Level data (C/F blocks, P = [I; W], Q = Q0 + P Z, Rhat = P^T):
    K_Z = K + H Z,  S = sum_{l<k} K_Z^l H,  P'_C = I + Z S,  P'_F = W + S + W Z S
    A_sigma(Z) = [I; W]^T A_cf [P'_C; P'_F]                       (next level)
    T_c(Z) = I - omega D_c^{-1} A_sigma(Z),  D_c = diag(A_sigma(Z))
Spectral-radius surrogate (n_c^{1/64} <= 1.11 overestimate for n_c <= 1000):
    s(Z) = (||T_c^32||_F^2 / n_c)^(1/64)
Objectives (Z on the sparsity pattern, L-BFGS with JAX gradients, start = LS):
    'Kk'        ||K_Z^k||_F^2 / n_f
    'Kk+jac'    ||K_Z^k||_F^2 / n_f + lam * max(0, s(Z) - s0)^2
    'LS+jac'    ||K_Z||_F^2 / n_f   + lam * max(0, s(Z) - s0)^2
Used per level through NZCycle(zmode=callable), so each level's objective sees
the operator actually produced by the finer levels.
"""

import numpy as np
import jax
import jax.numpy as jnp
import scipy.sparse as sp
from scipy.optimize import minimize
import ml_nz as mz

jax.config.update("jax_enable_x64", True)


def jac_surrogate(Asig, omega, p=5):
    nc = Asig.shape[0]
    T = jnp.eye(nc) - omega * Asig / jnp.diag(Asig)[:, None]
    for _ in range(p):                    # T^(2^p), rescaled to avoid overflow
        nrm = jnp.linalg.norm(T)
        T = (T / nrm) @ (T / nrm) * nrm ** 2
    return (jnp.sum(T ** 2) / nc) ** (1.0 / 2 ** (p + 1))


def make_zfunc(objective="Kk+jac", lam=10.0, s0=0.9, maxiter=150, log=None):
    def zfunc(d):
        K, H, W, patt, k, om = d["K"], d["H"], d["W"], d["patt"], d["k"], d["omega"]
        A, C, F = d["A"], d["C"], d["F"]
        nf, nc = K.shape[0], H.shape[1]
        if nf * nc > 4e6:                 # too large for the dense objective: LS
            return mz.z_ls(K, H, patt)
        perm = np.concatenate([C, F])
        Acf = jnp.array(sp.csr_matrix(A)[perm][:, perm].toarray())
        Kd, Hd, Wd = jnp.array(K.toarray()), jnp.array(H.toarray()), jnp.array(W.toarray())
        pc = sp.coo_matrix(patt)
        r_, c_ = pc.row, pc.col
        Pb = jnp.vstack([jnp.eye(nc), Wd])

        def Zof(z):
            return jnp.zeros((nc, nf)).at[r_, c_].set(z)

        def asig(Z):
            KZ = Kd + Hd @ Z
            S, V = Hd, Hd
            for _ in range(k - 1):
                V = KZ @ V
                S = S + V
            Pp = jnp.vstack([jnp.eye(nc) + Z @ S, Wd + S + Wd @ (Z @ S)])
            return KZ, Pb.T @ Acf @ Pp

        def f(z):
            Z = Zof(z)
            KZ, As = asig(Z)
            if objective.startswith("LS"):
                main = jnp.sum(KZ ** 2) / nf
            else:
                main = jnp.sum(jnp.linalg.matrix_power(KZ, k) ** 2) / nf
            if objective.endswith("+jac"):
                sj = jac_surrogate(As, om)
                main = main + lam * jnp.maximum(0.0, sj - s0) ** 2
            return main

        vg = jax.jit(jax.value_and_grad(f))
        Z0 = mz.z_ls(K, H, patt)
        z0 = np.asarray(sp.csr_matrix(Z0)[r_, c_]).ravel()

        def fun(x):
            v, g = vg(jnp.array(x))
            return float(v), np.array(g)
        res = minimize(fun, z0, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter))
        Z = sp.csr_matrix((res.x, (r_, c_)), shape=(nc, nf))
        if log is not None:
            _, As0 = asig(Zof(jnp.array(z0)))
            _, As1 = asig(Zof(jnp.array(res.x)))
            log.append(dict(level=d["level"], s_LS=float(jac_surrogate(As0, om)),
                            s_opt=float(jac_surrogate(As1, om)), f0=float(f(jnp.array(z0))), f1=float(res.fun)))
        return Z
    return zfunc
