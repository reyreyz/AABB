"""Translation-invariant (stencil) W and Z for 1D coarsening by 3: O(1) parameters.

Each F-point f has a type t (t = 0: right neighbour of a C-point, t = 1: left
neighbour of a C-point) and neighbouring C-points c_L < f < c_R.
    W[f, c_L] = a[t],  W[f, c_R] = b[t]       (interpolation, row f)
    Z[c_L, f] = c[t],  Z[c_R, f] = d[t]       (complement feedback, column f)
8 parameters, fitted once by min ||K_{W,Z}^k||_F^2 (multi-start L-BFGS) and
reused unchanged on every mesh: W, Z are fixed stencils, so the setup is O(1)
and every application O(N).
"""

import sys
import numpy as np
import jax
import jax.numpy as jnp
from scipy.optimize import minimize
import rbamg as rb
from joint_wz import Joint
from cd1d import cd1d, split3

jax.config.update("jax_enable_x64", True)


def layout(n):
    isF = split3(n)
    C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
    cidx = {c: i for i, c in enumerate(C)}
    entries = []                    # (fine index jf, coarse index ic, param slot) for W and for Z
    for jf, f in enumerate(F):
        t = 0 if (f - 1) in cidx else 1
        cl = f - 1 if t == 0 else f - 2
        cr = f + 2 if t == 0 else f + 1
        for side, cpt in ((0, cl), (1, cr)):
            if cpt in cidx:
                entries.append((jf, cidx[cpt], 2 * t + side))
    e = np.array(entries)
    return C, F, e


def build(n, pe, k, restrict="inj"):
    C, F, e = layout(n)
    perm = np.concatenate([C, F])
    A = cd1d(n, pe).toarray()[np.ix_(perm, perm)]
    patW = [e[e[:, 0] == j, 1] for j in range(len(F))]
    J = Joint(A, len(C), k, patW, patW, restrict=restrict)
    slots = np.concatenate([e[:, 2], 4 + e[:, 2]])   # W uses theta[0:4], Z uses theta[4:8]
    return J, slots


def fit(pe, k, n=95, starts=12, seed=0):
    J, slots = build(n, pe, k)
    Phi = np.zeros((len(slots), 8))
    Phi[np.arange(len(slots)), slots] = 1.0
    Phi = jnp.array(Phi)
    f = lambda th: jnp.sum(jnp.linalg.matrix_power(J.KH(Phi @ th)[2], k) ** 2)
    vg = jax.jit(jax.value_and_grad(f))

    def fun(th):
        v, g = vg(jnp.array(th))
        return float(v), np.array(g)
    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        th0 = rng.standard_normal(8) * (0.5 if s else 0.0)
        r = minimize(fun, th0, jac=True, method="L-BFGS-B", options=dict(maxiter=500))
        if best is None or r.fun < best.fun:
            best = r
    return best.x


def evaluate(theta, pe, k, n, restrict="inj"):
    J, slots = build(n, pe, k, restrict)
    x = theta[slots]
    r = J.report(x)
    s_q0 = J.report(np.zeros_like(x))
    return r, s_q0


def main(pes=(5, 8, 12, 20), ks=(3, 4, 6), ns=(191, 383, 767)):
    print("# Stencil (W, Z), 1D centered CD, coarsening by 3 — fitted at n = 95, reused on every n\n")
    print("| Pe | k | theta (W: a0 b0 a1 b1 / Z: c0 d0 c1 d1) | n | rho(K) Q0 | N (Q0) | rho(K_WZ) | ‖K^k‖₂ | **N_WZ** |")
    print("|---|---|---|---|---|---|---|---|---|")
    for pe in pes:
        for k in ks:
            th = fit(pe, k)
            for n in ns:
                r, q = evaluate(th, pe, k, n)
                f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                ths = " ".join(f"{t:.3f}" for t in th[:4]) + " / " + " ".join(f"{t:.3f}" for t in th[4:])
                print(f"| {pe} | {k} | {ths if n == ns[0] else ''} | {n} | {q['rhoK']:.3f} | {f(q['rhoE'])} | "
                      f"{f(r['rhoK'])} | {r['nKk']:.2g} | {f(r['rhoE'])} |")
                sys.stdout.flush()




# ---------------------------------------------------------------------------
# general stencils: w coarse neighbours per side, parameters periodic in p blocks
# ---------------------------------------------------------------------------
def layout_general(n, w, p):
    isF = split3(n)
    C, F = np.nonzero(~isF)[0], np.nonzero(isF)[0]
    cidx = {c: i for i, c in enumerate(C)}
    nslot = 2 * p * 2 * w          # type (2) x block phase (p) x neighbour (2w)
    eW = []
    for jf, f in enumerate(F):
        t = 0 if (f - 1) in cidx else 1
        cl = f - 1 if t == 0 else f - 2
        phase = (cl // 3) % p
        nbrs = [cl - 3 * s for s in range(w - 1, -1, -1)] + [cl + 3 * (s + 1) for s in range(w)]
        for q, cpt in enumerate(nbrs):
            if cpt in cidx:
                eW.append((jf, cidx[cpt], (t * p + phase) * 2 * w + q))
    return C, F, np.array(eW), nslot


def fit_general(pe, k, w, p, n=None, starts=16, seed=0):
    n = n or 3 * 16 * p - 1
    C, F, e, nslot = layout_general(n, w, p)
    perm = np.concatenate([C, F])
    A = cd1d(n, pe).toarray()[np.ix_(perm, perm)]
    pat = [e[e[:, 0] == j, 1] for j in range(len(F))]
    J = Joint(A, len(C), k, pat, pat)
    slots = np.concatenate([e[:, 2], nslot + e[:, 2]])
    Phi = np.zeros((len(slots), 2 * nslot))
    Phi[np.arange(len(slots)), slots] = 1.0
    Phi = jnp.array(Phi)
    f = lambda th: jnp.sum(jnp.linalg.matrix_power(J.KH(Phi @ th)[2], k) ** 2)
    vg = jax.jit(jax.value_and_grad(f))

    def fun(th):
        v, g = vg(jnp.array(th))
        return float(v), np.array(g)
    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        r = minimize(fun, rng.standard_normal(2 * nslot) * (0.5 if s else 0.0), jac=True,
                     method="L-BFGS-B", options=dict(maxiter=800))
        if best is None or r.fun < best.fun:
            best = r
    return best.x


def evaluate_general(theta, pe, k, w, p, n, restrict="inj"):
    C, F, e, nslot = layout_general(n, w, p)
    perm = np.concatenate([C, F])
    A = cd1d(n, pe).toarray()[np.ix_(perm, perm)]
    pat = [e[e[:, 0] == j, 1] for j in range(len(F))]
    J = Joint(A, len(C), k, pat, pat, restrict=restrict)
    x = theta[np.concatenate([e[:, 2], nslot + e[:, 2]])]
    return J.report(x), J.report(np.zeros_like(x)), A


def main_general(pes=(8, 12, 20), ks=(3, 4, 6), wps=((2, 1), (1, 2), (2, 2)), ns=(191, 383, 767)):
    print("# General stencils (w neighbours per side, period p), 1D centered CD, coarsening by 3\n")
    print("| Pe | k | w | p | #params | n | M | N (Q0) | rho(K_WZ) | **N_WZ** |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for pe in pes:
        for k in ks:
            for w, p in wps:
                th = fit_general(pe, k, w, p)
                for n in ns:
                    r, q, A = evaluate_general(th, pe, k, w, p, n)
                    nc = (n + 1) // 3
                    import rbamg as rb_
                    Pm = np.vstack([np.eye(nc), np.zeros((n - nc, nc))])
                    sM = rb_.Setup(A, np.diag(1 / np.diag(A)), Pm, Pm.T.copy(), k)
                    EM = rb_.rho(rb_.E_markov(sM))
                    f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                    print(f"| {pe} | {k} | {w} | {p} | {len(th)} | {n} | {f(EM)} | {f(q['rhoE'])} | {f(r['rhoK'])} | {f(r['rhoE'])} |")
                    sys.stdout.flush()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "general":
        main_general()
    else:
        main()
