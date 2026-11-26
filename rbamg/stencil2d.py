"""2D two-level example: centered CD with wind along x, semi-coarsening by 3 in x,
stencil decompositions (W, Z) acting along x-lines (y-invariant parameters).

Grid m_x x m_y (x fast), C-points where (j mod 3) == 2.  For each F-point the
coarse neighbours are the C-points of its own x-line on either side (w per side),
parameters periodic in p blocks along x, identical for all y: 16 numbers for
w = 1, p = 2, fitted once on a small grid by min ||K_{W,Z}^k||_F and reused.
"""

import sys
import numpy as np
import scipy.sparse as sp
import jax
import jax.numpy as jnp
from scipy.optimize import minimize
import rbamg as rb
from joint_wz import Joint

jax.config.update("jax_enable_x64", True)


def cd2d_x(mx, my, pe):
    """-eps Lap u + u_x, centered, h = 1/(mx+1) in both directions, cell Peclet h/eps."""
    h = 1.0 / (mx + 1)
    eps = h / pe
    Ix, Iy = sp.eye(mx), sp.eye(my)
    Lx = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(mx, mx))
    Ly = sp.diags([-1, 2, -1], [-1, 0, 1], shape=(my, my))
    Dx = sp.diags([-0.5, 0.5], [-1, 1], shape=(mx, mx))
    return (eps * (sp.kron(Iy, Lx) + sp.kron(Ly, Ix)) + h * sp.kron(Iy, Dx)).tocsr()


def layout2d(mx, my, w, p, wy=0):
    isC = np.zeros(mx * my, dtype=bool)
    for i in range(my):
        isC[i * mx + np.arange(2, mx, 3)] = True
    C, F = np.nonzero(isC)[0], np.nonzero(~isC)[0]
    cidx = {c: t for t, c in enumerate(C)}
    ny = 2 * wy + 1
    nslot = 2 * p * 2 * w * ny
    e = []
    for jf, f in enumerate(F):
        i, j = divmod(f, mx)
        t = 0 if (j - 1) % 3 == 2 and j >= 1 else 1
        cl = j - 1 if t == 0 else j - 2
        phase = (cl // 3) % p
        nbrs = [cl - 3 * s for s in range(w - 1, -1, -1)] + [cl + 3 * (s + 1) for s in range(w)]
        for dy in range(-wy, wy + 1):
            ii = i + dy
            if not (0 <= ii < my):
                continue
            for q, cj in enumerate(nbrs):
                g = ii * mx + cj
                if 0 <= cj < mx and g in cidx:
                    e.append((jf, cidx[g], ((t * p + phase) * ny + (dy + wy)) * 2 * w + q))
    return C, F, np.array(e), nslot


def joint2d(mx, my, pe, k, w, p, restrict="inj", wy=0):
    A = cd2d_x(mx, my, pe).toarray()
    C, F, e, nslot = layout2d(mx, my, w, p, wy)
    perm = np.concatenate([C, F])
    A = A[np.ix_(perm, perm)]
    pat = [e[e[:, 0] == j, 1] for j in range(len(F))]
    J = Joint(A, len(C), k, pat, pat, restrict=restrict)
    slots = np.concatenate([e[:, 2], nslot + e[:, 2]])
    return J, slots, nslot, A


def fit2d(pe, k, w=1, p=2, mx=23, my=12, starts=12, seed=0, wy=0):
    J, slots, nslot, _ = joint2d(mx, my, pe, k, w, p, wy=wy)
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


def evaluate2d(theta, pe, k, w, p, mx, my, wy=0):
    J, slots, nslot, A = joint2d(mx, my, pe, k, w, p, wy=wy)
    x = theta[slots]
    r = J.report(x)
    q = J.report(np.zeros_like(x))
    nc = J.nc
    Pm = np.vstack([np.eye(nc), np.zeros((J.nf, nc))])
    sM = rb.Setup(A, np.diag(1 / np.diag(A)), Pm, Pm.T.copy(), k)
    return r, q, rb.rho(rb.E_markov(sM))


def main(pes=(8, 12), ks=(4, 6), sizes=((23, 23), (47, 47), (71, 71)), wy=0):
    print(f"# 2D two-level: centered CD, wind along x, semi-coarsening by 3, stencils w = 1, p = 2, wy = {wy}\n")
    print("| Pe | k | m_x x m_y | n | M | N (Q0, P0) | rho(K_WZ) | **N_WZ** |")
    print("|---|---|---|---|---|---|---|---|")
    for pe in pes:
        for k in ks:
            th = fit2d(pe, k, wy=wy)
            for mx, my in sizes:
                r, q, EM = evaluate2d(th, pe, k, 1, 2, mx, my, wy=wy)
                f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                print(f"| {pe} | {k} | {mx} x {my} | {mx*my} | {f(EM)} | {f(q['rhoE'])} | {f(r['rhoK'])} | {f(r['rhoE'])} |")
                sys.stdout.flush()


if __name__ == "__main__":
    main(wy=int(sys.argv[1]) if len(sys.argv) > 1 else 0)
