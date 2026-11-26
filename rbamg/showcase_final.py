"""Showcase of the N_Z formalism: centered CD at high cell Peclet number.

Relaxation diverges (rho(T) > 1), so the Markovian scheme E^M = Pi_R T^k (Thm 4.2)
diverges, and so does the semi-Markovian one.  The non-Markovian scheme carries
T only through the CR term Q K_Z^k Q0sh (Thm 4.5): it converges as soon as the
complement Q = Q0 + P Z has a contracting CR propagator K_Z = K + H Z.  For the
local complement Q0, rho(K) > 1 as well.  Theorem 4.8 / Corollary 4.9: if (K, H)
is controllable (PBH > 0), Z can make K_Z nilpotent with index
kappa >= ceil(n_f / n_c); a k-step method then needs k >= kappa.

Z is built CR-centrically and simply: alternating least squares on ||K_Z^k||_F,
column j of the update minimizing ||K_{Z_t}^{k-1}(K e_j + H z)|| (dense pattern),
best iterate by ||K_Z^k||_F.  The schemes are run as algorithms on iterates
(Alg. 4.2 with complement Q, random b, x0 = 0).
"""

import sys
import numpy as np
import rbamg as rb
import nz_opt as no
from showcase_opt import make_setup
from showcase_scan import centered_cd
from reproduce_sm4 import pbh
from reproduce_sm10 import nz_alg


def best_als(pr, sweeps=30):
    s = pr.s
    z = pr.ls_start()
    best = (np.inf, z)
    for _ in range(sweeps):
        z = pr.als(z, sweeps=1)
        Z = np.array(pr.Z(z))
        val = np.linalg.norm(np.linalg.matrix_power(s.K + s.H @ Z, s.k))
        if val < best[0]:
            best = (val, z)
    return best[1]


def history(alg, s, ncyc=12):
    rng = np.random.default_rng(0)
    b = rng.standard_normal(s.n)
    x = np.zeros(s.n)
    nb = np.linalg.norm(b)
    h = [1.0]
    for _ in range(ncyc):
        x = alg(x, b)
        h.append(np.linalg.norm(b - s.A @ x) / nb)
        if not np.isfinite(h[-1]) or h[-1] > 1e30 or h[-1] < 1e-14:
            break
    return h


def main(ms=(15, 31), pes=(10, 20), ks=(3, 4)):
    print("# Showcase: only N_Z converges (centered CD, bilinear P, Galerkin, Jacobi 2/3)\n")
    print("| m | Pe | k | ceil(nf/nc) | PBH | rho(T) | rho(K) Q0 | rho(K_Z) | M | S | N (Q0) | **N_Z** |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    hists = {}
    for m in ms:
        for pe in pes:
            for k in ks:
                s = make_setup(centered_cd(m, pe), m, k)
                pr = no.Problem(s, no.pattern(s, "dense"))
                z = best_als(pr)
                r = pr.report(z)
                row = (m, pe, k, int(np.ceil(s.nf / s.nc)), pbh(s), rb.rho(s.T), rb.rho(s.K), r["rhoK"],
                       rb.rho(rb.E_markov(s)), rb.rho(rb.E_semi(s)), rb.rho(rb.E_nonmarkov(s)), r["rhoE"])
                print("| {} | {} | {} | {} | {:.1e} | {:.2f} | {:.2f} | {:.2e} | {:.3g} | {:.3g} | {:.3g} | **{:.2e}** |".format(*row))
                sys.stdout.flush()
                Z = np.array(pr.Z(z))
                hists[(m, pe, k)] = dict(
                    M=history(lambda x, b: rb.alg_markov(s, x[:, None], b[:, None])[:, 0], s),
                    N=history(lambda x, b: rb.alg_nonmarkov(s, x[:, None], b[:, None])[:, 0], s),
                    NZ=history(lambda x, b, a=nz_alg(s, Z): a(x[:, None], b[:, None])[:, 0], s))
    print("\n## Residual histories of the algorithms (relative residual per cycle)\n")
    for key, hv in hists.items():
        print(f"m = {key[0]}, Pe = {key[1]}, k = {key[2]}:")
        for lab in ("M", "N", "NZ"):
            print(f"- {lab:2s}: " + " ".join(f"{v:.1e}" for v in hv[lab]))
        print()
    return hists


if __name__ == "__main__":
    main()
