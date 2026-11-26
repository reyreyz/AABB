"""Centered CD at high cell Peclet: M and N(Q0) diverge; find a complement Z
with rho(K + H Z) < 1 (Theorem 4.8) so that N_Z converges.

CR-centric objectives only (nz_opt.Problem, JAX + L-BFGS, start from LS):
  LS, ALS (iterated k-step CR LS), Kk = ||K_Z^k||_F^2, K_rho = ||K_Z^m||_F^(2/m)
Patterns: interp, dist3, dense.  PBH distance min_w ||wP||/||w|| over left
eigenvectors of T measures controllability of (K, H) (Theorem 4.8(i)).
"""

import sys
import time
import numpy as np
import rbamg as rb
import problems as pb
import nz_opt as no
from showcase_scan import centered_cd
from reproduce_sm4 import pbh


def make_setup(A, m, k, omega=2 / 3):
    perm, P = pb.bilinear_cf(m)
    A = A[np.ix_(perm, perm)]
    return rb.Setup(A, omega * np.diag(1 / np.diag(A)), P, P.T, k)


def main(m=15, pes=(10, 20), ks=(2, 4), pats=("interp", "dist3", "dense"),
         objs=("LS", "ALS", "Kk", "K_rho")):
    print(f"# Showcase: centered CD at high Peclet, m = {m}\n")
    for pe in pes:
        for k in ks:
            s = make_setup(centered_cd(m, pe), m, k)
            print(f"\n## Pe = {pe}, k = {k}: rho(T) = {rb.rho(s.T):.3f}, rho(K) (Q0) = {rb.rho(s.K):.3f}, "
                  f"PBH = {pbh(s):.1e}\n")
            print(f"M: {rb.rho(rb.E_markov(s)):.3g}, S: {rb.rho(rb.E_semi(s)):.3g}, "
                  f"N (Q0): {rb.rho(rb.E_nonmarkov(s)):.3g}\n")
            print("| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |")
            print("|---|---|---|---|---|---|---|---|")
            for pk in pats:
                pr = no.Problem(s, no.pattern(s, pk))
                z0 = pr.ls_start()
                for ob in objs:
                    t = time.time()
                    z = z0 if ob == "LS" else (pr.als(z0, sweeps=10) if ob == "ALS" else pr.optimize(ob, z0))
                    r = pr.report(z)
                    f = lambda v: f"**{v:.3g}**" if v < 1 else f"{v:.3g}"
                    print(f"| {pk} | {ob} | {f(r['rhoK'])} | {r['nKk']:.2g} | {f(r['rhoE'])} | {r['nE']:.2g} | "
                          f"{r['nnz']:.1f} | {time.time()-t:.0f} |")
                    sys.stdout.flush()


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 15)
