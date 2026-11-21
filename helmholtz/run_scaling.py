"""Locality floor at the near-kernel: eps_loc(z=1; H) = sigma_min((I - S)[:, supp]) for a support of
2H+1 points centred mid-domain, versus H. Poisson (kh=0) vs Helmholtz."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from hh1d import helmholtz, OMEGA

n = 1023
Hs = [2, 4, 8, 16, 32, 64, 128]
print("eps_loc(z=1;H); in brackets the local slope d log eps / d log H")
print("kh     " + "".join(f"H={H:<9d}" for H in Hs))
for kh in (0.0, 0.1, 0.3, 0.6):
    A, S = helmholtz(n, kh)
    I = np.eye(n)
    c = n // 2
    vals = []
    for H in Hs:
        sup = np.arange(c - H, c + H + 1)
        vals.append(np.linalg.svd((I - S)[:, sup], compute_uv=False)[-1])
    sl = [np.nan] + [np.log(vals[i] / vals[i - 1]) / np.log(Hs[i] / Hs[i - 1]) for i in range(1, len(Hs))]
    print(f"{kh:<6} " + "".join(f"{v:.2e}({s:+.1f}) " for v, s in zip(vals, sl)))
    # prediction for Helmholtz: band of width |d lambda/d theta| * pi/(2H+2), lambda=4 sin^2(theta/2)-(kh)^2
    if kh:
        th0 = 2 * np.arcsin(kh / 2)
        pred = [OMEGA / (2 - kh**2) * 2 * np.sin(th0) * np.pi / (2 * H + 2) for H in Hs]
        print("  pred  " + "".join(f"{v:.2e}      " for v in pred))
