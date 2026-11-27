"""Best Markovian cycle on the same red-black/single-entry P~ hierarchy: Markov-Jacobi over (w_F, w_C, k) and the
twin over (w_F, k); multilevel, cheap coarse levels."""
import sys, itertools, numpy as np
from ml_twins import ML
from flows import cd
from twin_diverges import run
m = int(sys.argv[1]); A, xy = cd(m, "const", 6.0); n = A.shape[0]; b = np.random.default_rng(0).standard_normal(n)
I, J = np.divmod(np.arange(n), m); isF = (I + J) % 2 == 1
coarse = dict(k1=2, tauP_deep=1e-2, tau_filter=1e-2, tau_deep=1e-2, refine=0)
res = []
for wf, wc, k in itertools.product((0.2, 0.4, 0.6, 0.8, 1.0), (0.0, 0.2, 0.5, 1.0), (1, 2, 4, 8)):
    ml = ML(A, omega=wf, omega_c=wc * wf, wc_levels=1, k0=k, split0=isF, mode="lin", xy=xy, nkeep=1, lin_levels=1, **coarse)
    res.append((run(ml, b, "markovJ", 40)[0], "MJ", wf, wc * wf, k))
    if wc == 0.0: res.append((run(ml, b, "twin", 40)[0], "twin", wf, "-", k))
res.sort(key=lambda r: r[0])
print(f"m={m}: best Markovian cycles on the same hierarchy (rate, method, w_F, w_C, k), of {len(res)}:")
for r in res[:6]: print("  ", (round(float(r[0]), 3),) + r[1:])
