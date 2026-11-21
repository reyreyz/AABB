import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hh1d import *
n = 15; rng = np.random.default_rng(2)
for kh in (0.6, 0.9):
    A, S = helmholtz(n, kh); cs, _ = supports(n, 3); sups = [np.arange(n)] * len(cs)
    lam, V = np.linalg.eigh(S); idx = np.argsort(-np.abs(lam))[:len(cs)]
    starts = [pack(V[:, idx], sups)] + [rng.standard_normal(n * len(cs)) for _ in range(3)]
    v, P = true_opt_norm(A, S, sups, starts)
    print(f"n={n} kh={kh}: dense-P min ||E||2 = {v:.4f}; eigenvector space ||E||2 = {np.linalg.norm(twogrid(A,S,V[:,idx]),2):.4f}", flush=True)
