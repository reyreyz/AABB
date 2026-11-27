import sys, numpy as np
from ml_twins import ML, poisson, factor
from robust_n import recirc
m = int(sys.argv[1]) if len(sys.argv) > 1 else 127
probs = [("recirc Pe 20, Jacobi w=1 (rho(T)=4.5)", lambda: recirc(m), 1.0),
         ("recirc Pe 2, Jacobi w=0.8", lambda: recirc(m, 2.0), 0.8),
         ("Poisson, Jacobi w=0.8", lambda: poisson(m), 0.8),
         ("aniso Poisson eps=0.01, Jacobi w=0.8", lambda: poisson(m, 0.01), 0.8)]
print(f"m={m}\n\n| problem | truncation of S | P~ nnz/nnz(A) | op.cx | N factor | twin factor | Markov-Jacobi factor |")
print("|---|---|---|---|---|---|---|")
for pname, mk, w in probs:
    A = mk(); b = np.random.default_rng(0).standard_normal(A.shape[0])
    for tau, mode in ((0, "drop"), (1e-2, "drop"), (1e-1, "drop"), (3e-1, "drop"), (1e-1, "lump"), (3e-1, "lump")):
        ml = ML(A, omega=w, tauP=tau, mode=mode)
        oc, pc = ml.opcx()
        f = [factor(ml.solve(b, wh)) for wh in ("N", "twin", "markovJ")]
        print(f"| {pname} | {'none' if tau == 0 else f'{tau:g} {mode}'} | {pc:.2f} | {oc:.2f} | "
              + " | ".join("div." if x == float('inf') else f"{x:.3g}" for x in f) + " |")
        sys.stdout.flush()
