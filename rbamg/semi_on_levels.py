import sys, numpy as np, scipy.sparse as sp
from flows import cd as cdflow
from semi_on import S72
m = int(sys.argv[1]); A = sp.csr_matrix(cdflow(m, "const", 20)[0])
for tau in (3e-3, 1e-2):
    ml = S72(A, L=3, tau=tau)
    print(f"m={m} tau={tau}: n_l ratio / nnz(A_l)/nnz0 / nnz per row:")
    print("   " + "  ".join(f"{l.n}:{l.A.nnz/A.nnz:.2f}:{l.A.nnz/l.n:.0f}" for l in ml.levels()))
