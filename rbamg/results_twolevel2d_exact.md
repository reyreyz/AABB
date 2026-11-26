| m | coarse | factor | cycles | WU/cycle | stored/nnz(A) |
|---|---|---|---|---|---|
| 47 | LU (exact) | 0.224 | 18 | 115.3 | 88.3 |
| 47 | AIR+GMRES tol 1e-6 | 120 | 6 | 2981.0 | 88.3 |
| 47 | AIR+GMRES tol 1e-9 | 120 | 6 | 2981.0 | 88.3 |
| 95 | LU (exact) | 0.301 | 25 | 135.7 | 113.7 |
| 95 | AIR+GMRES tol 1e-6 | 139 | 4 | 3858.2 | 113.7 |
| 95 | AIR+GMRES tol 1e-9 | 136 | 4 | 4541.1 | 113.7 |
| 191 | LU (exact) | 0.465 | 30 | 146.4 | 127.1 |
/root/.local/lib/python3.13/site-packages/pyamg/krylov/_gmres_householder.py:319: LinAlgWarning: An ill-conditioned matrix detected: slice 0 has rcond = 1.2876656968822447e-19.
  y = sp.linalg.solve(H[0:(inner+1), 0:(inner+1)], g[0:(inner+1)])
| 191 | AIR+GMRES tol 1e-6 | 1.24e+16 | 1 | 5727.9 | 127.1 |
/root/.local/lib/python3.13/site-packages/pyamg/krylov/_gmres_householder.py:319: LinAlgWarning: An ill-conditioned matrix detected: slice 0 has rcond = 1.2876656968822447e-19.
  y = sp.linalg.solve(H[0:(inner+1), 0:(inner+1)], g[0:(inner+1)])
| 191 | AIR+GMRES tol 1e-9 | 1.24e+16 | 1 | 5727.9 | 127.1 |
