# Claims: N and its Markovian twin


## C. Gauss-Seidel: K = Qsh T Q is dense, N needs a matrix-free history

| problem | m | nnz(K)/n_f^2 (Jacobi) | nnz(K)/n_f^2 (GS) | ‖x_N - x_twin‖/‖x_N‖ | N: GS sweeps per cycle (matrix-free K, H) | twin: GS sweeps per cycle |
|---|---|---|---|---|---|---|
| CD eps=1e-2 | 15 | 0.155 | 0.547 | 1.5e-15 | 9 | 3 |
| CD eps=1e-2 | 31 | 0.041 | 0.390 | 1.7e-15 | 9 | 3 |
| CD eps=1e-3 | 15 | 0.155 | 0.547 | 1.4e-15 | 9 | 3 |
| CD eps=1e-3 | 31 | 0.041 | 0.391 | 1.9e-15 | 9 | 3 |
