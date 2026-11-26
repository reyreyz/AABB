# Multilevel 1D, centered CD Pe = 8, coarsening by 3, Jacobi omega = 1, injection

finest: k = 4, stencil w = 1, p = 2; coarse levels: k = 10, w = 1, p = 4; coarse operators capped to bandwidth 3

| n | cycle | levels | factor | cycles to 1e-10 | WU/cycle | op.cx | fit res per level | setup [s] |
|---|---|---|---|---|---|---|---|---|
| 689 | M (Markovian) | 4 | 378 | 5 (no conv.) | 4.6 | 1.16 | – | 0 |
| 689 | N (Q0, P0) | 4 | 3.4e+03 | 4 (no conv.) | 10.6 | 1.48 | – | 0 |
| 689 | N_WZ | 4 | 3.57e+14 | 1 (no conv.) | 55.3 | 2.10 | 4e-06,2e-04,1e-01 | 13 |
| 2066 | M (Markovian) | 5 | 378 | 5 (no conv.) | 4.6 | 1.16 | – | 0 |
| 2066 | N (Q0, P0) | 5 | 3.39e+03 | 4 (no conv.) | 10.9 | 1.49 | – | 0 |
| 2066 | N_WZ | 5 | 1.86e+31 | 1 (no conv.) | 59.4 | 2.14 | 4e-06,3e-04,6e-05,5e+00 | 15 |
| 6197 | M (Markovian) | 6 | 378 | 5 (no conv.) | 4.7 | 1.17 | – | 0 |
| 6197 | N (Q0, P0) | 6 | 3.39e+03 | 4 (no conv.) | 10.9 | 1.50 | – | 0 |
| 6197 | N_WZ | 6 | 1.12e+35 | 1 (no conv.) | 61.0 | 2.16 | 4e-06,3e-04,1e-04,5e-05,3e-04 | 20 |
