# O(N) N_Z showcase: centered CD, cell Peclet 6.0, bilinear hierarchy, Galerkin, Jacobi 2/3, tauA = 0.001, V-cycles

| m | n | k | cycle | theta per level | factor | cycles | WU/cycle | op.cx | max nnz(Z)/nf | setup [µs/unknown] |
|---|---|---|---|---|---|---|---|---|---|---|
| 31 | 961 | 5 | M (Markovian) | – | 2.29e+03 | 3 (no conv.) | 9.4 | 1.47 | 0.0 | 14 |
| 31 | 961 | 5 | N (Q0) | – | 3.19 | 17 (no conv.) | 57.4 | 3.20 | 0.0 | 18 |
| 31 | 961 | 5 | N_Z Kk dist3 | 1,1 | 1.14e+03 | 3 (no conv.) | 182.3 | 5.71 | 49.0 | 355942 |
