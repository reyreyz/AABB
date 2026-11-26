# Multilevel N_Z V-cycles on CD (bilinear hierarchy, Galerkin Rhat = P^T, Jacobi 2/3, filter tauA = 0.001)


## CD eps = 0.01

| m | k | cycle | levels | op.cx | cycles | factor | WU/cycle | WU total | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | 2 | M (Markovian) | 4 | 1.54 | 80 (no conv.) | 8.15e-01 | 5.7 | inf | 0.0 |
| 63 | 2 | N (Q0) | 4 | 2.61 | 27 | 4.54e-01 | 25.9 | 699.9 | 0.0 |
| 63 | 2 | N_Z LS, all levels | 4 | 5.01 | 24 | 4.29e-01 | 68.9 | 1652.8 | 1.3 |
| 63 | 2 | N_Z LS, finest only | 4 | 4.78 | 24 | 4.29e-01 | 59.4 | 1426.8 | 0.9 |
| 63 | 2 | N_Z ALS, finest only | 4 | 5.09 | 23 | 4.10e-01 | 61.2 | 1407.1 | 4.1 |
| 63 | 2 | N_Z Kk, finest only | 4 | 5.10 | 24 | 4.22e-01 | 61.6 | 1479.4 | 7.6 |
| 63 | 4 | M (Markovian) | 4 | 1.54 | 29 | 4.28e-01 | 8.7 | 252.9 | 0.0 |
| 63 | 4 | N (Q0) | 4 | 2.86 | 14 | 2.13e-01 | 51.2 | 717.3 | 0.1 |
| 63 | 4 | N_Z LS, all levels | 4 | 5.09 | 12 | 1.85e-01 | 138.7 | 1664.8 | 1.5 |
| 63 | 4 | N_Z LS, finest only | 4 | 4.84 | 12 | 1.85e-01 | 122.9 | 1475.1 | 1.0 |
| 63 | 4 | N_Z ALS, finest only | 4 | 4.77 | 17 | 2.89e-01 | 121.6 | 2067.2 | 6.0 |
