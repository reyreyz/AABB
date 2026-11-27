m=127, 5 smoothing levels + LU

| level | n | k | rho(T) | growth of C-error over k sweeps (random e) |
|---|---|---|---|---|
| 0 | 16129 | 4 | 4.78 | 359 |
| 1 | 8064 | 8 | 2.24 | 625 |
| 2 | 5099 | 8 | 1.43 | 15.3 |
| 3 | 3317 | 8 | 1.21 | 6.76 |
| 4 | 2304 | 8 | 1.16 | 1.77 |

| shift j=k applied on | factor | cycles |
|---|---|---|
| no level (Alg 7.3, j=0 everywhere) | 2.25e-06 | 2 |
| level 0 only | 5.89e-05 | 2 |
| levels >= 1 only | 3.61e-05 | 2 |
| all levels | 0.00323 | 4 |

m=255, 7 smoothing levels + LU

| level | n | k | rho(T) | growth of C-error over k sweeps (random e) |
|---|---|---|---|---|
| 0 | 65025 | 4 | 4.86 | 353 |
| 1 | 32535 | 8 | 2.29 | 723 |
| 2 | 20443 | 8 | 1.49 | 23.6 |
| 3 | 13118 | 8 | 1.33 | 15.1 |
| 4 | 9008 | 8 | 1.10 | 1.34 |
| 5 | 6472 | 8 | 1.11 | 0.817 |
| 6 | 4760 | 8 | 1.00 | 0.665 |

| shift j=k applied on | factor | cycles |
|---|---|---|
| no level (Alg 7.3, j=0 everywhere) | 0.00317 | 3 |
| level 0 only | 0.00314 | 3 |
| levels >= 1 only | 0.22 | 12 |
| all levels | inf | no conv. |

