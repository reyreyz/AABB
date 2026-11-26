# O(N) two-level example: centered CD, wind along x, Pe = 8, k = 4

stencil (W, Z): 48 numbers, fitted once on a 23 x 12 grid in 127 s; coarse solver: air (2 cycles)

| m | n | cycle | factor | cycles to 1e-10 | WU/cycle | stored nnz / nnz(A) | setup [µs/unknown] |
|---|---|---|---|---|---|---|---|
| 47 | 2209 | M | 24.7 | 10 (no conv.) | 5.4 | 0.80 | 11.5 |
| 47 | 2209 | N (Q0, P0) | 3.07e+21 | 1 (no conv.) | 21.5 | 3.14 | 24.8 |
| 47 | 2209 | N_WZ | 1.94e+15 | 1 (no conv.) | 230.0 | 88.34 | 134.3 |
| 95 | 9025 | M | 24.6 | 10 (no conv.) | 5.5 | 0.80 | 1.4 |
| 95 | 9025 | N (Q0, P0) | 2e+35 | 1 (no conv.) | 27.0 | 3.24 | 35.3 |
| 95 | 9025 | N_WZ | 2.8e+33 | 1 (no conv.) | 311.9 | 113.69 | 138.5 |
| 191 | 36481 | M | 24.7 | 9 (no conv.) | 5.6 | 0.80 | 1.1 |
| 191 | 36481 | N (Q0, P0) | 1.48e+52 | 1 (no conv.) | 30.3 | 3.29 | 24.4 |
| 191 | 36481 | N_WZ | 4.38e+52 | 1 (no conv.) | 369.6 | 127.09 | 222.8 |
| 383 | 146689 | M | 24.7 | 9 (no conv.) | 5.6 | 0.80 | 0.8 |
| 383 | 146689 | N (Q0, P0) | 6.71e+83 | 1 (no conv.) | 35.3 | 3.31 | 31.9 |
