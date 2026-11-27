# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | M | 3 | 463 | 4 (no conv.) | 6 | inf | 1.10 | 1.7 | 0.1 |
| 63 | N k=4, dc=0.3 | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.2 |
| 63 | N_Z k=4, dc=0.3, Z on A_cf | 13 | 3.55e+04 | 2 (no conv.) | 294 | inf | 12.98 | 213.7 | 3.5 |
| 63 | N_Z k=4, dc=0.3, Z on A_cf, theta 1 | 18 | 18.9 | 3 (no conv.) | 469 | inf | 22.16 | 335.4 | 2.4 |
| 63 | N k=3, dc=0.3 | 10 | 0.625 | 39 | 59 | 2296 | 8.42 | 42.0 | 0.2 |
| 63 | N_Z k=3, dc=0.3, Z on A_cf | 14 | 45.4 | 4 (no conv.) | 334 | inf | 18.60 | 257.0 | 7.0 |
| 63 | N k=4, dc=0.6 | 7 | 0.917 | 40 (no conv.) | 61 | inf | 6.29 | 36.0 | 0.1 |
| 63 | N_Z k=4, dc=0.6, Z on A_cf | 23 | 321 | 3 (no conv.) | 534 | inf | 44.30 | 286.4 | 20.3 |
| 63 | N k=4, dc=0.45 | 8 | 0.808 | 40 (no conv.) | 63 | inf | 6.76 | 38.9 | 0.1 |
| 63 | N_Z k=4, dc=0.45, Z on A_cf | 29 | 33.5 | 5 (no conv.) | 673 | inf | 56.18 | 379.1 | 9.9 |
