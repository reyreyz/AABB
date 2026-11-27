# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | M | 3 | 463 | 4 (no conv.) | 6 | inf | 1.10 | 1.7 | 0.1 |
| 63 | N k=4, dc=0.3 | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.2 |
| 63 | N_Z k=4, dc=0.3, full P'_C | 12 | 0.18 | 12 | 203 | 2441 | 12.95 | 149.5 | 3.2 |
| 63 | N_Z k=4, dc=0.3, P'_C mask acc | 10 | 5.51 | 14 (no conv.) | 112 | inf | 9.18 | 92.7 | 2.9 |
| 63 | N_Z k=4, dc=0.3, P'_C mask diag | 9 | 2.8 | 22 (no conv.) | 97 | inf | 7.97 | 80.7 | 3.0 |
| 63 | N k=3, dc=0.3 | 10 | 0.625 | 39 | 59 | 2296 | 8.42 | 42.0 | 0.1 |
| 63 | N_Z k=3, dc=0.3, P'_C mask diag | 9 | 0.946 | 40 (no conv.) | 72 | inf | 7.78 | 60.3 | 3.7 |
| 63 | N k=4, dc=0.6 | 7 | 0.917 | 40 (no conv.) | 61 | inf | 6.29 | 36.0 | 0.1 |
| 63 | N_Z k=4, dc=0.6, P'_C mask diag | 7 | 3.5 | 20 (no conv.) | 95 | inf | 7.02 | 68.9 | 39.6 |
| 63 | N_Z k=4, dc=0.6, P'_C mask acc | 8 | 12 | 10 (no conv.) | 97 | inf | 7.93 | 66.8 | 33.4 |
