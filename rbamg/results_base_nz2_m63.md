# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | M | 3 | 463 | 4 (no conv.) | 6 | inf | 1.10 | 1.7 | 0.1 |
| 63 | N k=4 | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.1 |
| 63 | N k=3 | 10 | 0.625 | 39 | 59 | 2296 | 8.42 | 42.0 | 0.1 |
| 63 | N k=2 | 10 | 0.88 | 40 (no conv.) | 42 | inf | 8.02 | 31.4 | 0.1 |
| 63 | N_Z k=4, finest | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.9 |
| 63 | N_Z k=3, finest | 10 | 0.625 | 39 | 59 | 2296 | 8.42 | 42.0 | 0.8 |
| 63 | N_Z k=2, finest | 10 | 0.88 | 40 (no conv.) | 42 | inf | 8.02 | 31.4 | 0.8 |
| 63 | N_Z k=4, finest, theta auto | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.9 |
| 63 | N_Z k=3, finest, theta auto | 10 | 0.625 | 39 | 59 | 2296 | 8.42 | 42.0 | 1.0 |
| 63 | N_Z k=4, all, theta auto | 12 | 0.18 | 12 | 203 | 2441 | 12.95 | 149.5 | 3.3 |
| 63 | N_Z k=3, all, theta auto | 11 | 3.41 | 18 (no conv.) | 145 | inf | 11.17 | 109.9 | 3.4 |
