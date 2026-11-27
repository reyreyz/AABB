# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | M | 3 | 463 | 4 (no conv.) | 6 | inf | 1.10 | 1.7 | 0.1 |
| 63 | N k=4 | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.2 |
| 63 | N_Z k=4, Z on level 1 | 11 | 1.85 | 39 (no conv.) | 115 | inf | 10.71 | 89.7 | 0.5 |
| 63 | N_Z k=4, Z on level 1, theta auto | 10 | 0.373 | 19 | 107 | 2038 | 9.87 | 81.2 | 0.6 |
| 63 | N_Z k=4, Z on levels 1-2, theta auto | 11 | 0.553 | 31 | 147 | 4559 | 11.16 | 115.3 | 1.1 |
| 63 | N k=3 | 10 | 0.625 | 39 | 59 | 2296 | 8.42 | 42.0 | 0.1 |
| 63 | N_Z k=3, Z on level 1, theta auto | 10 | 0.149 | 10 | 77 | 774 | 9.39 | 60.3 | 0.9 |
| 63 | N_Z k=3, Z on levels 1-2, theta auto | 10 | 0.191 | 11 | 95 | 1045 | 9.48 | 75.4 | 1.5 |
