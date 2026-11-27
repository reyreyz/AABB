# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 127 | M | 3 | 464 | 4 (no conv.) | 6 | inf | 1.10 | 1.9 | 0.2 |
| 127 | N k=4 | 16 | 0.623 | 34 | 94 | 3184 | 9.75 | 77.1 | 0.7 |
| 127 | N k=3 | 16 | 2.67 | 25 (no conv.) | 73 | inf | 9.89 | 60.0 | 0.6 |
| 127 | N_Z k=3, Z on level 1, theta auto | 16 | 0.57 | 30 | 91 | 2728 | 10.67 | 78.1 | 16.4 |
| 127 | N_Z k=4, Z on level 1, theta auto | 16 | 0.462 | 23 | 122 | 2813 | 10.57 | 104.9 | 2.8 |
| 127 | N_Z k=3, Z on levels 1-2, theta auto | 17 | 0.31 | 17 | 118 | 2001 | 11.97 | 106.4 | 22.1 |
| 255 | M | 3 | 954 | 4 (no conv.) | 6 | inf | 1.10 | 1.9 | 0.8 |
| 255 | N k=4 | 22 | 0.896 | 40 (no conv.) | 103 | inf | 10.32 | 95.5 | 3.3 |
| 255 | N k=3 | 23 | 29.1 | 8 (no conv.) | 78 | inf | 10.42 | 70.1 | 2.5 |
| 255 | N_Z k=3, Z on level 1, theta auto | 35 | 30.1 | 8 (no conv.) | 118 | inf | 14.46 | 112.6 | 80.4 |
| 255 | N_Z k=4, Z on level 1, theta auto | 25 | 68.8 | 6 (no conv.) | 136 | inf | 11.59 | 128.8 | 14.5 |
| 255 | N_Z k=3, Z on levels 1-2, theta auto | 51 | 72.7 | 6 (no conv.) | 159 | inf | 17.65 | 150.5 | 209.3 |
