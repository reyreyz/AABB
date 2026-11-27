# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | M | 3 | 463 | 4 (no conv.) | 6 | inf | 1.10 | 1.7 | 0.1 |
| 63 | N k=4, dc=0.3 | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.2 |
| 63 | N k=4, dc=0.6 | 7 | 0.917 | 40 (no conv.) | 61 | inf | 6.29 | 36.0 | 0.1 |
| 63 | N_Z k=4, dc=0.6 | 9 | 0.564 | 30 | 192 | 5757 | 10.91 | 126.3 | 18.5 |
| 63 | N k=4, dc=1.0 | 5 | 0.991 | 40 (no conv.) | 41 | inf | 3.95 | 18.6 | 0.1 |
| 63 | N_Z k=4, dc=1.0 | 9 | 0.976 | 40 (no conv.) | 165 | inf | 11.09 | 99.2 | 49.2 |
| 63 | N k=3, dc=0.6 | 8 | 0.792 | 40 (no conv.) | 48 | inf | 6.43 | 29.9 | 0.1 |
| 63 | N_Z k=3, dc=0.6 | 10 | 0.654 | 40 (no conv.) | 166 | inf | 12.62 | 119.9 | 21.1 |
