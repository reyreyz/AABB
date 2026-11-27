# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | N k=4, dc=0.3, tauA=1e-3 | 10 | 0.271 | 13 | 75 | 979 | 8.53 | 50.4 | 0.1 |
| 63 | N k=4, dc=0.3, tauA=1e-2 | 10 | 0.514 | 23 | 49 | 1120 | 5.55 | 31.3 | 0.1 |
| 63 | N_Z k=4, dc=0.3, tauA=1e-2 | 10 | 1.41 | 40 (no conv.) | 101 | inf | 6.15 | 74.0 | 2.4 |
| 63 | N_Z k=4, dc=0.3, tauA=1e-2, capP=12 | 10 | 0.915 | 40 (no conv.) | 55 | inf | 6.16 | 24.7 | 2.6 |
| 63 | N k=4, dc=0.6, tauA=1e-2 | 7 | 0.918 | 40 (no conv.) | 42 | inf | 4.35 | 24.3 | 0.1 |
| 63 | N_Z k=4, dc=0.6, tauA=1e-2 | 10 | 3.41 | 19 (no conv.) | 165 | inf | 7.84 | 121.2 | 18.0 |
| 63 | N_Z k=4, dc=0.6, tauA=3e-2 | 9 | 10.5 | 10 (no conv.) | 127 | inf | 5.19 | 96.1 | 15.7 |
| 63 | N_Z k=4, dc=0.6, tauA=1e-2, capP=12 | 9 | 2.08e+04 | 3 (no conv.) | 58 | inf | 6.52 | 25.2 | 17.2 |
