# 2D multilevel: convergent N, then optimized N_Z (recirculating flow, Pe 20)

| m | cycle | levels | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | stored | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 127 | M | 3 | 464 | 4 (no conv.) | 6 | inf | 1.10 | 1.9 | 0.2 |
| 127 | N k=4 | 16 | 0.623 | 34 | 94 | 3184 | 9.75 | 77.1 | 0.7 |
| 127 | N k=3 | 16 | 2.67 | 25 (no conv.) | 73 | inf | 9.89 | 60.0 | 0.6 |
