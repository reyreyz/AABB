# Convergent N, then N_Z by a safeguarded theta search (recirculating flow, Pe 20, refinement on levels 0-3)

| m | cycle | theta | test contraction | factor | cycles to 1e-8 | WU/cycle | WU total | op.cx | setup [s] |
|---|---|---|---|---|---|---|---|---|---|
| 63 | N k=4/6 | – | 5.03e-06 | 9.65e-06 | 2 | 1087 | 2174 | 8.48 | 0.2 |
| 63 | N_Z k=4/6, Z on levels 1-2 | 0 | 5.03e-06 | 9.65e-06 | 2 | 1087 | 2174 | 8.48 | 0.2 (+8 search) |
| 63 | N k=3/4 | – | 0.0315 | 0.0323 | 5 | 787 | 3934 | 8.53 | 0.1 |
| 63 | N_Z k=3/4, Z on levels 1-2 | 0 | 0.0315 | 0.0323 | 5 | 787 | 3934 | 8.53 | 0.1 (+5 search) |
| 127 | N k=4/6 | – | 0.000115 | 0.000187 | 2 | 1590 | 3181 | 9.60 | 0.8 |
| 127 | N_Z k=4/6, Z on levels 1-2 | 0 | 0.000115 | 0.000187 | 2 | 1590 | 3181 | 9.60 | 0.8 (+45 search) |
| 127 | N k=3/4 | – | 0.171 | 0.23 | 11 | 1071 | 11782 | 9.75 | 0.7 |
| 127 | N_Z k=3/4, Z on levels 1-2 | 0 | 0.171 | 0.23 | 11 | 1071 | 11782 | 9.75 | 0.7 (+22 search) |
