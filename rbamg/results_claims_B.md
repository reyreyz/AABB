# Claims: N and its Markovian twin


## B. Cost with Jacobi: sparse two-level N vs twin (recirculating flow, Pe 20, greedy split, P0, injection, exact coarse LU)

| m | complement | k | ‖x_N - x_twin‖/‖x_N‖ after 3 cycles | rel. residual | N: WU per cycle | twin: WU per cycle | N: ms per cycle | twin: ms per cycle |
|---|---|---|---|---|---|---|---|---|
| 63 | Q0 (twin = F-relaxation) | 4 | 1.5e-15 | 4.6e-15 | 15.39 | 10.89 | 1.7 | 0.3 |
| 63 | Q0 + P0 Z (bounded, CR-optimized) | 4 | 1.1e-15 | 1.3e-09 | 16.58 | 12.03 | 1.1 | 0.5 |
| 127 | Q0 (twin = F-relaxation) | 4 | 1.1e-14 | 7.8e-15 | 19.89 | 14.99 | 5.2 | 2.0 |
| 127 | Q0 + P0 Z (bounded, CR-optimized) | 4 | 4.8e-15 | 1.1e-05 | 107.33 | 98.92 | 21.1 | 13.1 |
| 255 | Q0 (twin = F-relaxation) | 4 | 6.7e-15 | 1.9e-08 | 24.77 | 19.87 | 21.2 | 11.9 |
| 255 | Q0 + P0 Z (bounded, CR-optimized) | 4 | 2.8e-15 | 2.9e-04 | 172.91 | 164.52 | 160.8 | 79.0 |
