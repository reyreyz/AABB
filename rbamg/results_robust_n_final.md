# Robust N (Algorithm 7.3 family): section 7.3 recirculating flow, Pe 20

Greedy F-points (delta = 0.3), Jacobi omega = 1, P0, injection, exact P', coarse
filter 1e-3 on levels < 8 and 3e-5 below (as in Algorithm 7.3), refinement step on
the listed levels. k = finest / coarse relaxation steps. Cycles to 1e-8 / total WU.

| configuration | m = 63 | m = 127 | m = 255 | m = 511 | stored nnz/nnz(A) at 511 |
|---|---|---|---|---|---|
| Alg. 7.3: k=8/12, refine 0-3 | 2 / 3811 | 2 / 7363 | 2 / 10938 | 2 / 14254 | 605 |
| **k=4/8, refine 0-1** | **2 / 952** | **2 / 1543** | **4 / 4092** | **4 / 4970** | 371 |
| k=4/12, refine 0-1 | 2 / 1312 | 2 / 2251 | 4 / 6345 | 4 / 7987 | 603 |
| k=4/6, refine 0-1 | 2 / 758 | 2 / 1167 | 4 / 3001 | 11 / 9749 | 259 |
| k=4/6, refine 0-2 | – | – | – | 10 / 16894 | 259 |

The Markovian cycle with the same hierarchy diverges at every size.
