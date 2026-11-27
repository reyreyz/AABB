# Claims: N and its Markovian twin


## F. Latency model: N (pipelined history) vs twin, recirculating flow Pe 20, P0, injection

| m | complement | k | pipelined N = N (iterates) | N: rounds on critical path | twin: rounds | N: flops / nnz(A) | twin: flops / nnz(A) |
|---|---|---|---|---|---|---|---|
| 127 | Q0 | 2 | 6.5e-16 | 5 | 4 | 4.90 | 2.00 |
| 127 | Q0 | 4 | 2.6e-14 | 7 | 6 | 7.90 | 3.00 |
| 127 | Q0 | 8 | 7.8e-13 | 11 | 10 | 13.90 | 5.00 |
| 127 | Q0 + P0 Z | 2 | 2.3e-15 | 7 | 6 | 12.77 | 6.35 |
| 127 | Q0 + P0 Z | 4 | 1.8e-14 | 9 | 10 | 32.90 | 18.24 |
| 127 | Q0 + P0 Z | 8 | 9.3e-13 | 13 | 18 | 98.21 | 58.69 |
| 255 | Q0 | 2 | 3.8e-15 | 5 | 4 | 4.90 | 2.00 |
| 255 | Q0 | 4 | 1.2e-13 | 7 | 6 | 7.90 | 3.00 |
| 255 | Q0 | 8 | 2.7e-11 | 11 | 10 | 13.90 | 5.00 |
| 255 | Q0 + P0 Z | 2 | 4.9e-15 | 7 | 6 | 12.79 | 6.37 |
| 255 | Q0 + P0 Z | 4 | 1.0e-14 | 9 | 10 | 33.27 | 18.53 |
| 255 | Q0 + P0 Z | 8 | 1.5e-11 | 13 | 18 | 101.81 | 61.17 |
