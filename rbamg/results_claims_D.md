# Claims: N and its Markovian twin


## D. Black-box smoother (iterate generator, not restartable): only N converges

| m | n | method | uses only the black-box iterates? | factor | cycles to 1e-8 | WU per cycle | note |
|---|---|---|---|---|---|---|---|
| 127 | 16129 | N (Algorithm 4.2 on the iterates) | yes | 2.22e-06 | 2 | 386 | black-box calls on level 0: 2 |
| 127 | 16129 | Markovian, non-intrusive (x^k + coarse correction) | yes | 5.23e+09 | 2 (no conv.) | 201 | black-box calls on level 0: 2 |
| 127 | 16129 | twin M_F(4,0) (needs k = 4 single-step restarts per cycle) | no | 1.43e-05 | 2 | 87 | not available under the black-box constraint |
| 255 | 65025 | N (Algorithm 4.2 on the iterates) | yes | 0.00317 | 3 | 555 | black-box calls on level 0: 3 |
| 255 | 65025 | Markovian, non-intrusive (x^k + coarse correction) | yes | 1.07e+12 | 1 (no conv.) | 271 | black-box calls on level 0: 1 |
| 255 | 65025 | twin M_F(4,0) (needs k = 4 single-step restarts per cycle) | no | 0.00317 | 3 | 105 | not available under the black-box constraint |
| 511 | 261121 | N (Algorithm 4.2 on the iterates) | yes | 0.00317 | 3 | 702 | black-box calls on level 0: 3 |
| 511 | 261121 | Markovian, non-intrusive (x^k + coarse correction) | yes | 2.48e+12 | 1 (no conv.) | 333 | black-box calls on level 0: 1 |
| 511 | 261121 | twin M_F(4,0) (needs k = 4 single-step restarts per cycle) | no | 0.00317 | 3 | 113 | not available under the black-box constraint |
