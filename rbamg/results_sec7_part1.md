# Section 7 reproduction

## 7.1 Markovian V-cycle, stagnation flow, m = 180

Algorithm 7.1: 5 cycles, 26.0 WU, op. complexity 3.72  (paper: 4 cycles, 19.3 WU, 3.7)

best PyAMG AIR, F/C-relaxation (24 configs): 49.6 WU, 7 cycles, {'theta': 0.1, 'degree': 2, 'f': 2, 'c': 1}
best PyAMG AIR, F-relaxation (24 configs): 44.8 WU, 6 cycles, {'theta': 0.3, 'degree': 2, 'f': 3, 'c': 0}
(paper: AIR F/C best of 41 configs 149.6 WU / 12 cycles; F-relaxation, 64 configs: 57.4 WU / 7 cycles; Alg. 7.1 cheaper by 7.8 and 3.0)
ratio vs best F/C: 1.9, vs best F: 1.7

WU of PyAMG solvers = ml.cycle_complexity() x cycles (PyAMG's own count).

## Table 4: Alg. 7.2 on the shear flow, Alg. 7.3 on the recirculating flow

| m | n | 7.2 rho | WU/cyc | WU/digit | Markov 7.2 | 7.3 cycles to 1e-8 | WU/cyc | WU/digit | Markov 7.3 (1 cycle) |
|---|---|---|---|---|---|---|---|---|---|
| 127 | 1.6e+04 | 0.063 | 131 | 109 | diverges | 2 | 3524 | 610 | 7.3e+08 |
| 255 | 6.5e+04 | 0.083 | 139 | 129 | diverges | 2 | 5287 | 1011 | 1.8e+10 |
Killed
