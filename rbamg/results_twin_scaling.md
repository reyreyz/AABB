stagnation: best pyAMG at m=180: AIR {'strength': ('classical', {'theta': 0.25, 'norm': 'min'}), 'restrict': ('air', {'theta': 0.05, 'degree': 1}), 'postsmoother': ('fc_jacobi', {'omega': 1.0, 'iterations': 1, 'withrho': False, 'f_iterations': 2, 'c_iterations': 1}), 'smoother': ('fc_jacobi', 1)} -> 149.6 WU, 12 cycles
| stagnation | 90 | 8100 | twin: 18.4 WU / 4 cyc / oc 3.48 | M direct: 51.6 WU / 24 cyc / oc 2.73 | pyAMG: 125.8 WU / 11 cyc / oc 3.17 |
| stagnation | 180 | 32400 | twin: 19.3 WU / 4 cyc / oc 3.72 | M direct: 91.0 WU / 40 cyc / oc 2.86 | pyAMG: 149.6 WU / 12 cyc / oc 3.49 |
| stagnation | 360 | 129600 | twin: 27.0 WU / 5 cyc / oc 4.07 | M direct: inf WU / 300 cyc / oc 2.95 | pyAMG: 185.2 WU / 14 cyc / oc 3.72 |
| stagnation | 720 | 518400 | twin: 28.0 WU / 5 cyc / oc 4.22 | M direct: inf WU / 12 cyc / oc 3.03 | pyAMG: 207.5 WU / 15 cyc / oc 3.90 |
curved: best pyAMG at m=180: RS {'strength': ('classical', {'theta': 0.5, 'norm': 'min'}), 'smoother': ('gauss_seidel', 2)} -> 159.1 WU, 14 cycles
| curved | 90 | 8100 | twin: 19.9 WU / 4 cyc / oc 3.75 | M direct: 57.3 WU / 27 cyc / oc 2.80 | pyAMG: 112.1 WU / 10 cyc / oc 1.99 |
| curved | 180 | 32400 | twin: 21.1 WU / 4 cyc / oc 4.09 | M direct: 89.9 WU / 41 cyc / oc 2.88 | pyAMG: 159.1 WU / 14 cyc / oc 2.01 |
| curved | 360 | 129600 | twin: 28.2 WU / 5 cyc / oc 4.22 | M direct: 169.0 WU / 75 cyc / oc 2.93 | pyAMG: 243.8 WU / 21 cyc / oc 2.06 |
| curved | 720 | 518400 | twin: 29.7 WU / 5 cyc / oc 4.42 | M direct: inf WU / 20 cyc / oc 2.97 | pyAMG: 407.5 WU / 35 cyc / oc 2.06 |
jump: best pyAMG at m=180: AIR {'strength': ('classical', {'theta': 0.3, 'norm': 'min'}), 'restrict': ('air', {'theta': 0.05, 'degree': 2}), 'postsmoother': ('fc_jacobi', {'omega': 1.0, 'iterations': 1, 'withrho': False, 'f_iterations': 2, 'c_iterations': 1}), 'smoother': ('fc_jacobi', 1)} -> 90.3 WU, 8 cycles
| jump | 90 | 8100 | twin: 12.8 WU / 4 cyc / oc 2.44 | M direct: 35.6 WU / 20 cyc / oc 2.20 | pyAMG: 85.5 WU / 8 cyc / oc 3.14 |
| jump | 180 | 32400 | twin: 12.9 WU / 4 cyc / oc 2.61 | M direct: 55.7 WU / 30 cyc / oc 2.28 | pyAMG: 90.3 WU / 8 cyc / oc 3.33 |
| jump | 360 | 129600 | twin: 15.7 WU / 4 cyc / oc 2.98 | M direct: 83.4 WU / 43 cyc / oc 2.61 | pyAMG: 105.5 WU / 9 cyc / oc 3.46 |
| jump | 720 | 518400 | twin: 21.0 WU / 5 cyc / oc 3.25 | M direct: 142.0 WU / 72 cyc / oc 2.68 | pyAMG: 105.9 WU / 9 cyc / oc 3.47 |

[exited with code 0]

Jump with the exact twin configuration (one-point W0, dP = nu = 3, filt 1e-4):
| jump | 90 | 13.0 WU / 4 cyc / oc 2.44 |
| jump | 180 | 13.0 WU / 4 cyc / oc 2.60 |
| jump | 360 | 16.2 WU / 4 cyc / oc 2.98 |
| jump | 720 | 17.2 WU / 4 cyc / oc 3.23 |
