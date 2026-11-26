# N_Z complement optimization (CD, m = 15)

### CD eps=1e-2, m = 15, k = 2 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 2.112e-01, N with Q0: rho = 2.881e-01, S (canonical): rho = 2.359e-01

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 3.475e-01 | 2.31e-01 | 1.406e-01 | 2.72e-01 | 2.6 | 3.3 |
| interp | ALS | 3.506e-01 | 2.63e-01 | 2.163e-01 | 5.73e-01 | 2.6 | 3.4 |
| interp | Kk | 3.452e-01 | 1.96e-01 | 1.285e-01 | 2.96e-01 | 2.6 | 3.9 |
| interp | noise | 4.670e-01 | 4.03e-01 | 2.192e-01 | 4.51e-01 | 2.6 | 3.5 |
| interp | EN_F | 3.506e-01 | 2.28e-01 | 1.231e-01 | 2.47e-01 | 2.6 | 10.4 |
| interp | EN_rho | 3.192e-01 | 7.28e-01 | 8.895e-02 | 1.06e+00 | 2.6 | 65.2 |
| interp | K_rho | 3.090e-01 | 1.29e+00 | 1.568e-01 | 2.05e+00 | 2.6 | 15.3 |
| interp (dup) | LS | 3.475e-01 | 2.31e-01 | 1.406e-01 | 2.72e-01 | 2.6 | 0.8 |
| interp (dup) | ALS | 3.506e-01 | 2.63e-01 | 2.163e-01 | 5.73e-01 | 2.6 | 1.2 |
| interp (dup) | Kk | 3.452e-01 | 1.96e-01 | 1.285e-01 | 2.96e-01 | 2.6 | 3.8 |
| interp (dup) | noise | 4.670e-01 | 4.03e-01 | 2.192e-01 | 4.51e-01 | 2.6 | 4.7 |
| interp (dup) | EN_F | 3.506e-01 | 2.28e-01 | 1.231e-01 | 2.47e-01 | 2.6 | 8.1 |
| interp (dup) | EN_rho | 3.192e-01 | 7.28e-01 | 8.895e-02 | 1.06e+00 | 2.6 | 67.0 |
| interp (dup) | K_rho | 3.090e-01 | 1.29e+00 | 1.568e-01 | 2.05e+00 | 2.6 | 15.5 |

### CD eps=1e-2, m = 15, k = 4 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 1.073e-01, N with Q0: rho = 8.872e-02, S (canonical): rho = 5.566e-02

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 3.475e-01 | 3.87e-02 | 1.662e-02 | 3.65e-02 | 2.6 | 3.0 |
| interp | ALS | 3.419e-01 | 5.33e-02 | 1.708e-02 | 9.12e-02 | 2.6 | 3.3 |
| interp | Kk | 3.433e-01 | 3.37e-02 | 1.395e-02 | 3.30e-02 | 2.6 | 1.6 |
| interp | noise | 3.487e-01 | 4.93e-02 | 1.504e-02 | 5.18e-02 | 2.6 | 5.5 |
| interp | EN_F | 3.453e-01 | 3.58e-02 | 1.428e-02 | 3.40e-02 | 2.6 | 3.3 |
| interp | EN_rho | 3.304e-01 | 5.25e-02 | 1.180e-02 | 9.82e-02 | 2.6 | 8.7 |
| interp | K_rho | 3.090e-01 | 1.54e-01 | 1.184e-02 | 2.06e-01 | 2.6 | 13.8 |
| interp (dup) | LS | 3.475e-01 | 3.87e-02 | 1.662e-02 | 3.65e-02 | 2.6 | 1.1 |
| interp (dup) | ALS | 3.419e-01 | 5.33e-02 | 1.708e-02 | 9.12e-02 | 2.6 | 1.8 |
| interp (dup) | Kk | 3.433e-01 | 3.37e-02 | 1.395e-02 | 3.30e-02 | 2.6 | 3.1 |
| interp (dup) | noise | 3.487e-01 | 4.93e-02 | 1.504e-02 | 5.18e-02 | 2.6 | 5.8 |
| interp (dup) | EN_F | 3.453e-01 | 3.58e-02 | 1.428e-02 | 3.40e-02 | 2.6 | 5.5 |
| interp (dup) | EN_rho | 3.304e-01 | 5.25e-02 | 1.180e-02 | 9.82e-02 | 2.6 | 7.4 |
| interp (dup) | K_rho | 3.090e-01 | 1.54e-01 | 1.184e-02 | 2.06e-01 | 2.6 | 16.8 |

### CD eps=1e-3, m = 15, k = 2 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 3.602e-01, N with Q0: rho = 1.501e-01, S (canonical): rho = 2.157e-01

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 4.446e-01 | 3.35e-01 | 2.303e-01 | 3.83e-01 | 2.6 | 2.9 |
| interp | ALS | 4.061e-01 | 4.06e-01 | 2.412e-01 | 5.58e-01 | 2.6 | 1.6 |
| interp | Kk | 4.083e-01 | 3.25e-01 | 1.953e-01 | 4.13e-01 | 2.6 | 3.7 |
| interp | noise | 4.227e-01 | 4.18e-01 | 1.763e-01 | 4.67e-01 | 2.6 | 3.2 |
| interp | EN_F | 3.892e-01 | 3.45e-01 | 1.637e-01 | 3.74e-01 | 2.6 | 8.6 |
| interp | EN_rho | 3.267e-01 | 1.73e+00 | 8.194e-02 | 2.33e+00 | 2.6 | 87.3 |
| interp | K_rho | 3.119e-01 | 2.42e+00 | 2.262e-01 | 4.67e+00 | 2.6 | 19.7 |
| interp (dup) | LS | 4.446e-01 | 3.35e-01 | 2.303e-01 | 3.83e-01 | 2.6 | 1.1 |
| interp (dup) | ALS | 4.061e-01 | 4.06e-01 | 2.412e-01 | 5.58e-01 | 2.6 | 1.6 |
| interp (dup) | Kk | 4.083e-01 | 3.25e-01 | 1.953e-01 | 4.13e-01 | 2.6 | 2.0 |
| interp (dup) | noise | 4.227e-01 | 4.18e-01 | 1.763e-01 | 4.67e-01 | 2.6 | 4.7 |
| interp (dup) | EN_F | 3.892e-01 | 3.45e-01 | 1.637e-01 | 3.74e-01 | 2.6 | 9.4 |
| interp (dup) | EN_rho | 3.267e-01 | 1.73e+00 | 8.194e-02 | 2.33e+00 | 2.6 | 85.1 |
| interp (dup) | K_rho | 3.119e-01 | 2.42e+00 | 2.262e-01 | 4.67e+00 | 2.6 | 18.1 |

### CD eps=1e-3, m = 15, k = 4 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 1.327e-01, N with Q0: rho = 2.466e-02, S (canonical): rho = 4.652e-02

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 4.446e-01 | 8.34e-02 | 4.498e-02 | 8.83e-02 | 2.6 | 3.1 |
| interp | ALS | 4.106e-01 | 1.70e-01 | 3.222e-02 | 3.08e-01 | 2.6 | 3.7 |
| interp | Kk | 3.264e-01 | 5.95e-02 | 1.162e-02 | 8.40e-02 | 2.6 | 6.0 |
| interp | noise | 3.341e-01 | 7.91e-02 | 1.289e-02 | 7.46e-02 | 2.6 | 7.9 |
| interp | EN_F | 3.365e-01 | 6.70e-02 | 1.289e-02 | 6.96e-02 | 2.6 | 20.9 |
| interp | EN_rho | 3.114e-01 | 3.27e-01 | 8.175e-03 | 3.83e-01 | 2.6 | 32.6 |
| interp | K_rho | 3.119e-01 | 6.00e-01 | 2.170e-02 | 1.04e+00 | 2.6 | 16.1 |
| interp (dup) | LS | 4.446e-01 | 8.34e-02 | 4.498e-02 | 8.83e-02 | 2.6 | 3.0 |
| interp (dup) | ALS | 4.106e-01 | 1.70e-01 | 3.222e-02 | 3.08e-01 | 2.6 | 3.6 |
| interp (dup) | Kk | 3.264e-01 | 5.95e-02 | 1.162e-02 | 8.40e-02 | 2.6 | 6.0 |
| interp (dup) | noise | 3.341e-01 | 7.91e-02 | 1.289e-02 | 7.46e-02 | 2.6 | 8.1 |
| interp (dup) | EN_F | 3.365e-01 | 6.70e-02 | 1.289e-02 | 6.96e-02 | 2.6 | 26.2 |
| interp (dup) | EN_rho | 3.114e-01 | 3.27e-01 | 8.175e-03 | 3.83e-01 | 2.6 | 29.5 |
| interp (dup) | K_rho | 3.119e-01 | 6.00e-01 | 2.170e-02 | 1.04e+00 | 2.6 | 17.9 |
