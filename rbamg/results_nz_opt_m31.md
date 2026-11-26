# N_Z complement optimization (CD, m = 31)

### CD eps=1e-2, m = 31, k = 2 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 3.339e-01, N with Q0: rho = 3.719e-01, S (canonical): rho = 3.470e-01

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 4.014e-01 | 2.33e-01 | 2.028e-01 | 2.98e-01 | 2.6 | 15.7 |
| interp | ALS | 3.698e-01 | 2.34e-01 | 2.012e-01 | 4.46e-01 | 2.6 | 16.5 |
| interp | Kk | 3.706e-01 | 2.00e-01 | 1.714e-01 | 3.44e-01 | 2.6 | 18.8 |
| interp | noise | 5.631e-01 | 4.38e-01 | 3.239e-01 | 4.72e-01 | 2.6 | 29.5 |
| interp | EN_F | 3.769e-01 | 2.22e-01 | 1.805e-01 | 2.40e-01 | 2.6 | 75.1 |
| interp | EN_rho | 3.575e-01 | 6.32e-01 | 1.248e-01 | 1.03e+00 | 2.6 | 454.7 |
| interp | K_rho | 3.601e-01 | 1.27e+00 | 2.344e-01 | 2.12e+00 | 2.6 | 62.3 |

### CD eps=1e-2, m = 31, k = 4 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 7.835e-02, N with Q0: rho = 1.430e-01, S (canonical): rho = 1.204e-01

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 4.014e-01 | 4.23e-02 | 3.369e-02 | 5.61e-02 | 2.6 | 15.9 |
| interp | ALS | 4.250e-01 | 1.18e-01 | 3.547e-02 | 1.62e-01 | 2.6 | 14.4 |
| interp | Kk | 3.691e-01 | 3.31e-02 | 1.841e-02 | 3.01e-02 | 2.6 | 16.7 |
| interp | noise | 3.658e-01 | 4.48e-02 | 1.802e-02 | 5.92e-02 | 2.6 | 31.5 |
| interp | EN_F | 3.693e-01 | 3.35e-02 | 1.843e-02 | 3.04e-02 | 2.6 | 32.2 |
| interp | EN_rho | 3.689e-01 | 3.40e-02 | 1.834e-02 | 3.06e-02 | 2.6 | 44.0 |
| interp | K_rho | 3.601e-01 | 2.55e-01 | 2.132e-02 | 3.85e-01 | 2.6 | 104.2 |

### CD eps=1e-3, m = 31, k = 2 (bilinear P, Galerkin, Jacobi 2/3)

References: Markovian rho(E_M) = 3.407e-01, N with Q0: rho = 1.819e-01, S (canonical): rho = 2.070e-01

| pattern | objective | rho(K_Z) | ‖K_Z^k‖₂ | rho(E_N) | ‖E_N‖₂ | nnz(Z)/nf | time [s] |
|---|---|---|---|---|---|---|---|
| interp | LS | 4.554e-01 | 3.30e-01 | 2.462e-01 | 3.79e-01 | 2.6 | 32.6 |
| interp | ALS | 4.125e-01 | 4.05e-01 | 3.066e-01 | 6.78e-01 | 2.6 | 34.1 |
| interp | Kk | 4.164e-01 | 3.22e-01 | 2.057e-01 | 4.40e-01 | 2.6 | 39.5 |
| interp | noise | 4.435e-01 | 4.21e-01 | 1.961e-01 | 4.89e-01 | 2.6 | 48.6 |
| interp | EN_F | 4.042e-01 | 3.42e-01 | 1.792e-01 | 3.82e-01 | 2.6 | 183.5 |
