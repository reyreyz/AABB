# Table SM4 (m = 15, k = 4, bilinear P): ours / paper


## Poisson: PBH 1.8e-04 / 1.8e-04, Markovian rho(E) 0.186 / 0.186

| complement | rho(K) | ‖K^4‖₂ | rho(E) | nnz(Z)/nf |
|---|---|---|---|---|
| Q0 | 0.648 / 0.648 | 0.176 / 0.18 | 0.176 / 0.176 | 0 |
| ker(Rhat A), S scheme | 0.647 / 0.647 | 0.177 / 0.18 | 0.176 / 0.176 | – |
| Q0 + PZ, LS | 0.552 / 0.552 | 0.0969 / 0.097 | 0.115 / 0.115 | 2.61 / 2.6 |
| Q0 + PZ, placement | failed (empty projection); paper: failed | | | |

## CD eps=1e-2: PBH 1.4e-01 / 1.4e-01, Markovian rho(E) 0.107 / 0.107

| complement | rho(K) | ‖K^4‖₂ | rho(E) | nnz(Z)/nf |
|---|---|---|---|---|
| Q0 | 0.563 / 0.563 | 0.261 / 0.26 | 0.0887 / 0.089 | 0 |
| ker(Rhat A), S scheme | 0.486 / 0.486 | 0.148 / 0.15 | 0.0557 / 0.056 | – |
| Q0 + PZ, LS | 0.347 / 0.347 | 0.0387 / 0.039 | 0.0166 / 0.017 | 2.61 / 2.6 |
| Q0 + PZ, placement (cond X 9.1e+03) | 0.05 / 0.05 | 0.00578 / 0.0021 | 0.000107 / 5.2e-05 | 64 / 64 |
| same, semi-Markovian | 0.05 | | 6.25e-06 / 6.2e-06 | |

## CD eps=1e-3: PBH 4.3e-01 / 6.6e-01, Markovian rho(E) 0.133 / 0.133

| complement | rho(K) | ‖K^4‖₂ | rho(E) | nnz(Z)/nf |
|---|---|---|---|---|
| Q0 | 0.388 / 0.388 | 0.306 / 0.31 | 0.0247 / 0.025 | 0 |
| ker(Rhat A), S scheme | 0.464 / 0.464 | 0.163 / 0.16 | 0.0465 / 0.047 | – |
| Q0 + PZ, LS | 0.445 / 0.445 | 0.0834 / 0.083 | 0.045 / 0.045 | 2.61 / 2.6 |
| Q0 + PZ, placement (cond X 3.0e+04) | 0.0855 / 0.05 | 0.0759 / 0.0038 | 0.000348 / 5.7e-05 | 64 / 64 |
| same, semi-Markovian | 0.0855 | | 5.34e-05 / 6.2e-06 | |
