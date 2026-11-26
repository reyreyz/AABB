Problems marked * use our own discretization (not specified in the paper).

## 1. Algorithms vs closed forms (m=15, k=3, bilinear P)

| problem | M | M' | S | N | exact: ||e_out||/||e0|| |
|---|---|---|---|---|---|
| Poisson | 4.8e-15 | 5.2e-15 | 4.9e-15 | 5.7e-15 | 1.3e-15 |
| Random diffusion* | 1.4e-14 | 3.5e-14 | 9.7e-15 | 2.3e-14 | 5.5e-14 |
| CD eps=1e-2* | 1.8e-15 | 2.4e-15 | 2.8e-15 | 2.0e-15 | 1.3e-15 |
| CD eps=1e-3* | 1.1e-15 | 1.3e-15 | 2.3e-15 | 1.7e-15 | 1.1e-15 |
| Shifted Laplacian* | 1.1e-14 | 1.7e-14 | 1.1e-14 | 2.2e-14 | 3.8e-15 |

## 2. Structural statements (m=15)

| problem | k | (4.4) split | rho(E^S) vs rho(K_can)^k | k=1: rho(E^S)-rho(E^M) | rho(E^S) <= rho(E^M) |
|---|---|---|---|---|---|
| Poisson | 1 | 1.6e-15 | 1.7e-15 | 1.4e-15 |  |
| Poisson | 2 | 1.7e-15 | 2.8e-15 |  | 0.4190 <= 0.4231: True |
| Poisson | 4 | 2.8e-15 | 1.6e-16 |  | 0.1756 <= 0.1861: True |
| Random diffusion* | 1 | 1.8e-15 | 2.9e-15 | -1.1e-16 |  |
| Random diffusion* | 2 | 3.2e-15 | 2.0e-15 |  |  |
| Random diffusion* | 4 | 4.1e-15 | 1.2e-14 |  |  |
| CD eps=1e-2* | 1 | 6.7e-16 | 1.3e-11 | 6.2e-12 |  |
| CD eps=1e-2* | 2 | 7.6e-16 | 4.4e-12 |  |  |
| CD eps=1e-2* | 4 | 7.5e-16 | 3.8e-12 |  |  |
| CD eps=1e-3* | 1 | 6.7e-16 | 6.4e-15 | 2.8e-15 |  |
| CD eps=1e-3* | 2 | 5.8e-16 | 1.6e-14 |  |  |
| CD eps=1e-3* | 4 | 5.9e-16 | 2.6e-14 |  |  |
| Shifted Laplacian* | 1 | 2.8e-15 | 9.9e-16 | 7.8e-16 |  |
| Shifted Laplacian* | 2 | 5.2e-15 | 2.4e-15 |  |  |
| Shifted Laplacian* | 4 | 8.8e-15 | 8.4e-15 |  |  |

## 3. Table SM1, Poisson, bilinear P, m=31

| k | | M | M' | S | N | rho(K)^k |
|---|---|---|---|---|---|---|
| 1 | ours | 0.662 | 0.661 | 0.662 | 0.662 | 0.662 |
| 1 | paper | 0.662 | 0.661 | 0.662 | 0.662 | 0.662 |
| 2 | ours | 0.439 | 0.437 | 0.438 | 0.438 | 0.438 |
| 2 | paper | 0.439 | 0.437 | 0.438 | 0.438 | 0.438 |
| 4 | ours | 0.195 | 0.192 | 0.192 | 0.192 | 0.192 |
| 4 | paper | 0.195 | 0.192 | 0.192 | 0.192 | 0.192 |

worst identity deviation (1, 2): 1.3e-11
worst |ours - paper| on Table SM1 Poisson: 0.0003
PASS
