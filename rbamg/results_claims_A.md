# Claims: N and its Markovian twin

## A. Equivalence of N and its Markovian twin (propagators and iterates)

| problem | complement Q | relaxation M | rho(T) | rho(E_N) | ‖E_N - E_twin‖ / max(‖E_N‖, ‖E_twin‖, 1) | ‖x_N - x_twin‖/‖x_N‖ (one cycle) |
|---|---|---|---|---|---|---|
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 (injective) | Jacobi 2/3 | 0.809 | 0.159 | 1.4e-15 | 3.7e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 (injective) | Jacobi 1 | 0.714 | 0.352 | 1.8e-15 | 6.1e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 (injective) | Gauss-Seidel | 0.523 | 0.0333 | 1.9e-15 | 3.3e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Jacobi 2/3 | 0.809 | 0.0486 | 3.7e-15 | 6.2e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Jacobi 1 | 0.714 | 0.83 | 5.4e-15 | 5.6e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Gauss-Seidel | 0.523 | 0.414 | 8.1e-15 | 1.3e-15 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 + P Z, random sparse | Jacobi 2/3 | 0.809 | 0.234 | 3.5e-15 | 8.5e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 + P Z, random sparse | Jacobi 1 | 0.714 | 0.619 | 3.4e-15 | 4.0e-16 |
| CD eps=1e-2 (bilinear P, R = P^T) | Q0 + P Z, random sparse | Gauss-Seidel | 0.523 | 0.1 | 2.8e-15 | 5.6e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 (injective) | Jacobi 2/3 | 0.545 | 0.0608 | 2.7e-15 | 4.3e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 (injective) | Jacobi 1 | 0.289 | 0.102 | 1.7e-15 | 4.2e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 (injective) | Gauss-Seidel | 0.11 | 0.172 | 1.3e-15 | 3.1e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Jacobi 2/3 | 0.545 | 0.105 | 2.6e-15 | 7.3e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Jacobi 1 | 0.289 | 0.876 | 1.4e-15 | 5.6e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Gauss-Seidel | 0.11 | 0.633 | 1.8e-15 | 1.1e-15 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 + P Z, random sparse | Jacobi 2/3 | 0.545 | 0.233 | 2.4e-15 | 5.2e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 + P Z, random sparse | Jacobi 1 | 0.289 | 0.539 | 2.8e-15 | 7.3e-16 |
| CD eps=1e-3 (bilinear P, R = P^T) | Q0 + P Z, random sparse | Gauss-Seidel | 0.11 | 0.487 | 2.0e-15 | 3.2e-16 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 (injective) | Jacobi 2/3 | 1.03 | 0.3 | 1.2e-14 | 1.5e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 (injective) | Jacobi 1 | 1.05 | 1.14 | 1.3e-14 | 1.5e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 (injective) | Gauss-Seidel | 1.09 | 0.0287 | 7.0e-15 | 2.6e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Jacobi 2/3 | 1.03 | 0.228 | 2.0e-14 | 1.9e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Jacobi 1 | 1.05 | 1.14 | 2.5e-14 | 2.1e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 + P Z, LS feedback | Gauss-Seidel | 1.09 | 0.105 | 2.6e-14 | 1.5e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 + P Z, random sparse | Jacobi 2/3 | 1.03 | 0.302 | 1.4e-14 | 2.0e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 + P Z, random sparse | Jacobi 1 | 1.05 | 1.14 | 1.1e-14 | 3.1e-15 |
| Shifted Laplacian (bilinear P, R = P^T) | Q0 + P Z, random sparse | Gauss-Seidel | 1.09 | 0.0315 | 1.4e-14 | 3.5e-15 |
| recirc. Pe 20 (greedy split, P0, injection) | Q0 (injective) | Jacobi 2/3 | 3.02 | 0.503 | 1.4e-14 | 2.3e-15 |
| recirc. Pe 20 (greedy split, P0, injection) | Q0 (injective) | Jacobi 1 | 4.5 | 0 | 3.0e-14 | 2.3e-15 |
| recirc. Pe 20 (greedy split, P0, injection) | Q0 (injective) | Gauss-Seidel | 20.2 | 8.27e+03 | 3.7e-18 | 5.2e-16 |
| recirc. Pe 20 (greedy split, P0, injection) | Q0 + P Z, random sparse | Jacobi 2/3 | 3.02 | 0.538 | 2.1e-14 | 1.1e-14 |
| recirc. Pe 20 (greedy split, P0, injection) | Q0 + P Z, random sparse | Jacobi 1 | 4.5 | 2.09 | 6.9e-15 | 1.1e-13 |
| recirc. Pe 20 (greedy split, P0, injection) | Q0 + P Z, random sparse | Gauss-Seidel | 20.2 | 8.27e+03 | 3.5e-16 | 6.0e-16 |
