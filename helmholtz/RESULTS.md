# 1D Helmholtz: pseudomodes and true sparse-interpolation optima

Small grids, exact dense linear algebra, pure NumPy (scipy was unavailable).

**Setup.** A = tridiag(-1, 2, -1) - (kh)^2 I (scaled by h^2, Dirichlet, n = 31, so k = 32 kh).
Smoother: damped Jacobi, omega = 2/3. Its diagonal is constant, so S = I - (omega/d) A is a polynomial in A and is symmetric.
Coarse points are every other point (n_c = 15). Column j of P is supported on w = 3, 5 or 7 points centred at the coarse point.
Two-grid operator: E = S (I - P (P^T A P)^{-1} P^T A) S, with an exact Galerkin coarse solve.
A is indefinite, so there is no energy norm. Two figures are reported:
- rho(E): the asymptotic rate.
- ||E||_2: the worst-case contraction of one cycle.

"True optimum" means a numerical minimisation over every entry of P inside the given sparsity pattern:
- BFGS with finite-difference gradients;
- a smooth surrogate (||E^m||_F^{1/m} for rho, a Schatten-p norm for ||E||_2), with m and p increased step by step;
- 4 to 8 starting points each, including every constructed P.

For width 3, the minimisation over translation-invariant stencils was also done by exhaustive grid search.

| kh (points per wavelength) | 0 | 0.3 (21) | 0.6 (10.5) | 0.9 (7) |
|---|---|---|---|---|
| negative eigenvalues of A | 0 | 3 | 6 | 9 |
| hat: z_p / eps_p | 0.78 / 0.25 | 0.80 / 0.26 | 0.88 / 0.30 | 1.08 / 0.42 |
| hat weight on amplified modes (s > 1) | 0 | 33 % | 47 % | 72 % |

## Main table: rho / ||E||_2

| interpolation | width | kh=0.0 | kh=0.3 | kh=0.6 | kh=0.9 |
|---|---|---|---|---|---|
| linear (hat) | 3 | 0.111 / 0.35 | 0.701 / 13.93 | 5.130 / 21.67 | 2.778 / 5.77 |
| local pseudomode opt (each column) | 3 | 0.842 / 2.65 | 2.962 / 12.31 | 0.580 / 2.98 | 0.554 / 2.40 |
| constr. min A^2, B = 1 near-kernel eigvec | 3 | 0.115 / 0.26 | 0.116 / 0.83 | sing. | 0.547 / 1.36 |
| constr. min A^2, B = 2 near-kernel eigvecs | 3 | 0.189 / 0.27 | 0.200 / 0.46 | 0.187 / 0.74 | 0.238 / 1.90 |
| constr. min A^2, B = one-sided sines (BC-aware) | 3 | 0.111 / 0.35 | 0.111 / 0.40 | 0.111 / 0.67 | 0.111 / 1.90 |
| **min rho, all P with this support** | 3 | 0.111 / 0.28 | 0.111 / 5.63 | 0.111 / 2.26 | 0.111 / 2.38 |
| **min ||E||2, all P with this support** | 3 | 0.142 / 0.23 | 0.170 / 0.33 | 0.193 / 0.48 | 0.514 / 0.91 |
| local pseudomode opt (each column) | 5 | 0.402 / 3.16 | 1.063 / 16.97 | 0.213 / 1.21 | 0.604 / 1.34 |
| constr. min A^2, B = 1 near-kernel eigvec | 5 | 0.112 / 0.19 | 0.133 / 1.10 | sing. | 4.858 / 8.24 |
| constr. min A^2, B = 2 near-kernel eigvecs | 5 | 0.174 / 0.21 | 0.171 / 0.35 | 0.198 / 0.40 | 0.606 / 0.64 |
| constr. min A^2, B = one-sided sines (BC-aware) | 5 | 0.111 / 0.26 | 0.120 / 0.30 | 0.191 / 0.40 | 0.605 / 0.64 |
| **min ||E||2, all P with this support** | 5 | 0.111 / 0.11 | 0.131 / 0.13 | 0.216 / 0.22 | 0.559 / 0.56 |
| **min rho, all P with this support** | 5 | 1/9 (SPD bound) | 0.111 / 5.12 | 0.111 / 2.01 | 0.111 / 4.81 |
| local pseudomode opt (each column) | 7 | 0.240 / 1.86 | 0.243 / 2.27 | 0.288 / 1.13 | 0.609 / 0.86 |
| constr. min A^2, B = 1 near-kernel eigvec | 7 | 0.281 / 0.29 | 0.264 / 0.64 | sing. | 3.981 / 6.48 |
| constr. min A^2, B = 2 near-kernel eigvecs | 7 | 0.221 / 0.22 | 0.226 / 0.27 | 0.226 / 0.23 | 0.509 / 1.19 |
| constr. min A^2, B = one-sided sines (BC-aware) | 7 | 0.281 / 0.30 | 0.305 / 0.34 | 0.226 / 0.24 | 0.873 / 1.66 |
| **min ||E||2, all P with this support** | 7 | 0.111 / 0.11 | 0.131 / 0.13 | 0.215 / 0.22 | 0.546 / 0.55 |
| **min rho, all P with this support** | 7 | 1/9 (SPD bound) | 0.111 / 10.28 | 0.111 / 5.44 | 0.111 / 2.40 |
eigenvector coarse space (dense), kh=0.0: 0.111 / 0.11
eigenvector coarse space (dense), kh=0.3: 0.131 / 0.13
eigenvector coarse space (dense), kh=0.6: 0.215 / 0.22
eigenvector coarse space (dense), kh=0.9: 0.493 / 0.49

## Findings

1. **Linear interpolation fails for Helmholtz:** rho = 0.70 at kh = 0.3, and the cycle diverges at kh >= 0.6.
   The hat stops being near the near-kernel. Its Rayleigh quotient z_p moves above 1, and 33 to 72 % of its weight sits on smooth modes that Jacobi *amplifies*.
2. **The rho-optimum is (1-omega)^2 = 1/9 for every kh, for every width, and even for a dense P** (checked at n = 15).
   The mode at theta = pi/2 has s = 1 - omega whatever kh is, so the two-grid is limited by the smoother. Any symmetric stencil leaves this mode untouched, because P^T v = 0 for it.
   The optimum is reached at width 3 by the symmetric stencil [a, 1, a] with **a = 1/(2 cos theta_0)**, where cos theta_0 = 1 - (kh)^2/2. The values are a = 0.50, 0.524, 0.610 and 0.840.
   This is phase-matched interpolation: it reproduces cos(theta_0 x) and sin(theta_0 x) exactly. The ρ = 1/9 plateau is a short interval of a whose lower end is exactly this a (for example [0.610, 0.648] at kh = 0.6).
3. **The right constraint is a two-dimensional local near-kernel that respects the boundary.**
   Constrained minimisation with B = {sin theta_0 (i+1), sin theta_0 (n-i)}, each imposed only away from the wall it violates, is exactly feasible at width 3. It fixes P completely and reproduces the optimum, rho = 1/9, for every kh.
   - One near-kernel vector is not enough: the coarse operator becomes singular at kh = 0.6, and rho = 0.55 at kh = 0.9.
   - Two eigenvectors give 0.19 to 0.24.
   - Unmasked plane waves solved by least squares give up to rho = 1.6, so the treatment at the boundary decides the outcome.
4. **Pseudomode quality of single columns is the wrong criterion here as well.** Columns that are individually optimal local pseudomodes give rho = 0.24 to 3.0.
5. **rho alone is misleading for indefinite A.**
   - The rho-optimal P's found by the search are strongly non-normal, with ||E||_2 between 2 and 10.
   - Even the clean phase-matched stencil has ||E||_2 = 0.40, 0.67 and 1.90 for kh = 0.3, 0.6 and 0.9.
   - The eigenvector coarse space ("spectral") is not a lower bound on rho for indefinite A: 0.13, 0.22 and 0.49 versus 1/9.
6. **The ||E||_2-optimum is the dense eigenvector space** (checked at n = 15 against a dense-P search). Sparse P reaches it at width 5 for kh <= 0.3, and at width 5 to 7 for kh = 0.6, but only with cond(P) of about 500 to 1400.
   At kh = 0.9, width 7 is still short (0.55 versus 0.49), with cond(P) up to 2*10^4.
   With width 3, the best norm is 0.23, 0.33, 0.48 and 0.91: wider support pays off in norm, unlike for Poisson.
7. **Locality floor (`scaling.txt`, n = 1023): eps_loc(z=1; H) ~ H^{-1} for Helmholtz, against H^{-2} for Poisson.**
   The near-kernel sits at theta_0, where the dispersion relation lambda(theta) = 4 sin^2(theta/2) - (kh)^2 has slope 2 sin theta_0 instead of 0. A support of length L = 2H+2 then smears the vector over a band of width |lambda'(theta_0)| pi / L, which gives
   eps_loc ~= (omega/d) * 2 sin(theta_0) * pi / (2H+2).
   This matches the computed values to 1 % for H >= 32. The crossover happens at H ~ 1/(kh).
   So local basis functions are intrinsically worse pseudomodes for Helmholtz than for Poisson, by a factor of about k H h.

## Files
- `hh1d.py`: operators, two-grid, constructions, exact constrained solves, and optimisers for rho and ||E||_2.
- `run_tables.py`: hat statistics, constructions, invariant-stencil scan, and the rho-optimum. Logs are in `log_n31_kh*.txt` and `res_*.json`.
- `run_norm.py`, `run_norm_restart.py`: the ||E||_2-optimum (`lognorm_*`).
- `run_dense_norm.py`: dense-P check at n = 15.
- `run_scaling.py`: locality floor (`scaling.txt`).
- `summarize.py`: the table above (`table.md`).
- `Popt_*.npy`, `Pnorm_*.npy`: the optimised P matrices.
