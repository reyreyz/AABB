# Relaxation-based AMG: independent checks

This is a dense NumPy reference implementation of the two-level schemes in
*A Theory of Relaxation-Based Algebraic Multigrid* (Moussa & Kahl). It is used
to check the paper's identities and to reproduce part of its tables.
Run `python3 run_checks.py`; it takes about 6 s and its output is in `results.md`.

## 1. Objects

| symbol | meaning |
|---|---|
| $T = I - MA$, $\widehat A = MA$ | error propagator of one relaxation step |
| $B=[P\;Q]$, $B^{-1}=[P^\sharp;Q^\sharp]$ | coarse/fine decomposition and its dual blocks |
| $K = Q^\sharp T Q$ | compatible-relaxation (CR) propagator: the orthogonal dynamics |
| $H = Q^\sharp T P$ | coupling from coarse to fine |
| $P^{(\ell)} = QK^{\ell-1}H$ | memory kernels (generalized prolongations) |
| $P' = \sum_{\ell=0}^{k} P^{(\ell)} = P + Q\sum_{\ell<k}K^\ell H$ | effective prolongation (4.11) |
| $\Pi_R = I - P(\widehat R A P)^{-1}\widehat R A$ | projection onto the canonical complement $\ker(\widehat R A)$ |

Setup used for the reproduction: 5-point Poisson on an $m\times m$ grid, coarse
points at even indices in both directions, bilinear $P=[I;W]$ in C/F order,
damped Jacobi $M=\tfrac23 D^{-1}$, Galerkin $\widehat R = P^T$, and the local
complement $Q_0=[0;I]$ with $Q_0^\sharp=[-W\;I]$.

## 2. Schemes as equations, and 3. their closed-form propagators

| scheme | algorithm | propagator | code |
|---|---|---|---|
| M | $\widehat RAP\,\epsilon = \widehat R r^{(k)}$, return $x^{(k)}+P\epsilon$ | $E^M = \Pi_R T^k$ (4.3) | `alg_markov`, `E_markov` |
| M′ | same with $P'$ in place of $P$ | $(I-P'(\widehat RAP')^{-1}\widehat RA)\,T^k$ | `E_markov_pp` |
| S | Alg. 4.1: exact coarse solve, shift relation $\epsilon^{(\ell)}=\epsilon^{(k)}+P^\sharp(x^{(k)}-x^{(\ell)})$, Horner replay $z\leftarrow P\epsilon^{(\ell)}+\mathcal Q T z$ | $E^S=(\Pi_R T\Pi_R)^k$ (4.8) | `alg_semi`, `E_semi` |
| N | Alg. 4.2: $y\leftarrow Ky+Hd^{(\ell)}$, $h=Qy$, $\widehat RAP'\epsilon=\widehat R r^{(k)}-\widehat RAh$ | $(I-P'(\widehat RAP')^{-1}\widehat RA)\,QK^kQ^\sharp$ (4.13) | `alg_nonmarkov`, `E_nonmarkov` |
| exact | Thm. 4.11: (4.15)–(4.18) with $(I-K^k)^{-1}$ | $0$ | `alg_exact` |

## 4. How the check works

Each algorithm is affine in $(x^{(0)},b)$. Set $b=0$, so the exact solution is
$x=0$ and $e^{(0)}=-x^{(0)}$. Feeding all unit vectors at once,
$x^{(0)}=-I$, the algorithm returns $x_{\text{out}}=-E_{\text{alg}}$.
The check then compares $E_{\text{alg}}$ with the closed form in the 2-norm.

## 5. Where this differs from the paper's runs

- Poisson is fully specified in the paper and is reproduced exactly.
- The paper does not give formulas for random diffusion, CD or the shifted
  Laplacian, so those problems (marked `*`) use our own choices, listed in
  `problems.py`. They enter only the identity checks, which hold for any $A$, $M$, $P$, $R$.
- Not reproduced: the $P_1$ interpolation, the Ruge–Stüben runs, the LS feedback
  complement $N_Z$, and the multilevel and PyAMG comparisons.

## 6. Results

- **Algorithms vs closed forms** (Table SM9 setup: $m=15$, $k=3$). The relative
  deviation is $10^{-15}$ to $4\cdot10^{-14}$ for M, M′, S and N. One cycle of
  the exact scheme leaves $\lVert E\rVert_2\le 6\cdot10^{-14}$.
- **(4.4) splitting** $E^M=\sum_\ell P^{(\ell)}P^\sharp T^{k-l}+QK^kQ^\sharp$ holds
  to $10^{-15}$.
- **$\rho(E^S)=\rho(K_{\text{can}})^k$** holds to $10^{-14}$. The exception is
  CD $\varepsilon=10^{-2}$, where it holds only to $10^{-11}$, because the
  eigenvalues of that non-normal $K$ are ill-conditioned.
- **Proposition 4.4** (Poisson): $\rho(E^S)\le\rho(E^M)$ holds, with
  $0.4190\le0.4231$ at $k=2$ and $0.1756\le0.1861$ at $k=4$.
- **Table SM1, Poisson, $m=31$**: every entry (M, M′, S, N, $\rho(K)^k$ for
  $k=1,2,4$) matches the paper to three digits.
