# Robust non-Markovian multilevel cycle with O(n) cost per cycle

Code: `robust_n.py` (one file, about 200 lines). Run `python3 robust_n.py 63 127 255 511`.

## 1. Objects

$Ax=b$ is the section 7.3 problem: $-\varepsilon\Delta u+\beta\cdot\nabla u$ on $(-1,1)^2$,
centered differences, $\beta=(2y(1-x^2),\,-2x(1-y^2))$, cell Péclet $|\beta|h/\varepsilon=20$.

On one level:

| symbol | meaning |
|---|---|
| $T=I-D^{-1}A$ | one undamped Jacobi step; $\rho(T)\approx4.5$, so relaxation diverges |
| $\mathcal C,\mathcal F$ | C/F split, greedy: $\sum_{j\in\mathcal F\setminus i}\lvert a_{ij}\rvert\le\delta\lvert a_{ii}\rvert$ for $i\in\mathcal F$, $\delta=0.3$ |
| $P_0=[I;0]$, $Q_0=[0;I]$ | local decomposition, so $P_0^\sharp=[I\;0]$ and $Q_0^\sharp=[0\;I]$ |
| $K=Q_0^\sharp TQ_0=I-D_{ff}^{-1}A_{ff}$ | CR propagator; $\lVert K\rVert_\infty\le\delta$ by the split |
| $H=Q_0^\sharp TP_0=-D_{ff}^{-1}A_{fc}$ | coarse-to-fine coupling |
| $S=\sum_{\ell<k}K^\ell H$ | memory sum (the history in Theorem 3.5) |
| $P'=P_0+Q_0S=[I;\,S]$ | effective prolongation (4.11) |
| $\widehat R=P'^T$ | restriction: Galerkin with $P'$ |
| $A_c=\widehat RAP'$ | coarse operator |

## 2. Algorithm as equations

**Setup on level $\ell$** ($k=4$ on level 0, $k=8$ below):

$$
S=\sum_{\ell'=0}^{k-1}K^{\ell'}H,\qquad
S\leftarrow\operatorname{drop}(S,\ 10^{-5}\cdot\text{row max})\ \ (\ell\ge1),
$$
$$
(AP')_C=A_{cc}+A_{cf}S,\qquad (AP')_F=A_{fc}+A_{ff}S,\qquad
A_c=(AP')_C+S^T(AP')_F .
$$

The next level is $A_{\ell+1}=\operatorname{filter}(A_c,\tau)$: off-diagonals below $\tau$ times the
row maximum are lumped into the diagonal, with $\tau=10^{-3}$ for $\ell+1<8$ and $3\cdot10^{-4}$ below.
Once $n_\ell<16\sqrt{n_0}$, that level is factored by sparse LU.

**Cycle on level $\ell$** (Algorithm 4.2 with the shift relation of Algorithm 7.3):

1. $x^{(j+1)}=x^{(j)}+D^{-1}(b-Ax^{(j)})$ for $j=0,\dots,k-1$
2. $d^{(j)}=(x^{(k)}-x^{(j)})_C$
3. $y\leftarrow Ky+Hd^{(j)}$ for $j=0,\dots,k-1$, starting from $y=0$; then $h=Q_0y$
4. $r=b-A(x^{(k)}+h)$
5. $r_c=\widehat Rr+A_cd^{(0)}$, with $A_c$ the exact (unfiltered) operator, applied matrix-free
6. $e=\mathrm{cycle}_{\ell+1}(r_c)$; on levels 0–1 one refinement step: $e\leftarrow e+\mathrm{cycle}_{\ell+1}(r_c-A_ce)$
7. $z=e-d^{(0)}$; return $x^{(k)}+h+P'z$

Why $z=e-d^{(0)}$: the coarse unknown is the coarse error at step 0,
$e_\sigma^{(0)}=e_\sigma^{(k)}+d^{(0)}$ (shift relation (3.6)). The coarse equation
$\widehat RAP'e_\sigma^{(k)}=\widehat R(r-\dots)$ then becomes $A_ce_\sigma^{(0)}=r_c$.

## 3. Matrix form

With an exact coarse solve, one cycle has the error propagator of Theorem 4.5:
$$
E_N=\bigl(I-P'(\widehat RAP')^{-1}\widehat RA\bigr)\,Q_0K^kQ_0^\sharp .
$$
The method has two factors, and each needs its own control:

- **The CR factor** $Q_0K^kQ_0^\sharp$ is small because the greedy split bounds $\lVert K\rVert$.
  Measured: $\rho(K)\approx0.23$ on every level and $\lVert K^8\rVert_2\approx5\cdot10^{-5}$, both independent of $m$.
- **The coarse-correction factor** $I-P'(\widehat RAP')^{-1}\widehat RA$ is an oblique projection.
  Its norm depends on $\widehat R$. With $\widehat R=P'^T$ it stays bounded as $m$ grows.

## 4. Code against the math

All lines are in `robust_n.py`:

| math | code |
|---|---|
| greedy $\mathcal F$, $\lVert I-D_{ff}^{-1}A_{ff}\rVert_\infty\le\delta$ | `greedy_f_points(A, delta)` |
| $K=I-D_{ff}^{-1}A_{ff}$, $H=-D_{ff}^{-1}A_{fc}$ | `K = (sp.eye(len(F)) - Dfi @ Aff)`, `H = (-Dfi @ Afc)` |
| $S=\sum_{\ell<k}K^\ell H$, truncated | the loop `V = K @ V; S = S + V`, then `drop_rel(S, tauP)` |
| $(AP')_C$, $(AP')_F$, $A_c=P'^TAP'$ | `APc`, `APf`, `Asig = APc + S.T @ APf` |
| $A_{\ell+1}=\operatorname{filter}(A_c,\tau)$ | `filter_lumped(Asig, tau)` |
| $A_cv$, matrix-free | `RAP(L, v) = APc @ v + S.T @ (APf @ v)` |
| steps 1–3 (relaxation and history) | the two loops in `cycle` |
| $r_c=\widehat Rr+A_cd^{(0)}$ | `rc = r[C] + ST @ r[F] + RAP(L, d0)` |
| refinement on levels 0–1 | `e = e + cycle(lev + 1, rc - RAP(L, e), ...)` |
| $x^{(k)}+h+P'z$, $P'z=[z;\,Sz]$ | `xh[C] += z`, `xh[F] += S @ z` |
| exact solve when $n_\ell<16\sqrt{n_0}$ | `sla.splu` |

## 5. What your Algorithm 7.3 computes instead

Your version takes injection, $\widehat R=[I\;0]$. Then
$$
A_c=(AP')_C=A_{cc}+A_{cf}S,\qquad r_c=r_C+(AP')_Cd^{(0)} .
$$
It also uses $k=8/12$ (finest/coarse), refinement on levels 0–3, no truncation of $S$, and
recursive cycles down to 150 unknowns.

The CR factor is the same in both versions. The difference is the projection
$I-P'((AP')_C)^{-1}[I\;0]A$, whose norm grows with $m$. Measured on level 1, the ratio of the
two-level factor to $\lVert K^8\rVert$ is 16, 32, 190 and 1000 for $m=63,127,255,511$, so the
per-level factors degrade (level 1: $4.8\cdot10^{-4}\to0.090$). Your algorithm compensates for
this with large $k$ and deep refinement, which is why it costs so much.

## 6. Check

Output of `robust_n.py` (the work count includes every sparse product):

| m | n | levels | factor | cycles to $10^{-8}$ | WU per cycle | op. complexity |
|---|---|---|---|---|---|---|
| 63 | 3,969 | 4 | 2.6e-6 | 2 | 173 | 4.23 |
| 127 | 16,129 | 6 | 2.2e-6 | 2 | 392 | 6.33 |
| 255 | 65,025 | 8 | 3.2e-3 | 3 | 563 | 7.76 |
| 511 | 261,121 | 11 | 3.2e-3 | 3 | 709 | 9.40 |

The development code (`ml_z73.py`, with `restrict="ppt"` and `coarsest=16√n`) also ran
m=1023 (n=1,046,529): factor 3.7e-3, 4 cycles.

For comparison, your Algorithm 7.3 needs 14,254 WU at m=511 (with the development code's work count).
