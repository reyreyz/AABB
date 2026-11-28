# Numerical examples

The three examples below each isolate one mechanism of the theory.

1. **The Markovian twin of N.** With injections for all operators and F-relaxation, the history terms of N vanish
   identically. N becomes a Markovian cycle M' whose interpolation is built by relaxation.
2. **Semi-Markovian S.** Memory enters only the interpolation of the coarse correction. Without the history shift,
   S reduces exactly to M', and M' diverges.
3. **True N.** The smoother moves the C-values, so the history is non-zero, and it enters both the coarse
   right-hand side and the interpolation.

Each example is compared with Markovian cycles that use the same smoother and hierarchy, and with PyAMG.

**Work and stopping rules.**
- One work unit (WU) is one product with the finest matrix. Every sparse product, including the coarsest LU solve,
  is charged with its number of nonzeros. Operator complexity (op.cx) is $\sum_\ell \mathrm{nnz}(A_\ell)/\mathrm{nnz}(A_0)$.
- Example 1 follows the tutorial protocol: $b=Ax_{\rm rand}$, $x^{(0)}=0$, stop at relative residual $10^{-10}$, at most
  300 cycles.
- Examples 2 and 3 use a random right-hand side, $x^{(0)}=0$, stop at $10^{-8}$, at most 60–80 cycles.
- "Diverges" means the residual grows; the number in brackets is the asymptotic growth factor per cycle.

---

## 1. Mathematical forms

**Notation.**
- Splitting $B=[P\;Q]$ with $B^{-1}=[P^\sharp;Q^\sharp]$. Relaxation $T=I-MA$.
- CR propagator $K=Q^\sharp TQ$ and coupling $H=Q^\sharp TP$; memory sum $S_k=\sum_{\ell<k}K^\ell H$.
- Effective interpolation $P'=P+QS_k$. Coarse correction $\Pi(X)=X(\widehat RAX)^{-1}\widehat RA$.
- In C/F ordering, injection is $P_0=[I;0]$, $Q_0=[0;I]$, $[I\;0]$, $[0\;I]$.

### 1.1 Example 1: N with injections and F-relaxation is its own Markovian twin

**Operators.**
- $P=[I;W_0]$, $Q=Q_0$, $P^\sharp=[I\;0]$, $Q^\sharp=[-W_0\;I]$, $\widehat R=[I\;0]$. The choice $W_0=0$ makes every
  operator an injection.
- F-relaxation $M=Q_0D_{ff}^{-1}Q_0^T$, so $T_F=I-Q_0D_{ff}^{-1}Q_0^TA$.

**The history vanishes.** F-relaxation never changes the C-values, so $P^\sharp T_F=P^\sharp$ and
$d^{(\ell)}=P^\sharp(x^{(k)}-x^{(\ell)})=0$. The history $y$, the shift $A_cd^{(0)}$ and the correction shift all vanish,
and Algorithm 4.2 reduces to the Markovian update
$$x^{\rm new}=x^{(k)}+P'\big(\widehat RAP'\big)^{-1}\widehat R\,r(x^{(k)}),\qquad
E_N=(I-\Pi(P'))\,T_F^k=(I-\Pi(P'))\,Q_0K^kQ^\sharp .$$

**The interpolation is k relaxation sweeps on $W$.** With $K=I-D_{ff}^{-1}A_{ff}$ and $H_0=-D_{ff}^{-1}A_{fc}$,
$$P'=\begin{bmatrix}I\\ W_k\end{bmatrix},\qquad W_{j+1}=KW_j+H_0\ \Big(=W_j-D_{ff}^{-1}(A_{ff}W_j+A_{fc})\Big),\qquad
W_k=K^kW_0+\sum_{\ell<k}K^\ell H_0 .$$
As $k\to\infty$, $W_k\to-A_{ff}^{-1}A_{fc}$ (ideal interpolation), and $\widehat RAP'=A_{cc}+A_{cf}W_k$ approaches the
Schur complement.

**Markovian reference M.** The same cycle with $W=W_0$ (direct interpolation): $E^M=(I-\Pi(P))\,T_F^k$.

### 1.2 Example 2: semi-Markovian S, memory in the interpolation only

**Operators.**
- Jacobi on all points, $T=I-\omega D^{-1}A$.
- $P=P_0$ (injection).
- $\widehat R=[I\;\;-A_{cf}(I+G+G^2)D_{ff}^{-1}]$ with $G=I-D_{ff}^{-1}A_{ff}$, giving $\Pi_R=I-\Pi(P_0)$.

**S** (Alg. 4.1 / 7.2): coarse solve at step $k$, shift $\varepsilon_\ell=\varepsilon_k+P^\sharp(x^{(k)}-x^{(\ell)})$, and
Horner replay $z\leftarrow P_0\varepsilon_\ell+\Pi_RTz$. Its error operator is
$$E^S=QK^kQ^\sharp=(\Pi_RT\Pi_R)^k=(I-\Pi_RMA)^k\,\Pi_R .$$
The last form shows that the Markovian equivalent of S needs $k+1$ coarse solves per cycle.

**Markovian references** (same $T$ and $\widehat R$):
$$E^M=\Pi_R\,T^k\ (\text{injection}),\qquad E^{M'}=(I-\Pi(P'))\,T^k\ (\text{static memory interpolation}).$$

**Identity:** S without the shift ($\varepsilon_\ell=\varepsilon_k$) has error operator
$\big(I-\sum_{j=0}^{k}(\Pi_RT)^jP_0(\widehat RAP_0)^{-1}\widehat RA\big)T^k=E^{M'}$. M carries the full $T^k$; S carries only
the CR operator $\Pi_RT\Pi_R$, with $\rho(\Pi_RT\Pi_R)=1-\omega$.

### 1.3 Example 3: true N, history in the coarse right-hand side and the interpolation

**Operators.**
- Jacobi with weights $\omega_F$ on F-points and $\omega_C$ on C-points, so $P^\sharp T\neq P^\sharp$ and $d^{(\ell)}\neq0$.
- $P_0$, $Q_0$ injections; truncated interpolation $\tilde P=[I;\tilde S]$; $\widehat R=\tilde P^T$; $\tilde\Pi=\Pi(\tilde P)$.

**Error operators** (two-level, exact coarse solve):
$$E_N=(I-\tilde\Pi)\big[\tilde PP^\sharp+Q(S_k-\tilde S)P^\sharp T^k+QK^kQ^\sharp\big],$$
$$E_{\rm twin}=(I-\tilde\Pi)\big[\tilde PP^\sharp+Q(S_k-\tilde S)P^\sharp+QK^kQ^\sharp\big],\qquad E_{\rm MJ}=(I-\tilde\Pi)T^k .$$
- The twin (k F-sweeps) puts the truncation defect on the raw C-error $c_0=P^\sharp e_0$.
- N puts it on $c_k=P^\sharp T^ke_0$, which the C-relaxation damps.
- MJ carries $T^k$.

**Identity (simplest Markovian rewrite).** With MJF = k sweeps of $T$, k F-sweeps, then correction:
$E_N-E_{\rm MJF}=(I-\tilde\Pi)QK^kQ^\sharp(I-T^k)$. This is zero when $K=0$.

---

## 2. Multilevel algorithms

### Algorithm A (Example 1): relaxation-built interpolation with injection restriction

Setup, on each level until fewer than 20 unknowns remain:
1. (optional) filter $A$: drop $a_{ij}$ with $|a_{ij}|<\tau|a_{ii}|$, add the dropped part to the diagonal;
2. classical strength (θ) and Ruge–Stüben C/F splitting (PyAMG);
3. starting interpolation $W_0$ (one-point or direct);
4. $d_P$ sweeps $W\leftarrow W-D_{ff}^{-1}(A_{ff}W+A_{fc})$;
5. injection restriction, $A_c=A_{cc}+A_{cf}W$.

Cycle:
1. $\nu$ F-relaxation sweeps;
2. $r_C=b_C-A_{cf}x_F-A_{cc}x_C$;
3. $e=$ cycle$(r_C)$;
4. $x_C\mathrel{+}=e$, $x_F\mathrel{+}=We$;
5. at the coarsest level, sparse LU.

With $d_P=\nu=k$ this is exactly Algorithm 4.2 with all-injection operators (§1.1). M is the same cycle with $d_P=0$.

### Algorithm B (Example 2): semi-Markovian V-cycle (parametrized Alg. 7.2)

Setup per level:
1. greedy F-points with $\|I-D_{ff}^{-1}A_{ff}\|_\infty\le0.3$, computed on a strength-filtered copy of $A$ (filter 0.3)
   on coarse levels;
2. $K=I-\omega D_{ff}^{-1}A_{ff}$ and $H=-\omega D_{ff}^{-1}A_{fc}$;
3. $R_f=-A_{cf}(I+G+G^2)D_{ff}^{-1}$, with entries below $10^{-2}\times$ the row maximum dropped;
4. $A_c=A_{cc}+R_fA_{fc}$, filtered (row-max $10^{-2}$, lumped).

Parameters: ω=0.9, k=2 on level 0 and k=4 below, LU at 150 unknowns.

Cycle:
1. $x^{(0)},\dots,x^{(k)}$ by Jacobi;
2. $d^{(\ell)}=x^{(k)}_C-x^{(\ell)}_C$, and $y\leftarrow Ky+Hd^{(\ell)}$ for $\ell<k$;
3. $r_c=r_C+R_fr_F+A_cd^{(0)}$;
4. $e=$ cycle$(r_c)$, plus two refinement cycles $e\mathrel{+}=$ cycle$(r_c-A_ce)$ on level 0;
5. $z=e-d^{(0)}$; $x_C\mathrel{+}=z$; $x_F\mathrel{+}=y+\sum_{\ell<k}K^\ell Hz$.

M is the same hierarchy with $x_C\mathrel{+}=e$ only. M' uses $P'$ statically.

### Algorithm C (Example 3): non-Markovian V-cycle with truncated interpolation

**Level 0:**
- red-black split (F = red), so $K=0$ and $S_k=H$;
- Jacobi with $\omega_F=1$, $\omega_C=0.2$, k=2;
- one-point $\tilde S$: the largest entry of each row of $S$, carrying the full row sum;
- $\widehat R=\tilde P^T$ and $A_c=\tilde P^TA\tilde P$.

**Coarse levels:** greedy split (δ=0.3), Jacobi ω=1, k=2, $S$ truncated at $10^{-2}$, $A_c$ filtered at $10^{-2}$; LU
below $16\sqrt n$.

**Cycle:**
1. $x^{(0)},\dots,x^{(k)}$;
2. $d^{(\ell)}=x^{(k)}_C-x^{(\ell)}_C$, and $y\leftarrow Ky+Hd^{(\ell)}$ for $\ell<k$; $h=Qy$;
3. $r_c=\widehat R\,r(x^{(k)}+h)+A_cd^{(0)}$;
4. $e=$ cycle$(r_c)$; $x\leftarrow x^{(k)}+h+\tilde P(e-d^{(0)})$.

The twin is the same hierarchy with k F-sweeps. MJ uses k Jacobi sweeps. MJF uses k Jacobi sweeps plus one F-sweep.

---

## 3. Numerical results

### 3.1 Example 1: the Markovian twin (upwind transport, ε = 0 / 0 / 10⁻⁶)

Problems:
- stagnation, $b=(x,-y)$;
- curved, $b=(1+y/2,\ \tfrac12\sin\pi x)$;
- jump, $b=(\mathrm{sign}\,y,\ 0.3)$.

All are first-order upwind on $(-1,1)^2$, with $m^2$ unknowns.

Configurations:
- **Twin:** θ=0.25, one-point $W_0$, $d_P=\nu=3$, τ=$10^{-4}$ for all three problems. With $d_P=\nu$ this is exactly N (§1.1).
- **M:** the best $d_P=0$ configuration (direct $W_0$, ν=1).
- **PyAMG:** the best of 41 configurations at m=180 (AIR for stagnation and jump, RS for curved), used at every size.

| problem | m | n | twin (= N): WU / cycles / op.cx | M (direct): WU / cycles | best PyAMG: WU / cycles |
|---|---|---|---|---|---|
| stagnation | 90 | 8,100 | **18.4** / 4 / 3.48 | 51.6 / 24 | 125.8 / 11 |
| | 180 | 32,400 | **19.3** / 4 / 3.72 | 91.0 / 40 | 149.6 / 12 |
| | 360 | 129,600 | **27.0** / 5 / 4.07 | no conv. (300) | 185.2 / 14 |
| | 720 | 518,400 | **28.0** / 5 / 4.22 | diverges | 207.5 / 15 |
| curved | 90 | 8,100 | **19.9** / 4 / 3.75 | 57.3 / 27 | 112.1 / 10 |
| | 180 | 32,400 | **21.1** / 4 / 4.09 | 89.9 / 41 | 159.1 / 14 |
| | 360 | 129,600 | **28.2** / 5 / 4.22 | 169.0 / 75 | 243.8 / 21 |
| | 720 | 518,400 | **29.7** / 5 / 4.42 | diverges | 407.5 / 35 |
| jump | 90 | 8,100 | **13.0** / 4 / 2.44 | 35.6 / 20 | 85.5 / 8 |
| | 180 | 32,400 | **13.0** / 4 / 2.60 | 55.7 / 30 | 90.3 / 8 |
| | 360 | 129,600 | **16.2** / 4 / 2.98 | 83.4 / 43 | 105.5 / 9 |
| | 720 | 518,400 | **17.2** / 4 / 3.23 | 142.0 / 72 | 105.9 / 9 |

- The twin converges in 4–5 cycles at every size, and its total work grows slowly (13.0–29.7 WU from 8k to 518k
  unknowns; operator complexity 2.4–4.4).
- M with direct interpolation needs 20–75 cycles and 2.7–8.3× more work where it converges. It fails on the two
  largest stagnation grids and the largest curved grid. One-point
  M fails on curved and jump, and needs 133.7 WU on stagnation at m=180.
- Against the best PyAMG configuration the twin is 5.6–13.7× cheaper. The gap widens with m on curved (5.6× → 13.7×)
  and stays around 6–8× on stagnation and jump.

### 3.2 Example 2: semi-Markovian S (centered convection, cell Pe 20)

Problems:
- shear flow, $b=(-y,0.3)$;
- constant flow, $b=(\cos0.4,\sin0.4)$.

Both are centered differences with $\varepsilon$ chosen so that the cell Péclet number is 20. Jacobi with ω=0.9 has
$\rho(T)\approx5$–$6$. Algorithm B.

| problem | m | n | S: factor / cycles | WU/cycle | WU total | op.cx | M (same hierarchy) |
|---|---|---|---|---|---|---|---|
| constant flow | 127 | 16,129 | 0.010 / 4 | 180 | 722 | 7.95 | |
| | 255 | 65,025 | 0.011 / 5 | 190 | 949 | 8.29 | diverges (316) |
| | 511 | 261,121 | 0.015 / 5 | 195 | 976 | 8.47 | |
| | 1023 | 1,046,529 | 0.036 / 6 | 198 | 1,188 | 8.57 | |
| shear flow | 127 | 16,129 | 0.010 / 5 | 147 | 736 | 6.44 | |
| | 255 | 65,025 | 0.011 / 5 | 154 | 771 | 6.71 | diverges (156) |
| | 511 | 261,121 | 0.017 / 5 | 160 | 799 | 6.87 | |
| | 1023 | 1,046,529 | 0.053 / 7 | 161 | 1,126 | 6.92 | |

**Markovian references.**
- **M' = S without shift.** Two-level (m=31): it diverges for every ω ∈ [0.02, 0.9] and k ∈ {2, 4, 16, 64}, e.g.
  161 (shear) and 295 (constant flow) at ω=0.7, k=4. S converges at $(1-\omega)^k$.
- **M' multilevel** (on the ω=0.7, k=4 hierarchy of Alg. 7.2): diverges (241–251 shear, 921–953 constant flow,
  m=127–511).
- **M with ideal interpolation (two-level):** diverges (165, 301).
- **PyAMG AIR and RS:** fail on both problems, stand-alone and as GMRES preconditioners (m=127, 255).

S converges in 4–7 cycles with bounded work per cycle (147–198 WU) and operator complexity (6.4–8.6) up to $10^6$
unknowns.

### 3.3 Example 3: true N (constant flow, cell Pe 6)

Problem: centered convection–diffusion, $b=(\cos0.4,\sin0.4)$, cell Pe 6. Jacobi with ω=1 has $\rho(T)=1.6$.
Algorithm C.

| m | n | levels | N: factor / cycles | WU/cycle | WU total | op.cx | twin | MJ | MJF: WU total |
|---|---|---|---|---|---|---|---|---|---|
| 127 | 16,129 | 6 | 0.365 / 19 | 44 | 840 | 4.16 | diverges (1.36) | diverges (2.06) | 666 |
| 255 | 65,025 | 7 | 0.364 / 19 | 49 | 925 | 4.48 | diverges (1.40) | diverges (2.07) | 709 |
| 511 | 261,121 | 9 | 0.363 / 19 | 53 | 1,013 | 4.82 | diverges (1.39) | diverges (2.07) | 758 |
| 1023 | 1,046,529 | 10 | 0.362 / 19 | 55 | 1,053 | 4.93 | diverges (1.38) | diverges (2.07) | 783 |

Two-level spectral radii ($\rho(E_N)$ / twin / MJ): 0.403 / 1.21 / 1.76 (m=63), 0.407 / 1.35 / 2.03 (m=127),
0.406 / 1.39 / 1.81 (m=255).

- N's factor is mesh-independent (0.36) and its cost O(n): 19 cycles, 44–55 WU per cycle, op.cx at most 4.9 up to
  $10^6$ unknowns.
- On the same hierarchy, the twin diverges (its truncation defect acts on the undamped $c_0$) and MJ diverges (it
  carries $T^k$).

**Remarks.**
- Because $K=0$ on level 0, MJF reproduces N exactly (same factor and cycles) and costs 20–26% less (identity of
  §1.3). The history can be replaced by one extra F-sweep.
- PyAMG AIR solves this problem in 14 cycles (132 WU at m=127/255); RS fails.

Example 3 therefore shows that N converges in O(n) where its own twin and the same-smoother Markovian cycle diverge.
It does not show N beating all Markovian methods.

### 3.4 Summary

| example | memory used at solve time | Markovian cycles on the same smoother and hierarchy | equivalent Markovian cycle | PyAMG |
|---|---|---|---|---|
| 1. twin (N with injections) | none (history vanishes) | M (direct): 2.7–8.3× more work, fails at large m | is itself Markovian (M') | 5.6–13.7× more work |
| 2. S | shift in the interpolation | M, M' (= S without shift), M(ideal) diverge | needs k+1 coarse solves per cycle | AIR, RS fail |
| 3. true N | history in coarse rhs and interpolation | twin, MJ diverge | MJF (one F-sweep), 20–26% cheaper | AIR converges (132 WU) |
