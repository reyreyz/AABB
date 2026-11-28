# N and its Markovian twin: claims and evidence

Script: `claims.py` (sections A–D). Outputs: `results_claims_{A,B,C,D}.md`.

**Notation.** $B=[P\;Q]$ with duals $P^\sharp,Q^\sharp$; relaxation $T=I-MA$; $K=Q^\sharp TQ$,
$H=Q^\sharp TP$; $P'=P+Q\sum_{\ell<k}K^\ell H$; $\Pi'=P'(\widehat RAP')^{-1}\widehat RA$.

- **N** is Algorithm 4.2 on the iterates.
- **The twin** takes $k$ steps $x\leftarrow x+Q\,Q^\sharp M(b-Ax)$, a relaxation of the form
  $I-\widetilde MA$ with $\widetilde M=QQ^\sharp M$ that uses only the residual, followed by the
  Petrov–Galerkin correction with $P'$.

## Claim 1. N and the twin are the same method, for any complement and any relaxation

$P^\sharp T_Q=P^\sharp$ and $Q^\sharp T_Q=Q^\sharp T$ give $T_Q^k=P'P^\sharp+QK^kQ^\sharp$. With
$(I-\Pi')P'=0$ this yields $E_{\rm twin}=(I-\Pi')QK^kQ^\sharp=E_N$.

**Evidence (A).** The propagators and the iterates after one cycle agree to $10^{-18}$–$10^{-13}$ in
all 33 cases:
- **4 problems:** CD $\varepsilon=10^{-2},10^{-3}$; shifted Laplacian ($\rho(T)>1$); recirculating
  flow at Pe 20 ($\rho(T)$ up to 20).
- **3 complements:** injective $Q_0$, LS feedback $Q_0+PZ$, random sparse $Q_0+PZ$.
- **3 relaxations:** Jacobi with $\omega=2/3$, Jacobi with $\omega=1$, Gauss–Seidel.

## Claim 2. With Jacobi and $Q=Q_0+PZ$, the twin is cheaper

Here $Q^\sharp=Q_0^\sharp$, so the twin needs only the F-rows of the residual. N needs full
sweeps, the history ($K$, $H$, $P^\sharp$), the extra product $Ah$, and $h=Qy$.

**Evidence (B).** Sparse two-level comparison on the recirculating flow, Pe 20 (greedy split, $P_0$,
injection, $k=4$, exact coarse LU, which both methods share):

| m | complement | iterates N vs twin | N: WU per cycle | twin: WU per cycle | N: ms | twin: ms |
|---|---|---|---|---|---|---|
| 63 | $Q_0$ | 1.5e-15 | 15.4 | 10.9 | 1.7 | 0.3 |
| 63 | $Q_0+P_0Z$ | 1.1e-15 | 16.6 | 12.0 | 1.1 | 0.5 |
| 127 | $Q_0$ | 1.1e-14 | 19.9 | 15.0 | 5.2 | 2.0 |
| 127 | $Q_0+P_0Z$ | 4.8e-15 | 107.3 | 98.9 | 21.1 | 13.1 |
| 255 | $Q_0$ | 6.7e-15 | 24.8 | 19.9 | 21.2 | 11.9 |
| 255 | $Q_0+P_0Z$ | 2.8e-15 | 172.9 | 164.5 | 160.8 | 79.0 |

The WU includes the shared coarse LU, which dominates for $Q_0+P_0Z$ at larger m. The difference,
N minus twin, is 4.5–8.4 WU per cycle at every size. In wall time the twin is 1.6–5.7× faster.

## Claim 3. With Gauss–Seidel the gap widens

$T=I-(D-L)^{-1}A$ makes $K=Q^\sharp TQ$ dense, so N must apply $K$ and $H$ matrix-free.

**Evidence (C).** On CD with $Q=Q_0+PZ$:

| | m=15 | m=31 |
|---|---|---|
| fill of $K$ with Jacobi ($\mathrm{nnz}/n_f^2$) | 0.155 | 0.041 |
| fill of $K$ with Gauss–Seidel ($\mathrm{nnz}/n_f^2$) | 0.547 | 0.39 |

The iterates agree to $2\cdot10^{-15}$. Per cycle, N needs 9 Gauss–Seidel applications ($k$ sweeps
plus 2 per history step); the twin needs 3.

## Claim 4. With a black-box smoother, only N converges

The smoother is an opaque iterate generator: given $(b,x_0)$ it returns $x^{(0)},\dots,x^{(k)}$ for
undamped Jacobi ($\rho(T)\approx4.5$) and cannot be modified or restarted mid-run.
- N post-processes the iterates.
- The non-intrusive Markovian cycle uses only $x^{(k)}$ before the coarse correction.
- The twin needs $x\leftarrow x+QQ^\sharp(S(x)-x)$ fed back after every step, i.e. $k$ restarts per
  cycle.

**Evidence (D).** Robust-N hierarchy (`robust_n.py`) on the recirculating flow at Pe 20:

| m | n | N on the black-box iterates | non-intrusive Markovian | twin (needs restarts) |
|---|---|---|---|---|
| 127 | 16,129 | 2 cycles, factor 2.2e-6, 386 WU/cycle | diverges (5e9) | 2 cycles, 87 WU/cycle |
| 255 | 65,025 | 3 cycles, factor 3.2e-3, 555 WU/cycle | diverges (1e12) | 3 cycles, 105 WU/cycle |
| 511 | 261,121 | 3 cycles, factor 3.2e-3, 702 WU/cycle | diverges (2.5e12) | 3 cycles, 113 WU/cycle |

Under the black-box constraint N is the only one of these methods that converges, with a
mesh-independent factor and cost per cycle close to O(n). If restarts are allowed, the twin does
the same work about 6× more cheaply.

## Claim 5. A dense $Q^\sharp$ cannot make N cheaper than its twin (`claims.py E`)

The twin uses $Q^\sharp$ only through the product $QQ^\sharp = I - PP^\sharp$. This is the oblique projector
onto $\mathrm{col}(Q)$ along $\mathrm{col}(P)$, and it doesn't depend on the basis chosen for $\mathrm{col}(Q)$.
N needs $P^\sharp$ anyway, to form the history. A basis change $Q\to QS$ (Prop. 3.6: $Q^\sharp\to S^{-1}Q^\sharp$,
$P^\sharp$ unchanged) can make $Q^\sharp$ dense, but N then works with $K_S=S^{-1}KS$ and $H_S=S^{-1}H$, which are
dense too. Example: $S=I+0.9\,\mathrm{subdiag}$, which is bidiagonal with a dense inverse.

| problem | m | fill $Q^\sharp$: $Q_0$ / $Q_0S$ | fill $K$ | fill $H$ | $\|QQ^\sharp$ change$\|$ | N iterates change | twin iterates change |
|---|---|---|---|---|---|---|---|
| CD 1e-2 | 15 | 0.016 / 0.522 | 0.057 / 0.614 | 0.092 / 0.672 | 9.6e-16 | 2.5e-16 | 3.4e-16 |
| CD 1e-2 | 31 | 0.004 / 0.304 | 0.014 / 0.375 | 0.025 / 0.393 | 1.9e-15 | 2.9e-16 | 3.7e-16 |
| CD 1e-3 | 15 | 0.016 / 0.522 | 0.057 / 0.614 | 0.092 / 0.672 | 9.6e-16 | 2.6e-16 | 3.5e-16 |
| CD 1e-3 | 31 | 0.004 / 0.304 | 0.014 / 0.374 | 0.025 / 0.392 | 1.9e-15 | 4.7e-16 | 5.4e-16 |

So the "dense $Q^\sharp$, sparse $K,H$" case is empty. If $K$ and $H$ are sparse in some basis, the twin uses the
same projector, and that projector does not depend on the basis.

## Claim 6. Latency: N has fewer critical-path rounds for a non-injective $Q$ and $k\ge4$, at higher flop cost (`claims.py F`)

The history can be pipelined: $y = P'_F\,P^\sharp x^{(k)} - u$, where $u\leftarrow Ku + HP^\sharp x^{(\ell)}$ is computed
alongside the sweeps. That takes the history off the critical path. The twin's projected update has to finish before
the next residual, so it stays on the critical path. The test problem is the recirculating flow at Pe 20 with $P_0$ and
injection; rounds are sparse matvecs on the critical path, counting one halo exchange each, with the coarse solve excluded.

| m | complement | k | pipelined N vs N | N rounds | twin rounds | N flops/nnz(A) | twin flops/nnz(A) |
|---|---|---|---|---|---|---|---|
| 127 | $Q_0$ | 4 | 2.6e-14 | 7 | 6 | 7.90 | 3.00 |
| 127 | $Q_0$ | 8 | 7.8e-13 | 11 | 10 | 13.90 | 5.00 |
| 127 | $Q_0+P_0Z$ | 2 | 2.3e-15 | 7 | 6 | 12.77 | 6.35 |
| 127 | $Q_0+P_0Z$ | 4 | 1.8e-14 | **9** | 10 | 32.90 | 18.24 |
| 127 | $Q_0+P_0Z$ | 8 | 9.3e-13 | **13** | 18 | 98.21 | 58.69 |
| 255 | $Q_0+P_0Z$ | 8 | 1.5e-11 | **13** | 18 | 101.81 | 61.17 |

The rounds are $k+5$ for N and $2k+2$ for the twin when $Q$ is non-injective, and $k+3$ against $k+2$ for $Q_0$. Full
table: `results_claims_F.md`. **Caveat:** these are counts from a latency model, not timings of a parallel run.

## Claim 7. In the multilevel, sparsified setting N and its twin differ only through truncation of $P'$ (`ml_twins.py`)

Setup: both cycles share the hierarchy $\tilde P=[I;\tilde S]$, $R=\tilde P^T$, $A_c=R A\tilde P$ (applied exactly), and a
filtered next-level operator. The coarse cycle $B_c$ is any linear map: filtering, recursion, refinement or LU.
N's history uses the exact $K,H$. With $c_\ell=P^\sharp e^{(\ell)}$ and $d^0=c_0-c_k=P^\sharp(I-T^k)e_0$:

$$E_N e_0-E_{\rm twin}e_0=-(I-\tilde P B_c R A)\,Q(S-\tilde S)\,d^0,$$

which is equivalent to

$$E_{\rm twin}e_0=(I-\tilde PB_cRA)\big[\tilde Pc_0+Q(S-\tilde S)\,c_0+QK^kq_0\big],\qquad
E_N e_0=(I-\tilde PB_cRA)\big[\tilde Pc_0+Q(S-\tilde S)\,c_k+QK^kq_0\big].$$

So an inexact coarse solve, filtering and recursion change both cycles in exactly the same way. The only difference is
which C-error the truncation defect $S-\tilde S$ acts on: the raw $c_0$ for the twin, and $c_k=P^\sharp T^ke_0$, the
C-error after $k$ global sweeps, for N.

Check (`check_twins.py`, m=63): with no truncation, N and the twin agree after one cycle to 1e-13 (Poisson) and 1e-7
(recirc Pe 20, where $\rho(T)^k$ amplifies roundoff), with filtering, recursion and refinement all active. With
truncation 1e-2 on two levels, the formula matches the measured difference to 1.2e-11.

What follows from this:
- If $T$ amplifies the C-error (recirc Pe 20, $\rho(T)=4.5$), N is **less** robust. With truncation 1e-2, N's factor
  is 0.42 against the twin's 0.010; at 1e-1, N diverges while the twin's factor is 0.90.
- If $T$ smooths the C-error and the truncation keeps row sums (smooth $c$ is then reproduced exactly), N is **more**
  robust. Recirc Pe 2, Jacobi $\omega=0.8$, $k$=4/8, truncation 0.5 lump ($\tilde P$ nnz $\approx0.5\,$nnz(A)):

| m | N: factor / cycles / WU | twin | Markov-Jacobi (same $\tilde P$) |
|---|---|---|---|
| 127 | 0.159 / 11 / 1307 | 0.45 / 23 / 810 | 0.30 / 16 / 1383 |
| 255 | 0.167 / 11 / 1625 | 0.66 / no conv. | 0.32 / 16 / 1673 |
| 511 | 0.175 / 11 / 1719 | 0.71 / no conv. | 0.32 / 17 / 1874 |

- **But with each method at its best truncation (m=511, Pe 2):** the twin at 0.1 lump reaches 0.114 / 9 cycles / **639 WU**,
  while N's best is 1719 WU and Markov-Jacobi's is 1756 WU. N's advantage is robustness over a wider range of
  truncation (its factor stays at or below 0.18 for $\tau\le0.5$), not cost. It wins only if a memory or communication
  budget forces $\tilde P$ below about 0.6 nnz(A). Tables: `results_twins_trunc_127.md`, `results_twins_mesh*.md`.

### Claim 7b. Shift variant: coarse unknown = step-$j$ error (`ml_twins.py`, `shift=j/k`)

N with coarse right-hand side $R\,r(x^{(k)}+h)+A_c d^j$, $d^j=c_j-c_k$, and correction $z=e-d^j$:

$$E_je_0=\tilde P(I-B_cA_c)\,c_j+(I-\tilde PB_cRA)\,Q\big[(S-\tilde S)\,c_k+K^kq_0\big].$$

The shift moves only the **coarse-cycle defect** $(I-B_cA_c)$, onto $c_j$. The truncation defect stays on $c_k$ for
every $j$. So the shift cannot interpolate between N and the twin in the truncation term. The twin has both defects
on $c_0$; $j=k$ puts both on $c_k$. (`check_shift.py`: the formula holds to 1e-13 for $j=0,k/2,k$ and $\tau=0,0.5$.)

Results (factor / cycles / WU, `results_shift.md`):
- Recirc Pe 20 ($T$ amplifies): with no truncation, $j=0$ and $j=k/2$ match the twin (0.0032 at m=255), while
  $j=k$ fails (0.0032 → 4 cycles at m=127, diverges at m=255). The coarse defect acting on the amplified $c_k$ hurts.
- Recirc Pe 2 and 5 ($T$ smooths): $j=0$, $k/2$ and $k$ agree to 3 digits in every row (m up to 511). The coarse
  cycle is accurate enough that its defect is negligible; the gap to the twin comes entirely from the truncation term.
- So $j=0$ (Alg 7.3) is the right choice: never worse, and much better when $T$ amplifies.

### Claim 7c. Diverging smoother: putting the coarse solve on the step-0 error is what makes N converge (`diverge_levels.py`)

With no truncation, $E_je_0=\tilde P(I-B_cA_c)\,c_j+(I-\tilde PB_cRA)QK^kq_0$, and $\|c_k\|\approx g\,\|c_0\|$, where $g$
is the growth of the C-error over $k$ sweeps. A coarse solve on the step-$k$ error multiplies the coarse defect by
$g$; a coarse solve on the step-0 error ($j=0$, Alg 7.3) never sees the divergence.
Recirc Pe 20, Jacobi $\omega=1$, m=255:

| level | k | $\rho(T_\ell)$ | growth $g_\ell$ of the C-error |
|---|---|---|---|
| 0 | 4 | 4.86 | 353 |
| 1 | 8 | 2.29 | 723 |
| 2 | 8 | 1.49 | 23.6 |
| 3 | 8 | 1.33 | 15.1 |
| 4–6 | 8 | ≈1.0–1.1 | 0.7–1.3 |

| coarse unknown = step-$k$ error on | factor | cycles |
|---|---|---|
| no level (Alg 7.3) | 0.00317 | 3 |
| level 0 only | 0.00314 | 3 |
| levels ≥ 1 only | 0.22 | 12 |
| all levels | diverges | – |

Level 0 tolerates $g=353$ because its coarse defect is about 5e-6 (refinement plus accurate nested cycles). On the
deeper levels the nested cycle is less accurate, so $g\times$defect > 1 there and the error compounds.
How the three methods treat the divergence:
- Markov-Jacobi has no history and no shift, so both the F- and C-divergence stay in the error (one cycle multiplies
  it by 5e7).
- N with $j=k$ removes the F-divergence through the history ($QK^kq_0$) but passes the C-divergence to the inexact
  coarse solve.
- N with $j=0$, like the twin, removes both.
The truncation defect of Claim 7 also acts on the amplified $c_k$, which is why N fails under truncation at Pe 20.

## Claim 8. Convergence: N with the same smoother is never less robust than the Markovian cycle (`conv_compare.py`)

Two-grid, same $T$, same $P'$, $R=P'^T$: $E_N=(I-\Pi)Q_0K^kQ_0^T$ and $E_M=(I-\Pi)T^k$. For injection,
$K=T_{ff}$ is a principal submatrix of $T$, so in the 1-, 2- and ∞-norms
$$\|E_N\|\le\|I-\Pi\|\,\|K\|^k\le\|I-\Pi\|\,\|T\|^k ,$$
while the best norm bound available for M is $\|I-\Pi\|\,\|T\|^k$. With Jacobi and the greedy split,
$\|K\|_\infty\le\delta$ for any $A$, so N's bound decays like $\delta^k$ whatever $\rho(T)$ is. If $\rho(T)>1$, then
$E_M$ grows like $\rho(T)^k$ unless the unstable modes lie in the kernel of $I-\Pi$.

Sweep (m=15): 432 cases, from 12 problems (Poisson, anisotropic, recirc centered Pe 2/20/100, constant-flow upwind
Pe 2/20/100) × 3 smoothers (Jacobi 2/3, Jacobi 1, Gauss–Seidel) × 6 splits (greedy δ=0.3/0.6/0.9, red-black, full,
aggressive) × $k=1,2,4$:

| outcome | cases |
|---|---|
| both converge | 324 |
| only N converges | 48 (all recirc Pe 20/100 with Jacobi) |
| only M converges | **0** |
| neither | 60 (Gauss–Seidel on centered Pe 20/100; Jacobi with non-greedy splits) |

When both converge, $\rho(E_N)\le\rho(E_M)$ in 242 of 324 cases, with median ratio 0.83. The largest advantage for M
is a factor of 1.20 (Pe 2, Jacobi 2/3, $k=1$). N's largest advantage is several orders of magnitude (Poisson,
Jacobi 1, red-black: 0 vs 0.96). The multilevel runs (Claims 7 and 7c) also contain no case where Markov-Jacobi
converges and N does not, provided the coarse unknown is the step-0 error.
Caveats:
- This compares N with Markovian cycles on the **same** smoother. The twin, which uses a modified smoother, converges
  exactly when N does (Claim 1).
- N's guarantee needs a split that controls $K$ for the smoother actually used. The greedy split does this for
  Jacobi, not for Gauss–Seidel (the 36 "neither" Gauss–Seidel cases).

## Claim 9. With a truncated $P'$, the twin can diverge while N converges (`twin_limit_scan.py`, `twin_div_k.py`, `twin_div_m63.py`)

From Claim 7, with a fixed sparse $\tilde S$ and $k\to\infty$:
$$E_{\rm twin}\to(I-\tilde\Pi)\,P'_\infty P^\sharp=(I-\tilde\Pi)\big[\tilde P+Q(S_\infty-\tilde S)\big]P^\sharp,\qquad
E_N=(I-\tilde\Pi)\big[P'P^\sharp T^k+QK^kQ^\sharp\big]\to0\ \text{ if }\rho(T)<1 .$$
- The twin's limit does not depend on the smoother, because $S_\infty=-A_{ff}^{-1}A_{fc}$ is the ideal prolongation.
- For SPD $A$, $\|E_{\rm twin,\infty}\|_A\le1$, so the twin cannot diverge. It needs a nonsymmetric $A$.

Example: centered recirc, Pe 10. Damped Jacobi with $\omega=0.1$ gives $\rho(T)=0.9995<1$. $\tilde P$ keeps one entry per
row carrying the full row sum ($\tau=1$ lump, 0.1 nnz(A)). Exact coarse solve.

| m | δ | k | $\rho(E_N)$ / cycles | $\rho(E_{\rm twin})$ | Markov-Jacobi | twin limit |
|---|---|---|---|---|---|---|
| 31 (eig) | 0.3 | 16 | 0.825 | 1.32 | 0.822 | 1.40 |
| 31 (eig) | 0.3 | 64 | 0.136 | 1.40 | 0.157 | 1.40 |
| 31 (eig) | 0.6 | 64 | 0.258 | 1.87 | 0.266 | 1.87 |
| 63 (runs) | 0.3 | 64 | 0.20 / 15 cycles | 1.86, diverges | 0.165 / 14 | |
| 63 (runs) | 0.6 | 64 | 0.864, no conv. in 60 | 5.43, diverges | 0.282 / 18 | |

Caveats:
- The plain Markovian cycle on the **same** smoother (k Jacobi sweeps, then $\tilde P$ correction) converges as well as
  N here. The twin fails because F-relaxation never touches the C-error, so the truncation defect hits the raw $c_0$.
  This example breaks the twin equivalence in N's favour; it does not show N beating every Markovian scheme.
- It needs $k=64$ heavily damped sweeps. At δ=0.6 and m=63, N has not converged within 60 cycles.

## Claim 10. Twin and Markov-Jacobi both diverge, N converges (two-level) (`search_missing.py`, `missing3*.py`, `missing4.py`, `missing5.py`)

With $P_0,Q_0$ and $(I-\tilde\Pi)\tilde P=0$:
$$E_N=(I-\tilde\Pi)Q\big[(S-\tilde S)P^\sharp T^k+K^kQ^\sharp\big],\quad
E_{\rm twin}=(I-\tilde\Pi)Q\big[(S-\tilde S)P^\sharp+K^kQ^\sharp\big],\quad
E_M=(I-\tilde\Pi)Q\,(Q^\sharp-\tilde SP^\sharp)\,T^k .$$
What each method needs:
- Markov-Jacobi diverges when $T^k$ grows away from range($P'$).
- The twin diverges when the truncation defect on the raw $c_0$ is large.
- N converges when $\rho(K)<1$ and $\tilde S$ is accurate on $c_k=P^\sharp T^ke_0$.

**Example.** Centered recirc, cell Pe 10. Jacobi with $\omega_F=1$ on F-points and $\omega_C=0.8$ on C-points (all
three methods use this same $T$). Greedy δ=0.3, so $\rho(K)\approx0$. $k=1$. $\tilde S$ keeps 2 entries per row,
minimally adjusted to match $S$ on $1,x,y$ (0.2 nnz(A)). $R=\tilde P^T$, exact coarse solve.

| m | n | $\rho(E_N)$ | $\rho(E_{\rm twin})$ | $\rho(E_M)$ | N from a random rhs |
|---|---|---|---|---|---|
| 127 | 16,129 | 0.564 | 1.01 | 5.94 | 40 cycles to 1e-10 (twin, MJ diverge) |
| 255 | 65,025 | 0.609 | 1.13 | 6.03 | 48 cycles (twin, MJ diverge) |
| 511 | 261,121 | 0.649 | 1.21 | 6.22 | |

N's increments per refinement shrink (0.10, 0.06, 0.05, 0.04 from m=31) while the twin's radius keeps rising. At
m≤63 the twin still converges (0.69, 0.82), so this is a hit only from m=127 on.
Limits:
- The best Markovian cycle over $\omega$, $\omega_C$ and $k$ on the same $\tilde P$ converges at 0.92 (m=127:
  $\omega=0.4$, $k=2$, F-only), against 0.58 for N, i.e. about 6.5× more cycles per digit. So Markovian cycles on
  this $\tilde P$ do not all diverge; they converge badly.
- Multilevel version not achieved: with $k_c=1$, N stagnates (0.96–0.98); with $k_c=2$ it diverges.
- Earlier attempts that failed: single-entry $\tilde P$ (N drifts to 0.92 at m=255), and the shifted Laplacian with
  injection $R$ (a hit only at m=15).

## Claim 11. Multilevel version of Claim 10 (`ml_missing_final.py`)

**Hierarchy, shared by all three methods** ($R=\tilde P^T$ on every level):
- Level 0: the Claim 10 setup (Jacobi $\omega_F=1$, $\omega_C=0.8$, $k=1$, $\tilde S$ with 2 entries per row exact on
  $1,x,y$).
- Levels ≥ 1: robust N (Jacobi $\omega=1$, $k=8$, $S$ truncated at 1e-5, filter 1e-3, or 3e-4 from level 8, sparse LU
  below $16\sqrt n$, V-cycle).

| m | n | levels | N: factor / cycles | WU/cycle | WU total | op.cx | twin | Markov-Jacobi |
|---|---|---|---|---|---|---|---|---|
| 127 | 16,129 | 5 | 0.551 / 31 | 117 | 3,627 | 4.38 | 1.02, diverges | 6e4, diverges |
| 255 | 65,025 | 7 | 0.614 / 37 | 167 | 6,192 | 5.52 | 1.15, diverges | 3e6, diverges |
| 511 | 261,121 | 9 | 0.645 / 41 | 199 | 8,153 | 6.33 | 1.24, diverges | 3e7, diverges |
| 1023 | 1,046,529 | 12 | 0.696 / 50 | 239 | 11,936 | 7.40 | 1.26, diverges | 2e8, diverges |

**Cost per cycle.** nnz(A_ℓ)/nnz(A) per level, at m=1023: 1.00, 0.80, 0.95, 0.85, 0.74, 0.62, 0.51, 0.42, 0.46, 0.40,
0.35, 0.30. This decays geometrically (ratio ≈ 0.82–0.87), so the operator complexity and the cost per cycle are
bounded by a convergent series, i.e. O(n) per cycle. The bound is not yet reached at 1M unknowns: WU per cycle still
rises 20% from m=511 to 1023, and the series limit is about 9.5 in operator complexity.
**Convergence.** The factor creeps from 0.55 to 0.70, inherited from the two-level factor (0.56, 0.61, 0.65), so the
cycle counts grow (31 → 50). Total work is therefore not yet mesh-independent.
**Two more points.**
- The coarse levels need the robust settings. Cheap coarse cycles ($k=2$, truncation 1e-2) stall at m ≥ 255, because
  N's error contains the coarse-cycle defect acting on $c_0$.
- The family has fixed cell Pe (10), as in §7.3, so the physical Pe grows with $m$.

### Claim 11b. Levers against the factor creep (`tv_test.py`, `wscan.py`, `wc_mesh.py`)

- **Test vectors for $\tilde S$** (adding $\psi$, quadratics, or N's slowest mode): no gain. With 2 entries per row N
  stays at 0.55–0.61 at m=127/255. With 3 entries N gains less than 0.03 and the twin starts to converge.
- **Split threshold δ:** no gain. **$\omega_F$:** must be exactly 1 (0.9 gives $\rho(K)>0$ and N stalls at 0.99).
- **C-weight $\omega_C$:** helps. Two-level $\rho(E_N)$ with $\omega_C=0.5$/0.6: 0.29/0.29 (m=63), 0.31/0.34 (127),
  0.33/0.38 (255), 0.49/0.48 (511). Twin 0.82–1.22, Markov-Jacobi 2.7–4.7.
- **Multilevel with $\omega_C=0.6$:**

| m | N: factor / cycles | WU total | twin |
|---|---|---|---|
| 127 | 0.301 / 16 | 1,872 | diverges |
| 255 | 0.343 / 18 | 3,012 | diverges |
| 511 | 0.435 / 22 | 4,375 | diverges |
| 1023 | 0.703 / 51 | 12,175 | diverges |

  Total work is halved up to m=511. At m=1023 the factor returns to 0.70. An extra coarse cycle (refinement) does not
  change it (0.703), so the limit is the finest level, not the coarse solve. The creep is reduced but not removed.

## Claim 12. Truncated multilevel N is always a stationary Markovian cycle with a composite smoother (`check_hat.py`)

N's update regroups as $\hat x=x^{(k)}+Qy-\tilde Pd^0$ followed by the standard correction, because
$R\,r(x^{(k)}+h)+A_cd^0=R\,r(\hat x)$. Hence $E_N=(I-\tilde PB_cRA)(I-\hat MA)$ with
$$\hat T=\tilde PP^\sharp+Q(S-\tilde S)P^\sharp T^k+QK^kQ^\sharp,\qquad
\hat M=\sum_{j<k}T_Q^jQQ^\sharp M+Q(S-\tilde S)P^\sharp\sum_{j<k}T^jM .$$
This holds recursively on every level, for any truncation, filtering and coarse cycle. $\hat M$ uses only the
residual: $k$ F-sweeps, $k$ global sweeps, and $S$ applied matrix-free through $k$ products with $K$ and $H$.
Check: the iterates of N and of the cycle with $\hat M$ agree to 4e-15 (Claim 11 hierarchy, m=63/127, converging) and
3e-14 (lump 0.5, k=4, both diverging to 1e22).
Consequences:
- $\hat M$ is the twin's smoother plus a correction that uses the **untruncated** $S$. Truncation in N is truncation
  of the coarse correction only; its effective smoother still applies the exact $S$ (through the history), which
  moves the defect from $c_0$ to $c_k$.
- $\hat T$ is not the $k$-th power of one sparse stationary step. It contains $k$ global sweeps and $k$ F-sweeps
  combined through $Q(S-\tilde S)P^\sharp$, and costs as much as N, or slightly more because it recomputes the
  sweeps from 0.

## Claim 13. O(n), mesh-independent multilevel N where the twin and Markov-Jacobi on the same smoother diverge (`ml_rb.py`)

**Problem.** Centered convection–diffusion on $(-1,1)^2$, constant flow $\beta=(\cos0.4,\sin0.4)$, cell Pe 6.

**Hierarchy, shared by all methods** ($R=\tilde P^T$):
- Level 0: red-black split (F = red, C = black), so $K=0$ and $G=0$ exactly. Jacobi with $\omega_F=1$,
  $\omega_C=0.2$, $k=2$. $\tilde S$ = one entry per F-row carrying the full row sum of $S$.
- Coarse levels: greedy split (δ=0.3), Jacobi $\omega=1$, $k=2$, $S$ truncated at 1e-2, filter 1e-2, V-cycle, sparse
  LU below $16\sqrt n$.

| m | n | levels | N: factor / cycles | WU/cycle | WU total | op.cx | twin | Markov-Jacobi |
|---|---|---|---|---|---|---|---|---|
| 127 | 16,129 | 6 | 0.365 / 19 | 44 | 840 | 4.16 | 1.36, diverges | 2.06, diverges |
| 255 | 65,025 | 7 | 0.364 / 19 | 49 | 925 | 4.48 | 1.40, diverges | 2.07, diverges |
| 511 | 261,121 | 9 | 0.363 / 19 | 53 | 1,013 | 4.82 | 1.39, diverges | 2.07, diverges |
| 1023 | 1,046,529 | 10 | 0.362 / 19 | 55 | 1,053 | 4.93 | 1.38, diverges | 2.07, diverges |

Two-level spectral radii ($\rho(E_N)$ / twin / MJ): 0.403 / 1.21 / 1.76 (m=63), 0.407 / 1.35 / 2.03 (127),
0.406 / 1.39 / 1.81 (255).

- **Cost.** nnz(A_ℓ)/nnz(A) at m=1023: 1.00, 0.70, 0.72, 0.82, 0.61, 0.42, 0.28, 0.18, 0.12, 0.08. This is a geometric
  decay (ratio ≈ 0.65), so the operator complexity is bounded (≈ 5). The cycle count and factor are constant.
  **N is O(n).**
- **Caveat: the best Markovian cycle on the same hierarchy converges.** Sweeping $\omega_F$, $\omega_C$ and
  $k\le16$ (`rb_competitors.py`, `rb_mjwork.py`), damped Markov-Jacobi with $\omega_F=0.8$, $\omega_C=0.16$, $k=16$
  reaches 0.19–0.20, 12 cycles, about 500 WU, which is cheaper than N (840–925 WU). So this example shows N
  converging in O(n) where the twin and Markov-Jacobi on the **same** smoother diverge. It does not show that no
  Markovian cycle converges well.
- Higher cell Pe (8–16), where damped Jacobi would be useless: with the single-entry $\tilde S$, N fails; with
  2 entries, the twin converges (`rbpe.py`).

## Claim 14. Dichotomy: N against two simple Markovian cycles on the same T, Q and $\tilde P$ (`dichotomy_check.py`, `first_order.py`)

Two-level, exact coarse solve, any relaxation $T$, $P_0/Q_0$, any truncation $\tilde S$. Let ${\rm MJF}_k$ be the
Markovian cycle "k sweeps of $T$, then k F-sweeps ($T_Q$), then the $\tilde P$ correction". Using
$T_Q^m=P'_mP^\sharp+QK^mQ^\sharp$ and $(I-\tilde\Pi)\tilde P=0$:
$$E_N=(I-\tilde\Pi)Q\big[(S_k-\tilde S)P^\sharp T^k+K^kQ^\sharp\big],\qquad
E_{{\rm MJF}_k}=(I-\tilde\Pi)Q\big[(S_k-\tilde S)P^\sharp T^k+K^kQ^\sharp T^k\big],$$
$$\boxed{E_N-E_{{\rm MJF}_k}=(I-\tilde\Pi)\,QK^kQ^\sharp\,(I-T^k)}$$
With $K=0$ (red-black and $\omega_F=1$, or a greedy split with an independent F-set), one F-sweep suffices and
**$E_N=E_{\rm MJF}$ exactly**.
Check (32 cases, recirc Pe 6/20, greedy δ 0.3/0.6, k 2/4, $\omega_C$ 0.3/1, 1–2 entries): the identity holds to
1e-15. N and ${\rm MJF}_k$ agree to 3 digits when $\rho(K)=0$; with $\rho(K)=0.43$, ${\rm MJF}_k$ is equal or slightly
better (0.274 vs 0.344). On the red-black candidates (m=127), N 0.405 / 0.359 / 0.689 against "MJ + 1 F-sweep"
0.407 / 0.359 / 0.689.

**Consequence (with Claim 7, twin ≡ N without truncation).**
- **Untruncated $P'$:** N ≡ twin, even multilevel.
- **Truncated, $T$ contracts the C-error:** N ≈ ${\rm MJF}_k$. The difference is the $K^k$ term, which N needs to be
  small anyway.
- **Truncated, $T$ amplifies the C-error:** the twin beats N. Its truncation defect acts on $c_0$, N's on the
  amplified $c_k$.

So in comparisons that allow these simple Markovian cycles on the same $T$, $Q$ and $\tilde P$ (k T-sweeps, F-sweeps,
correction), N cannot win by a wide margin. The examples where "only N converges" (Claims 10, 11, 13) are N against
the twin and plain Markov-Jacobi only; ${\rm MJF}_k$ reproduces them. The one structurally non-Markovian advantage left
is the step-0 coarse unknown under inexact multilevel coarse solves (Claim 7c), which matters only when $T$
amplifies, where the twin is already as good as N.

**Addendum to Claim 14, multilevel** (`results_mjf_ml.md`). On the Claim 11 and Claim 13 hierarchies, the Markovian
cycle "k sweeps of T + 1 F-sweep + $\tilde P$ correction" (MJF, applied recursively) has **exactly the same factor and
cycle count as N** and costs less:

| hierarchy | m | N: factor / cycles / WU | MJF: factor / cycles / WU |
|---|---|---|---|
| Claim 11 (recirc Pe 10, ω_C 0.6) | 127 | 0.301 / 16 / 1872 | 0.301 / 16 / 1263 |
| Claim 11 | 255 | 0.343 / 18 / 3012 | 0.343 / 18 / 1968 |
| Claim 13 (const Pe 6, red-black) | 127 | 0.365 / 19 / 840 | 0.365 / 19 / 666 |
| Claim 13 | 255 | 0.364 / 19 / 925 | 0.364 / 19 / 709 |

**Correction to Claims 10, 11 and 13.** "Twin and Markov-Jacobi diverge, N converges" holds, but a simple Markovian
cycle on the same smoother and hierarchy (one extra F-sweep before the correction) converges identically and more
cheaply. Those examples do not show an advantage of the non-Markovian formalism over Markovian cycles in general.

**Addendum 2 to Claim 14, black-box setting** (`bb_mjf.py`, `results_bb_mjf.md`). In Claim 4's setting (recirc Pe 20,
undamped Jacobi with $\rho(T)\approx4.5$, robust hierarchy), MJF needs only the black-box sweeps plus one *algebraic*
F-Jacobi step (no smoother restart). It converges exactly like N and more cheaply: m=127: N 2.2e-6, 2 cycles, 646 WU
vs MJF 2.2e-6, 2 cycles, 445 WU; m=255: N 3.2e-3, 3 cycles, 1359 WU vs MJF 3.2e-3, 3 cycles, 901 WU.
**Correction to Claim 4:** N is not the only convergent method under the black-box constraint; MJF also is.

## Claim 15. Natural examples, separately: (A) N against the twin, (B) N against plain Markov-Jacobi (`natural_scan.py`, `natural_scanA.py`, `natural_mesh.py`, `natural_ml.py`)

"Natural" means centered convection–diffusion (constant flow, recirculating flow, Poiseuille flow), the **same**
damped-Jacobi weight ω on all points, a standard split (red-black or greedy δ=0.3), and a standard truncation of $S$
(one-point: one entry carrying the row sum; drop 0.3; lump 0.5; exact).

From the error formulas (Claim 7 and Claim 14):
$$E_N=(I-\tilde\Pi)Q\big[(S-\tilde S)c_k+K^kq_0\big],\quad E_{\rm twin}=(I-\tilde\Pi)Q\big[(S-\tilde S)c_0+K^kq_0\big],\quad E_{\rm MJ}=(I-\tilde\Pi)T^k .$$
- **(A) N converges, twin diverges.** Needs a crude $\tilde S$ (one-point) and a smoother that *damps* the C-error, so
  that $c_k\ll c_0$. Damped Jacobi inside its stability range, $\omega<2/(1+b^2)$, does this, and then Markov-Jacobi
  converges too.
- **(B) N converges, MJ diverges.** Needs a smoother that *amplifies* ($\rho(T)>1$: centered convection at cell Pe ≳ 4
  with ω ≥ 0.6) and an accurate $\tilde S$. Then N and the twin coincide or nearly so, and both converge.

**Two-level scan** (m=31, dense, `results_natural_scan31.md`, `results_natural_scanA31.md`):
- Class (B) is everywhere: 200+ configurations, e.g. constant flow Pe 6, ω=0.8, k=2: N 0.04, twin 0.04, MJ 3.1;
  recirc Pe 20, ω=1: N = twin ≈ 0 (K=0, exact P'), MJ 20–500.
- Class (A) alone never occurs with ω ≥ 0.6; every twin-divergent case there also has MJ divergent. With damped
  Jacobi inside its stability range it does occur, e.g. constant flow Pe 8, ω=0.1, k=16, one-point: N 0.56, twin 1.96,
  MJ 0.43.

**Mesh check, two-level** (`results_natural_mesh.md`):

| case | m=31 | m=63 | m=127 (N / twin / MJ) |
|---|---|---|---|
| (A) const Pe 8, ω 0.1, k 16, one-point | 0.56 / 1.96 / 0.43 | 0.577 / 2.08 / 0.634 | 0.602 / 2.31 / 0.658 |
| (B) const Pe 6, ω 0.8, k 2, drop 0.3 | 0.042 / 0.055 / 3.18 | 0.063 / 0.056 / 3.38 | 0.065 / 0.056 / 3.54 |
| (B) recirc Pe 20, ω 1, k 2, exact | 0 / 0 / 20.4 | 1e-14 / 1e-14 / 21.8 | 6e-14 / 7e-14 / 22.9 |

**Multilevel** (`natural_ml.py`, `results_natural_ml*.md`):

(A2) constant flow Pe 8, damped Jacobi ω=0.1 and k=16 on **all** levels, greedy, one-point $\tilde S$ on level 0,
$S$ truncated at 1e-2 below, filter 1e-2:

| m | n | N: factor / cycles | WU/cycle | WU total | op.cx | twin | Markov-Jacobi |
|---|---|---|---|---|---|---|---|
| 127 | 16,129 | 0.573 / 33 | 163 | 5,391 | 4.02 | 2.37, diverges | 0.782, 74 cycles, 8,594 WU |
| 255 | 65,025 | 0.576 / 33 | 184 | 6,065 | 4.29 | 2.38, diverges | 0.821, not converged in 80 cycles |
| 511 | 261,121 | 0.59 / 35 | 199 | 6,978 | 4.45 | 2.37, diverges | 0.828, not converged in 80 cycles |

With undamped Jacobi on the coarse levels instead (A), MJ diverges as well (18–20), so the example becomes AB.

(B2) constant flow Pe 6, Jacobi ω=0.8 and k=2 on all levels, greedy, $S$ truncated at 0.1, filter 1e-3:

| m | n | N: factor / cycles | WU/cycle | WU total | op.cx | twin | Markov-Jacobi |
|---|---|---|---|---|---|---|---|
| 127 | 16,129 | 0.060 / 7 | 56 | 393 | 5.77 | 0.050 / 7 / 163 WU | 3.08, diverges |
| 255 | 65,025 | 0.066 / 7 | 65 | 454 | 6.64 | 0.076 / 8 / 189 WU | 3.10, diverges |
| 511 | 261,121 | 0.124 / 9 | 69 | 619 | 6.94 | 0.163 / 11 / 273 WU | 3.10, diverges |

**Caveats (Claim 14).**
- In (A) the smoother contracts, so the Markovian cycle "k sweeps + F-sweeps" (MJF) ≈ N.
- In (B) the twin ≈ N and is cheaper.

Each class shows N beating one specific Markovian cycle in a natural setting, not every Markovian cycle.

### Claim 15, addendum: m=1023 and other flows (`natural_ml.py`, `natural_ml2.py`, `results_natural_ml*.md`)

**Constant flow at m=1023** (n=1,046,529):

| case | N: factor / cycles | WU/cycle | WU total | op.cx | twin | Markov-Jacobi |
|---|---|---|---|---|---|---|
| (A2) Pe 8, ω 0.1, k 16 | 0.68 / 47 | 201 | 9,444 | 4.48 | 2.37, diverges | 0.842, not converged in 80 cycles |
| (B2) Pe 6, ω 0.8, k 2 | 0.283 / 14 | 71 | 991 | 7.21 | 0.355 / 17 / 411 WU | 3.11, diverges |

Both classes hold at 1M unknowns, and the cost per cycle is bounded (A2: 163, 184, 199, 201 WU; B2: 56, 65, 69, 71).
The factors creep at the largest size (A2: 0.57, 0.58, 0.59, 0.68; B2: 0.06, 0.07, 0.12, 0.28).

**Recirculating flow.**
- (B) is the robust-N problem (§7.3, Pe 20, undamped Jacobi, k = 4/8; ROBUST_N.md, Claim 4): N converges in 2–4
  cycles (factor 2.6e-6 to 3.7e-3) up to m=1023. Markov-Jacobi diverges (5e9–2.5e12). The twin (M_F) converges and is
  cheaper. Cost per cycle 173–709 WU (op.cx 4.2–9.4).
- (A) Pe 6, damped Jacobi ω=0.3 and k=16 on all levels, greedy, lump-0.5 $\tilde S$:

| m | N: factor / cycles / WU | twin | Markov-Jacobi |
|---|---|---|---|
| 127 | 0.15 / 10 / 1,225 | 1.26, diverges | 0.36 / 18 / 1,567 WU |
| 255 | 0.247 / 13 / 1,909 | 1.45, diverges | 0.465 / 24 / 2,498 WU |
| 511 | 0.597 / 33 / 5,088 | 1.54, diverges | 0.581 / 33 / 3,606 WU |

  The twin diverges at every size. N creeps, and at m=511 Markov-Jacobi catches up (the closed-streamline effect,
  iteration 1 of SEARCH_LOG).

**Poiseuille flow.**
- (B) Pe 10, Jacobi ω=0.8, k=2, S truncated at 0.01: N = twin 0.294 → 0.477 (m=127 → 255); Markov-Jacobi diverges
  (7.2). More accurate coarse levels (finer truncation and filtering, refinement, k_c=8) keep N at 0.28–0.29 for
  m=127 and 0.44–0.47 for m=255.
- (A) Pe 8, ω=0.1, k=16, red-black, one-point: twin diverges (1.43–1.48); N 0.52 → 0.77 and Markov-Jacobi 0.60 → 0.77
  (equal at m=255, Markov-Jacobi cheaper).

**Summary.** The class statements, (A) twin diverges while N converges and (B) Markov-Jacobi diverges while N
converges, hold for all three flows at every size tested. Only the constant-flow examples (and recirc (B), robust
N) are close to mesh-independent. With cell Pe held fixed, the physical Péclet number grows like $1/h$. The
recirculating (closed streamlines) and Poiseuille (wall-bounded) flows then degrade the coarse-grid correction, so
N's factor creeps, and in (A) Markov-Jacobi catches up by m=255–511.

## Claim 16. M vs S: memory in the interpolation only (`semi_scan.py`, `semi_ml.py`, `semi_twin.py`)

**Setup.** Same smoother (damped Jacobi, uniform ω), same greedy split, same restriction
$\widehat R=[I,\;-A_{cf}(I+G+G^2)D_{ff}^{-1}]$ with $G=I-D_{ff}^{-1}A_{ff}$, same coarse operator $\widehat RAP_0$
($P_0$ = injection).
- **M:** after k sweeps, $x\leftarrow x+P_0(\widehat RAP_0)^{-1}\widehat Rr$, so only the C-points are corrected.
- **S** (Alg. 4.1 / 7.2): the same coarse solve, then the correction is interpolated to the F-points from the
  stored iterates: $x_F\leftarrow x_F+y+\sum_{\ell<k}K^\ell Hz$ (shift relation and Horner replay).

$$E^M=\Pi_RT^k,\qquad E^S=QK^kQ^\sharp=(\Pi_RT\Pi_R)^k,\qquad \Pi_R=I-P_0(\widehat RAP_0)^{-1}\widehat RA .$$
M carries the full relaxation $T^k$, which diverges when $\rho(T)>1$. S carries only the CR operator
$\Pi_RT\Pi_R$ on range($\Pi_R$). S equals N with the complement induced by $\widehat R$ ($\widehat RAQ=0$).

**S's Markovian twin is expensive.** Since $\Pi_R^2=\Pi_R$, $(I-\Pi_RMA)\Pi_R=\Pi_RT\Pi_R$, and so
$$E^S=(I-\Pi_RMA)^k\,\Pi_R$$
(`semi_twin.py`: holds to 1e-14 to 3e-13). The equivalent Markovian cycle is "coarse correction, then k projected
relaxation steps $x\leftarrow x+\Pi_RM(b-Ax)$". Each step applies $\Pi_R$, i.e. one coarse solve, so it needs
**k+1 coarse solves per cycle against S's one**. Unlike N (Claims 1 and 14), S has no cheap Markovian equivalent.

**Two-level scan** (m=31, `results_semi_scan31.md`):
- **Rescue** (M diverges, S converges) for all centered CD at cell Pe 6–20 with ω ≥ 0.6 (M 1.5–1700, S 0.002–0.16).
- **Cost** (both converge, S ≥ 2× cheaper per digit): anisotropic Poisson with ω=1 (M 0.94–0.98, S 1e-4–0.26,
  30–170×); recirc Pe 2 with ω=1 (M 0.98–0.99, S 0.006–0.29, 54–130×); Poisson ω=1, L=3, k=4 (M 0.79, S 0.12, 4×).

**Multilevel, Alg. 7.2 vs its Markovian version** (`results_semi_ml*.md`):

| problem | smoother | m | S: factor / cycles / WU | S WU/cycle | op.cx | M |
|---|---|---|---|---|---|---|
| shear Pe 20 (paper Table 4) | ω 0.7, k 4 | 127 | 0.072 / 7 / 914 | 131 | 11.6 | diverges (123) |
| | | 255 | 0.12 / 9 / 1,254 | 139 | 11.8 | diverges (130) |
| | | 511 | 0.167 / 11 / 1,607 | 146 | 12.2 | diverges (133) |
| constant flow Pe 6 | ω 0.8, k 2 | 127 | 0.161 / 11 / 619 | 56 | 7.5 | diverges (2.1) |
| | | 255 | 0.227 / 13 / 841 | 65 | 8.1 | diverges (2.2) |
| | | 511 | 0.286 / 15 / 1,087 | 72 | 8.7 | diverges (2.5) |
| shear Pe 6 | ω 0.8, k 2 | 127 | 0.216 / 12 / 430 | 36 | 5.0 | diverges (1.43) |
| | | 255 | 0.322 / 17 / 706 | 42 | 5.5 | diverges (1.66) |
| | | 511 | 0.427 / 21 / 918 | 44 | 5.6 | diverges (1.75) |
| shear Pe 2 (cost) | ω 0.8, k 2 | 127 | 0.421 / 21 / 979 | 47 | 5.9 | 0.929, about 3,200 WU (est.) |
| | | 255 | 0.59 / 35 / 1,832 | 52 | 6.3 | 0.980, about 11,900 WU (est.) |
| | | 511 | 0.716 / 53 / 2,833 | 53 | 6.3 | 0.977, about 10,300 WU (est.) |

M's WU for shear Pe 2 is extrapolated from its factor over 80 cycles (13 WU/cycle).
**Weaker cases.** Anisotropic Poisson (ω 0.8 k 4: S 0.40, 0.74, 0.88; M 0.99 → diverges) and Poisson / recirc
Pe 2–6 creep or stall for S in this hierarchy.

**Reading.**
- **Rescue:** adding memory to the interpolation turns a divergent cycle (M, $\rho(T)$ up to 5) into a convergent
  one at modest extra cost per cycle (S/M WU per cycle ≈ 3–4).
- **Cost:** where M converges slowly (the injection correction leaves F-error to the smoother), S is 3–6× cheaper
  in total work.
- S's factor creeps with m in this hierarchy, as in the fixed-cell-Pe families of Claim 15.

### Claim 16, addendum: M with a serious interpolation (`semi_vsP.py`, `mlm_goodP.py`, `semi_ml.py`, `semi_pyamg.py`)

Injection makes M a straw man, so M is rerun with the same smoother, split and restriction
$\widehat R=[I,\,-A_{cf}(I+G+G^2)D_{ff}^{-1}]$ but a good interpolation $P=[I;W]$:
- $W_1=-D_{ff}^{-1}A_{fc}$ (direct / one-step Jacobi);
- $P'=P_0+Q_0\sum_{\ell<k}K^\ell H$, S's own interpolation made static (truncated at 1e-2 in multilevel);
- in two-level also the ideal $W=-A_{ff}^{-1}A_{fc}$.

The structural reason S survives: $E^M(P)=(I-\Pi_P)T^k$ contains the divergent relaxation for **every** $P$, while
$E^S=(\Pi_RT\Pi_R)^k$ contains only the CR operator.

**Two-level** (m=31, `results_semi_vsP31.md`):

| problem | ρ(T) | M(P0) | M(W1) | M(P') | M(ideal) | S(P0) |
|---|---|---|---|---|---|---|
| shear Pe 20 | 3.87 | 214 | 165 | 161 | 165 | 0.0081 |
| const Pe 6 | 1.29 | 1.59 | 0.77 | 0.81 | 0.77 | 0.04 |
| shear Pe 6 | 1.02 | 0.94 | 0.45 | 0.47 | 0.45 | 0.04 |
| shear Pe 2 | 0.97 | 0.64 | 0.34 | 0.35 | 0.34 | 0.04 |
| Poisson | 1.00 | 0.72 | 0.25 | 0.15 | 0.15 | 0.20 |

(S with $P=[I;W_1]$ diverges on the convection cases: S relies on the injection-induced complement.)

**Multilevel** (`results_semi_goodP*.md`):

| problem | m | S (injection + memory): factor / cycles / WU | M, $W_1$ | M, static $P'$ |
|---|---|---|---|---|
| shear Pe 20, ω 0.7, k 4 | 127 | 0.072 / 7 / 914 | diverges (270) | diverges (241) |
| | 255 | 0.12 / 9 / 1,254 | diverges (272) | diverges (246) |
| | 511 | 0.167 / 11 / 1,607 | diverges (276) | diverges (251) |
| const Pe 6, ω 0.8, k 2 | 127 | 0.161 / 11 / 619 | diverges (3.3) | diverges (2.9) |
| | 255 | 0.227 / 13 / 841 | diverges (3.8) | diverges (3.2) |
| | 511 | 0.286 / 15 / 1,087 | diverges (4.2) | diverges (3.4) |
| shear Pe 6, ω 0.8, k 2 | 127 | 0.216 / 12 / 430 | 0.789 / 79 / 1,731 | 0.799, not converged in 80 |
| | 255 | 0.322 / 17 / 706 | diverges (1.43) | diverges (1.15) |
| | 511 | 0.427 / 21 / 918 | diverges (2.4) | diverges (2.1) |
| shear Pe 2, ω 0.8, k 2 | 127 | 0.421 / 21 / 979 | 0.458 / 24 / 645 | 0.441 / 23 / 621 |

In multilevel, M diverges with a good interpolation even where its two-level version converges (const Pe 6, shear
Pe 6). The divergent relaxation acts on every level, and the coarse levels have larger effective Pe. At low Pe
(shear Pe 2), M with a good interpolation is cheaper than S.

**External Markovian references** (PyAMG defaults, `results_semi_pyamg.md`):
- shear Pe 20: AIR (4.6–5.3) and RS fail, also as GMRES preconditioners.
- const Pe 6 and shear Pe 6: AIR converges cheaply (90–143 WU, using its own F/C relaxation); RS fails.

**Statement that survives the fair comparison.**
- With the same relaxation and restriction, **no interpolation** (injection, direct, S's own $P'$ made static, or
  ideal in two-level) makes the Markovian cycle converge on centered convection at cell Pe 6–20 in multilevel.
  Adding memory to the interpolation (S, injection + history replay) converges at 0.07–0.43 up to m=511.
- At shear Pe 20, standard AMG (AIR, RS) fails as well.
- At moderate Pe, a Markovian AMG with a different relaxation (AIR's F/C relaxation) is far cheaper. S's advantage is
  over Markovian cycles on the **same** relaxation.
