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
