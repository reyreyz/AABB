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

## Not demonstrated

- **A dense $Q^\sharp$ with sparse $K$, $H$**, the only flop-count case where N could be cheaper than
  its twin: no natural example is known.
- **Parallel or latency effects:** not measured.
