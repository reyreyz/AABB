# Search log: O(n) multilevel example where only N converges

## Target
Multilevel, sparse and bounded throughout:
- **N:** converges with a mesh-independent factor, with cycles to 1e-8 and WU/cycle bounded up to m ≥ 511
  (ideally 1023).
- **Markovian twin** (k F-sweeps + $\tilde P$ correction) on the same smoother and hierarchy: diverges.
- **Markov-Jacobi** (k sweeps of the same T + $\tilde P$ correction) on the same hierarchy: diverges.
- Report the best Markovian variant on the same $\tilde P$ as well (any ω, $\omega_C$, k).

## Key formulas (from CLAIMS.md 7, 7b, 12)
$$E_N=(I-\tilde PB_cRA)\big[\tilde PP^\sharp+Q(S-\tilde S)P^\sharp T^k+QK^kQ^\sharp\big]$$
$$E_{\rm twin}:\ \text{same with } P^\sharp T^k\to P^\sharp,\qquad E_M=(I-\tilde\Pi)T^k .$$
What each needs:
- Twin diverges: heavy truncation, nonsymmetric $A$, defect on raw oscillatory $c_0$.
- Markov-Jacobi diverges: $T^k$ grows away from range($P'$).
- N converges: $\rho(K)<1$, $\tilde S$ accurate on $c_k$, and an accurate coarse cycle.

## State at start of the autonomous search (Claim 11 / 11b)
Best so far: centered recirc, cell Pe 10. Level 0: Jacobi $\omega_F=1$, $\omega_C=0.6$, $k=1$, δ=0.3, $\tilde S$
with 2 entries per row exact on 1,x,y. Coarse levels: robust N ($k=8$, $S$ at 1e-5, filter 1e-3/3e-4, LU below
$16\sqrt n$). Script: `ml_missing_final.py m refine methods wC`.

| m | N factor / cycles | WU/cycle | twin | MJ |
|---|---|---|---|---|
| 127 | 0.301 / 16 | 117 | div | div |
| 255 | 0.343 / 18 | 167 | div | div |
| 511 | 0.435 / 22 | 199 | div | div |
| 1023 | 0.703 / 51 | 239 | div | div |

Problems:
- The factor creeps; it is a finest-level effect (refinement does not help).
- WU/cycle is still rising (geometric bound only).

## Idea queue (mark [x] when tried, add results below)
- [x] 11. Flows without closed streamlines (constant, Poiseuille, hyperbolic) at moderate cell Pe (5-8): the recirc
          creep is a closed-streamline coarse-grid effect.
- [x] 1. Diagnose the slow mode of two-level N at m=511/1023 with ω_C=0.6: smooth or oscillatory, where it is
         localized (vortex centre, boundary, stagnation points).
- [ ] 2. Choose the kept entries of $\tilde S$ by geometry (nearest C-points spanning a simplex) rather than magnitude.
- [x] 3. Constant-flow centered CD instead of recirc (no vortex centre / stagnation), cell Pe 5/10/20.
- [ ] 4. Make ω_C row-dependent (local cell Pe), or optimize ω_C per mesh.
- [ ] 5. Smaller δ (0.1–0.2), or the F-set chosen to make $\|K\|_\infty$ ≈ 0 exactly.
- [ ] 6. k=2 with very small ω_C (0.1–0.3), with C-weights tuned.
- [ ] 7. Other smoothers on level 0: SOR/GS-type with C/F weights, two-stage C/F Jacobi.
- [ ] 8. Other problems: 3-D centered convection, convection with reaction, Helmholtz-type with convection.
- [ ] 9. Bound the coarse cost: per-level stencil caps on coarse levels while keeping the coarse factor ≤ 0.2.
- [ ] 10. Variants of the coarse-correction restriction R (P~^T, hybrid) on level 0 only.

- [x] 12. Red-black split on level 0 (K = 0 and G = 0 exactly), single-entry lumped S~, k = 1-3, w_C scan.
- [x] 13. Multilevel O(n) version of the red-black candidate, with robust coarse levels.

- [ ] 14. STRONG target: also beat every damped Markovian cycle on the same hierarchy (MJ over w_F, w_C, k <= 32 and the
          damped twin) by a clear margin in WU. Needs cell Pe where damped Jacobi is useless (Pe >= 10-12), with N still
          converging under a crude S~.
- [x] 15. At Pe >= 8: S~ with 2 entries per row exact on constants only (tv = empty), or lumping at tau 0.5, so the twin
          stays crude while N gets 2 entries.
- [ ] 16. Red-black at Pe 6 with k = 3-4 and small w_C: N cheaper per digit? Compare against the MJ k=16 competitor.
- [x] 17. Other problem classes where damping cannot help: indefinite (Helmholtz-like) plus convection, or reaction
          with negative sign, combined with the red-black mechanism.

- [ ] 18. Blended S~ = (1-theta) S~_1 + theta S~_2 (both exact on constants) at Pe 16-20, k 8-16, small w_C: widen the
          margin over damped Markov-Jacobi.

## Results of the autonomous search
(append: date/iteration, idea, setup, numbers, verdict)

### Iteration 1
- **Idea 1, slow modes of two-level N** (recirc, ω_C=0.6; `diag_slow.py`, `results_diag_slow.md`):
  - (a) an oscillatory complex pair localized at the boundary where |β| is largest (|λ|=0.37 at m=255, 0.40 at
    m=511; peak near (0.19,-0.97));
  - (b) from m=511 on, a **smooth** real mode on the closed streamlines (0.2<r<0.7), |λ|=0.48, which dominates at
    m=1023 (0.70).
  (b) is the classical recirculating-flow coarse-grid difficulty: with fixed cell Pe the physical Pe grows like 1/h.
  So the creep is a closed-streamline effect, not an N effect.
- **Idea 3, constant flow at uniform cell Pe 10/20** (`constflow.py`): N diverges, or the twin converges. Too hard
  everywhere.
- **Idea 11, other flows at m=127** (`flowscan.py`, `results_flowscan.md`, two-level, ω_C 0.4/0.6, 2 entries/row):
  - Poiseuille Pe 10: N 0.39, but the twin converges (0.53).
  - Hyperbolic Pe 10: N 0.79, twin 2.2, MJ 3.7.
  - **Constant flow Pe 6, ω_C=0.6: N 0.25, twin 1.18, MJ 2.44, a hit.**
  - Constant flow Pe 3/4: the twin converges.
- Mesh study of constant flow (`constmesh.py`): at m=63 the twin is 0.95–0.99 for Pe 6 and converges for Pe 5 and
  Pe 8. m=255 and m=511 still running.

### Iteration 2
- **Constant flow, greedy split, Claim 11b setup** (`constmesh.py`): Pe 5: N 0.22-0.28 at m=255, but the twin
  converges (0.64-0.72). Pe 6: twin 0.95-0.99 at m=63 and **0.985 at m=511** (N 0.434). Marginal, dropped.
- **Idea 12, red-black split at m=63** (`redblack.py`, `results_redblack.md`): several hits with the single-entry
  S~ (twin 1.2-1.55).
- **Red-black mesh study** (`rbmesh.py`, `results_rbmesh.md`; columns N / twin / MJ):

| case | m=63 | m=127 | m=255 |
|---|---|---|---|
| recirc Pe 10, k=1, wC 0.8 | 0.30 / 1.42 / 2.44 | 0.36 / 1.73 / 2.97 | 0.81 / 2.09 / 3.46 (creeps) |
| const Pe 6, k=1, wC 1.0 | 0.40 / 1.22 / 2.64 | 0.58 / 1.35 / 3.11 | 0.72 / 1.38 / 3.45 (creeps) |
| const Pe 6, k=1, wC 0.8 | | 0.56 / 1.35 / 2.48 | 0.56 / 1.41 / 2.99 |
| **const Pe 6, k=2, wC 0.2** | **0.403 / 1.21 / 1.76** | **0.407 / 1.35 / 2.03** | **0.406 / 1.39 / 1.81** |
| recirc Pe 20, k=1, wC 0.4, 2 entries | 0.78 / 2.65 / 8.76 | 0.85 / 3.47 / 9.17 | 1.02 (fails) |

  **Candidate:** constant flow (angle 0.4), cell Pe 6, red-black level 0, k=2, ω_F=1, ω_C=0.2, single-entry
  lumped S~ (lin, nkeep=1). Two-level N is mesh-independent at 0.41; twin and MJ diverge. m=511 is running.

### Iteration 3
- **Idea 13, multilevel red-black candidate** (`ml_rb.py`, `results_ml_rb*.md`): constant flow Pe 6, cheap coarse
  levels.

| m | N factor / cycles | WU/cycle | WU total | op.cx | twin | MJ |
|---|---|---|---|---|---|---|
| 127 | 0.365 / 19 | 44 | 840 | 4.16 | div | div |
| 255 | 0.364 / 19 | 49 | 925 | 4.48 | div | div |
| 511 | 0.363 / 19 | 53 | 1013 | 4.82 | div | div |
| 1023 | 0.362 / 19 | 55 | 1053 | 4.93 | div | div |

  **O(n) and mesh-independent, recorded as Claim 13.**
- **Competitors** (`rb_competitors.py`, `rb_mjwork.py`): damped Markov-Jacobi on the same hierarchy (wF 0.8,
  wC 0.16, k=16) converges at 0.19-0.20 for about 500 WU, cheaper than N. The **strong** target is not met.
- **Higher Pe with red-black** (`rbpe.py`, m=63): at Pe 8-16 only 2 entries per row give N < 0.8, and then the twin
  converges (0.57-0.94).

### Iteration 4
- **Idea 15, 2 entries per row exact on constants only** (`rbpe2.py`): at Pe 8-12 the twin always converges
  (0.50-0.73). Not crude enough.
- **Blended S~_θ = (1-θ)S~_1 + θS~_2** (`blend.py`, `blend2.py`): θ tunes how crude S~ is. At m=63, Pe 12,
  red-black, k=8, wC=0.05: θ=0.5: N 0.359, twin 2.18, MJ 1.43. θ=0.6: 0.371 / 1.67 / 1.26. θ=0.7: 0.393 / 1.22 / 1.14.
  At Pe 10, MJ on the same T often converges.
- **Mesh check** (`blendmesh.py`; θ=0.5, N / twin / MJ): 0.359 / 2.24 / 1.65 at m=127 and 0.357 / 2.24 / 1.71 at
  m=255. Two-level N is mesh-independent.
- **Strong test, m=127** (`strong_test.py`, multilevel, cheap coarse levels): N 0.397, 20 cycles, **1159 WU**. Best
  damped Markovian over (wF, wC, k ≤ 32): MJ wF 0.15, wC 0.03, k=32, 0.545, **1496 WU**. Only 1.3× more work; the
  margin needs to be larger.

### Iteration 5
- **Pe 16-20, blended S~** (`blend3.py`): for k ≥ 12, MJ on the same T converges about as well as N.
  **Strong test at Pe 16** (θ 0.6, k 8, wC 0.03; m=127): N 0.536, 30 cycles, 1673 WU. MJ with wF=1, wC=0.03, k=16:
  **1087 WU**, so the Markovian cycle wins.
  Reason: with red-black and ω_F=1, every sweep solves the F-points exactly (K=0). MJ with small ω_C is then damped
  Richardson on the red-black Schur complement, which is well-behaved for centered CD. In this family, a small ω_C
  helps N and that Markovian cycle equally.
- **Idea 17, indefinite A = L - sI (+ small convection)** with S~ fitted to the C-parts of the amplified modes
  (`indef_search.py`, dense, m=15): whenever the twin diverges, N diverges too (1.3-1.8). The Galerkin coarse
  correction of strongly indefinite A breaks down. Dropped.
- Running: strong test at Pe 14 (θ 0.5, k 8, wC 0.05).
