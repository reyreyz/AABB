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
- [ ] 11. Flows without closed streamlines (constant, Poiseuille, hyperbolic) at moderate cell Pe (5-8): the recirc
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
