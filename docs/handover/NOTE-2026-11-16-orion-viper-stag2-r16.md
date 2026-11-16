# NOTE 2026-11-16 Orion -> viper: stag2 after the fixes (rounds 15-16): recommended set reinstated, consolidated table

Follows the correction c20b2fa85. Local branch clean-cc-stag2 head c8afa65ec (on d6b6ed4e2): keys off = reg bitwise (job 230704).
Full table: Orion /orion/ptmp/jinma/clean_1115/stag2/cxx/RESULTS16.md (bench: stag2/bench/BENCH.md, r15, r16).

## Fixes since the correction
- **F-solve dissipation:** only in cells with an R2 face. This makes rw second order again.
- **h-free steady forms:** for the bracket Phi(J) and the energy-flux memory d. This gives well-balancing.
- **d in restart files:** restarts of the full set are bitwise (hcyl 2x4, sp3 2 blocks, flz q).
- **Key coupling:** fkin needs consistency + fwork.
- **Gas-coupled write-back:** E takes the row residual divided by its diagonal. This corrects 753059c16 in absorbing cells.

## Recommended set
clean_layout = staggered, clean_stag_split (default), clean_stag_fsolve, clean_stag_fkin, clean_stag_fwork, clean_consistency;
clean_rays = fv1 on sp.

## Full set vs cell (d6b6ed4e2), 0 CLEAN WARNING everywhere
- **Equal or better than cell:**

  | test | staggered | cell |
  |---|---|---|
  | rw N256 L1, P 1/3, tau 0.1 / 10 / 1000 | 2.39e-4 / 5.89e-5 / 1.43e-4 | 3.08e-4 / 8.15e-5 / 1.44e-4 |
  | rw N256 L1, P 3.3, tau 0.1 / 10 / 1000 | 6.09e-4 / 5.36e-5 / 2.62e-4 | 6.10e-4 / 6.62e-5 / 3.00e-4 |
  | rw dE_tot | <= 3.6e-16 | <= 3.2e-12 |
  | rhse_sp2 WB, max v/c_s | 5.6e-9 | 1.7e-7 |
  | nltek1 | 9.118e-5 | 9.120e-5 |
  | thin_r F | 0.109 | 0.142 |
  | sp3 1 vs 2 blocks | 1.3e-9 | 2.4e-9 |
  | flz q | 25 cycles | 41 cycles |

  On flz q, dE and dLz are ~1e-15 and the wall v_r is inward, like cell. Final E stays >= 0 on cb/shd/shdg, where cell goes negative.
- **Cost per cycle:** 1.2-4.6x cheaper on Cartesian, 1.25-1.6x on sp (bench timings).
- **Worse than cell:**
  - Beams with the positivity chain: cb/cbg/ws66 1.448/1.279/1.146 vs 0.46/1.10/0.921, because the chain strips the bracket.
  - shd umbra 27.8 vs 18.9.
  - hcyl 0.0129 vs 0.0121.
  - pulse0 closed-box energy 2.8e-14 vs 8.9e-16. This is the solve residual with no gas; the stopping rule is unchanged.
- **With the chain off:** beams and shadows equal cell exactly, with cell-like negatives, but flz q then stalls. Orion is looking at that next.
- **Cubed sphere:** part 1 only; paused for the USER.
