# shdfb: cell layout (d6b6ed4e2) Newton/Picard "divergence" on the shadow deck with gas coupling

Deck: `shdfb_cell.athinput` = cc/decks shd_FIN3_r12g2 at x0.5 (256 x 128, one MeshBlock) with
`gas_feedback = true` (only change) and `scheme = clean; clean_layout = cell` (all clean keys at their
defaults; tlim 1.027961e-3).
Binary: d6b6ed4e2, /orion/ptmp/jinma/clean_1115/gpufix_reg/build_new/src/athena (md5 e0899b8c).
Run: `OMP_NUM_THREADS=16 OMP_PROC_BIND=false OMPI_MCA_hwloc_base_binding_policy=none athena -i shdfb_cell.athinput`
(about 11 min for cycle 0 on 16 threads; the run finishes, 2 CLEAN WARNING, ~27 min).

## Symptom (log_excerpt_d6b6ed4e2.txt)
cycle 0, backward Euler (t2s 1), first attempt: `it -7 (sum 1566) pic 12 1.46e+295 ... neg 8 nreal 18 solves 250`
-> CLEAN WARNING, dt_retry; second attempt the same (pic 1.55e+293); third accepted with neg 96 nreal 1722.
Every later cycle prints pic ~1e287-1e291 and needs 180-250 solves per step.

## First bad quantity and location
`de`, the relative E change of CleanNewtonUpdate (rad_m1_clean.cpp, the cell layout calls it from
CleanCellNewtonUpdate):
    const Real de = fabs(e - e_old)/fmax(fmax(e, e_old), 1.0e-300);
For a cell whose iterate and previous iterate are both NEGATIVE (the cell layout keeps X < 0 for its positivity
limiter) the denominator is 1e-300, so de ~ 1e299: the Picard/Newton change never settles and the loop runs to
clean_maxnl (with clean_maxpic Picard passes) every positivity round: 250 solves, not converged.
Located with a debug build (one printf at that line, same cell-layout code as d6b6ed4e2: its log is identical,
pic 1.46e+295); first hits, cycle 0 stage t2s 1 (dbg_first40.txt, dbg/log.txt.gz):
    m 0 k 0 j 74 i 150  e -4.04e-03  e_old -1.65e-01  de 1.6e+299  T 0.267
    m 0 k 0 j 75 i 147  e -1.51e-03  e_old -3.12e-01  de 3.1e+299  T 0.196
    ... always the same ~16 cells, (i 145-150, j 55-56 and 74-76) = the rim of the absorber (shadow_ax 0.1,
    ay 0.06 at the origin), mirror pairs about the x axis, ~695 hits each.
Two things in it: (1) the denominator (a fabs, or max(|e|, |e_old|), would make de the true relative change);
(2) underneath, the gas-coupled cell layout drives E strongly negative (-0.3 against E ~ 1) at the absorber rim
(neg 8-96 cells after the chain); without gas coupling the same deck (shd x0.5) also ends stages with
neg > 0 (8 of 8 stages), but CleanCellNewtonUpdate calls CleanNewtonUpdate (and so this de) only with gas coupling.
Nothing was changed in the cell layout (viper's scheme).
