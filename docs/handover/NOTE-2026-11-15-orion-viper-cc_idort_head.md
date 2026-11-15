# NOTE 2026-11-15 Orion -> viper: handover acknowledged; cc_idort_head results; cc_fixedD rerun; restyle_m3 skipped

## Acknowledgement
Read HANDOVER-2026-11-15-clean-scheme-viper-to-orion.md, the package handover_1115_clean and NOTEs 9ad8fb5ed / f6f1bbdb6.
Orion now owns fork/clean-cc (head dc67dba7b). restyle_m3 (53ecd50fb) was NOT run: the user skipped it after f6f1bbdb6
(restyle checks on Raven).

## cc_idort_head (dc67dba7b, arm CC1 = clean, cell, extrap, maxpic 4, pc_sweeps 2)
Builds: build (MPI, serial) md5 900f957959d025b87b2a3fac7dad54da; build_omp (Athena_ENABLE_OPENMP) md5 5a934c3092e52dcb35fff23d70bd3c4b.
gcc 13.1.0, openmpi 4.1.8, cmake 3.30, -ffp-contract=off. Runs used the OMP binary (6-100 threads per run) for speed;
thread count changes results only at solver tolerance (cb x2 at 28 vs 64 threads: BiCGStab sums 669/664, score
identical to 16 digits).

### score.txt
```
cb_FIN3_r4       x0.5 1.686e+00  x1 1.875e+00  x2 1.812e+00 | order -0.15  0.05
cbg_FIN3_g8x32   x0.5 1.624e+00  x1 1.631e+00  x2 1.639e+00 | order -0.01  -0.01
hcyl_FIN3_g8x32  x0.5 4.684e-02  x1 2.088e-01 | order -2.16
nltek1_FIN3_d    x0.5 9.335e-05  x1 9.116e-05  x2 9.073e-05 | order 0.03  0.01
nltek4_FIN3_d    x0.5 1.330e-03  x1 1.327e-03  x2 1.327e-03 | order 0.00  0.00
shd_FIN3_r12g2   x0.25 2.643e-01 | order 
shdg_FIN3_g8x32  x0.25 2.558e-01  x0.5 2.652e-01 | order -0.05
ws66             x0.5 1.345e+00  x1 1.447e+00  x2 1.415e+00 | order -0.11  0.03
```

### Realizability / fallback per run (from the CLEAN lines; "BE" = CLEAN lines with all 1)
```
cb_FIN3_r4/CC1_x0.5      lines   13  BE   10  max_nreal   3148  max_neg   12  dt_retry   0  wall    150 s
cb_FIN3_r4/CC1_x1        lines   24  BE   20  max_nreal  13394  max_neg   28  dt_retry   0  wall    494 s
cb_FIN3_r4/CC1_x2        lines   46  BE   40  max_nreal  55156  max_neg  158  dt_retry   0  wall   2032 s
cbg_FIN3_g8x32/CC1_x0.5  lines   13  BE   10  max_nreal   3126  max_neg   10  dt_retry   0  wall    163 s
cbg_FIN3_g8x32/CC1_x1    lines   24  BE   20  max_nreal  13266  max_neg   28  dt_retry   0  wall    566 s
cbg_FIN3_g8x32/CC1_x2    lines   46  BE   40  max_nreal  53778  max_neg  950  dt_retry   0  wall   3181 s
hcyl_FIN3_g8x32/CC1_x0.5 lines   64  BE    0  max_nreal    156  max_neg    0  dt_retry   0  wall    148 s
hcyl_FIN3_g8x32/CC1_x1   lines   74  BE   62  max_nreal  21236  max_neg    0  dt_retry   0  wall   1653 s
nltek1_FIN3_d/CC1_x0.5   lines   59  BE    0  max_nreal      0  max_neg    0  dt_retry   0  wall     24 s
nltek1_FIN3_d/CC1_x1     lines  119  BE    0  max_nreal      0  max_neg    0  dt_retry   0  wall     70 s
nltek1_FIN3_d/CC1_x2     lines  239  BE    0  max_nreal      0  max_neg    0  dt_retry   0  wall    191 s
nltek4_FIN3_d/CC1_x0.5   lines   59  BE    0  max_nreal      0  max_neg    0  dt_retry   0  wall     59 s
nltek4_FIN3_d/CC1_x1     lines  119  BE    0  max_nreal      0  max_neg    0  dt_retry   0  wall    137 s
nltek4_FIN3_d/CC1_x2     lines  239  BE    0  max_nreal      0  max_neg    0  dt_retry   0  wall    273 s
shd_FIN3_r12g2/CC1_x0.25 lines    4  BE    3  max_nreal   6442  max_neg   30  dt_retry   0  wall    139 s
shdg_FIN3_g8x32/CC1_x0.25 lines    4  BE    3  max_nreal   6600  max_neg   74  dt_retry   0  wall    144 s
shdg_FIN3_g8x32/CC1_x0.5 lines    8  BE    6  max_nreal  27060  max_neg  178  dt_retry   0  wall   1872 s
ws66/CC1_x0.5            lines   13  BE   10  max_nreal   3074  max_neg   16  dt_retry   0  wall    156 s
ws66/CC1_x1              lines   24  BE   20  max_nreal  13160  max_neg   52  dt_retry   0  wall    504 s
ws66/CC1_x2              lines  113  BE   99  max_nreal  53100  max_neg  464  dt_retry   7  wall   4256 s
```

### Failed / not finished
- shd_FIN3_r12g2 x0.5 and x1, shdg_FIN3_g8x32 x1: the backward-Euler solve does not converge at cycle 1
  ("CLEAN WARNING cycle 1: backward-Euler solve not converged"), dt_retry halves dt 9 times (2e-4 -> 4e-7), dt then
  recovers ~10 %/step: shd x0.5 reached t = 4.4e-5 of 1.03e-3 after 50 min; x1 still in cycle 1. Stopped (user).
  x0.25 (your only previous resolution) is fine.
- hcyl_FIN3_g8x32 x2: still running (t = 0.32 of 0.8 after 75 min, 28 threads); result in a follow-up NOTE.

### Findings
1. hcyl degrades with resolution: x0.5 0.0468 (no BE step, max nreal 156) but x1 0.209 with BE on 62/74 steps and up to
   21236 non-realizable cells -> the realizability problem (open problem 1) reaches hcyl at x1 and drives open problem 2.
2. Beam decks: every high-order stage non-realizable at every resolution (max nreal 3.1e3 / 1.3e4 / 5.4e4 for x0.5/1/2),
   BE on ~85 % of the CLEAN lines; scores do not converge (cb 1.69/1.88/1.81, cbg 1.62/1.63/1.64, ws66 1.35/1.45/1.42; ws66 x2 also needed 7 dt_retry).
3. nltek1 vs exact 9.3e-5 / 9.1e-5 / 9.1e-5, nltek4 1.33e-3 (all three resolutions), 0 non-realizable cells.
4. cb x0.5 scores 1.686 here vs your 1.677 (CC17); cbg/ws66/hcyl/shd/shdg x0.5 or x0.25 match your numbers to 3-4 digits.

## cc_fixedD rerun on head (as asked, part of our own work)
Head dc67dba7b, staggered decks with clean_tol_nl 1.0e-12 (78 decks). 152/156 runs rc = 0 (0 failures among finished runs). The 4 staggered N4096 rows (P 3.3, tau 10, nt 512..4096)
are still running with ~1230 backward-Euler non-convergences each (dt_retry; t = 0.054 of 0.46): staggered needs more than
clean_tol_nl 1e-12 at N4096. The 4 cell N4096 rows finished clean. Tables when complete (not reported separately).
Also: extrap time order on the head confirmed (rw_cc19 lists): cell P 1/3 and 3.3 order 2.00, staggered (tol_nl 1e-12) 2.00;
space 1.8-2.4; energy drift round-off.

## Work in progress on Orion (clean-cc, local branches, not pushed yet)
- clean_closure = extrap1 (clean-cc-one 9df83267d): ONE formal solution per step, order 2.00/2.00 in time (P 1/3, 3.3), nltek1 x1
  9.12e-5; stage-1 first-iterate rays feed stage 1, stage 2 and the extrapolation history. IDORT-suite check running.
- Beam realizability: diagnosis = the low-order D = I flux is not realizability preserving in 2-D (|f +- D n| <= 1 +- f.n
  holds only for f_t = 0), and the moment solver's steady state differs from the rays' moments by ~6 % (2-D prototype).
  Consistency term F_ray,face - F_scheme(J, H) on the face fluxes (clean-cc-holo, v2): shd x0.25 0.264 -> 0.021, shdg 0.256 -> 0.030,
  cb x0.5 1.69 -> 0.55, hcyl x0.5 0.047 -> 0.027, but nltek1 x1 9.1e-5 -> 6.4e-4 and rw 3x worse, realizability not fixed (SC ray
  fluxes are not divergence-free per cell). Next round running: divergence-consistent ray fluxes, no bracket on low-order faces,
  C1-weighted bracket.
