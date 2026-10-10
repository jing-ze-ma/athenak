# NOTE gdthin-1011: the full thin switch (fse, hcs / hc) on the vet_gd path

Branch `gdthin-1011` (local, NOT pushed) in worktree /raven/u/jinma/ATHENAK/wt-gdthin, on thinsw-1010 61cab6c0.
Work dir /raven/ptmp/jinma/gdthin_1011.

Commits:
- 24ca5495: C1, eps and J of the vet_gd rays, plus the `m1_thsw` output.
- 257ae190: env `VGD_THIN_EPS`, a diagnostic only. Unset, it is a no-op, bitwise.
- (this NOTE: docs commit on top)

Binaries:
- CPU: bin/athena_cpu_24ca5495 (md5 7f1e745d...) and bin/athena_he_cpu_24ca5495.
- GPU: bin/athena_he_a100_24ca5495 (md5 2a8713ad...) and bin/athena_he_a100_257ae190 (md5 f16e98c4...).

## Verdict

- **`fse` works on vet_gd.**
  - C1 and eps are accumulated in the half-range post-pass of the vet_gd sweep. It is the same kernel that already
    reads every intensity, so there is no extra sweep.
  - Run-to-run GPU bitwise holds. The halo_exact lists, the twin_det reductions and the halo are untouched.
  - The default and `hj` stay bitwise unchanged (AG Car A GPU, BSG GPU, and 6 vet_sc CPU cases).
- **`hcs` maps cleanly and is implemented, but it is inert on both production models.**
  - On every cell above the vet_gd cut, eps = min(kappa_P/kappa_R, 1) = 1 exactly in both AG Car A and BSG.
  - So the scattering weight 1 - eps = 0, and `fse+hcs` = `fse` bitwise.
  - This is the right reading of the code's own closure source. With vet_scatter = ma, S = J + (kappa_P/kappa_R)(B - J)
    has no lagged-J (scattering) part where kappa_P >= kappa_R, and that lag is what hcs's gate is about.
- **`implicit_thin_corr = hc` (the same correction without the 1 - eps factor) is the useful variant on vet_gd.**
  - **AG Car A:** it is stable and run-to-run bitwise, and cheaper than `all` (Picard 4.00 vs 4.33).
  - **Clump shadow test:** Picard 7.7 vs 19.4 for `all`, 0.49 vs 0.71 s/cycle.
  - **Thin-envelope accuracy:** it removes the 2-3 % excess of E over the formal-solution J above the photosphere.
    Median E/J at 1.05-2 R_ph goes from 1.02-1.03 to 1.001.
  - **BSG:** fse+hc is also stable and bitwise across runs (Picard 2.63 vs 2.53 default). Median E/J above the
    photosphere goes from 1.013 to 1.000.
  - **Candidate for production:** `implicit_thin_switch = fse` + `implicit_thin_corr = hc` (with xthin 0). On vet_gd,
    `fse` = `fs` there, because eps = 1.
  - **Caveat:** the thinsw CPU battery found hc unstable on the Cartesian THERMAL beams (xb20/ba0/ba20, vet_sc). The
    vet_gd analogue tested here (a fixed cold, optically thick clump at 1.5 R_ph) is stable for 30 cycles, but this is
    not a beam proof. A vet_sc-style thermal beam cannot be run on vet_gd: vet_gd needs the sp wedge and vet_sc is
    Cartesian only (see "CPU consistency").

## Design: how vet_gd / vet_col data map to the switch

- **h+- and |H|/J** already came from vgd_hr (VetGdHalfRange, src 2 of M1HrCell): the 6 half-range ratios along r,
  theta and phi, twin-corrected.
- **C1 = sum_d w_d |I_d - S| / J** (new, component M1_VGD_HRC1 = 6):
  - **Where:** in the same stage-0 kernel (m1_vgd_hr), one more pass over the n directions of the cell. S is the
    sweep's own cell source exp(vgd_cs(1)); the fused twin uses its own vgd_cst (flag vgd_hr_tw).
  - **Same quantity as vet_sc's M1_VET_NL / J:** sum over the rays through the cell of w |I - S|, with the cell
    source.
  - **Twin:** the twin's C1 pattern is subtracted like r+-: C1 -= (C1_twin - <C1_twin>_shell), clamped >= 0. It uses
    the same deterministic per-shell sum and, under twin_det, the rank-ordered Allgather.
- **eps = min(kappa_P/kappa_T, 1), with 1 where kappa_T = 0** (M1_VGD_HREP = 7):
  - It is cell-local, so it gets no twin correction.
  - The expression is identical to vet_sc's M1_VET_EPS (opac at build time).
- **J** (M1_VGD_HRJ = 8, E units), needed by the corrections:
  - **Twin:** relative correction, J *= <J_twin>/J_twin.
  - **Formal-solution flux:** H_d = (r+_d + r-_d) J, with the twin-corrected ratios.
- **Storage:**
  - vgd_hr / vgd_hrt get 9 components only when `vgd_thx` = (thsw_c1 or any implicit_thin_corr). hj and the default
    keep 6, so memory and kernels are unchanged. The extra memory is +144 MiB per A100 on AG Car A.
  - There are no new halo exchanges. The lateral faces get C1 / eps / J through the existing iw RES/TA/TC/TR halo of
    the thin switch.
- **Implicit side:**
  - M1ThinCellGd (new, src 2) mirrors M1ThinCell: g = max(fb, q_c C1/1.5) with q_c = eps (fse), the fsq gate, or 1;
    hd = H_d; jo = J. fsa and fsj need K and S per cell, which are not stored for vet_gd, so they stay 0.
  - tcorr is now active for src 2 when vgd_hr has the extra components.
  - esc = 1 - eps comes from vgd_hr. M1ThinFace itself is untouched, so the thinfix-1011 merge sees no conflict there.
- **hcs / hc face flux on vet_gd:** DG = q w c (Hf - H+_L - H-_R).
  - Hf = (H_d,L + H_d,R)/2, clamped to [-J_R, J_L]; H+- = r+- J.
  - On x1 faces (r) DG enters the sp row through the existing `rr -= nup DG` (dt A_f/V). On x2/x3 faces it enters
    through `f2_ += (cl/ch) DG`.
  - The geometric factors are the same as for the HCL/HCR part of the half-range flux.
- **vet_col below the cut:** vet_col (the 1-D column solution) has no per-cell J, H+-, so cells below vet_gd's first
  shell stay flagged (r+ = -1), with central faces and w = 0, as before.
  - In both models these cells are optically thick (taucut 30): w = 0 there for every arm anyway.
- **New output `m1_thsw`** (bin, derived; works for vet_sc and vet_gd). Components:
  - m1_tfb = |H|/J;
  - m1_tc1 = C1;
  - m1_teps;
  - m1_tw = smoothstep((g - h0)/(h1 - h0)) for the run's switch;
  - m1_tj = J.
  - The value is -1 where there is no formal solution, or where C1, eps and J are not stored.
- **env `VGD_THIN_EPS=x`** (diagnostic, read once): eps -> min(eps, x). It was used to force the scattering limit and
  exercise hcs.

## Tests

### 1. CPU (vet_sc paths; consistency)
There is no case where both vet_sc and vet_gd run:
- vet_gd needs the 3-D spherical-polar wedge (VgdFatal "needs a 3-D wedge").
- vet_sc is refused on sp by the rad_m1 geometry guard.

So consistency was checked through the identities C1 obeys in both solvers, plus a bitwise check of the vet_sc paths.

**Bitwise, vet_sc** (cpu/verify.sbatch, job 31060454; old = thinsw 9d8fd8c9, new = 24ca5495). Bin, tab and hst data are
identical, and so are the Picard counts:

| case | arms | Picard mean |
|---|---|---|
| pulse_k0.128 hd2 X 128 | hr | 1.997 |
| ba0 (20 cycles) | hr / tEs | 4.65 / 4.00 |
| atm2d G1 | hr / tEs | 2.060 / 2.066 |
| xb20 | tEs | 5.40 |

**vet_sc C1 profile**: 1-D scattering atm, tEs, m1_thsw output.

| tau | C1/(1.5 \|H\|/J) |
|---|---|
| 97 | 1.006 |
| 3.5 - 56 | 1.000 |
| 0.37 | 1.036 |
| 0.11 | 1.079 |
| 5e-4 | 1.154 |

At the top C1 -> 0.997.

**vet_gd C1**:

| run | C1/(1.5 \|H\|/J), median over cells with \|H\|/J < 0.2 |
|---|---|
| CPU AG Car A, 480 x 32 x 32, 4 MPI x 18 OMP, 6 cycles | 0.999 |
| GPU AG Car A | 0.999 |
| BSG | C1/1.5 tracks fb (table below) |

- The diffusion-limit identity C1 = 1.5 |H|/J holds as in vet_sc.
- In the thin envelope C1/(1.5 fb) = 1.3 at 1.05-1.5 R_ph and 1.58 above 1.5 R_ph. C1 > 1 where S > J (the hot-gas
  source in the cool outer envelope), and w = 1 there anyway.

**vet_gd "beam"**: a beam or shadow on vet_gd needs he_star_m1. The analogue used is a fixed cold clump: he_clump at
1.5 R_ph, a = 20 Rsun, rho 1e-9, 8000 K, re-imposed every step. Results in section 3.

### 2. Default and hj bitwise
- **AG Car A GPU (4 A100, 30 cycles):** gdthin `all` = thinsw 9d8fd8c9 `all_1` (j.31055398), bitwise for bin and hst.
  This holds for 24ca5495 (all_1, all_2) and for 257ae190 (all_3).
- **hj:** gdthin `hj` = thinsw `hj` (j.31054015) = thinsw `fse` (which reduced to hj there).
- **BSG, reduced to 256 x 128 x 128, 4 blocks, 20 cycles** (gpu_bsg.sh, bsg.31061077): old 9d8fd8c9 = new, bitwise
  for 8 bin files and hst.
- **vet_sc:** see section 1.

### 3. GPU AG Car A, gpudev, 1 node x 4 A100, 30 cycles

Runs: gpu/runs/j.31060629 (all_1 hj fse_1 fse_2 fses_1 fses_2), j.31060630 (clump arms, all_2),
j.31061305 (all_3 fse_3 fseE0 fsesE0_1/2 cfsesE0), j.31061624 (fseh_1/2 cfseh). Analysis: gpu/ana_gpu.py,
gpu/ana_thsw.py, gpu/shell_prof.py.

| arm | s/cycle | Picard mean/max | NONCONV | inner | peak MiB | run-to-run |
|---|---|---|---|---|---|---|
| all (prod. hr, xthin 30) | 0.424 | 4.333 / 9 | 0 | 3.331 | 23735 | all_1 = all_2 = all_3 |
| hj | 0.423 | 4.433 / 10 | 0 | 3.211 | 23771 | = thinsw hj |
| fse | 0.420-0.422 | 3.900 / 8 | 0 | 3.872 | 23879 | fse_1 = fse_2 = fse_3 |
| fse+hcs | 0.416 | 3.900 / 8 | 0 | 3.872 | 23879 | = fse bitwise (eps = 1) |
| fse+hcs, VGD_THIN_EPS=0 (= hj weight + full hc) | 0.435 | 4.567 / 10 | 0 | 2.730 | 23879 | E0_1 = E0_2 |
| **fse+hc** | 0.418 | **4.000 / 7** | 0 | 3.167 | 23879 | fseh_1 = fseh_2 |
| clump: all | 0.706 | 19.40 / 82 | 0 | 3.381 | 23751 | |
| clump: hj | 0.749 | 17.97 / 23 | 0 | 4.699 | 23771 | |
| clump: fse (= fse+hcs, bitwise) | 0.591 | 12.60 / 23 | 0 | 2.775 | 23879 | |
| clump: fse+hcs, VGD_THIN_EPS=0 | 0.685 | 15.03 / 21 | 0 | 4.803 | 23879 | |
| **clump: fse+hc** | **0.492** | **7.67 / 10** | 0 | 3.170 | 23879 | |

Other checks:
- `fseE0` (fse with eps capped at 0, no correction) = `hj` bitwise in hydro, m1 and hst. Only m1_thsw differs. This
  confirms that the weight reduces exactly.
- fsesE0 differs from hj by the hc correction alone.

**w map** (m1_thsw at cycle 30; vet_gd domain = r >= 0.977 R_ph, i >= 248, cut at tau 30):
- Below 1.000 R_ph, w = 0 for every arm (fb <= 0.27). Above 1.009 R_ph, w = 1 for every arm.
- `fse` moves the switch-on 1-2 shells inward. The C1/1.5 term exceeds |H|/J in 77 % of the 0.95-1.05 cells and in
  all cells above:

| r/R_ph (shell) | hj mean w | fse mean w | fb median | C1/1.5 median |
|---|---|---|---|---|
| 1.001 (269) | 0.09 | 0.16 | 0.36 | 0.38 |
| 1.002 (270) | 0.44 | 0.62 | 0.44 | 0.47 |
| 1.004 (271) | 0.72 | 0.92 | 0.50 | 0.55 |
| 1.005 (272) | 0.86 | 1.00 | 0.53 | 0.60 |
| 1.006 (273) | 0.94 | 1.00 | 0.55 | 0.63 |

- Over the band 0.95-1.05 R_ph: w > 0 in 0.661 of the cells for both; w = 1 in 0.532 (hj) vs 0.597 (fse).

**Differences vs `all`** (max over (theta, phi) per shell of |arm - all| / max|all|, max over the band):

| r/R_ph | hj m1_e | fse m1_e | fse+hc m1_e | hj velx | fse velx | fse+hc velx |
|---|---|---|---|---|---|---|
| 0.8-0.95 | 3.0e-5 | 3.8e-5 | 1.2e-4 | 9.9e-4 | 1.2e-3 | 4.1e-3 |
| 0.95-1.05 | 6.1e-3 | 3.1e-2 | 4.7e-2 | 0.12 | 0.86 | 1.3 |
| > 1.05 | 1.1e-4 | 1.4e-4 | 2.0e-2 | 1.3e-2 | 1.3e-2 | 1.8e-2 |

- The fse E difference peaks at 1.000 R_ph: 3.1e-2, vs hj 2.8e-3.
- The velx differences sit in a thin layer at 0.995-1.005 R_ph (shell 266: 0.86 x the shell max of |velx_all|), which
  is small at t = 3.8e4 s. They come from the earlier switch-on, not from noise: fse_1 = fse_2 bitwise.
- Above 1.05 R_ph, fse = hj (w = 1 for both). The 2e-2 of fse+hc there is the correction, which moves E to J:

| median E/J (p5-p95) | 1.00-1.05 | 1.05-1.5 | 1.5-2 | 2-3 R_ph |
|---|---|---|---|---|
| fse (no correction) | 1.033 | 1.031 (1.026-1.033) | 1.022 | 1.025 |
| fse+hc | 1.003 | 1.001 (1.000-1.002) | 1.002 | 1.010 |
| clump fse | 1.033 | 1.031 | 1.020 (0.98-1.07) | 1.014 (0.96-1.12) |
| clump fse+hc | 1.003 | 1.002 | 1.001 (0.97-1.03) | 0.999 (0.97-1.07) |

The excess of the first-order upwind half-range flux over the formal solution in the thin envelope (2-3 %) is what hc
removes. `all`/hj have the same excess: their E above 1.05 R_ph equals fse's to 1.4e-4.

### 4. BSG (reduced lateral resolution, 4 A100, 40 cycles; bsg.31061178, bsg.31061657)

| arm | s/cycle | Picard mean / max | inner | notes |
|---|---|---|---|---|
| default (hr) | 0.234 | 2.525 / 5 | 15.8 | |
| hj | 0.231 | 2.750 / 6 | 14.7 | max\|dE\|/E vs default 4.1e-3 (52-55 Rsun) |
| fse | 0.226 | 2.575 / 5 | 14.4 | 6.3e-3; fse+hcs = fse bitwise; fses_1 = fses_2 |
| fse+hc | 0.224 | 2.625 / 5 | 14.5 | fseh_1 = fseh_2; max\|dE\|/E vs default 2.1e-2 (52-55 Rsun), 5e-3 above 55 Rsun |

- eps = 1 in every vet_gd cell (r >= 51.2 Rsun).
- w > 0 in 46 % of the 50-55 Rsun cells and in all cells above.
- w = 1 at 50-55 Rsun: 4.2 % (hj) vs 20.8 % (fse).
- C1/1.5 medians 0.24 / 0.86 / 1.02 / 1.11 at 50-55 / 55-60 / 60-70 / 70-80 Rsun; fb medians 0.24 / 0.70 / 0.80 / 0.86.
- Median E/J at 55-60 / 60-70 / 70-80 Rsun: fse 1.0135 / 1.0133 / 1.0146; fse+hc 1.0003 / 0.9995 / 1.0000. The
  same thin-envelope excess is removed, as on AG Car.

## Caveats / open
- hcs's 1 - eps gate is identically 0 on the production tables (kappa_P >= kappa_R above the vet_gd cut). If a
  scattering-fraction gate is wanted for vet_gd, the natural alternative is the electron-scattering part
  kappa_e/kappa_R (opac M1_OP_S) or the ma-form weight 1 - kappa_P/kappa_R (<= 0 here). Neither is implemented.
  fse+hc (no eps gate) is the tested option.
- hc has no proof yet on a vet_gd thermal BEAM. The cold-clump shadow is the only obstacle test. The vet_sc thermal
  beams (xb20/ba0/ba20) diverged with fse+hc in the thinsw battery. Before production, run a longer AG Car A fse+hc
  segment, for example a 200-cycle restart chain.
- With fse the switch-on moves about 0.003 R_ph inward. The photospheric E and velx change by 3 % and 0.86 x (thin
  layer) vs `all` after 30 cycles, deterministically. Whether this is more accurate there was not measured: there is
  no reference.
- fsa / fsj (K/J anisotropy, |J - S|) are 0 on vet_gd. They are not stored, and their thinsw defaults are 0 anyway.
- `VGD_THIN_EPS` is read with getenv once per process (static). It is a diagnostic only.
- The CPU AG Car small run (480 x 32 x 32) used the he CPU binary (bin/athena_he_cpu_24ca5495): fse and fses, 6 cycles,
  Picard 6.2, w-map and identities as on GPU.

## Files
- Code: src/rad_m1/rad_m1_vetgd.cpp (VetGdHalfRange, VetGdTwinFusedMoments), src/rad_m1/rad_m1_implicit.cpp
  (M1ThinCellGd, the two face kernels, vgd_thx), src/rad_m1/rad_m1.hpp (M1_VGD_HR*, vgd_thx, vgd_hr_tw),
  src/outputs/* (m1_thsw).
- Scripts: build_{cpu,he,gpu}.sh, cpu/verify.sbatch, cpu/agc.sbatch, gpu/gpu_arms_raven.sh (ARMS all hj fs fse fses
  fseh call chj cfse cfses cfseh fsesE0 fseE0 cfsesE0), gpu/gpu_bsg.sh (old new newt hj fse fses fseh),
  gpu/inp/agcar_A_gdt(.._clump).athinput, gpu/inp/bsg_gdt(_t).athinput, bitcmp.sh.
