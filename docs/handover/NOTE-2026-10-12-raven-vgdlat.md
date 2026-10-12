# NOTE vgdlat-1012: the lateral error of the vet_gd closure (T6a / T6b / P1 thin envelope) and `vet_gd_symref = col`

Worker V, 2026-10-12. Branch `vgdlat-1012` (local, not pushed), worktree /raven/u/jinma/ATHENAK/wt-vgdlat, from
prodfinal-1011 0197f020. Work dir /raven/ptmp/jinma/vgdlat_1012.

Commits:
- 550e9673: `vet_gd_symref = col` (opt-in).
- bf50c2a6: a per-cell weight of the symmetric shift (needed for P2).

Binaries in bin/:
- athena_cpu_sr4 and athena_he_a100_sr4 are bf50c2a6.
- sr3 is 550e9673.
- athena_cpu_0197f020 is a copy of the prodfinal binary.

(The worktree's `kokkos` is a symlink to wt-prodfinal/kokkos for building. It is not committed.)

## 1. Diagnosis

Harness: a copy of sptestsC ana/.
- `ana/diag.py` compares the vet_gd formal-solution quantities with a lateral-exact 1-D reference of the same discretisation.
- `ana/symref.py` is a Python p-ray solve with the SC's segment rules. It is exact on T6a to 1e-5.
- New arm suffixes: W (twin_full), Y1 (symref), K2 (conservative2, shift off), Rn (4x finer symref rays), C0 (thin_corr none), Sm (vet_seg mean).
- Env SPT_DIAG2=1 adds the m1_fs, m1_thsw, m1_vet and m1_gdhr tab outputs.

**T6a, thin sphere, 16x16, base (production fse+hc + hesdirk2 keys).** Band r 5-10, N 64 / 128 / 256:

| quantity (vs the lateral-exact reference) | base |
|---|---|
| radial half-range h+ = H+/J | -1.1 % |
| D_rr | -2.1 / -2.3 / -2.4 % |
| twin-corrected J (hc's J) | +41 / +48 / +51 % |
| lateral h+_theta | +91 / +95 / +96 % |
| E | +0.2 / +0.7 / +1.0 % |

- In the thin band the flux is the half-range flux (w = 1), so E is about F/(c h+). The E error is minus the h+ error.
- **Cause:** the twin subtraction keeps the twin's SHELL MEAN as the symmetric value. That mean carries two errors of the short-characteristics sweep:
  1. **Bilinear lateral interpolation of the fixed-direction intensities at every shell.** This is numerical diffusion of first order in the lateral cell. It saturates as N_r grows, so the error does not converge with N_r.
  2. **Direction sampling (120 directions) of a narrow cone.**
- Both shrink with more lateral cells: hp -1.5 % (L8), -1.1 % (L16), -0.9 % (L32). They also shrink with more directions: Q (432 directions) gives -0.6 %.
- rotate_every 0 changes nothing.
- vet_gd_twin_full is irrelevant. E is unchanged with it, because the thin band does not use D_rr, and vet_col's own f_K is wrong for an emitting inner sphere (0.35).

**T6b (IDORT sphere)** is mostly NOT a lateral problem.
- Its +4.2 % exterior E (16x16, N 256) is +2.95 % of M1 core luminosity plus about +1.2 % of closure error. The closure part is the same lateral effect.
- The luminosity is lateral-independent: L32 = L16 = 4x4 to 1e-4.

## 2. Fix: `rad_m1/vet_gd_symref = col` (opt-in; read only when named)

- **What it computes.** Every vet_gd build solves a 1-D spherical p-ray formal solution of the twin's shell-mean chi and S on the host (`VetGdSymRef` in rad_m1_vetgd.cpp).
  - It uses the SC sweep's own segment rules term for term: VetSegSplit at the face, the inner-face start (diffusion / emit), the turned ray, and the top segment and bath.
  - It differs from the SC sweep in two ways: there is no lateral interpolation, and the moments are integrated exactly over the piecewise-linear-in-mu intensity.
  - It produces f, h+-_r, h_lat, C1 and J per shell.
- **Where it is applied.** It replaces the twin's shell mean as the symmetric reference in VetGdTwin(1) (LAT0) and in VetGdHalfRange(2) (h+-, C1, J). The result is main - twin + ref.
- **Per-cell weight (bf50c2a6).** The shift (ref - <twin>) is weighted per cell by w = min(J_twin/J_main, 1).
  - w = 1 on a laterally uniform state.
  - w goes to 0 in a channel or beam much brighter than the mean field, where the mean's bias says nothing.
  - J is shifted by the factor (J_ref/<J_tw>)^w.
  - Without the weight, P2 went from 1.5 to 2.2x.
- **Rays.** vet_gd_symref_ncore (default 64) and vet_gd_symref_nsub (default 4). A 4x finer set changes P1 by <= 1e-4.
- **Cost.** The solve runs on a host thread (std::async), started at the twin's stage 0 and joined at stage 1, so it overlaps the sweep. Env VGD_SR_SYNC=1 runs it inline, with bitwise identical numbers.
- **Correctness across ranks.** The shell means are the twin's gathered ones (rank-ordered under twin_det), so the reference is identical on every rank. It holds no state, so restarts are exact.
- **Refusals.** It is refused with vet_gd_beam and with twin_full, and it needs vet_gd_twin.
- **Diagnostic.** VGD_SR_DUMP=<file> writes the reference of the latest build.
- **Default** (key unnamed): bitwise.
  - CPU tab: thin_L16 HFfM, 0197f020 vs sr1.
  - GPU AG Car A hst: 0197f020 = sr3 = sr4.

## 3. Results (CPU, production keys hesdirk2/predict/vimp/vstage, xthin 0, fse, hc, vet_bc_diff face, vet_seg emit; K2 = conservative2 + cs2_shift false where gas couples)

| test | base (0197f020 arm HFfM) | symref (Y1) |
|---|---|---|
| **T6a** E L1, 16x16, N 64/128/256 | 1.8e-3 / 6.2e-3 / 9.5e-3 (grows) | **8.1e-4 / 2.0e-4 / 4.8e-5, order 2.01 / 2.08** |
| T6a, 4x4 / 8x8 / 32x32 | 4x4 3.2 % at N 256, grows | identical to 16x16 (lateral-independent) |
| T6a max\|F r^2 - 1\| | 1.9e-2 / 5.0e-3 / 1.3e-3 | same (2nd order) |
| **T6b** exterior r>3, 16x16, N 64/128/256/512 | +4.9 / +0.8 / +4.2 / +5.4 % | +5.3 / +0.3 / +3.3 / +4.1 % (4x4 at 512) |
| T6b F r^2 / exact (luminosity) | 1.045 / 0.9995 / 1.030 / 1.045 | 1.061 / 1.007 / 1.034 / 1.045 |
| T6b closure part (E/ex minus L/ex), N 256 | +1.2 % | **-0.1 %** |
| **P1** patm L16 poly, E/J-1 at R_ph, N 64/128/256/512 | -0.47 / -0.19 / -0.014 / **+0.091 %** (crosses zero) | -0.47 / -0.23 / -0.10 / -0.038 % (monotone, order 1.0 / 1.2 / 1.4) |
| P1 thin envelope env (signed mean r > 1.2 R_ph) | +0.45 / +0.88 / +1.14 / +1.34 % (grows) | **-0.23 / -0.07 / -0.02 / -0.006 %** |
| P1 band tau<0.1, L1 | 5.7e-3 / 9.2e-3 / 1.2e-2 / 1.4e-2 | 2.1e-3 / 6.3e-4 / 1.9e-4 / 5.4e-5 (order 1.75) |
| P1 band 0.1-1 / 1-10, N 256 | 2.9e-3 / 6.9e-4 | 7.5e-4 / 7.3e-4 |
| P1 max\|F r^2 - 1\| | 3.8e-4 (N 256) | same |
| P1 Picard | 3.0 | 2.0 |
| L pulse r0/w 1000 vs Cartesian vet_sc (N 128) | 1.4e-3 (the "solver floor") | **3.8e-4** |
| L pulse r0/w 1, X self-convergence of E (64-128-256) | order 0.27 | order 1.23 |
| T4 dif r0 = 3 | 7.2e-3 / 1.8e-3 / 4.6e-4, order 2.0 | identical |
| T7 shadow 48x48 | umbra 0.049, lit L1 0.39 / 0.51 | identical (sr4); sr3 unweighted umbra 0.044 |
| rsw (K2) X order dens / v / E | 1.83 / 1.99 / 2.04 | identical (diff 1e-11) |
| P2 pinc, L/L_ex r 4.3/4.6/4.9 | 1.50 / 1.53 / 1.56, 1 NONCONV, Picard 6.6 | 1.51 / 1.55 / 1.58, 0 NONCONV, Picard 5.0 (sr3 unweighted: 2.2) |

- **T6a is fixed.** E is second order and does not depend on the lateral grid.
- **P1 R_ph.** The base's "second order" at R_ph was a cancellation of the first-order radial FS error against the growing lateral envelope excess. At N 512 the base overshoots, to +0.09 %.
  - With symref, R_ph converges monotonically at order 1.0-1.4.
  - The FS itself at R_ph (J_fs -1.6 / -0.9 / -0.4 %, D_rr +1.2 / +0.6 / +0.2 %) is the same in both arms. It is the radial first order of the SC rules, not lateral.
- **P2 is not this cause.** P2 is a strongly asymmetric channel, and symref leaves it as it was.

## 4. What remains of T6b (the coordinator's question)

- After the fix the T6b exterior error is almost entirely the M1 luminosity of the opaque core: +3.4 % at N 256 and +4.5 % at N 512, not converging.
- That luminosity exceeds the black-body limit of the code's own core radius (face 0.997).
- It is lateral-independent: 4x4 = 16x16 to 1e-4, and symref changes it by < 0.5 %.
- It depends on hc: without hc (C0) it is +10 %. vet_seg mean gives +4.8 %.
- The FS's own luminosity at that edge is -10 to -15 %. This is the SC's cell-centred source at an opaque -> vacuum face (the vet_seg emit limiter keeps the face S at 0.2-0.4 of the core S).
- **Verdict.** It is a separate RADIAL issue: the M1 half-range/hc flux and the SC source at a resolved-to-unresolved opaque surface (tau_cell 0.13-0.5). It is worth a follow-up with a 1-D (4x4) core ladder, the face-flux budget of the surface cell, and the hc correction's Hf there. Fraction remaining: closure part about 0 (from +1.2 %); the luminosity part, +3.4 % at N 256, is unchanged.

## 5. GPU cost and determinism (full AG Car A, prodfinal FRESH deck, 1 node x 4 A100, gpudev, 30 cycles)

gpu/gpu_t.sh. Jobs 31089262 (sr3) and 31089949 (sr4).

| arm | s/cycle (cycles 3-29) | Picard | hst |
|---|---|---|---|
| 0197f020 | 0.3118 | 2.761 | = sr3/sr4 default (bitwise) |
| sr3 default, rep 2 | 0.3107 / 0.3112 | 2.761 | bitwise |
| sr4 default | 0.3121 | 2.761 | bitwise to sr3 |
| symref async, 2 reps (sr3) | 0.3152 / 0.3155 (+1.4 %) | 2.818 | rep1 = rep2 bitwise |
| symref async, 2 reps (sr4) | 0.3156 / 0.3157 (+1.1 %) | 2.818 | rep1 = rep2 bitwise |
| symref inline (VGD_SR_SYNC=1) | 0.3362 (+8 %) | 2.818 | = async bitwise |

- Host cost:
  - Inline, the solve is 0.85 s over the run (about 13 ms per build at 480 shells).
  - The async start, including the copies of the two shell means, is 3.6 ms in total.
  - The join, with an upload of 6 c1 values, is 5.5 ms.
- About half of the +1.1 % comes from Picard (2.76 -> 2.82), which is a physics change, not overhead.
- **Restart-exact (CPU):** thin_L16 and core_L16 at N 64, 1 and 2 MPI ranks (halo_compact 2, twin_det), 40 cycles continuous vs 20 + restart 20. All of m1, m1_fs, m1_thsw, m1_vet and m1_gdhr are identical. 1 vs 2 ranks: <= 3e-12.

## 6. Recommendation / open

- `vet_gd_symref = col` fixes the "known lateral vet_gd error" for the laterally symmetric part.
  - That part is what dominates the thin envelope of AG Car / BSG: +0.3 % at 0.008 rad in P1b.
  - It costs about 1 % on A100.
- Before any production default:
  - an AG Car A restart segment of a few hundred cycles, comparing E/J above the photosphere and the light curve;
  - a BSG gate.
  This is not done here.
- Open (radial, not lateral):
  - T6b opaque-surface luminosity (section 4);
  - the first-order FS at R_ph (SC segment rules);
  - P2 (asymmetric channel; moment step, as in sptests P2b).
