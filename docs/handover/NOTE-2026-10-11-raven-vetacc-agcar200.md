# NOTE vetacc-1011: merge check and the 200-cycle AG Car A time-accuracy and cost comparison

Branch `vetacc-1011` 42aa6bbc (local, NOT pushed), worktree /raven/u/jinma/ATHENAK/wt-vetacc.
It is thinsw-1010 61cab6c0 + thinfix-1011 fe910d82 + gdthin-1011 a0f9300e (TR-slot conflict resolved by hand) + vcpred-1010 b6c06786.
Work dir: /raven/ptmp/jinma/vetacc_1011. No commit was needed (the merge is correct as it is).

Binaries (bin/), built by build.sh (job 31063532) with gcc/13, cuda/12.6, openmpi_gpu/5.0 and cmake/3.30:
- athena_cpu_42aa6bbc, md5 eb288f23...
- athena_he_a100_42aa6bbc, md5 cb50f5ea...

`git submodule update --init kokkos` was needed in wt-vetacc.

## 1. Merge check: every branch's feature is bitwise against that branch's own binary

| check | reference | result |
|---|---|---|
| thinfix tJe, atm 128 and 256 (tab + hst) | athena_cpu_cfb3258f, same input | bitwise (4 tab + hst each) |
| thinfix tJe, ba0 and xb20 (20 cycles, bin) | cfb3258f | bitwise (6 bin each) |
| thinfix tJe, pulse_k0.128 X32, be and hesdirk2 (tab) | cfb3258f | bitwise |
| default / tEs on vet_sc: ba0 hr, ba0 tEs, G1 atm2d hr, pulse hr hd2 X128 | base 9d8fd8c9 (= 61cab6c0 code) | bitwise |
| gdthin fse+hc, AG Car A GPU, 30 cycles (hydro, m1, m1_thsw bin + 2 hst) | gdthin j.31061624/fseh_1 (bin 257ae190) | bitwise |
| vcpred hesdirk2 + time2_vet_col=predict + vimp + vstage, AG Car A GPU, 30 cycles | vcpred j.31060000/hp | bitwise |
| production `all` (be), AG Car A GPU, 30 cycles, ths input | thinsw j.31055398/all_1 | bitwise |
| production be, vc input | vcpred j.31060000/be | bitwise |

Where the scripts and runs are:
- CPU: cpu/verify.sbatch (job 31063668) and cpu/verify/*.
- GPU: gpu/gpu_va.sh, arms vall, vfseh, vhp, vbe (job 31063797), in gpu/runs/j.31063797.

The hand-resolved TR slot behaves as intended:
- `tjse` is restricted to hsrc == 1 (vet_sc), so on vet_gd (hsrc == 2) the slot is set only for hcs, with eps from vgd_hr HREP.
- jse / bq = lin are inert on vet_gd. On vet_sc, hcs and jse both use 1 - eps from M1_VET_EPS.

## 2. 200-cycle AG Car A accuracy and cost (gpudev, 1 node x 4 A100)

### Setup

**Input.** gpu/inp/agcar_A_va.athinput is agcar_A_xtf plus the following:
- the command-line keys are named;
- hst output every cycle (dcycle 1);
- hydro_w / m1 / m1_thsw / m1_fs dumps every 50 cycles of dt.

**Fixed dt.** The cap `rad_m1/implicit_cfl` (dt_rad = implicit_cfl dxmin/c) is used:
- Arm values: implicit_cfl = 1251.3 for dt, /2 for dt/2, /4 for dt/4. These give dt = 1239.996 s, 619.998 s and 309.999 s on every cycle (spread 0).
- The hydro dt is 1252.56 s, so the cap is 1 % below production.
- The halving is exact (cap/2 in floating point). Matched times are matched cycles: 200 / 400 / 800 cycles all end at t = 2.479992e5 s.

**Rotation.** `vet_gd_rotate_every` is scaled with dt (10 / 20 / 40), so the ray sets have the same angles at the same simulated times.
- The `r` arms keep 10 at dt/2, to measure the rotation's own part.

**Arms** (on top of production):

| arm | keys |
|---|---|
| be | production |
| hp | `time_scheme=hesdirk2 time2_vet_col=predict implicit_vimp=true time2_vstage=true` |
| hpt | hp + `implicit_blend_xthin=0 implicit_thin_switch=fse implicit_thin_corr=hc` |

Suffixes: H = dt/2 (400 cycles), Q = dt/4 (800 cycles).

**Jobs:** 31064135 (be beH hp hpt), 31064136 (hpH hptH), 31064137 (hpQ beHr), 31064138 (beQ hpHr), and the hptQ follow-up (see 2.5).
- Runs are in gpu/runs/j.*, with links in gpu/acc/.
- Analysis: gpu/ana_va.py (perf / acc / ej / vsame), gpu/summ.py, gpu/hwin.py.
- Output tables: gpu/acc/{perf,acc_*,ej,vsame}.txt.

### 2.1 Robustness and cost

| arm | rc / FATAL | NON-CONV | Picard mean / max (per solve) | solves | s/cycle, no rotation / last | s/cycle, average | s/cycle, rotation step | cost: s per 1e4 s sim. (no-rot / avg) |
|---|---|---|---|---|---|---|---|---|
| be  (dt)   | 0 / 0 | 0 | 4.01 / 9 | 200 | 0.320 | 0.377 | 0.78 | 2.58 / 3.04 |
| beH (dt/2) | 0 / 0 | 0 | 3.36 / 8 | 400 | 0.295 | 0.324 | 0.77 | 4.76 / 5.22 |
| beQ (dt/4) | 0 / 0 | 0 | 2.93 / 7 | 800 | 0.283 | 0.298 | 0.77 | 9.12 / 9.60 |
| hp  (dt)   | 0 / 0 | 0 | 3.46 / 9 | 399 | 0.430 | 0.498 | 0.89 | 3.47 / 4.02 |
| hpH (dt/2) | 0 / 0 | 0 | 3.00 / 8 | 799 | 0.391 | 0.425 | 0.88 | 6.30 / 6.86 |
| hpQ (dt/4) | 0 / 0 | 0 | 2.72 / 8 | 1599 | 0.372 | 0.389 | 0.87 | 11.99 / 12.54 |
| hpt (dt)   | 0 / 0 | 0 | 3.34 / 9 | 399 | 0.439 | 0.502 | 0.90 | 3.54 / 4.05 |
| hptH (dt/2) | 0 / 0 | 0 | 3.00 / 8 | 799 | 0.406 | 0.439 | 0.87 | 6.56 / 7.08 |
| hptQ (dt/4) | 0 / 0 | 0 | 2.71 / 7 | 1599 | 0.384 | 0.401 | 0.86 | 12.38 / 12.95 |
| beHr (dt/2, rotate 10) | 0 / 0 | 0 | 3.43 / 8 | 400 | 0.298 | 0.348 | 0.75 | 4.80 / 5.61 |
| hpHr (dt/2, rotate 10) | 0 / 0 | 0 | 3.11 / 8 | 799 | 0.396 | 0.452 | 0.85 | 6.38 / 7.29 |

How the table was measured:
- s/cycle comes from the per-cycle `elapsed` lines (ndiag = 1), cycles >= 2.
- "No rotation" excludes the steps that start at a cycle with c % rotate_every == 0, and the last cycle (Finalize output).

General results:
- No arm has a NaN, a positivity fallback, a line_jacobi fallback or a stage fallback.
- hesdirk2 takes one backward-Euler start step and then 199 / 399 / 799 stage steps.
- The minimum E from the solve is 0.53 (hp) and 0.94 (hpt) at dt, against 1.43 for be. These are stage values, all positive.
- **hesdirk2 costs 1.34x be per step (no rotation; 1.32x on average) at the same dt.** Picard per step is 6.9 vs 4.0, but each solve is cheaper.
- **fse+hc costs +2 % over hp.**
- A rotation step costs about 0.45 s extra for every arm.

### 2.2 Time error at matched simulated time (Richardson)

**Error estimate:** err(dt) = 2^p/(2^p - 1) |u_dt - u_dt/2|, with p = 1 for be and p = 2 for hesdirk2.

**Observed order:** log2(|u_dt - u_dt/2| / |u_dt/2 - u_dt/4|).

**Metrics per r/R_ph band:**
- L1 = sum |a - b| / sum |b| over the cells of the band. This is the primary metric.
- Sh = max over the shells of the band of |<a> - <b>| / max |<b>|, with <> the (theta, phi) mean.
- Lmax = max over the shells of the band of max over (theta, phi) of |a - b| / max |b|. This is the vcpred metric.

**Cycle 200 (t = 2.48e5 s), err(dt), L1 (Sh):**

| var | arm | interior < 0.8 | 0.8-0.95 | photosphere 0.95-1.05 | envelope > 1.05 |
|---|---|---|---|---|---|
| dens | be  | 6.8e-5 (9.0e-6) | 9.8e-6 (2.0e-6) | 4.5e-5 (2.9e-5) | 8.2e-3 (8.6e-2) |
| dens | hp  | **2.7e-6** (1.1e-7) | 7.8e-6 (6.3e-7) | **1.0e-4** (1.2e-4) | **7.7e-5** (8.2e-4) |
| dens | hpt | 2.7e-6 (1.1e-7) | 7.8e-6 (7.9e-7) | 1.5e-4 (1.6e-4) | 1.3e-4 (4.5e-3) |
| E | be  | 5.5e-5 (9.3e-6) | 2.3e-5 (2.6e-6) | 3.2e-5 (1.8e-5) | 1.6e-5 (8.8e-6) |
| E | hp  | **2.2e-6** (8.5e-8) | 8.3e-6 (2.8e-7) | **6.6e-5** (3.9e-5) | 1.5e-5 (2.9e-5) |
| E | hpt | 2.2e-6 (8.5e-8) | 8.5e-6 (3.6e-7) | 1.0e-4 (5.8e-5) | 1.8e-5 (3.5e-5) |
| F1 | be  | 9.8e-3 (3.4e-3) | 3.4e-5 (3.2e-5) | 6.8e-6 (2.9e-6) | 1.5e-5 (6.2e-6) |
| F1 | hp  | **6.3e-4** (7.1e-5) | 7.4e-6 (2.4e-6) | 1.8e-6 (1.3e-6) | 4.0e-6 (7.1e-7) |
| F1 | hpt | 6.3e-4 (7.1e-5) | 7.6e-6 (3.0e-6) | 1.2e-6 (1.2e-6) | 6.6e-7 (4.3e-7) |
| v_r | be  | 5.5e-3 (1.0e-2) | 9.2e-4 (5.0e-3) | 5.9e-3 (1.2e-2) | 1.4e-2 (2.2e-1) |
| v_r | hp  | **4.3e-4** (3.2e-4) | 4.8e-4 (4.7e-3) | **2.4e-3** (6.9e-3) | **3.8e-4** (1.8e-3) |
| v_r | hpt | 4.3e-4 (3.2e-4) | 5.2e-4 (6.1e-3) | 3.0e-3 (1.1e-2) | 3.6e-4 (4.8e-3) |

**Lmax at the photosphere (cycle 200), dens / E / v_r:**
- be: 4.4e-2 / 5.2e-4 / 1.7e-1
- hp: 2.3e-2 / 9.1e-4 / 2.7e-1
- hpt: 3.3e-2 / 2.1e-3 / 1.2e-1

**Observed orders (dt, dt/2, dt/4; L1; cycles 50 / 100 / 150 / 200):**

| band | be (dens, E, F1, v_r) | hp (dens, E, F1, v_r) |
|---|---|---|
| interior | 1.0 for all four | dens 1.0-1.4 and E 1.0-1.2 (L1; Sh = 2.0); F1 1.3-1.9; v_r 0.8-1.4 |
| photosphere | dens 0.94-0.96, E 1.0-1.1, F1 0.9-1.1, v_r 0.2-0.8 | **dens 1.9-2.5, E 1.6-2.1, v_r 1.2-1.6 (Sh 2.1-2.4)**; F1 2.0-2.6 at cycles 50-150 |
| envelope | dens 1.0-1.4, E 0.9-1.1, F1 1.0 | dens -0.5 to 0.1 (non-convergent, but 100x below be); E 2.2-4.5; F1 0.9-2.4 |

**The hp interior's pointwise L1 is order ~1 with a tiny coefficient** (2-3e-6 for dens and E). The shell means are order 2.0.
- A cycle-keyed first-order part remains: either the rotation block jumps, or the single BE start step.

**Rotation (r arms, rotate_every 10 at dt/2 instead of 20).**
- Photosphere dens and v_r: identical within the digits shown.
- It raises the envelope E / F1 dt-differences about 7x (be E L1 1.6e-5 -> 1.1e-4; F1 1.5e-5 -> 1.0e-4) and the photosphere F1 Lmax (5.8e-5 -> 2.1e-4).
- So cycle-keyed rotation noise is NOT the photospheric dens / v_r time error. It does pollute the thin-envelope radiation field if dt is changed at fixed rotate_every.

### 2.3 History quantities

- **Mdot_top = 0 in every arm for all 200 cycles.** No mass reaches the outer boundary in 2.5e5 s, so there is no error to measure.
- **L_top and L_int, max |dt - dt/2| relative error by window** (Richardson):

| quantity | arm | cycles 1-30 | cycles 31-100 | cycles 101-200 |
|---|---|---|---|---|
| L_top | be | 2.1e-1 | 1.3e-5 | 4.4e-6 |
| L_top | hp | 2.3e-1 | 4.6e-5 | 8.6e-6 |
| L_top | hpt | 2.7e-1 | 5.8e-5 | 1.3e-5 |
| L_int | be | 3.3e-3 | 8.5e-6 | 2.1e-6 |
| L_int | hp | 1.2e-1 | 4.4e-5 | 4.5e-6 |
| L_int | hpt | 2.3e-1 | 5.9e-5 | 5.8e-6 |
| KE_int | be | 2.5e-3 | 5.1e-3 | 9.2e-3 (order 1.0, growing) |
| KE_int | hp | 8.7e-4 | 8.4e-4 | 6.4e-4 |
| KE_int | hpt | 8.7e-4 | 8.4e-4 | 6.4e-4 |

  - **The first ~30 cycles carry a 10-27 % L_top start-up transient in every scheme**, and its dt-dependence does not converge (order about 0). After cycle 30 L_top is constant to about 1e-5. Its dt-error is 1e-5 (be) to 4e-5 (hp, hpt) and falls to below 1e-5 by cycles 100-200.
  - L_top orders after cycle 100 are about 0, i.e. a floor at the 1e-6 level.
- **Photospheric v_r** (shell-mean velx at r = R_ph, from the dumps; |dt - dt/2| / |v|):

| cycle | be | hp | hpt |
|---|---|---|---|
| 50 | 1.4e-3 | 6.6e-3 | 7.3e-3 |
| 100 | 5.3e-4 | 7.1e-3 | 9.9e-3 |
| 150 | 1.0e-3 | 5.2e-3 | 1.2e-2 |
| 200 | 6.6e-4 | 1.2e-3 | 2.7e-3 |

  - Richardson err(dt), cycle 200: be 1.3e-3, hp 1.7e-3, hpt 3.6e-3.
  - hp, dt/2 vs dt/4: 5.5e-4 / 2.2e-3 / 6.9e-4 / 2e-6 relative at cycles 50 / 100 / 150 / 200. That is order 3.6 / 1.7 / 2.9 (at cycle 200 the dt/2 and dt/4 values coincide). So hp v_r(R_ph) converges at about second order or better.
  - be's dt/2 - dt/4: 3.2e-4 / 3.2e-4 / 5.3e-4 / 3.9e-4, order 2.2 / 0.7 / 0.95 / 0.8.

### 2.4 E/J in the thin envelope (cycle 200, median and p5-p95)

| arm | 1.00-1.05 R_ph | 1.05-1.5 | 1.5-2 | 2-3 |
|---|---|---|---|---|
| be / beH / beQ | 1.0182 / 1.0180 / 1.0179 | 1.0179 / 1.0177 / 1.0176 (0.987-1.047) | 1.0071 / 1.0069 / 1.0068 | 1.0095 / 1.0092 / 1.0091 |
| hp / hpH / hpQ | 1.0178 / 1.0178 / 1.0181 | 1.0175 / 1.0176 / 1.0178 | 1.0067 / 1.0067 / 1.0070 | 1.0091 / 1.0091 / 1.0093 |
| hpt / hptH | 1.0027 / 1.0027 | **1.0010 / 1.0011** (0.972-1.030) | 0.9994 | 1.0069 |
| hpt, twin-corrected J (m1_thsw) | 1.0026 (1.0021-1.0030) | **1.0008** (1.0002-1.0018) | 1.0013 | 1.0102 |

- J is the raw vet_gd formal solution (m1_fs, no twin correction) for every arm, so the p5-p95 widths are ray noise. The m1_thsw J (twin-corrected) exists only with the switch on.
- E/J is independent of dt and of the time scheme.
- **The production / hp thin envelope carries a 1.8 % excess of E over J, and fse+hc removes it** (1.001). This agrees with gdthin.

### 2.5 hptQ (follow-up) and the hpt orders

- gpudev was congested. A copy of job 31064811 was submitted on `gpu` (31065176, 20 min limit).
- 31064811 started first, so 31065176 was cancelled. hptQ ran: rc 0, 800 cycles, 6 min.

**hpt (fse+hc) observed orders, dt / dt/2 / dt/4, L1, cycles 50-200:**

| band | dens | E | F1 | v_r |
|---|---|---|---|---|
| photosphere | **1.9-2.0** | **1.9-2.1** | 0.6-1.7 | 1.4-1.6 (Sh 1.8-1.9) |
| envelope | 0.3-1.4 (Sh 2.4-2.5; hp: non-convergent) | 1.9-2.1 | | |
| interior | as hp | as hp | as hp | as hp |

- v_r(R_ph), dt/2 - dt/4: 2.5e-3 / 1.4e-3 / 8.5e-4 / 1.7e-3 relative, order 1.5 / 2.8 / 3.8 / 0.7.
- L_int after cycle 30: order 1.8.
- So **fse+hc is also second order at the photosphere**. Its larger photospheric time error (about 1.5x hp's) is a coefficient, not an order loss.

## 3. Answers

**1. Does hesdirk2's photosphere error shrink with more cycles, i.e. was the 30-cycle result a start-up transient?**
- **No.** The photospheric dt-error of hp does not shrink between cycles 50 and 200. L1 err(dt) for dens, cycles 50 / 100 / 150 / 200: 1.1 / 1.3 / 1.0 / 1.0 e-4. Lmax grows from 1e-2 to 2.3e-2.
- At dt it stays about 2.2x be's for dens (1.0e-4 vs 4.5e-5) and 2x for E (6.6e-5 vs 3.2e-5).
- **But it is a genuine, converging second-order error, not a defect.**
  - The dt / dt/2 / dt/4 ladder gives order 1.9-2.5 for photospheric dens, 1.6-2.1 for E, and 2.1-2.4 for the shell-mean v_r (Sh), at every dump. The v_r(R_ph) orders are 1.7-3.6.
  - be is order 1.0 there.
  - hp - be at the photosphere drops 3.9x from dt to dt/2 (dens L1 1.1e-4 -> 2.9e-5). So the difference is hp's own error, and be's first-order error is small at the photosphere.
- At dt/2 the two are equal (dens 2.5e-5 vs 2.2e-5). At dt/4 hesdirk2 wins.
- **The real start-up transient is the first ~30 cycles.** The L_top dt-differences there are 10-27 % for every scheme (order 0). They fall by about 4 orders of magnitude after cycle 30. The vcpred 30-cycle numbers sit inside that transient.
- **What hesdirk2 buys at the production dt:**
  - interior: 25x lower dens / E error (2.7e-6 vs 6.8e-5) and 16x lower F1 error;
  - envelope: 100x lower dens error, 40x lower v_r error;
  - photosphere: 2.5x lower v_r error (L1 2.4e-3 vs 5.9e-3);
  - KE_int: 14x lower error (6.4e-4 vs 9.2e-3, and be's grows linearly in time).
  - It costs 1.34x per step, still cheaper than be at dt/2 (cost 3.5 vs 4.8 s per 1e4 s).
  - Only photospheric dens / E are about 2x worse at the same dt, and those errors are small in absolute terms (1e-4 L1).
- **The photospheric dens / v_r error is not from the cycle-keyed rotation** (the r arms give the same numbers).

**2. What does fse+hc change at the photosphere?**
- **Time error:** slightly larger, about 1.5x hp's, and still second order (photospheric dens / E orders 1.9-2.1 from hptQ).
  - dens L1 1.5e-4 vs 1.0e-4; E 1.0e-4 vs 6.6e-5; v_r L1 3.0e-3 vs 2.4e-3; v_r(R_ph) 3.6e-3 vs 1.7e-3.
  - The interior is unchanged (bitwise to 1e-7).
  - The envelope dens Sh error rises from 8e-4 to 4.5e-3.
- **Solution: a structural change about 30x larger than any time error.** hpt - hp at cycle 200, photosphere band:
  - L1: dens 6.3e-4, E 7.1e-4, v_r 1.3e-2.
  - Shell means: E +4-6 % at 1.000-1.006 R_ph (the switch-on layer, where w goes from 0 to 1 over 1.003-1.008).
  - dens -3 % at 1.001 R_ph and +20 % at 1.008-1.011 (a density shell displaced outward).
  - v_r +0.6 km/s at 1.008 and +1.2 km/s at 1.011.
  - Above 1.02 R_ph: E +0.24 %, J +1.9 %, E/J 1.018 -> 1.001.
  - The difference is the same at dt/2 (hptH - hpH: dens L1 6.3e-4, v_r 1.3e-2), so it is a model change, not a time error.
- **Robustness and cost:** NON-CONVERGED 0, Picard 3.34 vs 3.46, cost +2 %, min E 0.94.
- L_top / L_int after cycle 30 move by < 1e-5.

## 4. Recommendation

1. **The merge vetacc-1011 is good:** every branch's feature is bitwise against its own binary, and the defaults are bitwise against the base. It can be the base for the next work.
2. **hesdirk2 + predict + vimp + vstage is ready to be the AG Car A production time scheme.**
   - 200 cycles with no failures or fallbacks, at 1.34x be per step.
   - It is second order where it matters. At the production dt it beats be at dt/2 in the interior, envelope, v_r and KE, at lower cost.
   - Its photospheric dens / E error at dt is about 2x be's but second order, ~1e-4 L1. If that matters, dt/2 hesdirk2 (6.3 per 1e4 s) still costs less than be at dt/4 (9.1) and is more accurate everywhere.
   - The 30-cycle "photosphere does not improve" finding was the start-up transient plus a larger second-order coefficient, not a first-order defect.
3. **fse+hc (xthin 0):** adopt it for the thin-envelope physics (E/J 1.018 -> 1.001, a dt-independent flux model; the production xthin = 30 weight is dt-dependent).
   - It changes the photospheric layer structure (E +6 %, a displaced density shell) much more than any time error does. The accuracy question there is spatial / model, not temporal.
   - Before production, the earlier caveat still applies: there is no vet_gd thermal-beam proof yet. Validate the photospheric structure against a reference (e.g. resolution, or a vet_sc 1-D analogue).
4. **When comparing runs at different dt, scale `vet_gd_rotate_every` with 1/dt.** Otherwise the ray-noise pattern changes and the envelope E / F differences grow about 7x.
5. **The first ~30 cycles from the IC are a transient with O(10-25 %) L_top swings.** Discard them in any accuracy or science diagnostic.

GPU time used: about 45 node-minutes of 4 A100 (validation 4.5 min, accuracy 34 min, hptQ ~7 min).
