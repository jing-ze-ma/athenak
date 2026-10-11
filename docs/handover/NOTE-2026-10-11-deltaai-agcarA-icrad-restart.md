# NOTE DeltaAI -> raven: AG Car A follow-up (TASK-2026-10-11-deltaai-agcarA-icrad-restart) -- DONE

Picked up 2026-10-10 ~22:05 CDT; result ~22:25 CDT. Job **3357234** (ghx4-interactive gh095, 1 node x 4 GH200, 1 rank/GPU,
MPICH_GPU_SUPPORT_ENABLED=1), wall 2 min 42 s. Same binaries, files and deck handling as NOTE-2026-10-11-deltaai-agcarA-face-icrad
(B = icrad-1011 7045bc14, cell boundary; A = deepatm-1011 e84d2513). Fixed dt 1239.996 s. **All 5 runs rc 0, FATAL 0, NON-CONVERGED 0.**

## Per run
Picard = the code's end-of-run line (includes the ic_rad relax solves); s/cycle = median (mean) of per-cycle wall, excluding rotation cycles and the last.

| run | binary | what | cycles | rc | FATAL | NON-CONV | Picard mean / max | s/cycle |
|---|---|---|---|---|---|---|---|---|
| berx | B | ic_rad=relax, rst every 10 cycles | 0-20 | 0 | 0 | 0 | 4.09 / 9 (34 solves = 20 + 14 relax) | 0.193 (0.230) |
| rs10 | B | restart from berx cycle 10 | 10-20 | 0 | 0 | 0 | 2.90 / 4 | 0.196 |
| rs0 | B | restart from berx cycle 0 (rst written after the IC step) | 0-10 | 0 | 0 | 0 | 3.10 / 5 | 0.194 |
| hprx | B | relax + hesdirk2, time2_vet_col predict, implicit_vimp, time2_vstage | 60 | 0 | 0 | 0 | 4.10 / 13 (142 solves, 23 relax) | 0.263 (0.326) |
| hc1 | A | face + vet_gd_halo_compact=1 | 50 | 0 | 0 | 0 | 3.90 / 9 | 0.214 (0.300) |

## Restart (berx / rs10 / rs0)
- berx wrote rst at cycles 0 (after the IC step, 1.60 GB: no solver-history blocks yet), 10, 20 (and the final 20 again, 1.98 GB each).
  Cycle read from each rst header (time, dt, ncycle).
- **IC step re-run on restart: NO**, both rs10 and rs0. The only ic_rad line in either log is
  `<rad_m1> ic_rad = relax (fresh starts only) ... [RESTART: the IC passes are skipped]`; no relax/fs pass lines.
- **rs10 is bitwise vs berx**: user.hst and hydro.hst rows for cycles 11-20 identical text (10/10 each); cycle-20 bin
  hydro_w, m1, m1_fs, m1_thsw arrays identical; the cycle-20 rst identical after the parameter text (1 983 659 334 bytes).
- **rs0 is bitwise vs berx**: hst rows cycles 1-10 identical (10/10 each); its cycle-10 rst identical to berx's cycle-10 rst
  after the parameter text. (No cycle-10 bin in berx to compare; the rst covers the full state.)

## hprx (hesdirk2 + relax) vs B relax (be), L_top (erg/s)
| cycle | B relax (be) | hprx (hesdirk2) | rel diff |
|---|---|---|---|
| 0 | 2.844903e38 | 2.848614e38 | +1.3e-3 |
| 1 | 2.844750e38 | 2.848449e38 | +1.3e-3 |
| 2 | 2.844658e38 | 2.808015e38 | -1.3e-2 |
| 3 | 2.844589e38 | 2.635723e38 | -7.3e-2 |
| 4 | 2.844526e38 | **2.533218e38** | **-11.0 %** |
| 5 | 2.844466e38 | 2.789574e38 | -1.9e-2 |
| 6 | 2.844413e38 | 2.811953e38 | -1.1e-2 |
| 7 | 2.844367e38 | 2.765527e38 | -2.8e-2 |
| 8 | 2.844327e38 | 2.819960e38 | -8.6e-3 |
| 9 | 2.844293e38 | 2.837506e38 | -2.4e-3 |
| 10 | 2.844264e38 | 2.834742e38 | -3.3e-3 |
| 30 | 2.844063e38 | 2.842599e38 | -5.1e-4 |
| 60 | 2.844066e38 | 2.843561e38 | -1.8e-4 |

- With hesdirk2, relax still removes the +54 % step-1 spike (default be: 4.366e38 at cycle 1), but **L_top dips by up to 11 %
  (cycle 4) and rings over cycles 2-8**, then converges to the be relax curve (1.8e-4 by cycle 60). be + relax is flat to 6e-4.
- The relax step under hesdirk2 takes 23 steps (5.9 s) vs 14 (3.9 s) under be, and ends at L/L_IC 1.0019 vs 1.0006.
- hesdirk2 costs 0.263 vs 0.203 s/cycle (+30 %), Picard max 13.

## hc1 (vet_gd_halo_compact = 1) vs A face (= 2)
**Bitwise**: user.hst and hydro.hst md5-identical; all 12 bin files (hydro_w, m1, m1_fs, m1_thsw at cycles 0, 50, final)
array-identical. (The bin file md5s differ only because the bin header carries the input deck, which names halo_compact.)
Cost 0.214 vs 0.199 s/cycle median (single runs; within run-to-run noise of the mean, 0.300 vs 0.270, not resolved here).

## Addendum (item 5): does relax perturb the envelope gas? (job 3357130 runs, no new runs)
**Cycle 0 (dumps written after the IC step), relax and fs_relax vs default:**
- dens, velx, vely, velz: **bitwise identical**. Only the gas internal energy changes (gas_feedback in the relax solves).
- Max relative change of the shell-mean eint per band (relax = fs_relax to the digits shown):
  < 0.8 R_ph 1.1e-5 | 0.95-1.05 R_ph 0.39 | 1.05-1.5 R_ph 0.47 | > 1.5 R_ph 0.97.
- Sign: relax **cools** the gas above the photosphere. Shell-mean eint relax/default at cycle 0: r/R_ph 0.98 0.996, 1.00 0.953,
  1.02 0.897, 1.05 0.606, 1.2 0.597, 1.5 0.524, 2.0 0.127, 2.5 0.039, 2.9 0.029.
- **The default run gets to the same gas state by itself**: the same ratio at cycle 50 / 60 is 0.997-1.009 at every radius above.
  So relax does at t = 0 what the default IC does in its first cycles (with the L_top spike).
- Gas pressure: no P or T output exists in these binaries (general EOS, tabulated). At fixed rho, P_gas changes with eint by a
  similar fraction (exactly so only where the ionization state does not change); I did not evaluate the EOS table.
**Velocities over cycles 0-60:**
- max |velx| per dump: cycle 0 is 0 in all three; cycle 50: default 7.32e5, relax 7.58e5, fs_relax 7.58e5 cm/s; cycle 60: 8.84e5 / 9.11e5 / 9.11e5.
- KE_int (hst): relax is 9 % **below** the default at cycle 1 (1.430e39 vs 1.576e39), 1.4 % below at cycle 10, 0.5 % below at
  cycle 60 (1.731e42 vs 1.741e42). KEr_int follows the same pattern (0.6 % below at cycle 60). Both grow smoothly, monotonically.
- **Relax creates no velocity transient.** Neither does the default IC: in both, KE grows smoothly from the seed; the default's
  L_top spike is radiative, and its KE is slightly higher than relax's throughout (by 9 % at cycle 1, 0.5 % at cycle 60).
  max |velx| is 3 % higher with relax at cycles 50/60 (single-dump maxima, not a transient).

## GPU-h
3357234: 2.7 min x 4 GH200 = **0.18 GPU-h**. (Addendum: 0.) Running total for the two AG Car A tasks: 0.84 GPU-h.

## Files
In git: docs/handover/deltaai-agcarA-1011/ (job2.sub, rstcyc.py, rstcmp.py, addendum.py, rstcmp.txt, addendum.txt).
On DeltaAI: /work/nvme/bivj/jma20/agcarA_1011/runs/3357234/<run>/ and logs/agcA_icrad_rst.3357234.out.
