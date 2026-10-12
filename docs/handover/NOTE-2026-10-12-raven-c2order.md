# NOTE c2order: hesdirk2 time order of the conservative gas source (rlat, rsw, rswv, P = 10 waves)

Worker T, 2026-10-12, MPCDF Raven.

- Work dir: /raven/ptmp/jinma/c2order_1012.
- Branch: c2order-1012 (worktree /raven/u/jinma/ATHENAK/wt-c2order), not pushed.
  - 15385870: revert of a90314fe, so implicit_rad_inertia_solve acts on all faces again.
  - 1f453257: merge of prodfinal-1011 4253af30 (implicit_cs2_hll = drag).
- Binaries in bin/:
  - athena_cpu_a5: L's a90314fe build.
  - athena_cpu_c1 and athena_he_cpu_c1: 1f453257.
  - athena_*_c2: 1f453257 plus a probe change that was later reverted. It is only reached with time2_start = slope. The AG Car and P = 10 runs below use c2 without the probe.

## 1. Verdict

1. **G and L disagreed; both were partly right.**
   - L: the after-solve split of implicit_rad_inertia (the old ImplicitRadInertia rule) is not consistent, and it is NOT second order.
   - G: the after-solve split gave rlat 2.0 against 1.8 for the in-solve split. That 2.0 only holds at levels 1-8 and is a cancellation.
   - Over levels 1-64 on rlat_T (fresh), dens orders:

     | arm | 1-2 | 2-4 | 4-8 | 8-16 | 16-32 | 32-64 |
     |---|---|---|---|---|---|---|
     | D: matched, after-solve | 1.94 | 2.16 | 1.10 | 0.12 | 0.58 | |
     | conservative2 dbg 1: E part after the solve | same as D: floor of about 5e-7 | | | | | |
     | S5: x1 in-solve, x2/x3 after-solve (a90314fe) | 1.80 | 1.72 | 1.42 | 1.16 | 1.00 | |
     | conservative2 | 1.81 | 1.73 | 1.62 | 1.48 | 1.29 | |
     | no inertia | 1.78 | 1.76 | 1.58 | 1.46 | | |

   - D - S5 is a theta-dependent model difference. It is identical under time2_tableau = trbdf2, and its size is 6e-6 down to 3e-6 at levels 1-8, decreasing with refinement.
     - The cause: the after-solve rule leaves the E rows without -theta dA. In the stiff regime that term is -a dE/dt/(chat rho kappa), O(1) in dt.
     - In the inner (wall) cells the after-solve x2/x3 split is plainly first order. Errors against level 8 at levels 1/2/4:
       - D and S5: 4.4e-6, 1.9e-6, 0.86e-6;
       - in-solve on all faces, conservative2 or no inertia: 3.8e-6, 0.9e-6, 0.18e-6.
   - **Code: a90314fe (x1 faces only) is reverted (15385870).** implicit_rad_inertia_solve is in-solve on all faces again, opt-in and default false.
2. **The remaining rlat deficit is start-up and the vet_gd rotation, not the gas source.**
   - **Pre-relaxed restart ladder** (rlat_Tr: pre-run to t of about 1 at level 16, then levels 1-16 from that state; errors against level 16; dens / vely / E):

     | arm | dens | vely | E | orders |
     |---|---|---|---|---|
     | **conservative2, rotate 0** | 5.8e-7, 1.35e-7, 3.5e-8, 8.9e-9 | 1.2e-5 ... 1.4e-7 | 3.8e-7 ... 5.4e-9 | **2.0** |
     | matched + rinr_solve (all faces), rotate 0 | 6.7e-7 ... 1.1e-8 | | | 2.0-2.1 |
     | matched, no inertia, rotate 0 | | | | 2.0-2.1 |
     | conservative2, rotate 10 | | | | 1.86, 1.48, 1.14 |

     The rotation is per cycle, so each level samples different quadrature sets at a given time.
   - **Fresh, kappa_s 100, conservative2, rotate 0:** dens 2.13, 2.11, 2.12; E 1.98, 1.95, 1.97 (levels 1-16). With rotate 10 every arm flattened to 0.5-0.8.
   - **Fresh, kappa_s 1000:** conservative2 1.73, 1.74, 1.64, 1.45. It is the same with rotate 0.
     - This is the stiff initial layer: the IC has F = 0 with a lateral gradient of E in a thick medium.
     - time2_start_frac 1e-3 makes it worse (a floor at about 2e-7).
     - rlat's fresh ladder therefore measures the IC, not the scheme.
   - **The restart is bitwise exact.** A continuous run and a restart from its own t = 1 rst give 0.0 difference at t = 1.5 and t = 2, in D and in S5.
     - L's rlat_Tr "floor" is the rotation plus the after-solve inner-cell error, not lost restart state.
3. **rw_p10t1** (Cartesian, P = 10, tau = 1, amplitude 1e-5; tol 1e-14 / lin 1e-15):
   - **matched: the E and F1 first-order term is per step.**
     - Pre-relaxed restart ladder: E 1.81 1.74 1.47, F1 1.87 1.54 1.39; the gas is 2.0. With eddington: E 1.99 1.79 1.14.
     - Unchanged by: implicit_rad_inertia=false; the thin switch and thin corr; closure_lag pass; vet_x1_npass 40; vet_sc rebuild; bmom_half; f0_init; pred_cap 0; nofb_work; enthalpy central; trbdf2; time2_start_frac.
     - Removed by: hydro reconstruct = dc (E 1.99 1.99 1.99), and by conservative2.
     - Worse with: enth_vel old/start, enthalpy upwind, vimp off, vstage off (F1 exactly 1.00).
     - It is the matched F0 model's bookkeeping of the velocity part of d(aE): conservative2's lab-F face equation removes it. It is not an inertia-split term.
   - **conservative2: second order per step.** Pre-relaxed restart ladder, every variable 1.86-2.11 (eddington: E 2.11 2.22 2.43, F1 2.0-2.3).
     - Fresh: E 1.81 1.92 2.00 and F1 1.85 1.67 1.56, with dens 1.91 1.84 1.70. The gas part is start-up.
     - implicit_cs2_dbg=2 (no Cs2Vimp prediction) makes the fresh gas 2.0 but F1 1.66. Not adopted.
4. **No new default; no production change.**
   - The only kept code change is the revert of a90314fe, inside the opt-in implicit_rad_inertia_solve.
   - conservative2 needs no change for second order per step.
   - The production deck (conservative2 + cs2_hll drag + vimp_bc) is reproduced bitwise; see 5.4.

## 2. Orders with c1 (fresh T ladders, slope + vimp_bc, orders 1-2 ... 16-32)

**rsw_T**

| arm | E | F1 | dens | v |
|---|---|---|---|---|
| conservative2 | 2.00 1.99 1.99 1.98 | 2.00 2.00 2.00 1.99 | 2.03 2.06 2.11 2.15 | 2.00 2.00 2.00 1.99 |
| matched + rinr_solve (all faces) | 1.99 1.97 1.94 1.90 | 2.00 1.99 1.97 1.92 | 2.03 2.05 2.08 2.09 | 2.00 1.99 1.98 1.93 |
| matched D | 1.82 1.47 1.05 0.86 | | 2.03 1.76 0.69 0.62 | |

**rswv_T**

| arm | E | F1 | dens | v |
|---|---|---|---|---|
| matched + rinr all faces + f0_init zero + rotate 0 | 1.96 1.93 1.88 1.79 | 1.97 1.94 1.88 1.78 | **2.02 2.02 2.02 1.95** | 1.98 1.95 1.90 1.80 |
| conservative2 + f0z + rotate 0 + time2_start be (X's arm) | 1.91 1.77 1.42 1.23 | 1.90 1.78 1.56 1.28 | 2.01 2.01 1.96 1.67 | 1.90 1.78 1.56 1.28 |
| conservative2 + f0 auto | 1.90 1.77 1.72 1.54 | | 1.3 (start-up) | |
| D | about 0.5 | | | |

- **Caveat:** time2_start = slope with an IC whose u0 F differs from F0 + aE (f0_init zero) is wrong under conservative2.
  - Probe as it is: rswv order 1.00 everywhere.
  - Restoring gcs_first after the probe: a floor. Reverted.
  - The probe's K1 of the cell F carries (aE)/dt_p. Use time2_start = be for such ICs.

**rw_p10t100** (tol 1e-14): every arm 2.0 down to level 16. At 16-32 the tolerance / round-off floor is reached (about 2e-7 of the amplitude), in matched as well.

**rw_t10a** (F1 orders):
- conservative2: 2.25 2.32 2.16 1.90;
- matched + rinr: 1.95 1.96 1.99 1.82;
- every other variable 1.9-2.1.

**rlat_T:** see 1.2 (the scheme is 2.0 on the pre-relaxed ladder).

**rw_p10t1:** see 1.3.

## 3. Conservation (c1)

| check | conservative2, hesdirk2 | conservative2, BE | matched (with or without rinr_solve) |
|---|---|---|---|
| eq3f P_tot drift | 1.3e-15 | 8.9e-16 | 6.2e-4 (non-conservative source) |
| eq3f E / v | 1.13260 / 2.95538 | 1.13260 / 2.95538 | |
| rlat_om5 L_z drift / max\|L_rad\| at t = 4 | +5.26e-6 | +5.44e-6 | -4.15e-2 |

## 5. Regression

### 5.4 Reduced CPU AG Car A

- Setup: 480x32x32, 60 cycles, agc/runs/j.31092388. Base athena_he_cpu_4253af30 against athena_he_cpu_c2.
- FRESH A0 deck (matched): every hydro_w / m1 bin dump is **bitwise**. Picard mean 3.02.
- A0 + conservative2 + cs2_shift false + vimp_bc true (the v2 production keys): **bitwise**, Picard 3.02.
- A0 + implicit_rad_inertia_solve (all faces) + vimp_bc: rc 0, Picard mean **4.18** (max 7) against 3.02.
  - The cost is the stop without a contraction estimate (L saw the same).
  - So matched + rinr_solve stays opt-in. conservative2 has no extra cost.

### 5.5 P = 10 waves, conservative2, hesdirk2, against vetqs_d (rw/, logs/rw_c2.txt)

| tau | L1(rho) N 16 / 64 / 256 | omega(256) | model root |
|---|---|---|---|
| 0.01 | 6.1e-2 / 5.4e-3 / 2.7e-4 | 7.99235 - 0.51941i | 7.99183 - 0.51937i |
| 1 | 3.7e-2 / 2.7e-3 / 1.6e-4 | 6.28409 - 0.85814i | 6.28356 - 0.85844i |
| 100 | 6.2e-2 / 3.8e-3 / 2.3e-4 | 14.38446 - 0.37942i | 14.38553 - 0.37926i |

These are the same as G's 27c35877 table: the dispersion is intact, including cs2_hll = drag.

## 4. Files

- ana/:
  - h.py: L's harness, W = here.
  - dif.py: arm difference per level.
  - ref.py, reft.py: errors against a reference level.
  - eqtot.py: eq3f totals.
- sp/ana: G's sp harness, with arm suffix S = implicit_rad_inertia_solve.
- rw/: G's radwave harness.
- reg/run_cs.sh: eq3f runs.
- agc/agc.sh: reduced AG Car A.
- rst.sh: restart ladder.
- jobs_1 ... jobs_10, job*.sbatch, runs/<case>/<arm>_<lev>.
