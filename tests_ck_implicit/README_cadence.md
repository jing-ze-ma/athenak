# ck-cadence: implicit ck RT every N steps, with a linearised step in between (2026-09-24)

Branch `ck-cadence` (worktree /viper/ptmp2/jinma/wt_ckcad, from rt-integration 90b01f2f).
Runs and binaries are in /viper/ptmp2/jinma/ckcad_0924. Scripts are in
`tests_ck_implicit/cadence/`.

The goal was to bring the correlated-k RT from about 2x hydro (c2) toward 1x hydro.

**Result.** e4 (every = 4, no guard) costs **RT = 0.94x hydro** (28.64 against
14.78 ms/cycle; c2 is 2.16x in the same job). Its accuracy against c2: day rms 2.0e-5,
kinks in the < 1e-7 bar band 733 (c2: 716, +2.4%; t4: 737), 0 non-converged calls.
Without a guard, N = 8 fails (166 non-converged calls, 851 kinks). The guard (thr) does
not save cost: a partial call costs almost as much as a full one (section 5).

## 1. Design

All switches are in `<problem>` and default off. `ck_impl_every = 1` is bitwise the
code without them (section 3).

| switch | meaning |
|---|---|
| `ck_impl_every = N` (1) | Full implicit call (c2 levers, `ck_impl_once`) on the cycles with `ncycle % N == 0`, and on the first call of a run or a restart. Linearised step on the other cycles. |
| `ck_impl_every_thr = thr` (0) | Per-column refresh guard. A column gets a full solve on this step if, since its last full solve, any RT cell has moved by more than thr in \|T/T0 - 1\| or \|rho/rho0 - 1\|. |

Requirements: ck_implicit, ck_impl_once, ck_impl_fuse, ck_impl_colskip, debug <= 0,
glob none, warm off. A missing requirement is a FATAL error.

**What is stored.** After each full call, for every cell of every solved column:
- `Q0 = (e - e*)/dt`, the heating rate the gas actually received. This is `ck_dep`, the
  call's applied increment. It is exact even when ck_impl_pred ended the Newton.
- `D = dQ/dT` at fixed rho. For a thick cell this is the Newton diagonal `jac(m,1)`. For
  a thin cell it is `-4 E/T`, which is what CkThinSolve linearises. D is clamped to
  <= 0, and NaN becomes 0.
- `T0` (EOS, K) and `rho0`, both at the final state.

Cells below the cut store Q0 = D = 0.

**The linearised step.** It is in `CkCadLin`, one team per column. Each cell takes

    dT = (Q0 + D (T^n - T0)) dt / (c_v - D dt),   de = c_v dT

where c_v is the EOS heat capacity, `rho SpecificHeatCv(rho, e, T)/temp_cgs` per kelvin.
The step is implicit in T. With D <= 0 the denominator is >= c_v for any dt, so the step
is stable in the thin top. It relaxes T toward T0 - Q0/D. |de| is capped at
min(ck_impl_dtmax, ck_impl_demax) e, the caps of the Newton step. On the production A/B
the cap never fired except in e16 (5090 cells).

**Guard.** CkCadLin first measures each column's deviation. A flagged column gets mask 0
and is not touched. Every other column gets mask 2. When any column is flagged,
`picket_fence_two_stream_RT` then runs with `ck_cad_partial`:
- `ck_done` starts from the mask.
- The pass kernels skip done > 0 columns through the existing colskip test.
- The fused step skips done > 1.5 columns before its residual.

A partial call always stores its operator. The operator is then dropped
(`xs_cyc = -1`) because it is valid for the flagged columns only.

**Interval: each step covers its own dt.** Neither "sum of dt over the elapsed interval"
nor "lead over the next N" is used. Reasons:
- With the linearised step in between, a full call over the elapsed interval would count
  that interval twice, unless it first removes the linear increments.
- That removal is per cell, not per fluid parcel. Over N steps the gas has advected, so
  the removal can take energy out of cells that never received it, including toward
  e < 0.
- A per-column guard would also need a separate interval, and a separate bdt, for every
  column.
- Leading over N steps front-loads N dt of heating from a predicted dt.

With each step covering its own dt:
- Every step deposits exactly one RT increment over exactly its own dt, so nothing is
  lagged, lumped or double counted.
- A full-call step is c2's own backward-Euler balance, and its ckdesum closes to the
  Newton tolerance.
- The only energy error is on the linear steps. There, Q0 dt is the flux divergence of
  the last full solve, and D (T - T0) dt is a local, non-flux-form correction.

This error is measured on every step (verbose line `### ck_cadence`):
- `edef = (sum dx de - sum dx Q0 dt)/sum dx |Q0| dt` is the part of the linear deposit
  that is not the conservative Q0 part.
- `linerr = sum dx |Q_new - (Q0 + D (T - T0))_old| / sum dx |Q_new|` is the error of the
  linear model, measured at the next full call.

The step counter is `ncycle`, which is in the restart. The stored arrays are **not**
written to the restart, so the first call after a restart is a full call (section 3).

Code (+375 lines):
- `src/utils/two_stream_column_ck.hpp`: the switches and state, and the mask test in the
  fused step.
- `src/utils/two_stream_rt.hpp`: the ck_done mask, partial calls force a store, and
  CkCadState / CkCadStore / CkCadLin / CkCadStep at the end of the file.
- `src/pgen/deep_hot_jupiter_rt.cpp`: the keys, the requirement check, and the dispatch
  in DhjCkRtSplit.

## 2. Binaries

| binary | built from |
|---|---|
| CPU base: `athena.cpu.base` | 90b01f2f |
| CPU cad: `athena.cpu.cad1` | working tree = the ck-cadence commit |
| GPU: `athena.gpu.cad` | same source (git stash 224aaf64); md5 9b1dda83929059ab8098ddefe732a629 |

All are in /viper/ptmp2/jinma/ckcad_0924. The GPU binary serves every arm, including
c2 and h.

## 3. Gate A (CPU)

`cadence/cpugate.sh` writes `cpugate.out`.

| comparison | result |
|---|---|
| base vs cad, every = 1, well-posed A2 t4 and semi (rec.txt + rst) | BITWISE |
| base vs cad, full hydro 6 cycles from the wp IC, T4 and c2 (rst + hst) | BITWISE |

**Restart (not bitwise).** `cadence/rstgate3.sh` and `rstcmp.py` write
`rstgate3_cmp.txt`. The runs are the wp full-hydro IC, 8 cycles straight against 4 +
restart + 4. The restart lands on a scheduled full cycle, 4 % 4 = 0. The table gives the
max relative difference of e (and rho) in the dumps after the restart.

| arm | cycle 5 | cycle 8 |
|---|---|---|
| e1 (= c2) | 8.6e-5 | 3.1e-5 (rho 2.0e-5) |
| e4 | 3.3e-4 | 1.8e-3 (rho 1.4e-4) |
| e4g2 | 8.7e-5 | 2.6e-3 (rho 7.2e-5) |

c2 itself is not a bitwise restart. The Newton carries state from call to call that the
restart does not hold, for example the pass-0 thin/thick classification `ck_thk`, so
its first post-restart call already differs. This was not known before; it predates
this branch.

The cadence adds nothing new at the restart cycle, because the forced full call is the
scheduled one. It does amplify c2's own difference through the stored Q0/D over the
linear steps. This wp case is a violent transient from the IC: c2 itself has
NOT-CONVERGED calls on cycles 9-11.

Writing the stored arrays to the restart would not make the restart bitwise until c2's
own cross-call state is also saved. This is stated, not fixed.

## 4. Gate B: production A/B (GPU)

Setup:
- Jobs 11956028, 11956029 and 11956379 (apudev, 2 ranks, HSA_XNACK=1,
  HSA_NO_SCRATCH_RECLAIM=1).
- Restart: bench/cs_hyd4_prod/rst/dhj.00138.rst, read in place: t = 2.1045e7 s,
  **rotation 69.0**, the same state as job 11955220.
- Each arm runs to t0 + 3.05e4 s (0.1 rotation, about 1660 cycles).
- Reference: c2, rerun with the same binary.

The analysis writes `ab_ana.txt` (`fast/ab_ana.py <root> c2`) and `ab_cad.txt`
(`cadence/cad_ana.py`). Columns:
- "T rel" is relative to c2.
- The night side is chaotic: t4x against t4 had a night rms of 7.8e-3 (README_fast
  section 4).
- "kinks" is the number of columns with osc > 0.1 in the < 1e-7 bar band. c2 = 716 in
  this run; README_fast has t4 = 737.
- "nc" is the number of non-converged calls.
- "fref" is the mean fraction of columns refreshed per linearised step.
- edef and linerr are defined in section 1. |gap| is the mean |ckdesum| of the full
  calls.

| arm | N | thr | max day | rms day | rms night | <1e-7 day/night rms | kinks | nc | pass | fref | edef | linerr | \|gap\| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| c2 | 1 | - | 0 | 0 | 0 | 0 | 716 | 0 | 3.07 | - | - | - | 4.1e-5 |
| e2 | 2 | off | 5.2e-4 | 6.2e-6 | 6.1e-3 | 4.8e-6/6.0e-3 | 705 | 0 | 3.08 | 0 | 1.5e-5 | 2.3e-3 | 4.1e-5 |
| **e4** | 4 | off | 1.8e-3 | 2.0e-5 | 8.3e-3 | 3.5e-6/1.0e-2 | **733** | 0 | 4.05 | 0 | 3.5e-5 | 4.5e-3 | 4.1e-5 |
| e4g5 | 4 | 0.05 | 1.6e-3 | 1.9e-5 | 6.2e-3 | 5.2e-6/6.2e-3 | 707 | 0 | 3.06 | 0.22 | -1.8e-6 | 2.3e-3 | 1.9e-5 |
| e4g10 | 4 | 0.10 | 1.9e-3 | 2.0e-5 | 9.4e-3 | 1.3e-5/9.2e-3 | 724 | 0 | 3.02 | 0.012 | 3.5e-5 | 4.9e-3 | 2.2e-5 |
| e8 | 8 | off | 4.6e-3 | 4.7e-5 | 1.9e-2 | 2.5e-5/3.3e-2 | **851** | **166** | 7.91 | 0 | 8.2e-5 | 9.0e-3 | 4.0e-5 |
| e8g2 | 8 | 0.02 | 4.4e-3 | 4.4e-5 | 9.0e-3 | 4.1e-6/7.5e-3 | 719 | 0 | 3.08 | 0.42 | -1.7e-4 | 2.1e-3 | 1.7e-5 |
| e8g5 | 8 | 0.05 | 4.1e-3 | 4.5e-5 | 4.4e-3 | 1.1e-5/6.0e-3 | 711 | 0 | 3.05 | 0.24 | -2.5e-5 | 2.7e-3 | 1.5e-5 |
| e8g10 | 8 | 0.10 | 4.5e-3 | 4.7e-5 | 7.9e-3 | 1.2e-5/8.2e-3 | 725 | 0 | 3.01 | 0.030 | 5.2e-5 | 7.9e-3 | 1.8e-5 |
| e8g20 | 8 | 0.20 | 4.5e-3 | 4.7e-5 | 1.1e-2 | 1.6e-5/1.6e-2 | 772 | 61 | 3.85 | 0.010 | 7.1e-5 | 9.8e-3 | 2.1e-5 |
| e16 | 16 | off | 1.0e-2 | 1.0e-4 | 4.0e-2 | 7.5e-5/7.2e-2 | 710 (131 > 0.3) | 90 | 7.82 | 0 | 1.9e-4 | 1.8e-2 | 4.1e-5 |
| e16g2 | 16 | 0.02 | 7.3e-3 | 9.3e-5 | 7.4e-3 | 9.2e-6/7.3e-3 | 717 | 0 | 3.08 | 0.42 | -3.4e-4 | 2.2e-3 | 1.5e-5 |
| e16g5 | 16 | 0.05 | 9.1e-3 | 9.9e-5 | 7.9e-3 | 1.2e-5/6.8e-3 | 713 | 0 | 3.06 | 0.25 | -1.0e-4 | 3.2e-3 | 1.4e-5 |
| e16g20 | 16 | 0.20 | 9.8e-3 | 1.0e-4 | 1.7e-2 | 2.7e-5/2.4e-2 | 793 | 1 | 3.24 | 0.014 | 1.2e-4 | 1.5e-2 | 2.0e-5 |

The arm set was widened from the brief's thr in {0.02, 0.05} to 0.10 and 0.20, after
the cost job showed that 0.02 and 0.05 refresh 22-42% of the columns on every step.

Energy:
- The linear steps deposit Q0 dt to within |edef| <= 3.5e-5 (e4) to 3.4e-4 (e16g2) of
  sum |Q0| dt.
- The linear model's error against the true source at the next full call (linerr) grows
  with the time since the refresh: 2.3e-3 (N = 2), 4.5e-3 (4), 9.0e-3 (8) and 1.8e-2
  (16).
- The full calls close to the Newton tolerance, as in c2 (|gap| 4.1e-5).
- The capped count is 0 in every arm except e16 (5090 cells).

## 5. Gate C: cost (GPU, apudev, 2 ranks, interleaved fwd/rev/fwd, median of 3)

Jobs 11956036 (t1) and 11956380 (t2), same binary. Each arm runs t0 + 9000 s (480
cycles, 780 for h). The analysis is `cadence/ana_time.py time/<tag>`, written to
`time_t1.txt` and `time_t2.txt`. RT = arm - h.

| arm | ms/cycle (t1 / t2) | RT ms/cycle | RT x hydro |
|---|---|---|---|
| h | 14.82 / 14.78 | - | - |
| c2 | 46.64 / 46.67 | 31.9 | 2.16 |
| e2 | - / 39.08 | 24.3 | 1.64 |
| **e4** | 28.83 / 28.64 | **13.9** | **0.94** |
| e8 (fails accuracy) | 24.71 / - | 9.9 | 0.67 |
| e4g10 | - / 38.07 | 23.3 | 1.58 |
| e8g10 | - / 39.05 | 24.3 | 1.64 |
| e8g5 | 50.29 / - | 35.5 | 2.40 |
| e16g5 | 49.52 / - | 34.7 | 2.35 |

**Why the guard does not pay.** A guard call costs about as much as a full call, even
when it refreshes 1-3% of the columns:
- e8g10 has guard calls on 1319 of 1453 linearised steps and costs 39 ms/cycle, against
  24.7 ms for e8.
- The fixed per-call cost does not scale with the number of live columns: rt_pre_tp /
  geom / cut over all cells, kernel launches, and the store sweep, which every partial
  call does.

## 6. Verdict

Lag does not matter at N = 4:
- e4 keeps the kinks within c2 + 2.4% (733 against 716; t4 737) with 0 non-converged
  calls.
- The day-side error is 2.0e-5 rms, far below lever 1's 3e-4 relative to t4.
- The RT cost drops from 2.16x to **0.94x hydro**.

N = 8 without a guard fails: 166 non-converged calls and kinks +19%. Each refresh after
8 stale steps jolts the Newton to 8 passes.

The guard fixes N = 8 and N = 16 for accuracy (thr 0.05-0.10). It costs more than it
saves, because a partial call is almost a full call.

**Recommended: N = 4, thr off** (`ck_impl_every = 4` on top of c2).

## 7. HANDOVER

Commit on ck-cadence: code, scripts and this README, **not merged**.

Before merging:
- Verify the diff.
- Note that the restart is not bitwise, for c2 as well as for the cadence (section 3).

Next steps, in order:
1. Make a guard call cheap:
   - Let a partial call re-apply the stored operator instead of storing.
   - Restrict rt_pre_tp / geom / cut and the launches to the flagged columns, or batch
     the flagged columns into a compact index list.
   - Then N = 8 + thr 0.10 (accurate, 0 non-converged) could reach about 0.7x hydro.
2. Save c2's cross-call state (ck_thk / ck_thu, the xstep operator) together with the
   cadence arrays in the restart, if a bitwise restart is wanted.
3. A production A/B longer than 0.1 rotation for e4 (a few rotations) before switching
   production to it.

Cleanup: the build directories (src_*) are deleted and the binaries are kept in
/viper/ptmp2/jinma/ckcad_0924.
