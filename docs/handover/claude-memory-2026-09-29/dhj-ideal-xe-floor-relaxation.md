---
name: dhj-ideal-xe-floor-relaxation
description: "CLOSED 2026-08-17, committed as 6003c1ce: floors, max_eta and STS all settled by measurement for the ideal-gas resistive run; docs/ideal_gas_resistive.md + inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput carry it."
metadata: 
  node_type: memory
  type: project
  originSessionId: adaab265-737e-41c3-9d7d-e43d74682019
  modified: 2026-08-17T02:41:27.045Z
---

Answered 2026-08-17. The user asked: run the **ideal gas with the new EOS resistivity**,
and can the floors be relaxed even at **bbot = 10 G**? Four runs, one rotation period each.

## SHIPPED as `6003c1ce` (pushed)

`docs/ideal_gas_resistive.md` and `inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput` now
carry every number below, so this note is the working record and those two files are the
handover channel. Recommended config: `ohmic_resistivity = eos`, **`max_eta = 1e13`**,
**`use_rkg_sts = false`**, floors unchanged at 5.44e-12 / 1e-1 / 9.186e9.

Whole-rotation wall times, 16x7, t = 3.05e5, all finite, mass and tot-E to 5-6 digits:
1e14+STS 6347.7 s | pfloor 1e-3x 6257.5 s (0.986) | 1e14 no STS 5877.8 s (0.926) |
**1e13 no STS 4129.4 s (0.651)**.

Two things that overturned earlier beliefs:
- **STS is not a standing win at max_eta = 1e14.** The sign FLIPS during the run as dt_hydro
  decays: 1.12x ahead at t = 7.5e4, level at 1.3e5, 1.34x BEHIND at 2.5e5. The
  "REQUIRED at 1e14" line in the general-EOS input came from a constant-eta test at
  bbot = 3 G where tau was 29. Left that input alone -- the general EOS runs COLDER and
  needs its own measurement (its `perna`+1e14+STS run stalled outright).
- **`max_eta` 1e14 -> 1e13 is nearly free in accuracy.** Measured against the uncapped Saha
  eta over a real snapshot, not argued from x_e tables: 1e13 caps 6.1% of cells holding
  0.87% of ME, 1e14 caps 3.3% holding 0.26%, 1e12 caps 9.5% holding 2.68%. The older claim
  in [[resistivity-perna-uhj]] that 1e12 "overrides most of the domain" overstates it.

## THE ANSWER

| run | dfloor | outcome | wall / simulated second |
|---|---|---|---|
| `base` | 5.44e-12 | completed 1 rotation, healthy | 0.0208 |
| `d3` | 5.44e-15 | healthy, physics unchanged | 0.1225 (**5.9x**) |
| `d6` | 5.44e-18 | **NaN at t ~ 7.4e4**, whole domain at once | -- |
| `all3` | all three floors 1e-3x | unusable, ~500x slower; cancelled | -- |

**dfloor x1000 is SAFE and changes nothing that matters.** Matched-time structure compare
(t = 1.22e5, base vs d3): rho relative difference 1.5e-5..3.3e-4 for shells i <= 32,
3e-3..4.6e-2 at i = 40-48, 2.3e-2 at i >= 56. Shell magnetic energy ratio 0.997-1.003 below
i = 48 and **0.94 above i = 56**. Mass and tot-E agree to 6 digits over the whole run;
integrated MEs to <= 3%. So it moves the outer ~8 radial shells of 64 by ~6% in magnetic
energy and nothing below. Velocity fields decorrelate chaotically -- relative RMS is
meaningless there because the mean flow is ~0.

**The floor lives only in the top 12% of the domain in radius.** No cell below i = 48 ever
touches it; 37% of cells in the outermost shell sit on it. Relaxing it lets those fall to
~1.7e-14 instead of 5.44e-12, and max v_A there goes 2.5e6 -> 1.8e7 cm/s. THAT is the 5.9x,
and it is an honest Alfven CFL, not an instability.

**Why x1e6 dies.** At 5.44e-18 the floor is BELOW the density the gas actually reaches
(rho_min ~2e-14), so it fires only 3.6e-6 of the time -- it has stopped being a bulk floor
and become a rare-event catcher. Those rare events are the whole point: one cell can then
reach v_A ~1e9 cm/s inside a single step. d6 and d3 were identical to 3-4 digits at
t = 6.1e4; by the next dump every cell in d6 was non-finite. The failure is stochastic, not
a clean function of the floor value.

**pfloor alone is harmless AND pointless** (measured, run `p3`): full rotation in 0.986x the
wall time, same mass and tot-E, and the flooring simply MOVES -- efloor 7.97% -> 2.12%,
tfloor 0.44% -> 6.95%. Because `dfloor*Rgas*tfloor = 0.050` already bounds p and pfloor is
0.1, it only has teeth in that factor-2 window. **tfloor is the one that must not move.**
Relaxing pfloor+tfloor together with dfloor at BASELINE (run `t3`) went NaN at t = 1.2e4 in
95 s -- so the thermodynamic floors kill it on their own, not via dfloor.

**Do NOT relax pfloor/tfloor alongside.** This is the trap specific to
`ohmic_resistivity = eos`: cooler gas is MORE resistive, eta pins at max_eta, dt_diff
collapses, and the RKG stage count explodes (258 sub-stages). `all3` advanced 3053 s of
simulated time in 7 HOURS against base 3.05e5 in 1.8 h, with dt itself still 5.3 s -- the
cost is entirely per-cycle stage count.

**RECOMMENDATION given to the user:** keep dfloor = 5.44e-12 for production. If the
upper-atmosphere field structure is ever the science target, 5.44e-15 is safe and validated,
but only with the BASELINE pfloor and tfloor.

## Where the runs are

`/orion/u/jinma/ATHENAK/idflr/` (NOT in [[run-directory-untouchable]]):

- `dhj_ideal_xe.athinput` — the ideal-gas variant of the tracked
  `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`. `eos = ideal`, `ohmic_resistivity = eos`,
  `max_eta = 1e14`, `use_rkg_sts = true`, `bbot = 1.0e1`, mesh 64 x 56 x 128.
  tlim = 3.05e5 (**one rotation**, 2pi/omega), hst+log every 3.05e3, **bin every 3.05e4 =
  11 snapshots**. The user asked whether the cadence was the reference input's 2.16e6: it is
  NOT, and cannot be -- dt has fallen to ~1.5 s (base) and ~0.28 s (relaxed), so 2.16e6
  would take days per snapshot.
- `run_one.sh` — `sbatch -o <dir>/slurm.out -e <dir>/slurm.err run_one.sh <tag> <dfloor>
  <pfloor> <tfloor>`, p.shared, 16 x 7, 12 h.
- Four dirs, jobs 190247-190250: `base` (5.44e-12 / 1e-1 / 9.186e9), `d3` (dfloor 1e-3x),
  `d6` (dfloor 1e-6x), `all3` (all three floors 1e-3x).
- `short_{base,relax3,relax6}` — the first, too-short 3000-cycle attempt. Superseded; the
  user correctly said to run long enough for >= 10 snapshots.

**tfloor is stated in CODE units here** (p/d = T_K * Rgas = 200 K * 4.593e7 = 9.186e9), NOT
via `tfloor_kelvin`. `tfloor_kelvin` converts with `<units>/mu`, and this input has mu = 1
while the run's Rgas = 4.593e7 implies mu = 1.797 -- so it would be 1.8x off. Either set
`<units>/mu = 1.797` or use `tfloor` as done here.

The x_e table needed `eos_logd_min = -20` (was -14) to bracket the relaxed dfloor.
Table 401 x 451, 22 MB.

## The mechanism, for whoever asks "why" again

**The floors are NOT inert at 10 G.** Baseline fires dfloor in ~7% of all C2P calls
(4.1e8 hits over 221210 cycles) and efloor in ~5%. They are not a mere nan-guard here --
but see above: the cells they touch carry no mass, energy or field structure.

**Two different constraints are tangled up in the word "floors":**

1. *dfloor vs the Alfven CFL* -- the B-coupled one. v_A = B/sqrt(4 pi rho) at 10 G:
   rho 4.0e-10 (real outer shell) -> 1.4e5; 5.44e-12 (dfloor) -> 1.2e6;
   5.44e-15 -> 3.8e7; 5.44e-18 -> 1.2e9 cm/s. Sound speed is only 4.1e5 at 2500 K and the
   run's total signal speed ~6.8e5. **So the baseline dfloor already sits where v_A ~ 2 c_s**
   -- the floor is what bounds the Alfven speed, and relaxing it converts straight into dt.
   The old STABILITY reason for a high dfloor is gone ([[dhj-highB-outer-bc]] fixed the
   outer-BC Maxwell term); this second, genuinely physical reason is untouched by it.
2. *pfloor/tfloor vs thermal collapse* -- nothing to do with B, and fatal to the runtime
   rather than to the solution. See the `all3` row above.

**A wrong intermediate conclusion, recorded so it is not re-derived:** at 21 minutes d3 and
d6 had the same dt (0.29 vs 0.28 s) and it looked like the dfloor penalty SATURATED, so that
extra decades were free. That is false -- d6 went NaN later. Equal dt does not mean equal
robustness: the two differ only in the worst-case single-step Alfven speed a rare floor
activation permits, which does not show up in the mean timestep at all.

## Tools and gotchas worth keeping

- The `log` output (`file_type = log`, `EventLogOutput`) IS the floor diagnostic: columns
  eos_dfloor / eos_efloor / eos_tfloor / eos_vceil / eos_fail / c2p_it / fofc. Cumulative,
  and written only when some counter is non-zero. Set its `dt` small or the only row you get
  is the last cycle.
- Relaxing dfloor MOVES the flooring rather than removing it: dfloor hits fall ~400x
  (4.1e8 -> 1.0e6) while efloor hits rise ~4x (2.3e8 -> 1.0e9). Something still catches the
  tenuous gas.
- d3 hit `nlim = 600000` at t = 1.39e5, i.e. 46% of a rotation. Budget ~1.3e6 cycles for a
  full rotation at that floor.
- **A banner comparison that is NOT evidence of anything** (recorded because it looks like
  it is): this ideal run reports n_e/n_tot = 3.6e-3 "at the grid centre" against the general
  EOS reference's 1.5e-4. That is NOT the mu/T mismatch -- it is because this input sets
  `eos_logd_min = -20` instead of -14, so the TABLE GRID CENTRE moved from log rho = -7 to
  -10. x_e goes as n^-1/2, and 10^1.5 = 32x accounts for the whole factor. The banner
  evaluates the table at its own grid centre, not at any state the run occupies.
  The real mu/T mismatch is separate and does hold: Rgas = 4.593e7 -> mu = 1.810 against the
  table's ~1.26, so the ideal T is ~1.44x too high, and with d ln x_e/d ln T ~ 13 that is
  ~2 orders of magnitude in x_e. See [[resistivity-xe-table-for-ideal]].
- Comparison scripts were in the session scratchpad (gone). They are ~30 lines:
  `sys.path.insert(0,'athenak/run'); import bin_convert; d=bin_convert.read_binary(f)`,
  then `d['mb_data']` is a dict of (nmb, nz, ny, nx) arrays keyed
  dens/velx/vely/velz/eint/bcc1/bcc2/bcc3. nx1_mb = 64 = the full radial extent, so the
  last index IS the radial shell and no block bookkeeping is needed.

## STILL OPEN on this thread

- Nothing forced. Candidates: restore a real `min_xe` in the code (the max_eta cap is a
  sqrt(T)-dependent x_e floor, see [[resistivity-perna-uhj]]); run the GPU port on Viper
  (d52f3610 is compile-time + CPU-bitwise only, no kernel has executed on an accelerator);
  Stage 4 (rho,e) table, still deferred by the user.

Related: [[resistivity-xe-table-for-ideal]], [[dhj-highB-outer-bc]], [[resistivity-perna-uhj]],
[[freya-job-submission]].
