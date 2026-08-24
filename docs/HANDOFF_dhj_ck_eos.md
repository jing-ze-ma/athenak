# Handoff: deep_hot_jupiter_rt, tabulated EOS + correlated-k

Written 2026-08-24 from the viper session, for whoever picks this up elsewhere
(e.g. a session on orion looking at a tabulated-EOS + correlated-k run).
Branch `polar-average-perf`, HEAD **c1d513c2**, pushed to `fork`
(`git@github.com:jing-ze-ma/athenak.git`). `src/` is clean; nothing is running.

**Read this before proposing an explanation for a blow-up.** Nine hypotheses have
already been tested and refuted with evidence. Six of them were the assistant's own
and three died the same day they were proposed.

---

## 1. What shipped recently

| commit | what |
|---|---|
| `ea823e3a` | RT source limiter, `problem/rt_de_max` — caps the explicit operator-split radiative energy change per cell per step |
| `c1d513c2` | two diagnostics: `problem/ad_dump_file` (initial p, T, nabla, grad_ad, 0.9*grad_ad, before and after `adjust_ad_pT_arr`) and Gamma_1 + grad_ad appended as columns 6–7 of the existing `ck_dump_file` column dump (existing parsers unaffected) |
| `ba2f0943` | general-EOS runs with `tfloor_kelvin` can be restarted at all (they previously could not) |
| `0d1f6f6a` | fixed an out-of-bounds read: `wtemp` is general-EOS-only and the RT read it unconditionally |
| `3882e37f` | polar MPI host-mirror bug — every multi-rank polar run died on cycle 1 from 2026-08-19 until this landed |

Earlier: `38311a8a` is the correlated-k chain-parallel split kernel, still the live
RT architecture. Correlated-k itself is finished and documented in
`docs/correlated_k_rt.md`, with regression tests in `tst/test_suite/rad/`.

## 2. There are TWO blow-ups, not one

### Failure one — ideal EOS + correlated-k — **CLOSED, do not reopen**

NaN at t ≈ 5.2e3 s (0.017 rotations). Cause: the RT source is explicit and
operator-split with **no radiative timestep constraint**, and `pfloor = 1.0` pinned
the internal energy at 2.115 erg/cm^3 at the top so a single step could overshoot.
Fixed by the `rt_de_max` limiter (`ea823e3a`) and A/B verified.

### Failure two — tabulated/general EOS — **OPEN**

Dies under any RT: general+grey at 2.509 rot, general+correlated-k at 0.37–0.67 rot
(1.07 with RKG). Rotation = 3.05e5 s. `eos_fail = 0`, `c2p_it = 0`, `vceil = 0`,
`fofc = 0` in every run — **floors are the symptom, not the disease.**

The closing 2x2, all at the `xe_long` grid and floors:

| EOS | RT | outcome |
|---|---|---|
| ideal | grey | **never fails** — clean to tlim 3.934 rot, mass to 0.15 % |
| general | grey | dies 2.51 rot |
| general | correlated-k | dies 0.37–0.58 rot |
| ideal | correlated-k | dies 0.017 rot — this one is now FIXED, see failure one |

## 3. Refuted — do NOT re-propose without new evidence

1. **Temperature floor.** Backwards: lowering tfloor 200 → 35 K failed 4x sooner and
   as a hard NaN. The floor was protecting the run.
2. **Resistivity cap.** Both `max_eta` 1e13 and 1e14 fail.
3. **RKG super-time-stepping.** A ~2x delay, not a cure.
4. **Radial resolution.** `nx1` 64 → 128 failed at 0.395 rot. Kills the
   "marginally resolved instability" reading.
5. **Metal condensation.** `eos_metal_condensation = false` failed at 0.932 rot.
6. **Floor choice.** Shipped floors (pfloor 1e-1, dfloor 5.44e-12) failed at 0.470 rot.
7. **Convectively unstable initial condition.** The new `ad_dump_file` shows ZERO
   genuinely super-adiabatic levels after `adjust_ad_pT_arr` (1 %-of-grad_ad
   tolerance), H2 on or off. It does not start unstable; it develops the front.
8. **Cells per scale height.** Measured per column from the bin dumps: the run that
   NEVER fails (ideal+grey) is the *most* sub-scale-height of the four (2.06 % of
   (col,level) pairs vs 1.87 % for the fastest-dying one). Not the discriminator.
9. **RT source limiter as a cure.** It survives to tlim but loses 67 % of its mass at
   0.240 rot and the density floor refills it — not a physical state.

`eos_h2 = false` runs clean to tlim under BOTH grey and correlated-k — but that is a
much bigger change than "H2 off": grad_ad's minimum goes 0.0914 → 0.2687 and the
initial interior goes 6728 K → 14098 K at 300 bar. A factor-of-two hotter, differently
stratified planet, and also the only perfectly scale-height-resolved case. Strongest
signal we have, weak as evidence that dissociation is the *mechanism*.

## 4. The one lead left

Gamma_1. It is identically 1.4728 for the ideal gas; for the table it spans 1.084–1.666
and jumps by up to **0.58 between ADJACENT cells** (T 1963 → 4006 K), already in the
healthy state at t = 2e5. Measured with the `c1d513c2` dump at r/r0 = 1.26–1.33,
p ~ 1e-6..1e-5 bar. A strong correlate, **not a proven cause** — it is entangled with
the floored top, so it does not cleanly separate "steep Gamma_1 breaks the
reconstruction" from "the floored top is chaotic and Gamma_1 merely reports it".

The only clean test, which **holds diffusivity fixed**:

> reconstruct pressure and evaluate Gamma_1 AT THE INTERFACE from the reconstructed
> state, instead of reconstructing the stored cell-centred value. Same order, same
> dissipation, so a change in outcome can only come from the Gamma_1 handling.

**Donor-cell reconstruction is not a valid test** — it changes the scheme's
diffusivity, and dissipation alone is already known to just delay the failure
(failure times 0.370 < 0.485 < 0.584 < 1.072 track increasing dissipation). The user
rejected it on exactly that ground.

Related, recorded but **fenced off**: `src/mhd/mhd.hpp:588-596` documents that under a
general EOS the Riemann solvers do not recompute pressure from the reconstructed (d,e)
— they read it from `wder` — and states *"Gamma_1 is left alone: it is a smooth O(1)
quantity that is not stratified over many scale heights, so plain PLM is appropriate."*
The measurement contradicts that assumption. **The user said "let's not go that way"
on 2026-08-23; do not open it without them raising it first.**

## 5. Bars to clear, and traps

- Fast reproducer: general + correlated-k, input `blowup.athinput`, dies at
  0.37–0.67 rot ≈ 36 min of apu1 wall.
- **Set `tlim` past 2.51 rot or the result is unreadable.** General+grey dies at
  2.509; a short tlim already produced one wrong call this campaign. Clean means
  reaching `tlim = 1.2e6` = 3.934 rot.
- On a restart, do **not** override `time/nlim` — a restart resumes at a large cycle
  number, so a small nlim runs zero cycles and the one-shot dumps never fire.
- GPU runs need `export HSA_XNACK=1` or they die. On viper, GPU scratch was
  `/viper/ptmp/jinma/claude_eos_gpu` (session /tmp is node-local and invisible to
  compute nodes).
- Don't design an experiment whose result can't discriminate between the live
  hypotheses — the user vetoes those, and calls a halt when hypotheses churn without
  converging.

## 6. The exact input

Start from the shipped **`inputs/mhd/deep_hot_jupiter_rt_eos.athinput`** (already in this
repo at `c1d513c2`) and apply the deltas below. That reconstructs the fast reproducer
exactly; on viper the live copy was `/viper/ptmp/jinma/claude_eos_gpu/blowup.athinput`.

| block / param | shipped | reproducer | why |
|---|---|---|---|
| `mesh/x1max` | `12.54e9` | `13.04e9` | the `xe_long` radial extent |
| `meshblock/nx2`, `nx3` | `8`, `8` | `16`, `16` | 32 blocks on one GPU — see the decomposition note below |
| `mhd/max_eta` | `1.0e14` | `1.0e13` | the long-run resistivity cap |
| `mhd/use_rkg_sts` | `true` | `false` | RKG is required at 1e14, optional at 1e13; leaving it off is what makes this the *fast* reproducer |
| `mhd/pfloor` | `1.0e-1` | `1.0e0` | the `xe_long` floors |
| `mhd/dfloor` | `5.44e-12` | `7.26e-12` | " |
| `problem/bbot` | `5.0e0` | `1.0e1` | " |
| `problem/rt_ck` | `false` | `true` | correlated-k on |
| `problem/ck_table`, `ck_data_dir` | relative to `data/exo_fms_ck` | absolute | only because the job ran from a scratch dir; make them resolve, however |
| `output1/dt` | `8.64e6` | `2.0e2` | fine cadence — it dies inside 0.67 rot |
| `output2/dt` | `2.16e6` | `2.0e3` | " |
| — | — | add `<output3> file_type=rst, dt=1.0e4` and `<output4> file_type=log, dt=2.0e2` | restarts + the event log |

**Watch `rt_de_max`.** The live `blowup.athinput` carries `rt_de_max = 0.5`, i.e. the
limiter ON — that is the `lt.limon` variant, which survives to tlim but only by losing
67 % of its mass and having the density floor refill it. **The baseline that dies at
0.674 rot has the limiter off** (`rt_de_max` absent, or `<= 0`). Set it deliberately;
do not inherit it by accident.

The other variants were all one-line edits of that same file, which is worth reproducing
as a habit — it is what made the results comparable:

| variant | the single change | outcome |
|---|---|---|
| `grey` | `rt_ck = false` | dies 2.509 rot |
| `v_noh2` | `eos_h2 = false` | clean to tlim |
| `v_noh2grey` | `eos_h2 = false` + `rt_ck = false` | clean to tlim |
| `v_hr` | `mesh/nx1` and `meshblock/nx1` 64 → 128 | dies 0.395 rot |
| `v_nocond` | `eos_metal_condensation = false` | dies 0.932 rot |
| `v_shipfl` | `pfloor 1.0e-1`, `dfloor 5.44e-12` | dies 0.470 rot |
| `tfloor` | `tfloor_kelvin` 200 → 35 | fails 4x SOONER, hard NaN |
| `ideal_ck` / `ideal_grey` | separate ideal-gas files, differing from each other in `rt_ck` ALONE | see failure one |

The ideal-gas pair is derived from `inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput`.

Also relevant to sizing a job: 32 meshblocks on one GPU beats 2 by 1.35–1.45x, and a
second APU buys only ~1 %; `nx1` is capped at 264 by MHD LDS (production uses 256), and
x1 cannot be split because the RT is a column solve.

## 7. Repo conventions

- Push to `fork`, never `origin` (origin is upstream IAS-Astrophysics and is denied).
- `run/` is ~165 GB of output: read-only, and never `git add -A` in this repo.
