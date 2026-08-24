# Handoff: deep_hot_jupiter_rt, tabulated EOS + correlated-k

Branch `polar-average-perf`. **Updated 2026-08-24 (evening).** An earlier version of this
file, written the same day at `397b49b1`, is superseded in three important ways — read
section 0 before anything else. Push to `fork`
(`git@github.com:jing-ze-ma/athenak.git`), never `origin`.

---

## 0. READ FIRST — three things that invalidate single-run reasoning

**(a) The baseline is a chaotic ensemble, not a repeatable experiment.** Death times on
the *identical* setup now read **0.370, 0.485, 0.584, 0.685, 1.464 rot — plus at least two
survivals past 0.69.** Members differing only in `output2/dt` (a physics-inert knob that
merely shifts the dt sequence) land all over that range. **No single run of this baseline
is evidence of anything**, in either direction. Every "arm X survived" and every "arm X
failed" claim below — including the refuted list in section 3 — was called from one or two
runs and needs ensemble treatment before it can be trusted.

**(b) The orion run DID blow up.** `or.gpu`, orion's literal input replayed on a viper
GPU, passed 0.690 rot clean and then **died at t = 4.464e5 = 1.464 rotations** — mass
1.0002 → 0.105 inside ONE 200 s history interval, from dt = 1.34 s, with no precursor.
Orion stopped watching at 0.713 rot, which is the only reason it looked like a
disagreement. **The commit message on `cd2e8815` ("which did NOT blow up") is therefore
wrong**, and CPU-vs-GPU is NOT the explanation — it reproduces on viper hardware.

**(c) The code is nondeterministic run to run, and that is now the prime suspect.**
Two byte-identical invocations of the same binary on one GPU rank diverge. `c39b794c`
fixed one race (the outer-x1 `bcc0` cross-thread read) but did not fix this. **Fix the
race before drawing any more physics conclusions** — see section 4.

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

Dies under any RT: general+grey at 2.509 rot, general+correlated-k anywhere in
0.37–1.46 rot, and sometimes not at all (see section 0a). Rotation = 3.05e5 s. `eos_fail = 0`, `c2p_it = 0`, `vceil = 0`,
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
9. **RT source limiter as a cure — REVISED.** `rt_de_max = 0.5` does clear tlim 3.934 rot
   with general EOS + correlated-k, against 0.685 rot with it off, and `eos_h2 = false`
   cures it independently under BOTH grey and c-k. But the limiter arm takes a 33 % mass
   transient getting there, and both act on the LAST link of the chain (the marginal
   explicit RT source), not the first — pure hydro survives with the limiter off entirely.
   Treat them as real but downstream cures.

`eos_h2 = false` runs clean to tlim under BOTH grey and correlated-k — but that is a
much bigger change than "H2 off": grad_ad's minimum goes 0.0914 → 0.2687 and the
initial interior goes 6728 K → 14098 K at 300 bar. A factor-of-two hotter, differently
stratified planet, and also the only perfectly scale-height-resolved case. Strongest
signal we have, weak as evidence that dissociation is the *mechanism*.

## 4. Where the real lead is now: races, and the fact that MHD is required

### The nondeterminism (REOPENED 2026-08-24)

Intermittent; correlated-k runs so far (2 of 3 pairs; grey 0 of 2 over 32 cycles). Bisected
with byte-identical invocation pairs compared bitwise on every dump and the history:

| config | cycles | result |
|---|---|---|
| MHD + c-k, polar ON (baseline) | 600 | diverges, first at cycle 40 |
| MHD + grey / RT off / outflow BCs / `bbot=1e-30`, polar ON | 600 | all diverge |
| **HYDRO, polar ON, `user_srcs=false`** | 600 | **diverges at cycle 20 — FASTEST REPRODUCER** |
| **MHD, polar OFF** (`use_polar_boundary=false`) | **3000** | **BITWISE IDENTICAL** |
| **HYDRO, polar ON, built `-DDHJ_POLE_PACK=0`** | 600 | **BITWISE IDENTICAL** |
| MHD, polar ON, srcs ON, `-DDHJ_POLE_PACK=0` | 600 | diverges → a SECOND source |
| HYDRO, polar OFF, srcs ON | 600 | diverges → the second source again |

So **at least TWO independent races**:

1. In the `do_pole` branch of `MeshBoundaryValuesCC::PackAndSendCC` (`bvals_cc.cpp`
   ~104–197, and presumably the `bvals_fc` twin). Disabling that branch alone makes an
   otherwise-identical, equally active hydro+polar run bitwise reproducible.
2. A second, **not yet isolated** — shows up as MHD+polar with the pack branch disabled,
   and as hydro+`user_srcs` with polar off. Candidates: `PolarAzimuthalAverageEr`
   (`mhd_corner_e.cpp:460+`) and `polar_local_sum_b` (`bfield_bcs.cpp:339`, not yet read).

**Already excluded, do not redo.** The pack kernel is internally deterministic (double-pack
comparison: 0 differing buffer entries over 30 cycles) and the unpack writes no ghost cell
twice (atomic write-counter, 4 cycles). Physical BCs never touch polar ghosts
(`hydro_bcs.cpp` has `default: break;` for x2). AthenaK uses a single execution space, so
kernels cannot overlap. Not an out-of-bounds access — a full Debug build (Kokkos bounds
checking) runs the fastest reproducer 40 cycles clean. The amplification confound is
excluded: RMS d(rho)/rho between dumps 0 and 30 is 5.149e+04 in the diverging run and in
BOTH reproducible ones, identical to four digits — the clean runs are exactly as active.

Both (1) and the double-pack result can only be true together if the race lives in a kernel
**downstream** that consumes polar ghost data, and is therefore worth hunting without
reference to the pole at all.

**Next instrument, designed but NOT built:** checksum `u0` after every task in the stage,
run twice, report the first task whose checksum differs. That names the guilty kernel
directly instead of bisecting configurations. Everything needed is in `driver.cpp`'s task
loop.

### The failure needs MHD — the strongest physics signal

| run | config | outcome |
|---|---|---|
| `hy.hoff` | **pure hydro**, c-k, limiter OFF | **cleared tlim 3.934 rot, mass +0.068 %** |
| `hy.hon` | pure hydro, c-k, limiter ON | cleared tlim, same mass to three digits |
| `mv.weakb` | `bbot` 10 G → **1 mG** | **cleared tlim, mass +0.034 %** |
| `mv.noeta` | resistivity REMOVED (ideal MHD, 10 G) | past 1.03 rot, healthy |
| `lt.limoff` | full MHD baseline | died 0.685 rot |

Removing the field cures it outright with the RT limiter OFF, and the hydro mass history is
flat, not merely survivable. **And it is not a dt effect**: median dt is 6.97 s for hydro
against 1.28 s for the baseline, so hydro takes ~5x LONGER steps and still lives. Something
magnetic manufactures the low-e jagged state that makes the RT source marginal. Subject to
caveat (a) — these are one run each.

### Gamma_1 — demoted, not dead

Gamma_1 is identically 1.4728 for the ideal gas but spans 1.084–1.666 for the table, jumping
up to **0.58 between ADJACENT cells**, already in the healthy state. Still the only measured
EOS-side discriminator, but it is entangled with the floored top and is now behind the race
and the MHD requirement in priority. If tested, the only clean test **holds diffusivity
fixed**: evaluate Gamma_1 at the interface from the reconstructed state rather than
reconstructing the stored cell-centred value. **Donor-cell is not a valid test** — it changes
the scheme's diffusivity, and dissipation alone is known to just delay the failure; the user
rejected it on exactly that ground.

Related, **fenced off**: `src/mhd/mhd.hpp:588-596` documents that under a general EOS the
Riemann solvers read pressure from `wder` rather than recomputing it, and asserts Gamma_1 is
"a smooth O(1) quantity". The measurement contradicts that. **The user said "let's not go
that way" on 2026-08-23; do not open it without them raising it first.**

## 5. Bars to clear, and traps

- Fast reproducer: general + correlated-k, now IN THE REPO as
  **`inputs/mhd/deep_hot_jupiter_rt_blowup.athinput`** (`cd2e8815`, pushed from orion).
  Dies anywhere in 0.37–1.46 rot, or not at all — budget for the long tail, and run an
  ensemble, not one job.
- **Set `tlim` past 2.51 rot or the result is unreadable**, and do not stop watching
  early — a short tlim and an early stop have now EACH produced a wrong call this
  campaign (the second was orion's, see section 0b). General+grey dies at 2.509. Clean means
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

**Easiest path: use `inputs/mhd/deep_hot_jupiter_rt_blowup.athinput`** (`cd2e8815`), which
is the literal reproducer, and just repoint `ck_table`/`ck_data_dir`. It differs from
viper's `/viper/ptmp/jinma/claude_eos_gpu/blowup.athinput` only in `tlim` (8.64e7 vs 1.2e6),
`rt_de_max` (-1.0 vs 0.5), the two table paths, and one commented-out line — i.e. they are
the same physics, and `cd2e8815` adds no source change.

The delta table below reconstructs the same thing from the shipped
`inputs/mhd/deep_hot_jupiter_rt_eos.athinput`, and is what documents WHY each knob is set.

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
