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
runs and needs re-running on a fixed binary before it can be trusted.

**That scatter was the race, not chaos.** The spread in (a) was produced by the bug in (c):
identical inputs gave different answers. With `6e600f12` in place, identical inputs give
bitwise identical answers again, so single runs are meaningful once more — subject only to
whatever genuine physical chaos the problem has, which is now measurable instead of being
confounded with a bug. Every arm in sections 3 and 4 is being re-run on the fixed binary.

**(b) The orion run DID blow up.** `or.gpu`, orion's literal input replayed on a viper
GPU, passed 0.690 rot clean and then **died at t = 4.464e5 = 1.464 rotations** — mass
1.0002 → 0.105 inside ONE 200 s history interval, from dt = 1.34 s, with no precursor.
Orion stopped watching at 0.713 rot, which is the only reason it looked like a
disagreement. **The commit message on `cd2e8815` ("which did NOT blow up") is therefore
wrong**, and CPU-vs-GPU is NOT the explanation — it reproduces on viper hardware.

**(c) The nondeterminism is FOUND AND FIXED — `6e600f12`.** It was a missing
`member.team_barrier()` in the general-EOS x1 flux kernel, core AthenaK code, not the
pgen. **Rebuild before running anything from this document**; every result below that
predates `6e600f12` was produced by a raced binary. See section 4.

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

### The nondeterminism — SOLVED 2026-08-24, commit `6e600f12`

**The bug.** `src/hydro/hydro_fluxes.cpp` and `src/mhd/mhd_fluxes.cpp`, general-EOS branch
of the **x1** direction:

```
PiecewiseLinearX1(member, m, k, j, il-1, iu, wder_, dl, dr);   // thread i writes dl(n,i+1)
par_for_inner(member, il, iu, [&](const int i) {
  dl(IDPR,i) = fmax(dl(IDPR,i), eos_.pfloor);                  // thread i reads dl(n,i)
  ...
});
...
member.team_barrier();      // the only barrier — and it came AFTER the floor
```

The reconstruction writes `dl(n,i+1)` from the thread that owns `i`, so the floor reads a
scratch slot **another thread wrote**, and successive `par_for_inner` loops are not
implicitly synchronised. Whether the floor saw the reconstructed value or a stale one was
scheduling-dependent. x2 and x3 are unaffected — those reconstructions write index `i` from
thread `i`, which is exactly why only the x1 flux ever differed.

**Scope:** general EOS only (the floor lives inside `if (nder > 0)`, and `nder` is 0 for an
ideal gas), and device only (a serial `par_for_inner` runs in order). The fix is one
`member.team_barrier()` in each file; it costs **1.0 %** (2.363e7 → 2.339e7 zone-cycles/s).

**Consequence beyond reproducibility:** a few interfaces per column silently kept an
UNFLOORED pressure. That is the mechanism behind the scatter in section 0(a).

**How it was found**, after configuration-bisection stalled: per-task checksums. A hook in
`TaskList::DoAvailable` printing the wrapping integer sum of the raw bit patterns of `u0`,
`w0`, `wder`, `wtemp` and the three flux arrays after every completed task — exact and
order-independent, unlike a floating-point sum which can cancel a difference away. Ten
replicate runs grouped by the hash of their checksum stream: 8 agreed, 2 did not, and the
first differing line named it — `task=1` (`Hydro::Fluxes`) with `u`, `w`, `wder`, `wtemp`,
`f2`, `f3` all bitwise identical and **only the x1 flux different**.

**Verification** (clean Release build, 3 replicates, 600 cycles, bitwise on every dump and
every history row): hydro+polar with no source terms, MHD+correlated-k+polar (production),
and hydro with polar off and sources on — all three diverged before the fix (at cycles 20,
40 and 120) and are all bitwise identical after it.

**Regression test:** `tst/test_suite/nr/test_nr_geneos_repro_gpu.py`. Runs the same binary
three times and requires identical histories; validated FAILING without the fix and passing
with it. **Do not shrink its grid** — detection is very sensitive: the production grid gave
3 of 3 replicate pairs differing, while 32x32, 16x16 and a general-EOS linear wave gave
0 of 3. The race needs the pressure floor to bind on SOME interfaces but not all, and
enough resident teams for the two inner loops to overlap.

**CORRECTION to the earlier version of this section.** It claimed "at least TWO independent
races", because polar-off and `-DDHJ_POLE_PACK=0` builds looked clean. They were not clean
— the same single bug was present and merely took longer to grow into the observable. One
fix makes every configuration reproducible, and the `do_pole` A/B was a false lead:
disabling that branch changed how fast a perturbation amplified, not whether one was
created. The exclusions listed there (pack kernel internally deterministic, no double
writes in the unpack, not an out-of-bounds access, single execution space) all still stand
and are still worth not redoing.

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
