# Handoff: deep_hot_jupiter_rt, tabulated EOS + correlated-k

Branch `polar-average-perf`. **Updated 2026-08-25 — the campaign is CLOSED.** Push to
`fork` (`git@github.com:jing-ze-ma/athenak.git`), never `origin`.

---

## 0. READ FIRST — the "general-EOS blow-up" was the data race

There is no thermodynamic blow-up. The failure that this whole document was written to
chase was a **missing `member.team_barrier()` in the general-EOS x1 flux kernel**
(`6e600f12`, core AthenaK, not the pgen — section 4 has the anatomy). With that barrier in
place the unmodified baseline runs to `tlim` and conserves mass to 7e-4.

Six arms, each 8 h on one apu1 GPU, all on `athena_fixrel` = `6e600f12`, all from
`inputs/mhd/deep_hot_jupiter_rt_blowup.athinput` (= viper's `blowup.athinput`), all
`tlim = 1.2e6 s = 3.934 rotations` (rotation = 3.05e5 s):

| job | tag | the single change | pre-fix | **post-fix** |
|---|---|---|---|---|
| 11006656 | `rf.limoff` | **none — the control**, `rt_de_max = -1` | died 0.370–0.685 rot | **cleared 3.934 rot**, mass 1.00030 |
| 11006657 | `rf.limon` | `rt_de_max = 0.5` | cleared, via a 33 % mass transient | **cleared**, mass 1.00030, no transient |
| 11006658 | `rf.noh2` | `mhd/eos_h2 = false` | cleared | **cleared**, mass 1.00238 |
| 11006659 | `rf.hydro` | pure hydro (`hydro_ck.athinput`) | cleared | **cleared**, mass 1.00070 |
| 11006660 | `rf.weakb` | `problem/bbot = 1.0e-3` | cleared | **cleared**, mass 1.00030 |
| 11006661 | `rf.noeta` | resistivity removed | past 1.86 rot, unfinished | **cleared**, mass 1.00030 |

All six exited 0 at `time=1.200000e+06`. No NaN, no dt collapse anywhere (`limoff`:
dt min/median/final = 0.921 / 1.254 / 0.987 s), mass in [0.99992, 1.00068] for the whole
of `limoff`, and the only warning printed is the benign `Rgas`-vs-general-EOS notice at
`deep_hot_jupiter_rt.cpp:1401`. `limoff` and `limon` are now near-identical — the limiter
clipped **4 cells** over the entire run — so the limiter is effectively inert here.

**What this retracts.** Everything below that was measured before `6e600f12`:

- **The 0.370 / 0.485 / 0.584 / 0.685 / 1.464 rotation spread on byte-identical input was
  the bug, not chaos.** There is no chaotic ensemble to sample; identical inputs now give
  identical answers.
- **`rt_de_max` does NOT cure a general-EOS blow-up** — there was nothing to cure. The
  limiter remains the correct and verified fix for the *ideal* EOS + correlated-k NaN
  (section 2), which is a genuinely separate failure with a different mechanism.
- **`eos_h2 = false` is not a cure either**, and neither is weak field, no resistivity, or
  pure hydro. Those arms "survived" pre-fix because the race, which only bites when the
  general-EOS pressure floor binds on some x1 interfaces and not others, happened to grow
  slowly in them. Their apparent significance was an artefact.
- **"The failure needs MHD"** — retracted. Hydro and MHD both run clean.
- **The refuted list (section 3) is now moot**, not because those hypotheses were wrong but
  because there is no longer a phenomenon for them to explain. Do not resurrect any of them.
- **The commit message on `cd2e8815`** ("which did NOT blow up") turns out to have been
  right after all, for the wrong reason: that run's binary raced too, it just got lucky.

**What survives.** The race analysis and its regression test (section 4); the ideal+c-k
limiter fix (section 2); the run recipe and traps (sections 5–6).

**Open, and the only thing left of the physics thread:** `Gamma_1` spans 1.084–1.666 across
the table and jumps by up to 0.58 between adjacent cells, while `src/mhd/mhd.hpp:588-596`
asserts it is "a smooth O(1) quantity" when it lets the Riemann solvers read pressure from
`wder` instead of recomputing it. That is a latent accuracy question, no longer a
suspected-crash question, and **the user said "let's not go that way" on 2026-08-23 — do
not open it unless they raise it.**

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

## 2. The one real blow-up: ideal EOS + correlated-k — CLOSED

NaN at t ~ 5.2e3 s (0.017 rotations). Cause: the RT source is explicit and operator-split
with **no radiative timestep constraint**, and `pfloor = 1.0` pinned the internal energy at
2.115 erg/cm^3 at the top, so a single step could overshoot. Fixed by the `rt_de_max`
limiter (`ea823e3a`) and A/B verified. This is the only genuine blow-up the campaign found,
and it is unrelated to the EOS table.

The other three quadrants of the closing 2x2 all run clean to `tlim = 3.934 rot` on the
fixed binary: ideal+grey (never failed at any point), general+correlated-k (`rf.limoff`
above), and general+grey. The pre-fix numbers for those — "general+grey dies at 2.509 rot",
"general+c-k dies at 0.37-1.46 rot" — were the race.

## 3. Retracted hypotheses — do NOT resurrect

These were each proposed and killed with evidence while chasing the general-EOS "blow-up".
They are listed only so nobody re-derives them; **since the phenomenon itself was a bug,
none of them is a live question any more.**

Temperature floor (lowering `tfloor` 200 -> 35 K made it worse, so the floor was
protecting the run); the resistivity cap (`max_eta` 1e13 and 1e14 both "failed"); RKG
super-time-stepping (a ~2x delay); radial resolution (`nx1` 64 -> 128 "failed" sooner);
metal condensation; the choice of floors; a convectively unstable initial condition (the
`ad_dump_file` diagnostic shows zero genuinely super-adiabatic levels after
`adjust_ad_pT_arr`, H2 on or off — this measurement is still valid and still useful);
cells per scale height (the never-failing configuration was the *most* sub-scale-height of
the four).

Two claims from the pre-fix era were themselves wrong and are retracted outright: that the
RT source limiter cures the general-EOS case, and that the failure "needs MHD". Both were
one-run-per-arm calls on a raced binary. See section 0.

Also retracted: "denser bin output changed the dt sequence and hence the failure time".
`Mesh::NewTimeStep` (`mesh.cpp:637-716`) clips dt for `tlim` and nothing else, and outputs
fire after the step and are read-only. Output cadence cannot move the dt sequence; what
moved the failure time was the race.

## 4. The race: what it was, how it was found, and what it retracts

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

### MHD is NOT required — retracted

The pre-fix table here read "pure hydro, weak field and no-resistivity all clear `tlim`,
full MHD dies at 0.685 rot", and concluded that something magnetic manufactured the state
that made the RT source marginal. On the fixed binary **full MHD clears `tlim` too**
(`rf.limoff`), so there is no contrast left to explain. The arms differ now only in how
fast they run — median dt 1.25 s for full MHD against 6.97 s for hydro and 7.61 s for the
1 mG field, which is ordinary Alfvenic/resistive timestep behaviour, not a stability signal.

### Gamma_1 — still measured, no longer a suspect

Gamma_1 is identically 1.4728 for the ideal gas but spans 1.084-1.666 for the table, jumping
by up to **0.58 between ADJACENT cells**. That measurement stands. What has changed is what
it might explain: nothing crashes any more, so it is at most an accuracy question about
`src/mhd/mhd.hpp:588-596`, which under a general EOS lets the Riemann solvers read pressure
from `wder` rather than recomputing it, on the stated assumption that Gamma_1 is "a smooth
O(1) quantity". The measurement contradicts that assumption.

**The user said "let's not go that way" on 2026-08-23; do not open it without them raising
it first.** If it is ever opened, the only clean test **holds diffusivity fixed**: evaluate
Gamma_1 at the interface from the reconstructed state rather than reconstructing the stored
cell-centred value. Donor-cell is not a valid test — it changes the scheme's diffusivity —
and the user rejected it on exactly that ground.

## 5. Bars to clear, and traps

- **`inputs/mhd/deep_hot_jupiter_rt_blowup.athinput`** (`cd2e8815`) is kept as the
  *reference general-EOS + correlated-k run*, not as a reproducer of anything — nothing to
  reproduce survives. Repoint `ck_table`/`ck_data_dir` and it runs; 8 h on one apu1 GPU
  reaches `tlim = 1.2e6` = 3.934 rot (`rf.limoff`: 2.17e4 s of GPU time, 2.32e7
  zone-cycles/s).
- Clean means reaching `tlim = 1.2e6` = 3.934 rot with mass flat. Do not stop watching
  early: a short `tlim` and an early stop each produced a wrong call during this campaign.
- **Rebuild before comparing anything to a pre-`6e600f12` number.** Every result recorded
  before that commit came off a raced binary, and none of them is comparable with a fixed
  one — that is what invalidated an entire campaign's worth of arms.
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
is the reference run, and just repoint `ck_table`/`ck_data_dir`. It differs from
viper's `/viper/ptmp/jinma/claude_eos_gpu/blowup.athinput` only in `tlim` (8.64e7 vs 1.2e6),
`rt_de_max` (-1.0 vs 0.5), the two table paths, and one commented-out line — i.e. they are
the same physics, and `cd2e8815` adds no source change.

The delta table below reconstructs the same thing from the shipped
`inputs/mhd/deep_hot_jupiter_rt_eos.athinput`, and is what documents WHY each knob is set.

| block / param | shipped | reference run | why |
|---|---|---|---|
| `mesh/x1max` | `12.54e9` | `13.04e9` | the `xe_long` radial extent |
| `meshblock/nx2`, `nx3` | `8`, `8` | `16`, `16` | 32 blocks on one GPU — see the decomposition note below |
| `mhd/max_eta` | `1.0e14` | `1.0e13` | the long-run resistivity cap |
| `mhd/use_rkg_sts` | `true` | `false` | RKG is required at 1e14, optional at 1e13; leaving it off is what makes this the cheap configuration |
| `mhd/pfloor` | `1.0e-1` | `1.0e0` | the `xe_long` floors |
| `mhd/dfloor` | `5.44e-12` | `7.26e-12` | " |
| `problem/bbot` | `5.0e0` | `1.0e1` | " |
| `problem/rt_ck` | `false` | `true` | correlated-k on |
| `problem/ck_table`, `ck_data_dir` | relative to `data/exo_fms_ck` | absolute | only because the job ran from a scratch dir; make them resolve, however |
| `output1/dt` | `8.64e6` | `2.0e2` | fine cadence, from when this was thought to fail early |
| `output2/dt` | `2.16e6` | `2.0e3` | " |
| — | — | add `<output3> file_type=rst, dt=1.0e4` and `<output4> file_type=log, dt=2.0e2` | restarts + the event log |

**`rt_de_max` in that file is `0.5`, i.e. the limiter ON.** On the fixed binary it barely
matters — over the whole 3.934 rot it clipped **4 cells**, and `rt_de_max = -1` gives the
same mass to five digits. It is still REQUIRED for the ideal-EOS + correlated-k run
(section 2), which NaNs without it. Set it deliberately rather than inheriting it.

The variants below were all one-line edits of that same file. Their pre-fix outcomes are
listed only to mark them as retracted; **all six clear `tlim` on the fixed binary** (section
0), so none of these knobs is a stability control.

| variant | the single change | pre-fix (RACED — not valid) | post-fix |
|---|---|---|---|
| `grey` | `rt_ck = false` | "dies 2.509 rot" | clean |
| `noh2` | `eos_h2 = false` | clean | clean |
| `hr` | `mesh/nx1` and `meshblock/nx1` 64 -> 128 | "dies 0.395 rot" | not re-run |
| `nocond` | `eos_metal_condensation = false` | "dies 0.932 rot" | not re-run |
| `shipfl` | `pfloor 1.0e-1`, `dfloor 5.44e-12` | "dies 0.470 rot" | not re-run |
| `tfloor` | `tfloor_kelvin` 200 -> 35 | "fails 4x sooner, hard NaN" | not re-run |
| `weakb` | `bbot = 1.0e-3` | clean | clean |
| `noeta` | resistivity removed | past 1.86 rot | clean |
| `hydro` | `hydro_ck.athinput` — pure hydro | clean | clean |
| `ideal_ck` / `ideal_grey` | separate ideal-gas files, differing in `rt_ck` ALONE | see section 2 | `ideal_ck` needs `rt_de_max` |

The ideal-gas pair is derived from `inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput`.
The four "not re-run" rows would be worth 8 h each only if someone wants the record
complete; there is no open question they answer.

Also relevant to sizing a job: 32 meshblocks on one GPU beats 2 by 1.35-1.45x, and a
second APU buys only ~1 %; `nx1` is capped at 264 by MHD LDS (production uses 256), and
x1 cannot be split because the RT is a column solve.

## 7. Repo conventions

- Push to `fork`, never `origin` (origin is upstream IAS-Astrophysics and is denied).
- `run/` is ~165 GB of output: read-only, and never `git add -A` in this repo.
