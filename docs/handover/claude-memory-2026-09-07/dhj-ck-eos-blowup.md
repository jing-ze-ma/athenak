---
name: dhj-ck-eos-blowup
description: "CLOSED 2026-08-25. There was no general-EOS blow-up: it was the GPU data race (6e600f12). All six re-run arms, INCLUDING the unmodified baseline, clear tlim 3.934 rot. Every pre-fix arm comparison in this file is retracted."
metadata:
  node_type: memory
  type: project
---

## 2026-08-25 (RESULT): THE BASELINE CLEARS tlim. THE BLOW-UP WAS THE RACE.

Six arms, 8 h each on apu1, all on `athena_fixrel` = **6e600f12** (the missing
`member.team_barrier()` in the general-EOS x1 flux kernel, [[dhj-run-to-run-nondeterminism]]),
all `tlim = 1.2e6 s = 3.934 rot`. **All six exited 0 at t = 1.200000e+06.**

| tag | change | pre-fix | post-fix mass |
|---|---|---|---|
| `rf.limoff` | **none, the control** (`rt_de_max=-1`) | died 0.370-0.685 rot | **1.00030** |
| `rf.limon` | `rt_de_max = 0.5` | cleared via a 33% mass transient | 1.00030, no transient |
| `rf.noh2` | `mhd/eos_h2 = false` | cleared | 1.00238 |
| `rf.hydro` | pure hydro | cleared | 1.00070 |
| `rf.weakb` | `bbot = 1.0e-3` | cleared | 1.00030 |
| `rf.noeta` | resistivity removed | past 1.86 rot | 1.00030 |

No NaN, no dt collapse (`limoff` dt min/med/final 0.921/1.254/0.987 s), `limoff` mass stays
in [0.99992, 1.00068] throughout. The limiter clipped **4 cells** over the whole run, so
`limon` and `limoff` are now the same run.

**What is retracted** (all of it was measured on a raced binary):
* the 0.370/0.485/0.584/0.685/1.464 rot spread was the bug, NOT chaos -- there is no
  chaotic ensemble to sample;
* `rt_de_max` does NOT cure a general-EOS blow-up -- it is still the correct fix for the
  **ideal**-EOS + correlated-k NaN, a genuinely separate failure;
* `eos_h2 = false` is not a cure; nor weak field, nor no-resistivity, nor pure hydro;
* **"the failure needs MHD" -- retracted**, hydro and MHD both run clean;
* the nine "refuted hypotheses" are moot, not wrong: there is no phenomenon left to explain.

**Still open, and fenced off by the user:** Gamma_1 spans 1.084-1.666 and jumps 0.58 between
adjacent cells while `src/mhd/mhd.hpp:588-596` assumes it smooth. That is now an ACCURACY
question, not a crash question. The user said "let's not go that way" on 2026-08-23 -- do
not open it unless they raise it.

Doc rewritten to match: `docs/HANDOFF_dhj_ck_eos.md`. **Everything below this line is the
pre-fix record, kept only to show what a raced binary did. Do not quote it as physics.**

## 2026-08-25: EVERYTHING BELOW PREDATES THE RACE FIX -- RE-RUNS IN FLIGHT

The GPU race is found and fixed ([[dhj-run-to-run-nondeterminism]], commit **6e600f12**,
pushed): a missing `member.team_barrier()` in the general-EOS x1 flux kernel. **The
0.370 / 0.485 / 0.584 / 0.685 / 1.464 rotation spread on byte-identical input WAS that
bug, not chaos**, and every arm recorded in this file was run on a raced binary, so no
comparison between arms below is currently trustworthy.

Six arms re-launched 2026-08-25 on `athena_fixrel` (= 6e600f12), 8 h each on apu1,
tlim 1.2e6 = 3.934 rot, all from `/viper/ptmp/jinma/claude_eos_gpu`:

| job | tag | change from the baseline |
|---|---|---|
| 11006656 | `rf.limoff` | none -- the control, `rt_de_max = -1` |
| 11006657 | `rf.limon` | `rt_de_max = 0.5` |
| 11006658 | `rf.noh2` | `mhd/eos_h2 = false` |
| 11006659 | `rf.hydro` | `hydro_ck.athinput` -- pure hydro |
| 11006660 | `rf.weakb` | `problem/bbot = 1.0e-3` |
| 11006661 | `rf.noeta` | `blowup_noeta.athinput` -- resistivity removed |

Read them with the summary in this file's older tables as the PRE-FIX comparison, not as
truth. `or.cpu` (job 11002742, CPU build, 16 ranks x 6 threads, orion's configuration) was
left running and was healthy at 0.740 rot; a CPU build is far less exposed to the race
because a serial `par_for_inner` executes in order.

Cancelled as superseded: `or.gpu`, `en.e1`, `en.e2` -- all three had already blown up, and
they existed only to measure the spread that turned out to be the bug.

### PICK UP HERE — how to read the six re-runs

```bash
squeue -u $USER
cd /viper/ptmp/jinma/claude_eos_gpu
python3 -c "
import numpy as np, glob
for d in sorted(glob.glob('rf.*/')):
    g = glob.glob(d + 'dhj.*.hst')
    if not g: continue
    a = np.loadtxt(g[0]); t = a[:,0]; m = a[:,2]/a[0,2]
    i = np.argmax(np.abs(np.diff(m)) > 0.05) if len(t) > 1 else 0
    print('%-12s t=%.3e (%.3f rot) mass=%.4f  first >5%% jump: %s'
          % (d, t[-1], t[-1]/3.05e5, m[-1],
             ('%.3f rot' % (t[i]/3.05e5)) if i else 'none'))
"
```

Failure = a >5 % mass excursion or a NaN (dt collapsing to ~1e-33). Rotation = 3.05e5 s,
tlim = 1.2e6 = 3.934 rot.

**Snapshot 25 min in (all six healthy, mass to 4 decimals):** `limoff` 0.361 rot,
`limon` 0.359, `noeta` 0.374, `weakb` 1.734, `noh2` 1.938, `hydro` 2.184.

`limoff` is the one to watch: it is the UNCHANGED baseline, and its pre-fix repeats died at
0.370 and 0.485 rot. It was already past 0.361 and healthy. **If `limoff` now clears
3.934 rot, the "general-EOS blow-up" was substantially the race itself** and most of this
file collapses to a single cause. If it still dies, the thermodynamic problem is real and
the other five arms become meaningful one-variable tests again — for the first time, since
every earlier comparison was confounded by the bug.

Do NOT compare these against the old tables as if they were the same experiment: quote the
new numbers, and treat the old ones only as "what a raced binary did".


## 2026-08-24: THE QUEUED JOBS LANDED. TWO INDEPENDENT CURES, BOTH CONFIRMED

Read the histories of the jobs left running on 2026-08-23. Rotation = 3.05e5 s, tlim = 1.2e6
= 3.934 rot. All six are general EOS + the same grid/floors.

| run | job | config | outcome |
|---|---|---|---|
| `lt.limon` | 10991231 | c-k, **rt_de_max = 0.5** | **CLEARED tlim, 3.934 rot**, mass -0.12 %, warned "clipped in 4 cell(s)" |
| `lt.limoff` | 10991232 | c-k, rt_de_max = -1 | **died 0.685 rot**, mass to 0.6 % |
| `sw.noh2` | 10990782 | c-k, **eos_h2 = false** | **CLEARED tlim, 3.934 rot**, mass +0.24 % |
| `gl.noh2grey` | 10991190 | grey, **eos_h2 = false** | **CLEARED tlim, 3.934 rot**, mass +0.24 % (a genuinely distinct run from `sw.noh2` -- energies differ by 1e14, only mass agrees) |
| `sw.nocond` | 10990783 | c-k, no metal condensation | died 0.953 rot |
| `sw.hr` | 10990781 | c-k, nx1 = 128 | died 0.399 rot -- **resolution does NOT cure it** |

So **the RT source limiter (ea823e3a) DOES fix the general-EOS + correlated-k blow-up**,
not just the ideal+c-k NaN: `limon` vs `limoff` is a one-parameter A/B on the same binary,
inf vs 0.685 rot. That was the open question queued on 2026-08-23 and it is answered.

**And `eos_h2 = false` cures it independently, in BOTH grey and c-k** -- which kills the
reading that H2-off merely moved the run out of correlated-k's reach. The H2 dissociation
plateau (the recombination front at r/r0 ~ 1.2, documented above) is the thermodynamic half
of the story: it is what makes e small and jagged enough for the unlimited explicit RT
source to cross tau < dt. Two halves, one failure -- the limiter fixes the numerics, H2-off
removes the trigger. `nocond` and higher resolution do neither.

Caveat on `limon`: at t = 7.3e4 it takes a real transient -- mass falls to 33 % over ~4
history samples and then re-fills to 99.9 % and stays smooth for the remaining 3.7 rot. The
run is healthy afterwards but that excursion is not understood.

## 2026-08-24 (III): IT IS NOT THE TIMESTEP -- THE CODE IS NONDETERMINISTIC RUN TO RUN

Two things were wrong in the older notes and are now corrected by measurement.

**1. Output cadence cannot change the dt sequence.** `Mesh::NewTimeStep` (mesh.cpp:637-716)
clips dt for `tlim` and nothing else; `Driver` fires outputs AFTER the step on
`time_32 >= next_32`. Outputs are read-only. So "denser output changed the dt sequence,
hence the different failure time" -- repeated in these notes from earlier sessions -- is
FALSE and should not be repeated.

**2. The real cause: the same binary, same input, same GPU, one rank gives different
answers run to run.** Measured, not inferred (`det.sh`, job 11003903): arms `a` and `b`
were byte-identical invocations of `athena_head` on `orion_repro.athinput`. Their
histories differ from line 5 (t = 4.0e2) and **every bin dump from the first written one
differs**, so the STATE diverges, not just the diagnostic. It shows up first in the
history's `2-mom` column because theta-momentum is a ~1e-9 cancellation of ~1e9 terms and
so amplifies a roundoff-level difference into the printed digits.

**Localisation** (`det2.sh`, job 11003990, bin output every cycle, two identical runs):
dumps 0-21 identical, first difference at **cycle 22**, and it is confined to

* `bcc3` (B_phi) ONLY -- dens, velx-z, eint, bcc1, bcc2 all bitwise identical;
* i = 61, 62, 63 -- the outermost ACTIVE radial cells;
* j = 0, 14, 15 and k = 0, 1, 14, 15 -- cells touching meshblock edges (nx2 = nx3 = 16);
* 4 of the 32 blocks, all at lx2 in {1,2}, lx3 in {4,5} -- NOT the polar blocks.

**Which module** (`det3.sh`/`det4.sh`, jobs 11004037/11004040, identical pairs, ~32 cycles):

| pair | config | differing dumps |
|---|---|---|
| grey1/grey2 | `rt_ck = false` (grey, monolithic) | **0** |
| gs1/gs2 | `rt_ck = false`, `rt_split = true` (1 chain block) | **0** |
| nomx1/nomx2 | correlated-k, `bc_outer_maxwell = false` | **25**, first at cycle 10 |
| ck1/ck2 | correlated-k, default | **0** |
| p/q (det2) | correlated-k, default | differ, first at cycle 22 |

So it is **intermittent** -- the same correlated-k configuration reproduced once and
diverged twice -- and it is NOT the outer-x1 Maxwell term (turning that off made it worse,
cycle 10). Correlated-k has diverged in 2 of 3 pairs; grey in 0 of 2, but those were only
32 cycles and grey is NOT yet proven clean. Note c-k runs 22 chain blocks where grey runs
1, so "c-k" and "nblk > 1" are still confounded; `rt_split` with grey only ever gives
nblk = 1, so that A/B could not separate them.

The split-kernel design is race-free ON PAPER: per-block storage `rt_Fb(m,blk,...)` /
`rt_Qb`, summed by kernel C in block order, no atomics (`deep_hot_jupiter_rt.cpp:4096-4112`).
`grep -rn atomic src/` finds none in the evolution path. So the race is somewhere else and
is NOT yet found. **This is the top open item.**

Context: the previous session found and fixed exactly this class of bug in the outer-x1 BC
(c39b794c -- a cross-thread read of `bcc0` inside the kernel that also wrote it), and the
long comment at deep_hot_jupiter_rt.cpp:2196 describes it. That fix was real but did not
make the problem reproducible.

**Consequence for every result in this file:** the blow-up times are draws from a
distribution and no single run distinguishes anything. Before trusting any arm, it needs
repeating -- or the race needs fixing so runs are reproducible.

## 2026-08-24 (IV): THE ORION RUN BLOWS UP TOO -- it was stopped too early

`or.gpu` (orion's exact input, viper GPU) passed 0.690 rot clean, as orion reported at
0.713, and then **died at t = 4.464e5 = 1.464 rotations**: mass 1.0002 -> 0.105 in ONE
200 s history interval, from dt = 1.34 s with no precursor. Orion stopped watching at
0.713, which is why it looked like a disagreement.

Failure times on this identical setup now read **0.370, 0.485, 0.584, 0.685, 1.464 rot**.
Nothing survives; the spread is the race plus chaos. `en.e1` is at 0.767 and still healthy,
`mv.noeta` at 1.856, `or.cpu` (CPU build, 16 ranks) at 0.257 -- all still running.

## 2026-08-24 (later): ORION REPRODUCED **ON VIPER**. It is not a machine difference --
## the baseline is a chaotic ensemble and some members simply live.

The orion session pushed **cd2e8815**, one file: `inputs/mhd/deep_hot_jupiter_rt_blowup.athinput`,
the literal input its job 190900 read, which passed 0.674 rot and was still rising. Diffed
against viper's `blowup.athinput` ignoring comments: identical except `tlim` (1.2e6 vs
8.64e7), `rt_de_max` (-1.0 vs 0.5), the two table paths, and one commented-out line. Both
were run at the same physics -- `lt.limoff` set `tlim=1.2e6` and `rt_de_max=-1` on the
command line. **Same problem, and the code is the same too: cd2e8815 adds no source.**

Re-ran that literal file on viper (`or.gpu`, job 11002708, path-repointed only) on a GPU
binary freshly built at 397b49b1: **passed 0.690 rot with mass = 1.00029 and dt steady at
1.19 s**, where `lt.limoff` was already down to 0.6 % of its mass at 0.685. So orion's
result reproduces on viper hardware. **CPU-vs-GPU is NOT the explanation.**

What differs between `or.gpu` and `lt.limoff` is (a) the bin-output cadence, 2.0e3 vs the
1.0e4 that `lt.limoff` passed on the command line, and (b) the build (ea823e3a vs
397b49b1 -- and c1d513c2's 79 lines are all gated on `ad_dump_file`/`ck_dump_file` being
non-empty, verified by reading the diff, so it is inert). Both are physics-inert knobs that
shift the dt sequence. This is exactly the sensitivity recorded earlier, where two
byte-identical setups died 30 % apart in time purely from a different output cadence.

**Death times on this identical setup now read 0.370, 0.485, 0.584, 0.685 rot, plus at
least two survivals past 0.69 (orion, `or.gpu`).** No single run of this baseline is
evidence of anything. Ensemble members `en.e1` (11003130) and `en.e2` (11003131) differ
from `or.gpu` ONLY in `output2/dt` (2.5e3, 3.0e3) and are running to measure the survival
rate. Until that lands, treat every "arm X survived" claim below -- including the hydro,
weak-field and no-resistivity results -- as needing the same ensemble treatment.

A `or.cpu` arm (job 11002742) is also running: the CPU MPI+OpenMP build at 397b49b1, 16
ranks x 6 threads on an APU node's host cores, same file. NOTE: viper's `general`/`small`
CPU partitions cannot be reached from the `viper13` login node (their `AllocNodes` is
`viper01-08`, and ssh there needs Kerberos), so a Sapphire-Rapids-class CPU run must be
submitted by the user from `viper01-08`.

## 2026-08-24: the MHD ladder -- weak field is clean, no-resistivity is at least better

Both rt_de_max = -1, same grid/floors/tlim, `athena_srclim`.

| run | job | config | outcome |
|---|---|---|---|
| `mv.weakb` | 11002375 | `bbot` 10 G -> **1 mG** | **cleared tlim, 3.934 rot, mass +0.034 %**, dt ~6-7 s |
| `mv.noeta` | 11002374 | `ohmic_resistivity` REMOVED (ideal MHD, 10 G) | past **1.03 rot** and healthy, dt ~1.5 s and rising |
| `lt.limoff` | 10991232 | full baseline | died 0.685 rot |

Both point at the field, and `mv.weakb` going the full distance is the strongest arm --
but see the ensemble caveat directly above before calling either one a cure.

## 2026-08-24 RESULT: PURE HYDRO SURVIVES. THE FAILURE NEEDS MHD.

The user's test. `hydro_ck.athinput` = `blowup.athinput` with `<mhd>` -> `<hydro>`,
hlld -> hllc, the three resistivity lines commented out, output2 variable -> `hydro_w`.
Everything else byte-identical, same binary `athena_srclim`, same tlim 1.2e6 = 3.934 rot.
The pgen already branches phydro/pmhd throughout and needed no code change.

| run | job | config | outcome |
|---|---|---|---|
| `hy.hoff` | 11001828 | hydro, c-k, **rt_de_max = -1** | **cleared tlim, 3.934 rot, mass +0.068 %** |
| `hy.hon` | 11001829 | hydro, c-k, rt_de_max = 0.5 | cleared tlim, mass +0.068 % (clipped 10 cells) |
| `lt.limoff` | 10991232 | **MHD**, c-k, rt_de_max = -1 | died 0.685 rot |

**Removing the magnetic field cures it outright, with the RT limiter OFF**, and the hydro
mass history is not merely survivable but flat -- +0.068 % over four rotations, no
excursion anywhere (contrast `lt.limon`, which survives but takes a 33 % mass transient).
The limiter is irrelevant in hydro: both arms give the same mass to three digits.

**And it is not a dt effect.** Median dt: hydro 6.97 s, `lt.limoff` 1.28 s (min 0.0156 s).
Hydro takes ~5x LONGER steps and still lives, so "explicit RT source vs dt" cannot by
itself be the story -- something magnetic manufactures the low-e jagged state that makes
the source term marginal. The RT limiter and `eos_h2 = false` are still real cures, but
they act on the last link of the chain, not the first.

## THE 2x2, CLOSED (2026-08-23) -- BOTH factors destabilise, correlated-k more

| EOS | RT | fails at |
|---|---|---|
| ideal | grey | **never** -- clean to 3.93 rot, hit `tlim`, mass to 0.15 % |
| general/table | grey | 2.51 rot (t = 7.652e5) |
| general/table | correlated-k | 0.37 / 0.49 / 0.58 rot; 1.07 with RKG |
| ideal | correlated-k | **0.017 rot** (NaN at t = 5.2e3) |

Two clean single-variable comparisons, each isolating one factor:

* **ideal+grey vs general+grey** -> the tabulated EOS destabilises (inf -> 2.51 rot).
* **ideal+grey vs ideal+c-k** -> correlated-k destabilises (inf -> 0.017 rot). The two
  inputs `ideal_grey.athinput` and `ideal_ck.athinput` differ in `rt_ck` ALONE.

**RETRACTED: "ideal + c-k is an unsound pairing, disregard."** I could find no evidence
behind that label and it is wrong to lean on it. Re-run today with the current binary
`athena_after2` (job 10990853, `ck.ick2`): NaN reproduces at t = 5.14e3, matching the
original to 1 %. It is NOT the `wtemp` out-of-bounds bug either -- `athena_after` was built
15:16, after 0d1f6f6a landed at 14:47. Also RETRACTED: "the RT is innocent, chase the EOS."
Both factors are real; c-k is the stronger by ~150x in time-to-failure.

The one residual physics argument for treating the ideal+c-k cell as weaker evidence -- NOT
yet tested -- is that the c-k opacities are keyed on FastChem equilibrium composition at
(T,p), while the ideal EOS's T = p/(Rgas rho) at fixed mu is inconsistent with that
composition. That is a hypothesis, not a reason to discard the cell.

## MEASURED: correlated-k cools the tenuous top far harder than the grey picket fence

Same ideal EOS, same floors, only `rt_ck` differs. Minimum values over the whole domain
(dfloor = 7.26e-12; the 200 K tfloor is eps = Rgas*200/(gamma-1) = 1.943e10):

| t | ideal+c-k min eps | ideal+grey min eps |
|---|---|---|
| 0 | 3.22e11 | 3.22e11 |
| 2008 | 1.44e11 | 1.88e11 |
| 4001 | **1.943e10 = ON the tfloor** | 9.01e10 (4.6x above it) |

Both reach the density floor by t = 4000, but only the c-k run reaches the temperature
floor -- and it NaNs ~1200 s later. Event counters (`ck.ick2/dhj.log`) show `eos_dfloor`
climbing to 2.9e6 firings per ~20-cycle interval on 524288 cells, i.e. **~28 % of the
domain on the density floor every cycle**, with `eos_fail = c2p_it = vceil = fofc = 0`
throughout. The run is being evacuated by radiative cooling, not failing an inversion.

This is the natural mechanism for c-k accelerating the general-EOS failure ~7x as well:
11 resolved bands cool the optically thin upper atmosphere much faster than a grey
picket fence, and the general EOS's chemical reservoir is what buys the extra ~20x of
survival time.

## What happens

Running `xe_long/b10_e13`'s grid and floors (64x64x128, 32 blocks, bbot = 10 G,
x1max = 13.04e9, pfloor 1.0e0, dfloor 7.26e-12) but with **`eos = general/table` +
`rt_ck = true`**, the run dies within one rotation (rotation = 3.05e5 s):

| config | fails at | rotations |
|---|---|---|
| max_eta 1e13, no RKG | t = 1.48e5 | 0.49 |
| max_eta 1e13, no RKG (rerun, denser output) | t = 1.13e5 | 0.37 |
| max_eta 1e14, no RKG | t = 1.78e5 | 0.58 |
| max_eta 1e14, **RKG on** | still healthy at 2.05e5 | >0.67 |

**The resistivity cap is NOT the cause** -- both caps fail. The two 1e13 runs are physically
identical and fail 30 % apart in time (denser output changes the dt sequence, hence rounding),
so this is a marginal/chaotic instability, not a deterministic bug at a fixed time.

## Signature (from the instrumented run, `/viper/ptmp/jinma/claude_eos_gpu/as.a13n_dump`)

Angular-mean radial profiles, last good dump (t=1.12e5) vs first bad (t=1.14e5):

* the base is EVACUATED to ~45 % of its density;
* **above r/r0 ~ 1.2 the stratification is replaced by a UNIFORM state** -- rho ~ 3.7e-7
  everywhere (3500x too high at the top) and eint ~ 1.4e7 (65,000x too hot at the top);
* total mass falls to 7.6 % in the history, 37 % in the angular mean;
* it is NOT a density-floor cascade: afterwards NO cell is on the density floor and the
  minimum density is 45,000x ABOVE it.

Event counters (`file_type = log`) through the failure:

* `eos_fail = 0`, `c2p_it = 0`, `vceil = 0`, `fofc = 0` -- **the EOS inversion never fails**;
* before: dfloor ~1.5e7, efloor ~1.0e7, tfloor ~1.0e5 firings per ~150-cycle interval;
* at the failure **tfloor jumps 12x to 1.26e6** while dfloor collapses to 2.9e4.

## TEST 1 DONE -- the temperature floor is REFUTED as the cause; it was PROTECTING the run

`tfloor_kelvin` 200 -> 35 K (35 K is just inside the EOS table, whose floor is
`eos_logt_min = 1.5` = 31.6 K; the long runs' `tfloor = 10` CODE units is ~1e-5 K and would
be off-table). Everything else identical.

| tfloor_kelvin | outcome | when |
|---|---|---|
| 200 K | mass loss, limps on | t = 1.13e5 (0.37 rot) |
| **35 K** | **outright NaN** | **t = 2.66e4 (0.09 rot)** |

**4x SOONER, and a hard NaN instead of a mass-loss event.** dt declined normally to 2.4 then
froze at 4.681 -- the NaN-timestep signature. Energy-floor firings per cell per interval:
2-10 at 200 K, but **141-145 at 35 K, i.e. every cell every cycle**. The whole domain lands
on the energy floor once the temperature floor stops intercepting it first.

So the floors are the SYMPTOM, not the disease: something upstream drives the gas to
energies the floors must rescue. Do not chase the floors again without new evidence.

## Superseded hypothesis (kept so it is not re-proposed)

The long runs use `tfloor = 10.0` in CODE units, which is **effectively no floor at all**
(~1e-5 K; cf. the note in the shipped input about tfloor = 400 being ~1e-5 K). These runs
used the shipped general-EOS `tfloor_kelvin = 200.0`, **a real 200 K floor**, and it fires
1e5-1e6 times per interval in the coldest, most tenuous cells and spikes exactly at the
failure. A temperature floor RAISES e, so it injects energy into the upper atmosphere
continuously -- consistent with the outer region ending up uniformly hot and dense.

Caveat on provenance: this floor combination (long-run pfloor/dfloor + shipped
tfloor_kelvin) was assembled for a timing comparison, not by the user, and may simply be
wrong for a general-EOS run.

## apudev CANNOT see this -- and that invalidates "it ran clean" as evidence

Wall time actually needed to reach the failure, measured:

| run | wall to blow-up | where a 12-min apudev slot stops |
|---|---|---|
| 1e13 original | **36.5 min** | t = 6.7e4, 45 % of the way |
| 1e13 instrumented (dense I/O) | 26.4 min | t = 6.6e4, 58 % |
| 1e14 no RKG | **73.4 min** | t = 5.3e4, 30 % |

`apudev` caps at 15 min, so it reaches t ~ 5-7e4 and the earliest failure is 1.13e5.
Every cost and profiling number in [[general-eos-optimization]] and
[[dhj-overall-gpu-profile]] was taken on apudev at t < 4e4, where this setup is still
perfectly healthy -- mass and energy conserved to 0.07 %, dt declining smoothly. **Those
numbers remain valid as COST measurements, but "the run completed cleanly on apudev" is
worth nothing as evidence of stability.** Use `apu1` (24 h) for anything about robustness.

## RESULTS MATRIX (all at the same grid/floors, incl. an effective 200 K tfloor)

Rotation = 3.05e5 s. Failure = >5 % mass excursion or NaN.

| EOS | RT | variant | outcome |
|---|---|---|---|
| general | correlated-k | 1e13 no RKG | **fails 0.370 rot** |
| general | correlated-k | 1e13 no RKG (repeat) | **fails 0.485 rot** |
| general | correlated-k | 1e14 no RKG | **fails 0.584 rot** |
| general | correlated-k | 1e14 **RKG on** | **fails 1.072 rot** |
| general | correlated-k | 1e13, tfloor 35 K | **NaN 0.087 rot** |
| general | **grey** | 1e13 no RKG | **fails 2.510 rot** (t=7.652e5) |
| ideal | **grey** | 1e13 no RKG | **CLEAN to 3.93 rot** (hit tlim 1.2e6) |
| ideal | correlated-k | 1e13 no RKG | **NaN 0.017 rot** -- REAL, reproduced 2026-08-23 |

**Only ideal+grey survives. Every other cell of the 2x2 fails.**

**RKG only DELAYS it (0.584 -> 1.072), it does not cure it** -- I wrongly called this "the
strongest lead" before a14r ran on. The ordering 0.370 < 0.485 < 0.584 < 1.072 tracks
increasing dissipation, so this behaves like a **marginally resolved instability being
damped, not cured**. Even if c-k is the trigger, "grey is fine" may only mean "grey is
further from the margin". **No resolution test has been done -- that is the obvious gap.**

**The failures are NOT deterministic:** two byte-identical 1e13 setups failed 30 % apart in
time (0.370 vs 0.485) purely because a different output cadence changed the dt sequence.
Any single-run comparison near the margin is weak evidence.

## Hypotheses tested and REFUTED (do not re-propose without new evidence)

1. **Temperature floor** -- refuted, and backwards: lowering tfloor 200 -> 35 K made it fail
   4x SOONER and as a hard NaN. The floor was PROTECTING the run. At 35 K the energy floor
   fires on every cell every cycle (141-145 per cell per interval vs 2-10 at 200 K).
2. **Resistivity cap** -- refuted, both 1e13 and 1e14 fail.
3. **RKG coupling** -- refuted as a cure; a ~2x delay only.

Floors are the SYMPTOM, not the disease: `eos_fail = 0`, `c2p_it = 0`, `vceil = 0`,
`fofc = 0` in every run.

## SUPERSEDED lead: RKG super-time-stepping (delay only, see above)

`a14r` -- general EOS + correlated-k, max_eta 1e14, **RKG ON** -- reached t = 2.61e5
(0.86 rotations) with mass conserved to 0.02 %, while every no-RKG variant died at
0.37-0.58 rotations. RKG applies the resistive update via Lie splitting AFTER the full RK2
step, instead of putting the resistive fluxes inside the RK2 stages. If that difference is
what decides stability, the bug is in **how the Ohmic term is coupled**, not in the EOS or
the RT. Worth testing directly: 1e13 WITH RKG.

## Isolation experiments

1. ~~tfloor off~~ **DONE, refuted -- see above.**
2. **grey RT** at otherwise identical settings -- isolates correlated-k. RUNNING.
3. **ideal EOS + correlated-k** -- pairs with (2) for a clean 2x2.
4. **shipped general-EOS floors** (pfloor 1e-1, dfloor 5.44e-12) instead of the long-run
   ones, which were imposed for a timing comparison and may be wrong here -- `pfloor = 1.0`
   in particular.
4. Why does RKG survive longest? Its Lie splitting applies resistivity after the full RK2
   step -- worth knowing whether that is genuinely stabilising or just delay.

## WHERE THE ENERGY IS: the atmosphere is parked ON the H2 dissociation plateau

Per-cell analysis of `iso.grey` (general+grey) vs `iso.igr` (ideal+grey) at the SAME time
t = 7.64e5, one dump (1.2e3 s) before the general run detonates. Specific internal energy
eps = eint/dens, angular mean, cgs:

| r/r0 | eps general | eps ideal |
|---|---|---|
| 1.00 | 1.70e12 | 1.18e12 |
| 1.10 | 1.14e12 | 7.9e11 |
| 1.20 | **5.7e11 (minimum)** | 2.4e11 |
| 1.25 | **1.90e12** | 2.1e11 |
| 1.37 | 1.79e12 | 3.0e11 |

Full H2 dissociation costs `0.5*chi_d` per H atom = **1.58e12 erg/g of gas** at X_H = 0.7381
(`chi_d = 4.478 eV`, `eos_composition.hpp:65`, zero point = all H in ground-state H2,
`e_href` at line 517). So:

* the general-EOS upper atmosphere sits at eps ~ 1.9e12, i.e. **just past complete
  dissociation, with the chemical term 5x the thermal term**;
* between r/r0 = 1.19 and 1.25 eps swings by 3.4x -- a **recombination front inside the
  domain**, only ~9 radial cells wide on a 64-cell grid;
* the front SHARPENS monotonically with time (eps at r/r0 = 1.20: 1.53e12 at t=0, 7.6e11 at
  2e5, 5.8e11 at 6e5) although the inversion itself is present in the INITIAL CONDITION;
* laterally at r/r0 = 1.25 the general run spans eps = 1.6e10 (the 200 K tfloor, nightside)
  to 1.8e13 (dayside, ~3e4 K and fully ionized) -- a 1000x contrast in one shell. The ideal
  run at the same radius spans only 1.9e10 to 8.9e11.

The ideal gas has NO chemical reservoir, so its upper atmosphere never leaves 2-3e11 and
never develops the front. **That is the structural difference between the run that lives
and the runs that die.** Only 0.06 % of cells are on the tfloor, confirming again that the
floors are a symptom.

## THE ideal+c-k NaN, LOCALISED (2026-08-23) -- floors build an artificial 3000 K top,
## and the c-k column solve returns a flux that cannot be physical

Dense-output rerun `dn.ickd` (input `v_ickdense.athinput`, hst+log every 5 s, bin every
25 s, job 10990905, 27 s of compute). Sequence:

* t = 5127.4: dt falls 5.00 -> 0.134 in one interval, mass still exact;
* next cycle: dt = 1.07e-33 and everything is NaN. So a FINITE, absurd sound speed comes
  first -- the NaN is downstream of it.

**The runaway is confined to ONE column, (k=99, j=14), cells i=50..63.** Every cell with
eps > 1e13 in the whole domain is in that single column. mu0 = -0.112, i.e. just past the
terminator on the NIGHT side, so `lit` is false and no stellar beam is involved -- this is
the longwave. The cell (99,14,54) sits at eps = 2.9e11 for 400 s with density creeping onto
the floor, then jumps to **2.25e16 in one 25 s interval, a factor of 77,000**, with no
precursor. The RT chain kernel is the only per-(k,j) column solve in the code.

**The implied flux is impossible.** rho = 7.29e-12, Delta eps = 2.25e16, bdt ~ 5 s, dr =
5.625e7 give a source of 3.3e4 erg/cm^3/s and a flux divergence of **1.8e12 erg/cm^2/s**.
The hottest gas anywhere in the domain is 3000 K, whose sigma T^4 is 4.59e9. A two-stream
solve cannot return a flux 400x the Planck flux of the hottest available source. That is
the strongest single piece of evidence that the c-k path has a real defect, or at minimum
no guard, in this regime.

**The regime is one the floors manufacture.** At t = 5102 the top of EVERY column is pinned
at `pfloor = 1.0` dyn/cm^2 AND `dfloor = 7.26e-12`, which forces exactly

    T = pfloor/(Rgas*dfloor) = 1.0/(4.593e7*7.26e-12) = 2998.93 K

and the dumps read 2998.93 K to five digits over i = 55..63. So the "upper atmosphere" is an
artificial isothermal 3000 K slab, sitting directly on cells that are at the 200 K `tfloor`
in neighbouring columns -- adjacent cells differing 15x in T, 5e4 in B. Density drops 38x
in one cell (i=52 -> 53). `ideal_ck.athinput:159` even says pfloor is redundant because
`dfloor*Rgas*tfloor = 0.050`; at 1.0 it is 20x that and it is the binding constraint.

**Reading: both.** The floors put the solver in a state no atmosphere reaches, and the c-k
column solve does not degrade gracefully there. The grey run reaches the density floor at
the same time but cools 4.6x less far, never builds the discontinuity, and survives.

NOT yet done, and the obvious next step: dump F(i) and Q(i) for column (m,k=99,j=14) a few
cycles before the jump, using the pgen's built-in `ck_dump_file` / `ck_dump_j` /
`ck_dump_k` one-shot column dump. That needs a restart saved within ~25 s of t = 5127
(the dense run had `rst` disabled). `shipfl` (job 10990784) already tests the floor half:
pfloor 1e-1 + dfloor 5.44e-12 gives a 400 K top instead of 3000 K.

## ANSWERED (2026-08-23): NOT a correlated-k bug. The RT source has NO radiative
## timestep constraint, and the floors push the top where that bites.

Used the pgen's own one-shot column dump (`problem/ck_dump_file` / `ck_dump_m` /
`ck_dump_j` / `ck_dump_k`) on the doomed column, restarting `ck.ick2`'s saved restarts
(job 10990962, `coldump/`). Global (k=99, j=14) maps to meshblock **m = 24**, and with
`nghost = 2` (plm) the dump wants absolute **k = 5, j = 16**; the header echoed
mu0 = -0.1436, confirming the column. Restart file n holds t = 20n.

`src = -(F(i+1)-F(i))/dr + Q_sw` is applied as `u0(IEN) += src*bdt` (`rt_apply`,
deep_hot_jupiter_rt.cpp:4487) with the HYDRO timestep and no limiter. Comparing the local
radiative time tau = e/|src| against dt, e = p/(gamma-1):

| | worst tau/dt in the column | verdict |
|---|---|---|
| t = 2000, healthy | **86** (range 86-570) | stable, F constant to 5 digits in the thin top |
| t = 5120, 22 s before NaN | **0.014** (five cells below 1) | explicit source, 70x over the limit |

At t = 2000 the c-k solve is demonstrably CORRECT: in the optically thin top F_lw is
3.1253e9 to five significant digits over i = 59..66, i.e. zero flux divergence, exactly as
it should be. There is no table-lookup or two-stream defect to find. The startup self-tests
(isothermal invariance 1.9e-17, transparent slab exactly 1) also pass.

At t = 5120 the column reads, bottom to top:

    i=54  p=1.12e-5 bar  T= 650 K
    i=55  p=1.000e-6 bar (EXACTLY pfloor)  T=2081 K   tau/dt = 0.014
    i=56  p=3.29e-6 bar  T=9806 K                     tau/dt = 0.25
    i=57  p=1.27e-6 bar  T=3804 K                     tau/dt = 0.29
    ...   p -> 1.000e-6 bar, T -> 2998.93 K

Note p is NON-MONOTONIC (1.00e-6 at i=55 below 3.29e-6 at i=56): the pressure floor has
inverted the gradient. `pfloor = 1.0` dyn/cm^2 pins e at 1.0/(gamma-1) = **2.115 erg/cm^3**,
so even a flux divergence of 29 erg/cm^3/s -- unremarkable in absolute terms -- gives
tau = 0.072 s against dt = 5 s. The scheme then overshoots to negative energy, the floor
rescues it, and it overshoots again: that oscillation IS the 650 / 2081 / 9806 K jaggedness,
and the 9806 K cell is a previous overshoot.

**Why correlated-k and not grey.** Same `rt_apply` structure, same exposure. c-k drives the
top 4.6x colder (onto the tfloor, measured above), and once the profile is jagged the 11
band fluxes respond far more sharply than one grey mean. Grey never crosses tau = dt.
Likewise the general EOS survives ~20x longer because its chemical reservoir keeps e large,
so tau/dt stays above 1.

**Two independent fixes, both wanted:**

1. **Setup.** `ideal_ck.athinput:159` already says pfloor is redundant --
   `dfloor*Rgas*tfloor = 0.050` -- yet pfloor = 1.0 is 20x that and binding. It manufactures
   the artificial isothermal 2998.93 K top. `shipfl` (job 10990784) tests pfloor 1e-1 +
   dfloor 5.44e-12, a 400 K top.
2. **Code.** The RT source needs a guard. A hard `dt <= C*min(e/|src|)` would be correct but
   would crush dt to ~0.07 s. Better: clamp the per-step energy change to a fraction of e,
   or apply the source semi-implicitly toward the local radiative equilibrium. This would
   protect the general-EOS runs too and is the reusable fix. NOT YET IMPLEMENTED.

## IMPLEMENTED 2026-08-23: the RT source limiter (`problem/rt_de_max`)

Committed as **ea823e3a** (5 files, PUSHED to `fork/polar-average-perf`): `src/pgen/deep_hot_jupiter_rt.cpp`,
`docs/correlated_k_rt.md`, both shipped `inputs/mhd/deep_hot_jupiter_rt_*.athinput`, and a
new `tst/test_suite/rad/test_rad_dhj_srclim_cpu.py`.

`LimitRTSource(de, eint, de_max)` caps every EXPLICIT radiative update at
`|de| <= rt_de_max*e_int`, default 0.5, `<= 0` disables. A **hard clamp**, deliberately:
it is the identity wherever `|de| < rt_de_max*e`, so it cannot move an answer in any
regime where the explicit step was legitimate -- the Exo-FMS validation survives untouched.
Wired into all three explicit sites: `rt_apply` (split path, grey + c-k), the monolithic
`2stream_rt`, and `double_gray_two_stream_RT`. `double_gray_two_stream_RT_source` already
solves its source implicitly and was left alone. Counting is free -- the three kernels
became `par_reduce_clip3/4`, local flattened `parallel_reduce` wrappers mirroring
athena.hpp's `par_for` -- and `RTSourceLimiterWarn` warns ONCE per run.

**A/B on the exact case that NaN'd** (`v_ickdense.athinput`, ideal EOS + c-k, tlim 2e4):

| | result |
|---|---|
| `rt_de_max = -1` (disabled) | NaN at **t = 5.12754e3** -- byte-identical failure time to the unpatched binary |
| `rt_de_max = 0.5` | **no NaN**, ran to tlim 2.0e4, mass to 0.16 %, warned "clipped in **10 cell(s)**" |

Ten cells is exactly the count found in the dump analysis above. The disabled arm is a
perfect control: the only thing that changed is the clamp.

CPU verification of the new test's assertions, on the real binary at nx = 64x8x8, 20 cycles:
c-k history bitwise IDENTICAL on/off; grey history bitwise IDENTICAL on/off; zero warnings
in both healthy runs; `rt_de_max = 1e-12` clips 2816 cells and warns exactly once.

**QUEUED 2026-08-23: does the limiter also fix the GENERAL-EOS failures?** Untested until
now, and the open question. Jobs **10991231 (`lt.limon`, rt_de_max = 0.5)** and **10991232
(`lt.limoff`, rt_de_max = -1)**, 8 h each on apu1, both on `athena_srclim` so the clamp is
the only difference. Input `blowup.athinput` = the baseline that dies at 0.37-0.58 rot.

* `limoff` must reproduce 0.37-0.58 rot, else something other than the clamp changed.
* `limon` clearing **2.51 rot** (where general+grey died) would say the same missing dt
  constraint drives the general-EOS failures too.
* `limon` dying at 0.4-0.6 anyway would say the general-EOS blow-up is a separate,
  thermodynamic problem -- which is what `noh2` and `gl.noh2grey` are chasing.

Note `iso.grey` also died in ONE cycle (dt 0.97 -> 0.0108, mass 34x in one interval), so
the same mechanism is plausible -- but plausible is not measured.

**Do not oversell it.** t = 2e4 is 0.066 rotations. This proves the limiter stops THIS NaN;
it does not make the setup sound. The floors still manufacture a 2998.93 K top, and the
general-EOS failure at 0.4-2.5 rot is a separate, thermodynamic problem (see the sweep).

## SWEEP RESULT (jobs 10990781-4, ~50 min in, all still running)

| tag | change | status |
|---|---|---|
| `noh2` | `eos_h2 = false` | **3.37 rot, mass to 0.24 %** -- past general+grey (2.51) and still going. Its dt is ~8 s vs the baseline's 0.015 s |
| `nocond` | `eos_metal_condensation = false` | 0.63 rot healthy -- past the baseline band (0.37-0.58) but not decisively |
| `hr` | `nx1` 64 -> 128 | 0.27 rot healthy, too early |
| `shipfl` | shipped floors | **FAILED ~0.478 rot**, mass to 1.6 % |

**Every sweep arm is general EOS + CORRELATED-K** (`rt_ck = true`, verified in the inputs
and in each run.log's startup banner). c-k is the fast reproducer -- it fails at 0.37-0.58
rot versus 2.51 for general+grey -- so it is the only affordable bed for varying one EOS
switch at a time. All four ran on `athena_after2`, i.e. BEFORE the RT source limiter
(ea823e3a); the comparison between arms is still clean because they share that.

**QUEUED 2026-08-23: `gl.noh2grey`, job 10991190**, 12 h on apu1, `v_noh2grey.athinput`,
also on `athena_after2` so it matches its twins. It is general EOS + `eos_h2 = false` +
**GREY** RT, and it differs from `grey.athinput` (which produced `iso.grey`, dead at 2.51
rot) in `eos_h2` ALONE. It separates "H2 off cures it" from "H2 off merely moves it out of
correlated-k's reach": clearing 2.51 rot means the cure is real, failing near 2.5 means
`noh2`'s c-k survival was only about c-k.

**`shipfl` failing matters:** the shipped floors alone do NOT cure the general-EOS blow-up,
so the floors are not the whole story there -- they are the story for the ideal+c-k NaN.
**H2 dissociation is the strong lead for the general-EOS failure**, consistent with the
dissociation-plateau structure measured above.

## RESULTS 2026-08-23, second half: what is ruled OUT for the H2/general-EOS failure

All at the xe_long grid/floors, general EOS + correlated-k unless stated. Rotation = 3.05e5 s.

| test | job | outcome |
|---|---|---|
| radial resolution `nx1` 64 -> 128 | `sw.hr` | **FAILED 0.395 rot** -- resolution is NOT the fix |
| `eos_metal_condensation = false` | `sw.nocond` | **FAILED 0.932 rot** -- delays only |
| shipped floors (pfloor 1e-1, dfloor 5.44e-12) | `sw.shipfl` | **FAILED 0.470 rot** |
| `eos_h2 = false`, correlated-k | `sw.noh2` | **CLEAN to tlim 3.934 rot**, mass 0.24 % |
| `eos_h2 = false`, GREY | `gl.noh2grey` | **CLEAN to tlim 3.934 rot**, mass 0.24 %, dt flat at 8.09 |
| RT source limiter ON | `lt.limon` | survives to 3.934 rot BUT loses 67 % of its mass at 0.240 rot and the density floor refills it to 0.9988 -- not a physical state |
| RT source limiter OFF (control) | `lt.limoff` | **FAILED 0.674 rot**, i.e. the usual baseline |

`hr` failing kills the "marginally resolved instability" reading that the dissipation-ordered
failure times (0.370 < 0.485 < 0.584 < 1.072) had suggested. Do not re-propose resolution.

## The initial condition is NOT convectively unstable -- HYPOTHESIS REFUTED

New diagnostic `problem/ad_dump_file` (**c1d513c2**, pushed) writes the initial (p,T) profile with
nabla, grad_ad and the 0.9*grad_ad the code actually enforces, before and after
`adjust_ad_pT_arr`. I had suspected that routine of fixing only ONE convective zone,
because it scans from the bottom, `break`s at the first crossing, and a general EOS has a
grad_ad dip at H2 dissociation AND at H ionization.

**Wrong.** After the adjustment there are ZERO genuinely super-adiabatic levels (tolerance
1 % of grad_ad), H2 on or off. The single unstable band beforehand (p > 3.17 bar, the deep
interior) is fully fixed. In the dissociation band itself (p = 0.02-0.08 bar, T ~ 3300-3500 K)
the profile runs nabla = -0.03..+0.075 against 0.9*grad_ad = 0.083-0.095: comfortably
STABLE, with a mild inversion. The run does not start unstable; it develops the front.
A first count of "1169 unstable levels after" was an artefact of testing `>` against a
profile the adjustment had set to exactly the threshold. The shipped dump now applies a
1 %-of-grad_ad tolerance itself and reports "2416 before, 0 after".

**But the same dump shows `noh2` is a much bigger change than "H2 off".** grad_ad's minimum
is 0.0914 with H2 and 0.2687 without, so the adiabat is far steeper and the initial interior
reaches **14098 K at 300 bar without H2 versus 6728 K with it**. `noh2` is a factor-of-two
hotter, differently stratified planet -- still the strongest signal we have, but much weaker
as evidence that dissociation is the MECHANISM than the clean 2x2 makes it look.

## MEASURED: Gamma_1 swings across its whole range between ADJACENT CELLS

`ck_dump_file` now also writes Gamma_1 and grad_ad of the current state (**c1d513c2**, pushed; appended as columns 6-7 so existing parsers are unaffected).
Restarting `iso.grey` and `as.a13n_dump` with `rt_ck=true` purely to make the one-shot dump
fire (it reports the restarted state, so the RT choice does not affect what is measured;
note `time/nlim` must NOT be overridden -- a restart resumes at a large cycle number and a
small nlim runs zero cycles and never calls the RT):

| state | Gamma_1 range | steepest |dGamma_1| across ONE cell |
|---|---|---|
| grey, healthy (t = 2e5) | 1.104 - 1.666 | **0.442** (1.657 -> 1.215, T 3020 -> 2296 K) |
| grey, pre-failure (t = 7.6e5) | 1.084 - 1.661 | **0.576** (1.084 -> 1.661, T 1963 -> 4006 K) |
| correlated-k, pre-failure | 1.087 - 1.666 | **0.579** |

So Gamma_1 changes by 30-40 % across a single cell, at r/r0 = 1.26-1.33, p ~ 1e-6..1e-5 bar
-- and it does so **already in the healthy state**. T there swings 1711 -> 5529 -> 2296 K
between neighbouring cells while p sits pinned at pfloor = 1.000e-06 bar and runs
NON-MONOTONIC (3.78e-6 at i=46 below 1.11e-5 at i=50). The ideal-gas run has
Gamma_1 = 1.4728 identically and can have none of this.

That is a concrete, measured difference between the run that lives and the runs that die,
independent of the RT. It is **a strong correlate, NOT a proven cause** -- and it is entangled
with the floors yet again, so it does not cleanly separate "steep Gamma_1 breaks the
reconstruction" from "the floored top is chaotic and Gamma_1 merely reports it".

Next candidates, untested: PLM reconstruction of Gamma_1 / the HLLD wave speeds built from
it (the general path reconstructs stored p and Gamma_1 to interfaces); a semi-implicit RT
source in place of the clamp (the pgen already Newton-solves one in
`double_gray_two_stream_RT_source`); Ledoux vs Schwarzschild in the adjustment.

## ALSO REFUTED: cells per scale height does NOT discriminate

I over-read the nightside column dump (9-15 levels with H_p < dr) into a claim that the
general-EOS top is unresolved. Checked per column across the whole grid, from the bin
dumps, using the density scale height:

| run | unresolved (col,level) pairs | columns affected | worst column |
|---|---|---|---|
| **ideal + grey -- NEVER FAILS** | **2.06 %** | **52.9 %** | 12 levels |
| general + grey -- dies 2.51 | 1.99 % | 48.5 % | 16 |
| general + c-k -- dies 0.37 | 1.87 % | 36.7 % | 10 |
| general nx1=128 -- dies 0.395 | 1.04 % | 31.9 % | 18 |
| general `noh2` -- NEVER FAILS | **0.00 %** | 0.0 % | 0 |

The run that never fails is the MOST sub-scale-height of the four. The angular-mean
stratification is fine in every run (min 1.8-4.7 cells per H). So this is not the
discriminator -- do not re-propose it. It does add a third independent sign that `noh2` is
simply an easier problem: perfectly resolved, grad_ad 0.269, interior 14098 K.

## PICK UP HERE (paused 2026-08-23, fourth session)

Everything launched this session has finished and been read; the cluster is idle. Nothing
is mid-flight. Two commits shipped and pushed: **ea823e3a** (RT source limiter) and
**c1d513c2** (the two diagnostics). Write-up artifact:
https://claude.ai/code/artifact/c41bb913-2883-4932-b199-79170287aa20

**Failure one is CLOSED.** ideal + correlated-k: no radiative timestep constraint on an
explicit operator-split source, tripped by `pfloor = 1.0` pinning e at 2.115 erg/cm^3.
Fixed, A/B verified. Do not reopen.

**Failure two is OPEN and is where to resume.** The tabulated EOS blows up under any RT.
Nine hypotheses tested, nine refuted -- the full list with evidence is in the sections
above. **Read those before proposing anything; six of them were mine and three died on the
same day I proposed them.**

**The code states the assumption my measurement contradicts.** `src/mhd/mhd.hpp:588-596`,
documenting `WbPiecewiseLinearDerX1`: under a general EOS the Riemann solvers do NOT
recompute pressure from the reconstructed (d,e) -- they read it from `wder`, reconstructed
to interfaces independently -- and *"Gamma_1 is left alone: it is a smooth O(1) quantity
that is not stratified over many scale heights, so plain PLM is appropriate."* Gamma_1
jumping 0.44-0.58 across one cell is exactly that assumption failing. Recorded as a fact
about the code, next to the measurement; **the user said "let's not go that way"
2026-08-23, so do NOT open this without them raising it first.**

**The one untested lead, and the only clean way to test it.** The sole surviving
discriminator is Gamma_1: identically 1.4728 for the ideal gas, 1.084-1.666 for the table
with jumps up to 0.58 between ADJACENT cells, present already in the healthy state. The
general-EOS path reconstructs stored cell-centred p and Gamma_1 to interfaces and builds
HLLD wave speeds from them. The diffusivity-neutral test:

> reconstruct pressure and evaluate Gamma_1 AT THE INTERFACE from the reconstructed state,
> instead of reconstructing the stored cell-centred value. Same order, same dissipation, so
> a change in outcome can only come from the Gamma_1 handling.

**Donor-cell reconstruction is NOT a valid test of this** -- it changes the scheme's
diffusivity, and dissipation alone is already known to delay the failure
(0.370 < 0.485 < 0.584 < 1.072 tracks increasing dissipation). The user rejected it on
exactly that ground. Any test of the Gamma_1 hypothesis has to hold diffusivity fixed.

**Bars to clear.** The general + c-k baseline dies at 0.37-0.67 rot (~36 min of apu1 wall,
the fast reproducer, input `blowup.athinput`). General + grey dies at 2.509 rot. Clean means
reaching `tlim = 1.2e6` = 3.934 rot, which is what ideal+grey and both `eos_h2=false` runs
did. **Set tlim past 2.51 rot or the result is unreadable** -- a short tlim already produced
one wrong call this campaign.

**Binaries in `/viper/ptmp/jinma/claude_eos_gpu`:** `athena_srclim` = ea823e3a,
`athena_g1` = ea823e3a + the Gamma_1 dump. Anything new should be rebuilt from c1d513c2.
`export HSA_XNACK=1` or every GPU run dies. Job scripts: `sweep.sh`, `limtest.sh`,
`greylong.sh`, `g1dump.sh`, `coldump.sh`.


## HANDOFF DOC (2026-08-24): `docs/HANDOFF_dhj_ck_eos.md` -- UPDATED to c628293e, pushed

**Current version is c628293e**, not 397b49b1. The first version was written before the
evening's results and reasoned from single runs; the update leads with a "READ FIRST"
section carrying (a) the chaotic-ensemble caveat with the 0.370/0.485/0.584/0.685/1.464
spread, (b) orion's run actually dying at 1.464 rot after being stopped at 0.713 -- noted
against `cd2e8815`, whose commit message says the opposite -- and (c) the reopened
nondeterminism as the prime suspect. Section 4 is now the two races, the bisection table,
what the double-pack and write-count instruments exclude, and the checksum-per-task
instrument to build next; plus the hydro/weak-field cures. Gamma_1 is demoted, not dropped.

**The git channel to orion WORKS and is bidirectional**: orion pulled 397b49b1 and pushed
`cd2e8815` on top of it. Use it.

## HANDOFF DOC, first version (2026-08-24): commit 397b49b1

Written so a session on another machine (orion) can pick this up: what shipped, the two
failures and which is closed, all nine refuted hypotheses, the Gamma_1 lead and its only
diffusivity-neutral test, and **the exact input deltas** from the shipped
`inputs/mhd/deep_hot_jupiter_rt_eos.athinput` that rebuild the fast reproducer
(`x1max` 13.04e9, meshblock nx2/nx3 16, `max_eta` 1e13, `use_rkg_sts` false,
pfloor 1.0, dfloor 7.26e-12, `bbot` 1e1, `rt_ck` true) plus the one-line delta for every
single-variable variant. **Trap recorded there:** the live
`/viper/ptmp/jinma/claude_eos_gpu/blowup.athinput` carries `rt_de_max = 0.5`, which is the
limiter-ON variant, NOT the 0.674-rot baseline.

Cross-machine messaging does not work: `ListAgents` sees no peers and `orion` does not
resolve from viper. Git via the fork is the channel.
