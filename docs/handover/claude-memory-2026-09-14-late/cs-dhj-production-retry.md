---
name: cs-dhj-production-retry
description: The cubed-sphere dhj PRODUCTION retry -- the gate FAILED and the fix WAS active (beta COLLAPSES from 26.3 to 2e-3 as the flow develops; 3.8% of cells trigger it). The well-balanced source and the GS07 EMF each delay the death 24-50% and neither cures it
metadata:
  type: project
---

**2026-09-03: THE GATE FAILED -- AND IT TESTED NOTHING.**

## The beta story: inert at t = 0, ACTIVE by half-way

**First reading, then corrected.** On the production INITIAL state the minimum plasma beta
is **26.3**, 53x above the 0.5 threshold, 0 of 786,432 cells below it -- which looked like
proof the switch could never fire. It is not: **beta collapses as the flow develops.**
Measured on `cs_arm_wb`'s dumps (beta = 2p/B^2 from the raw face-frame components, the same
quantity the switch tests):

    t [s]        0     1e4     2e4     3e4      4e4      5e4      6e4      7e4
    min beta   26.3    7.26    1.88   0.754    0.103   0.0401   0.0020   0.0029
    % b<0.5     0       0       0       0      0.67%    3.03%    3.81%    3.75%

So the fallback is **inert for the first half of the run and then fires on a few percent of
the domain** -- and the run still dies. The gate DID test the fix. **My same-day retraction
("the gate tested nothing") is itself RETRACTED**; it was drawn from the t = 0 snapshot
alone. The lesson stands in a different form: a switch's trigger condition must be checked
over the RUN, not the initial state, in both directions.

Note the production run reaches beta 2e-3 -- squarely inside the regime the reproducers
blow up in ([[cs-mhd-instability-characterized]]) -- so the low-beta instability is
plausibly PART of this after all, just not the whole of it.

## THE ARM MATRIX (2026-09-03, apudev, 15 min each)

Baseline death **t = 6.4407e4 s = 0.211 rot**; run-to-run noise ~8.5%. All arms are the
production input on the 1d7c4e5a binary, dumping every 1e4 s.

    arm                          flag                          dies at      vs baseline
    baseline                     (cs_lowbeta_fallback = 0.5)   6.4407e4        --
    wb    well-balanced source   cs_wellbalanced_src = true    7.9602e4      +23.6%
    emf   GS07 corner EMF        cs_gs07_emf = true            9.6391e4      +49.7%
    both  wb + emf                                             NEVER DIED   >1.2517e5
    diss  wb + emf + HLLE EVERYWHERE (threshold 1e9)           NEVER DIED   >1.3061e5
    ppm4  wb + emf + reconstruct=ppm4, nghost=4                2.8133e4     -56%  (see below)

**`both` IS THE FIRST CONFIGURATION THAT DOES NOT DIE.** It and `diss` each ran their whole
13.5 min slot -- 0.41 and 0.43 rotations, twice the baseline -- with dt stable at ~4.0 and
RISING, and ZERO NaN rows. Each flag ALONE still dies. The two therefore act on largely
INDEPENDENT error channels (+24% and +50% alone, no death together).

**DISSIPATION IS NOT THE LEVER AT PRODUCTION BETA.** `diss` is `both` plus HLLE on every
face and it bought 4% over `both` -- noise. This INVERTS the conclusion of
[[cs-mhd-lowbeta-fix]], which was measured on the reproducer: there dissipation was the
only thing that worked; here the geometric source and the corner EMF are, and forcing
dissipation on top is wasted smearing.

## THE LONG RUNS: wb+emf SURVIVES ROTATIONS AT nx=32, AND FAILS UNDER REFINEMENT

**`bench/cs_both_chain` (nx2=nx3=32): FINISHED at 5.50 rotations, clean.** 212,200 cycles
over 12 chained debug slots, NO dt collapse anywhere, ZERO NaN, final dt 12.49 -- ABOVE the
sp_dhj_ctl control's 11.6, so the "stable but throttled" worry is retired: dt fell to ~4 in
the early transient, plateaued, then recovered fully. mass/mass0 = 0.99328, i.e. 5e-4 drift
over 5.5 rotations after the initial 0.63% settling.

**THE HEADLINE NUMBER: a factor of 20 in survival between two resolutions.**
nx=32 alive at 5.50 rot; nx=64 dead at 0.281 and 0.324 rot. Same configuration, same
binary, 4x the angular cells.

(historical, superseded by the line above:)
dt settled at ~4 and then ROSE to ~6.1; zero NaN, zero c2p failures, mass flat to 1e-5
after the initial 0.63% settling transient (the baseline shows the SAME 0.99375 at
t=3.05e4, so that transient predates the fixes). Reproduced from a COLD START, so the
debug-slot survival was not a fluke.

**`bench/cs_both_n64` (nx2=nx3=64, 96 blocks): DIED at t = 8.5845e4 = 0.281 rot.**
Refinement made the ORIGINAL failure worse, so this was the experiment that mattered, and
it says **wb+emf IS NOT A FIX -- it is RESOLUTION-DEPENDENT**:

    nx2=nx3=32   survives > 1.5 rot   (7x the baseline), dt rising
    nx2=nx3=64   DIES at 0.281 rot    (1.33x the baseline)

Same death shape as always: dt 5.1e-17 -> 5.8e-73, simulation time frozen at 8.584494e4,
then the post-mortem dt frozen at 10.44661 with NaN in the history. The floor counter had
gone NEGATIVE (integer overflow) in the last interval before the death.

### CONFIRMED by a cold-start repeat, and the death time is NONDETERMINISTIC

    trajectory                              dies at            rot
    nx=64 cold start #1                     8.584494e4        0.281
    nx=64 cold start #2 (repeat)            9.903212e4        0.324
    nx=64 restarted from its own 6e4 rst    alive at 1.016e5  0.333  (slot ended)
    nx=64 restarted again (ulp control)     alive at 1.004e5  0.329  (slot ended)
    nx=32                                   alive at 1.679e6  5.50   (12 slots, finished)

**The nx=64 failure is PROBABILISTIC, not a threshold**: two trajectories died at 8.58e4 and
9.90e4, two others ran past both points alive. The growth time is comparable to the run
length, so whether a given trajectory dies is decided by round-off.

Both cold starts DIE. The two runs are **bitwise identical through cycle 10300** and then
differ by **1 ULP in dt at cycle 10400** (4.597414 vs 4.597413) -- after which they separate
and die 15% apart. A one-ulp seed, almost certainly a non-associative GPU reduction (dt is a
global min), amplifying exponentially IS the instability: round-off does not grow like that
in a stable scheme. The restarted run is just a third perturbed trajectory that had not died
yet when its slot ended -- NOT evidence of stability.

**Two retractions of my own, both from single trajectories.** I wrote "SURVIVES REFINEMENT"
at 0.25 rot (it died at 0.281), then wrote "the claim is unsupported" when the restart ran
past the death point (it was simply a different trajectory). **On a marginally stable
problem, one run is never a result -- and neither is one restart.** Same error class as
reading a growth rate off a truncated window ([[cs-mhd-instability-characterized]]).

So the honest summary of the whole arm matrix: the well-balanced source and the GS07 EMF
are real, independent improvements that between them buy a factor of ~7 in survival time at
nx=32 -- and they do not cure the instability, which returns under refinement. The
character of the bug is unchanged from [[cs-mhd-dhj-blowup]]: **refining makes it worse.**

### How to run rotations without a long queue: CHAIN DEBUG SLOTS

The surviving config does **155 sim-s per wall-s** on 2 GPUs, so 2 rotations is ~66 min --
NOT the 4 h the original slot was sized for. `apudev` 15-min slots chained with
`--dependency=afterany` deliver that with no apu1 queue wait. Two requirements:
  * set the RESTART cadence short enough that a slot writes one (`output4/dt = 6.0e4`,
    ~0.2 rot). The production 3.05e6 (10 rot) writes NONE in a 15-min slot.
  * each slot picks the newest `rst/dhj.NNNNN.rst` and runs `-r <rst> -i <input>`; that
    combination is legal and the input is applied AFTER the restart.
Verified continuous across the handoff in both cycle count and time.

**dt is 2.6x BELOW the spherical-polar control.** `bench/sp_dhj_ctl` ran 101 rotations at
mean dt 11.6 with zero NaN. cs+wb+emf sustains ~4-6. Stable, but throttled -- the open
question, and the reason this is not yet called a fix.

**`bench/cs_dhj_e12` is NOT a counter-example.** It "reaches" t=8.64e7 (283 rot) but has
2833 NaN rows of 2834 and mass/mass0 = 0.0000: it died at once and ran to tlim as garbage
with dt frozen at 104. Third instance this session of "reached tlim is not ran".

## THE 1-ULP PERTURBATION TEST: nx=32 is CHAOTIC BUT BOUNDED; nx=64 arm was INERT

`bench/ulp/` -- control and a twin differing by ONE FLIPPED BIT in the restart file, from the
SAME state (t = 6e4) at each resolution, `ndiag = 1` so dt is logged every cycle. Divergence
measured as fit-free CROSSING CYCLES of |d(dt)|/dt, never a fitted growth rate.

    nx=32:  first difference  +1148 cycles
            > 1e-6            +1724
            > 1e-3            +5532
            > 1e-1            NEVER (19,280 cycles);  final 1.04e-2

**So nx=32 amplifies round-off to ~1% and STOPS there while the run continues** -- chaotic
but bounded, which is what a stable turbulent simulation looks like. This RETIRES the simple
story that nx=32 is insensitive and nx=64 is unstable: both are sensitive.

**The nx=64 arm is INCONCLUSIVE -- the perturbation was inert.** Identical method, but the
two runs stayed BITWISE IDENTICAL for 4700+ cycles. A single flipped bit at an arbitrary
offset is a LOTTERY: it landed in something that does not influence the solution (quiescent
interior, or a buffer overwritten from a neighbour). `perturb_rst.py` now has an
**`everywhere` mode** -- flip the low bit of every 1000th double across the data region,
still ~1e-16 relative, guaranteed to seed the active region at any resolution, and a better
model of round-off anyway. Re-run both arms with it before comparing resolutions.

    python3 perturb_rst.py <ctl.rst> <pert.rst> everywhere

**Traps this test walked into, both worth remembering.** (1) dt is a global MIN reduction, so
as a divergence proxy it is a STEP FUNCTION -- decades crossed "in 0 cycles" is quantization,
not instant growth; only the wide gaps carry rate information. (2) The restart's parameter
header is not a multiple of 8 bytes, so the doubles are NOT aligned to the file's 8-byte
grid; reading the wrong phase yields plausible IEEE garbage (1e+264, 1e-285). Scan all eight
phases and pick the one whose values vary smoothly.

## PPM4 IS RULED OUT HERE -- and it is the FLOORS, not the scheme

`reconstruct = ppm4` + `nghost = 4` on top of `both` died at t = 2.8133e4, less than HALF
the baseline, and the collapse was a SINGLE 100-cycle window from a healthy dt = 12.4 to
1.4e-8 -- the shape of a state excursion, not of the baseline's slow squeeze. The event-log
counters say why:

    arm     cycles   eos_efloor fires   per CELL per CYCLE
    ppm4     2153     1,763,097,492          1.04
    both     1539           340,336          0.00028

**The internal-energy floor fires about ONCE PER CELL PER CYCLE** -- 3700x `both` -- with
dfloor and tfloor similarly inflated. PPM4's overshoots leave the TABULATED EOS's valid
range and the floor then rewrites the state wholesale; on cs the C2P floor also overwrites
u.e ([[cs-mhd-c2p-floor-corrupts-ue]]). This is why ppm4 behaved sensibly on the cs TEST
problems: they are ideal gas, with no table edge to fall off. Higher-order reconstruction
is not refuted -- it is untestable here until the floor/table interaction is handled.

**The floor counters are a good cheap instrument** we were not using: `both` sits at
2.8e-4 efloor hits per cell per cycle. Compare any new arm against that.

**Both fixes buy real time and neither cures it.** That is a genuinely different result
from the test problems, where [[cs-mhd-instability-characterized]] measured both as
changing NOTHING -- so the production failure is not simply the reproducer's instability.

**Reading the death correctly:** dt collapses (1e-29 or smaller) and the run then CONTINUES
with dt frozen at exactly **10.4466** while the history fills with NaN. Do not read the
frozen-dt tail as "still alive" -- that mistake was made once already this session. The
death time is the first cycle with dt < 1e-3.

## The superseded first reading

`cs_lowbeta_fallback` switches to HLLE on faces with plasma beta below the threshold,
**0.5**. Measured on the production initial state
(`bench/cs_prod_mhd_gate/bin/dhj.mhd_w_bcc.00000.bin`, 786,432 cells, beta = 2p/B^2 from
the raw face-frame components -- the same quantity the switch tests):

    min beta      26.3          <- 53x ABOVE the threshold
    0.1 pct       31.0
    1 pct         48.6
    median        9.6e4
    beta < 0.5    0 cells       beta < 10: 0 cells;  beta < 100: 4.4%

**Zero faces can ever trigger it.** The gate therefore ran plain HLLD with an inert switch,
and reproduced the original failure because it WAS the original scheme. The "0.194 ->
0.211 rotations, +8.5%" comparison below is **RETRACTED**: that is run-to-run noise between
two identical schemes, not a partial cure.

**What this means physically.** The dhj blow-up is NOT the beta < 0.05 instability the
reproducers exhibit -- production beta never goes near it. By the force-budget table in
[[cs-mhd-instability-characterized]], at beta ~ 26 the scheme's spurious Lorentz force is
orders below the physical one, so a low-beta dissipation fix was never going to touch this.
[[cs-mhd-lowbeta-fix]] remains true for the test problems and is IRRELEVANT here.

**The missing instrument.** `cs_lowbeta_fallback.hpp` has **no firing counter**, so nothing
in the code could report that the switch was dead. A whole gate cycle was spent on it.
Cf. [[validate-the-instrument]]: check that the switch you are testing can engage BEFORE
reading its result.

## The original (now mis-framed) result

Measured on `apudev` (job 11353535, `bench/cs_prod_mhd_gate`, same binary 1d7c4e5a, same
production input), because the gate never needed the 4 h slot -- see "How to measure this
in 15 minutes" below.

    quantity                        old cs_prod_mhd     WITH the fallback
    -------------------------------------------------------------------------
    t where dt collapses            5.933084e4 s        6.440731e4 s
    in rotations                    0.194               0.211        (+8.5%)
    dt at collapse                  2.08e-68 -> frozen  2.57e-59 -> 0
    cycle of collapse               7600                7600
    history at t = 6.1e4            all NaN             FINITE

The one genuine improvement: the history row at t = 6.10e4 s, which is *all NaN* in the old
run, is finite here (dt 5.93, mass 3.43864e26), and dt tracked consistently ABOVE the old
run's for 6000 cycles (e.g. cycle 4800: 6.23 vs 3.64). So the low-beta fallback IS doing
real work. It is not enough: dt still falls off a cliff at t = 6.44e4 and the simulation
time freezes exactly as before. **This is the same pathology, delayed by 8.5%, not a fix.**

**Do not launch the 283-rotation campaign.** `bench/cs_prod_mhd_long` stays unlaunched and
its auto-launch watcher is stopped. [[cs-mhd-lowbeta-fix]] remains true for the TEST
problems and is now known to be false for the production case; the blow-up
([[cs-mhd-dhj-blowup]]) is **still OPEN**.

## How to measure this in 15 minutes, not 4 hours

The 4 h arm sat `PD (Priority)` for nine hours behind 17 higher-priority jobs. It never
needed four hours: **the failure point is reached at elapsed 332 s and the gate point at
539 s** -- read that straight off the original `cs_prod_mhd/log.out.*`. So run it on
`apudev`: 15 min limit, `debug` QOS at priority 1e6, 2 nodes that are usually free, and the
answer arrives in nine minutes. Use a SEPARATE directory so it cannot touch a watched
history file.

Levers that do NOT work, all probed with `sbatch --test-only`: shortening the wall limit
(4 h, 3 h, 2 h, 1.5 h, 1 h, 30 min all give an identical start estimate) and switching QOS
(a0001/a0004/a0008 identical). Priority here is ~97% fairshare; age contributed 12,677 of
3,556,788 in nine hours, so cancel-and-resubmit costs essentially nothing -- and gains
nothing either.

**Slurm will not let you swap a queued job's batch script** (`Update of this parameter is
not supported: Command=`), though it WILL let you raise TimeLimit. So "make the 4 h job a
24 h job" is a trap: the stored script still runs `athena -t 3:55:00`, the code exits
gracefully at 3 h 55, and with restarts every 10 rotations (4 h buys ~2) it writes none.
To lengthen a queued run you must cancel and resubmit a script with the right internal -t.

### The ORIGINAL plan, kept for the record

**STATUS 2026-09-03 18:40 -- STILL PENDING, THE GATE IS UNMEASURED.** The job sat `PD
(Priority)` on apu1 for nine hours and has never started; Slurm now estimates
**2026-09-04T01:11** on vipa1283. The first watcher's "gave up after 8 h without meeting
the gate" line is therefore NOT a failure of the fix -- it is a watcher that timed out
while its job was still in the queue. Nothing about the fallback has been tested on the
production problem yet.

    dir     /viper/u2/jinma/ATHENAK/bench/cs_prod_mhd_fix
    binary  built from build_dhj_gpu at HEAD 1d7c4e5a (the later commits are comments only)
    slot    4 h, apu1, 2 GPUs, athena -t 3:55:00;  input tlim is the full 8.64e7
    check   tail dhj.mhd.hst  -- NaN anywhere means it failed;  1 rotation = 3.053e5 s

The input is the production one, unchanged except for an explicit (default-valued)
`cs_lowbeta_fallback = 0.5` line, so the run IS cs_prod_mhd with the fix in.  Restarts are
every 10 rotations, so this 4 h slot (~2 rotations) writes none.

## THE CAMPAIGN SLOT LAUNCHES ITSELF -- read watch.log FIRST

`bench/cs_prod_mhd_fix/watch_and_launch.sh`, started with `setsid` so it outlives the
session, polls the history every 2 minutes and `sbatch`es
`bench/cs_prod_mhd_long/submit.sh` (23:55, same binary, same input, from t = 0) **only if
all three hold**: history reached t >= 1.2e5 s (0.4 rot, TWICE past the 0.2 rot failure),
no NaN anywhere in the history, and dt >= 15 s (hydro control holds 22-25; the failure was
squeezed to 9.15 then FROZE at 10.4466). NaN or a collapsed dt aborts it without launching;
so does 8 hours without the gate being met. Every decision, with its reason, is appended to

    /viper/u2/jinma/ATHENAK/bench/cs_prod_mhd_fix/watch.log

**So the first thing to read next session is that log**, not the queue: it says whether the
fix held on the production problem and whether the campaign is running. The gate thresholds
were dry-tested before launch (1.30e5/25.0 -> LAUNCH, 1.30e5/9.1 -> abort, 1.19e5/25.0 and
6.10e4/25.0 -> keep waiting).

### Round 2 (the live watcher): `watch_and_launch2.sh`

The round-1 watcher **counted wall clock, not job state**, so a queue wait it could not
control burned its whole budget and it reported "gave up" for a job that had produced no
history at all. `watch_and_launch2.sh` (started 2026-09-03T18:38, `setsid`, same three gate
thresholds, 5-minute poll) instead polls `squeue -h -j 11334526 -o %T` and only quits when
the job has **left the queue**, with a 30 h backstop that covers the 01:11 start plus the
4 h slot. Same log file, same auto-`sbatch` of the campaign slot on PASS.

**The general lesson: a watcher on a batch job must key on the JOB, not on a wall clock**
-- otherwise a scheduler delay is indistinguishable from a physics failure in its log.
Same family as "a job that COMPLETES with exit 0 is not a job that ran"
([[cs-dhj-long-run]]) and [[validate-the-instrument]].

## Why

[[cs-mhd-lowbeta-fix]] fixed the instability on the test problems, but the run that started
the whole thread -- the cubed-sphere deep-hot-Jupiter production case, which blows up at
**~0.2 rotations** ([[cs-mhd-dhj-blowup]]) -- has never been retried with the fallback.
That is the only gate that matters for the science.

## The recipe

Old run: `/viper/u2/jinma/ATHENAK/bench/cs_prod_mhd` -- `deep_hot_jupiter.athinput`,
`submit.sh`, and `dhj.mhd.hst` whose last rows are all NaN. Grid nx1 = 128, nx2 = nx3 = 32
(cubed sphere), meshblock 128x16x16 = 24 blocks, general EOS + table, plm, hlld,
`ohmic_resistivity = eos`, tlim 8.64e7 (283 rotations).

1. Build a GPU binary at HEAD with `-D PROBLEM=deep_hot_jupiter_rt` and the flags in
   [[viper-hip-build-recipe]] (`build_cs_gpu` already exists in the repo and is the
   natural one to refresh). The login node viper13 has no GPU: it compiles, it cannot run.
2. New bench dir, copy the input and submit.sh unchanged. **The input needs no edit** --
   `<mhd>/cs_lowbeta_fallback` defaults to 0.5 on a cubed-sphere mesh, so the fallback is
   already on; record it explicitly in the input anyway so the run is self-documenting.
3. Run a SHORT arm first, `time/tlim=1.5e5` (0.5 rotation), 2 GPUs on `apu1`, a couple of
   hours. `apudev`'s 15-minute limit only reaches ~0.12 rot, which is BEFORE the failure
   and would prove nothing -- cf. [[validate-the-instrument]] and the partial-arm trap in
   [[cs-mhd-dhj-blowup]].

    1 rotation = 3.053e5 s.  The failure was at ~0.2 rot = 6.1e4 s.

## The gate

* **passes t = 6.1e4 s (0.2 rot)** at all, which the old run did not;
* **dt does not collapse**: the hydro control holds 22-25 s, the dying MHD run was squeezed
  28.98 -> 9.15 and then FROZE at 10.4466 while grinding on with a NaN state;
* `eos_efloor` does not run away (it integer-overflowed to -1135976192 in the failure);
* mass drift comparable to the hydro control's -0.71 % over 11.2 rotations.

A job that completes with exit 0 is not a job that ran -- check the .hst for NaN
([[cs-dhj-long-run]]).

## Then

* If it survives, consider also turning on `<mhd>/cs_wellbalanced_src` (c3992145): it buys
  10x on the residual once the instability is gone, but costs 1.68x wall clock until its
  face triads are cached per (m,k,j). Cache them first.
* `<mhd>/cs_gs07_emf` (2566a9f7) is the more correct CT scheme and worth 30x on its own,
  but adds nothing once the fallback is on. Leave OFF until it has a convergence gate.
* The threshold 0.5 has ~10x margin over the measured onset (beta ~ 0.05). The dhj
  atmosphere reaches beta ~ 1e-4 at the top, so most of the magnetised region uses HLLE
  either way; lowering it toward 0.1 would recover some HLLD accuracy if that matters.

## 2026-09-04 (later session): PRODUCTION RERUNS LAUNCHED WITH THE ROTATION FIX

Results in: **cs_rotfix (nx=64, full rotation) is ALIVE at 1.18 rot**, dt recovered 1.60 ->
2.73 s, KE tracking cs_rot32; cs_rot32 clean at 2.68 rot. The user asked whether cs is
ready for production from scratch; answer YES for hydro (it never failed -- cs_prod_hyd
158 rot, cs_dhj_hyd 283 rot -- but both ran the BETA-PLANE rotation and are physically
wrong), cautiously yes for MHD at nx=32 (12x past the old death, the run IS the test).
Submitted on the aa6ddb9b HIP binary (build_hip_diag), apu1, chained 24 h slots:
* `bench/cs_prod_hyd_rot` (jobs 11409652-53): cs_prod_hyd's input + cs_full_rotation=true.
* `bench/cs_prod_mhd_rot` (jobs 11409654-56): cs_rot32's validated input (wb src, GS07
  EMF, lowbeta 0.5, full rotation, EOS resistivity) at production cadence (bin 2 rot,
  rst 10 rot, hst/log 10 per rot).
Success = dt recovering toward ~12 s and no NaN; compare against cs_rot32 to 2.68 rot.
