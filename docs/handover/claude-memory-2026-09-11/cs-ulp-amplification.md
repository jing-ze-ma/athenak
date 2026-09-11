---
name: cs-ulp-amplification
description: The 1-ULP instrument, everywhere mode -- round-off grows to O(1) and KILLS nx=64 while it saturates near 1% at nx=32; nx=64 is dynamically unstable, and the restart-perturbation script must skip the tree metadata
metadata:
  type: project
---

**STALE AS OF 2026-09-06: measured on a pre-979edada binary, BEFORE 13a97399 (stretched source), c5c85e3b (stretched resistivity was ANTI-DIFFUSIVE near the top -- exactly the low-beta shell where this grew), 979edada, 248f1b77/551cd78a (bottom ghost drift, IC frame). RETEST: bench/cs_n64_fixed, job 11517116 (production input + nx2=nx3=64, 1 rot, bins 0.1 rot, binary 551cd78a); gate = alive past 0.35 rot with monotone dt. RESULT (22:10): ALIVE at rot 0.50, dt settled at 2.3 s (10.4 at 0.1, 3.9 at 0.3, 2.15 at 0.4, 2.28 at 0.5 -- falling then flat, NOT the 4.28 -> 0 collapse), mass flat at 3.455e26, KE1 flat 4.2-4.5e32, ME growing smoothly 4.7e31 -> 6.6e32. The dfloor counter wraps past 2^31 (1.4e9 by cycle 25k: ~2 % of the 3M cells per cycle sit at the density floor at the top, routine at this size, not a death). So on the fixed code nx=64 does NOT die at 0.3 rot: the old instability was (at least in part) the pre-fix bugs. Old claim STALE. 3.5 h slot ends ~0.6 rot.**

**Measured 2026-09-04.** Instrument: `bench/perturb_rst.py <src> <dst> everywhere`, which flips
the low mantissa bit of every 1000th double in a restart's field data (10,159 doubles at
nx=32, 40,635 at nx=64 -- a ~1e-16 relative perturbation, and a better model of round-off than
one bit in one cell). Runs live in `bench/ulp/{n32,n64}_{ctl,pert}`, compared with
`compare.py`, which reports FIT-FREE crossing cycles rather than a growth rate.

Both arms are the **wb+emf** configuration of [[cs-dhj-production-retry]] --
`cs_wellbalanced_src=true`, `cs_gs07_emf=true`, plm + hlld -- and the two inputs differ ONLY
in nx2=nx3 (32 vs 64). 1 rotation = 3.050e5 s; the restarts start at t = 6.0e4 = 0.197 rot.

## The answer: amplification IS what distinguishes the resolutions

    |d(dt)|/dt        nx=32                    nx=64
    first nonzero     cycle 6325               cycle 7482
    > 1e-6            6624                     8435
    > 1e-3            10921                    10245
    > 1e-1            19167 (a TRANSIENT        13316  -- and this is the DEATH
                      excursion, max 1.014e-1)
    end of slot       1.167e-2, BOTH ALIVE     dt -> 0.0, run destroyed
                      to 0.482 rot             at t = 9.4243e4 = 0.309 rot

At nx=64 the seeded noise reaches O(1) and the run dies: dt falls 4.28 -> 0 in 140 cycles
(4.28 -> 0.268 in ONE cycle at 13316), and `eos_efloor` integer-overflows to -1.28e9, the same
signature as every other death in this thread. Its own control was still healthy at those
cycles (dt = 4.27, monotone) and ran on to 0.335 rot, stopped by the wall clock, not by a
death. At nx=32 the identical perturbation saturates near 1% and both trajectories finish the
slot alive.

**So nx=64 is dynamically unstable in the sense that matters**: round-off-level noise is
amplified to O(1) within a third of a rotation, while at nx=32 it is bounded. The death at
0.309 rot lands inside the 0.281-0.324 rot band where the un-perturbed nx=64 arms died, which
is what makes the failure PROBABILISTIC rather than a threshold -- the arms that "ran past
both points alive" were simply luckier draws from the same distribution.

Caveat kept honestly: ONE perturbed trajectory per resolution ([[validate-the-instrument]]).
What is not N=1 is the band -- 0.281, 0.309, 0.324 are now three independent deaths.

**Correction to the earlier note**: nx=32 DOES cross 1e-1, once, at cycle 19167, then decays
back to 1.2e-2. The previous "never reaches 1e-1" came from the single-site arm, which is a
weaker seed. The claim that survives is *bounded and recovering*, not *never*.

## The trap that cost a submission: perturb_rst.py corrupted the TREE

The first `everywhere` run died on cycle 0 in both arms:

    n32: meshblock_tree.cpp:450  Neighbor search failed; MeshBlockTree broken.
    n64: build_tree.cpp:579      Tree reconstruction failed. tree=168, file=96

The sweep started at `<par_end>` and walked the whole file, so it flipped bits inside the
**LogicalLocation array**. `find_field_offset`'s "does this decode as a plausible double" test
CANNOT distinguish packed integer metadata from a field array -- the single-site version
escaped only because it landed 55% into the file. Fixed by `find_data_start()`, which does not
assume the header layout but IDENTIFIES it: `data_size` is the unique 8-byte word w in the
header with `filesize - (offset+8) == nmb_total * w`, and the function refuses if it is not
unique. Verified with `cmp -l`: the first differing byte is exactly the first byte of the data
region, and the differing-byte count equals the perturbed-double count (one low byte each).

Two runs exiting in 3 s would have read as "the perturbation is inert" if only `compare.py`
had been consulted -- it reports "bitwise identical", not "the job never started". Always
check that the perturbed run's FIRST LOGGED CYCLE matches the control's.

Also superseded: the single-site arm's nx=64 result ("bitwise identical for 9148 cycles") was
a dud landing, not a physical statement. Its logs are archived in `*_pert/single_site_arm/`.

## What this rules in and out

The next suspect named before this test -- the EOS floors -- is NOT exonerated but is now
downstream: amplification is the ROUTE, the floors are WHERE IT LANDS. A useful next question
is where in the domain the perturbation grows (cube vertex? the low-beta top?), which needs a
state-difference norm by region, not the dt scalar this instrument uses.

## WHERE it grows: a low-beta radial SHELL, not the cube vertex

`bench/ulp/where.py` bins the per-cell relative state difference between the matched ctl/pert
binary dumps by radial index, by Chebyshev distance to the nearest panel corner (a cube
vertex), and by panel. It refuses any dump pair whose CYCLES differ -- past the first
divergence the two dumps are no longer the same instant, and by then the difference has
saturated at O(1) anyway (at n32 dump 00003, t=9.15e4, max rel = 24 and the p99.9 is 6.8:
fully decorrelated turbulence). **Only the earliest matched dump can localise anything**;
everything here is dump 00002, t = 6.1001e4, ~150 (n32) / ~190 (n64) cycles after the seed.

By radial index (128 cells, 0 = bottom), max rel diff:

    n32   5e-10  3e-10  5e-10  1.7e-05  8.6e-06  1.9e-03  1.6e-04  4.8e-05
    n64   9e-10  2e-09  1.1e-07 1.9e-04  5.0e-04  3.7e-02  4.1e-03  7.1e-03

The bottom three bins sit at round-off in both. Growth switches on at bin 3-4 and peaks at
bin 5 (i = 80-96), NOT at the top of the domain. That shell is exactly where the field starts
to win: horizontally, median beta falls 1.2e8 (bottom) -> 2.3e3 there, and MIN beta crosses 1
at bin 3 (1.02) and is ~1e-2 from bin 4 up. So the amplification lives in the low-beta region
-- the same region [[cs-mhd-low-beta-divergent]] identified -- and it grows there even though
`cs_lowbeta_llf` (threshold 0.5) is active throughout it.

By distance to a cube vertex (n64, mean over the bin): 1.9e-6 at 0-2 cells, 2.3e-6 at 2-4,
then 4e-7 to 1e-6 further out -- a factor 2-5 enhancement in the mean at the corner, while
the GLOBAL max sits 30 cells away from any corner. **The vertex is mildly hot, not the
source.** For this route the blow-up is not a cube-vertex phenomenon, which is a real
distinction from [[cs-mhd-dhj-blowup]].

Panel 0 stays at round-off (1.5e-8 at n64, 2.5e-10 at n32) while the other five reach
4e-3..3.7e-2, so one panel is quiet even though the seed is uniform over the whole file.

The resolution comparison, in the growth shell at the same physical time and the same seed
FRACTION (1 double in 1000, 1 ULP each, at both resolutions):

                       mean       p99.9      frac of cells > 1e-6
        nx=32        1.87e-07    1.80e-05           0.56 %
        nx=64        5.43e-06    8.92e-04           7.2  %
        ratio           29x         49x             13x

That is the amplification difference measured on the STATE, ~190 cycles in, independent of
the dt proxy -- and it is already there long before dt notices (dt first differs at cycle
6325, a thousand cycles after this dump).

**Trap for the next reader**: the dumps are single precision, so any relative difference
below ~1e-7 is quantisation, not signal. The n32 percentiles above are close to that floor;
the n64 ones are not. Do not read the bottom radial bins as a measurement of anything but
"below the output's resolution".

## CORRECTION: the SEAM, not the vertex -- and the metric bug that hid it

The section above binned by `max(d2,d3)`, the Chebyshev distance to the nearest panel
CORNER. Distance to the nearest panel EDGE -- the SEAM -- is `min(d2,d3)`, and it was never
computed. "The vertex is mildly hot, not the source" was therefore a statement about the
vertex with NO seam measurement behind it. Do not bin a cubed sphere by one metric and
conclude about the other.

Re-done with both metrics, and with a RATE (fraction of cells above a fixed threshold)
rather than a max, because the vertex bins hold ~30x fewer cells than the seam bins and a
max grows with sample size all by itself. In the growth shell i=80-96:

    n32 (mean rel 1.9e-7, LEAST evolved)      >1e-6    >1e-5    >1e-4
      within 2 cells of a VERTEX              0.0007   0.0000   0.0000
      within 2 of a SEAM, >8 from a vertex    0.0133   0.0042   0.0011
      4-8 cells from a seam                   0.0027   0.0002   0.0000

    n64 (mean rel 5.4e-6, 30x further along)
      within 2 cells of a VERTEX              0.1120   0.0632   0.0267
      within 2 of a SEAM, >8 from a vertex    0.0698   0.0307   0.0082
      4-8 cells from a seam                   0.0587   0.0246   0.0063
      deep interior (>16 from a seam)         0.0661   0.0253   0.0056

**At the least-evolved sample the SEAM is the hot region (5x at 1e-6, 20x at 1e-5) and the
VERTEX is the QUIETEST thing measured -- zero cells above 1e-5.** At n64 the disturbance has
filled the domain (seam ~ interior) and the vertex leads by only 1.6-4x. Read together: a
seam-seeded disturbance that spreads, with the vertex hot LATE rather than first. Caveat:
n32 and n64 are different resolutions, not two times in one run, so "early vs late" is an
inference from the amplitude, not a controlled statement.

This is the SECOND independent reason to stop treating the cube vertex as the established
location of the dhj blow-up; the first was the retracted vertex order deficit
([[cs-mhd-instability-characterized]]). Supersedes the vertex framing in
[[cs-mhd-dhj-blowup]] for the ROUND-OFF-GROWTH route specifically.

Why it matters for the sp-vs-cs puzzle: the instability is NOT gnomonic-specific on the
idealised tests (sp_lowbeta reproduces it), but the PRODUCTION failure is -- sp_dhj_ctl runs
clean where cs dies. A seam-seeded mechanism is a difference spherical polar structurally
does not have: no seams, and its one singular structure is the polar axis, handled by a
different mechanism. That makes the seam halo -- the gnomonic basis transform plus the
along-seam resample -- a better suspect for the production route than the vertex.

## And it is NOT where beta is lowest

Growth vs per-cell beta at the same dump (n64), mean relative difference:

    log10 beta   -3..-1   -1..-0.3  -0.3..0.3   1..2     2..3     3..4    >6
    mean rel     3.5e-16   1.2e-09   1.3e-07   1.5e-06  3.2e-06  3.0e-06  5.8e-10

It PEAKS at beta ~1e2-1e4 and is dead below beta 0.1 and above beta 1e6 -- a band, not a
monotone trend, and the opposite of what "a low-beta instability" predicts. Confound to
keep in mind: beta < 0.5 is exactly where `cs_lowbeta_llf` acts, so the quiet low-beta
cells may mean "the fallback is working", not "low beta is safe". Second caveat: round-off
amplification in a turbulent shear flow is ordinary Lyapunov growth, and beta 1e2-1e4 at
i=80-96 is also where the zonal flow is fastest, so this maps where NOISE grows, which need
not be where the run DIES.

Localising the death directly FAILED: the only post-death dump is 100 % NaN in all 96
blocks and all 128 radial cells, written ~2400 cycles after dt hit zero, and the run's one
checkpoint is from the same moment. That would need a replay with dense output around cycle
13300, and replays here are only conditionally deterministic.

### The control that stops "seam" from being over-claimed

A panel seam is a block face too, so the seam finding needs the ordinary-block-face control:
same halo exchange, NO gnomonic transform and no along-seam resample. Same shell, same
rates (n32 is the least-evolved and therefore the one that localises):

    n32                                        frac>1e-6  frac>1e-5
      within 2 of a SEAM (transform+resample)    0.0133     0.0042
      within 2 of a same-panel BLOCK face        0.0056     0.0007
      block interior, >4 cells from any face     0.0017     0.0002

    n64
      within 2 of a SEAM                         0.0698     0.0307
      within 2 of a same-panel BLOCK face        0.0893     0.0411
      block interior                             0.0603     0.0229

So at n32 the growth is BOUNDARY-ASSOCIATED, not seam-exclusive: an ordinary block face is
already 3.3x the block interior, and the seam is a further 2.4x (1e-6) to 6x (1e-5) on top
of that. The right statement is "hottest at the seam, elevated at every halo boundary",
NOT "the seam machinery is the source". At n64 the contrast has washed out entirely and the
ordinary block face even edges ahead -- another reason to localise only on the least-evolved
sample. A generic halo/reconstruction effect would also exist on spherical polar, which does
NOT fail in production, so the seam excess is the part of this that could explain sp vs cs --
but 2.4-6x on one snapshot at one resolution is a lead, not a mechanism.
