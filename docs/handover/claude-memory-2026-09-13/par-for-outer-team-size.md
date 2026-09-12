---
name: par-for-outer-team-size
description: "AthenaK's par_for_outer asks Kokkos::AUTO for the team size, which picks 256 without knowing the inner loop is only ncells1 long. Optimal is ceil(ncells1/64)*64. But on the PRODUCTION geometry (nx1=128) it is worth only 3% of wall time. DECIDED 2026-08-21: not worth doing. Do not re-propose."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T23:30:00.000Z
---

**Nothing was changed.** `src/athena.hpp` is untouched. This is a recorded finding for a
possible future upstream PR, not a pending edit.

## The finding

`par_for_outer` in `src/athena.hpp` builds every team policy as
`Kokkos::TeamPolicy<>(exec_space, nkj, Kokkos::AUTO)`. AUTO sizes the team from scratch and
register pressure -- it cannot know that the `par_for_inner` inside will only have
`ncells1 = nx1 + 2*nghost` items to hand out. On MI300A it picks **256**.

With `ncells1 = 68` that means wave 0 gets 64 items, wave 1 gets 4, and waves 2 and 3 get
none and never issue. VALUUtilization over issuing waves is (64+4)/128 = **53.1 %**, which
is exactly what the counters read for all three MHD flux kernels (53.46 / 53.56).

**The optimal team size tracks ncells1: `team = ceil(ncells1/64)*64`.** Measured, per-cell
`CalculateFluxes` cost (job 10981149):

| ncells1 | AUTO(256) | team 128 | team 64 | best |
|---|---|---|---|---|
| 64 (nx1=60) | 3.7516 | 2.9669 | **2.5722** | 64 |
| 128 (nx1=124) | 2.9168 | **2.4972** | 2.9121 | 128 |
| 256 (nx1=252) | **2.2193** | 2.6377 | 3.5837 | 256 = AUTO |

Note 68 rounds UP to 128, not down to 64: 68 items over 64 threads is two ragged passes.
That is why team 128 beat team 64 at the production grid (373 vs 405 ms) while team 64 wins
at ncells1 = 64 exactly.

## DECISION: NOT WORTH DOING (2026-08-21)

Measured in WALL time, not kernel time, and on the production geometry:

| geometry | AUTO(256) | team 192 | team 128 |
|---|---|---|---|
| nx1=64 (ncells1=68), the A/B config | 5.147 s | -- | **4.791 s (1.074x)** |
| **nx1=128 (ncells1=132), PRODUCTION** | **7.486 s** | **7.260 s (1.031x)** | 7.629 s (1.9 % SLOWER) |

**About 3 % on the real grid.** The 1.11x figure quoted earlier was GPU KERNEL time on the
small A/B config and overstated it by nearly 4x, for three reasons: kernel time is only
~71 % of wall; production runs nx1=128 not 64, so ncells1=132 already gives AUTO 69 %
utilisation instead of 53 %; and the gain scales inversely with nx1.

3 % does not justify a change to `src/athena.hpp` that touches every problem generator and
needs the full test suite. **The user decided against it. Do not re-propose it** unless
someone is already working in `par_for_outer` upstream, or a run appears with much smaller
blocks -- at ncells1 = 64 it was 1.46x on the flux kernels, so heavy decomposition or AMR
would see far more.

Note team 128 is WORSE than doing nothing at ncells1 = 132. Any fixed constant is wrong;
only the ceil(ncells1/64)*64 rule holds.

## What it is worth, on the small A/B grid only

nx1 = 64, so ncells1 = 68, team 128 instead of AUTO(256), job 10981105:

| | AUTO(256) | team 128 |
|---|---|---|
| CalculateFluxes | 461.1 ms | **372.9** (1.24x) |
| x1 kernel (10368 teams) | 137 | **90** (1.52x) |
| x2/x3 kernels (576 teams) | 324 | 283 (1.14x) |
| all GPU | 912.8 ms | **820.1** (1.11x) |

**Bitwise identical output** -- team size does not change the answer.

The gain scales INVERSELY with nx1 and vanishes by nx1 ~ 250, so it is worth most to
small-block and heavily decomposed runs. It is grid-specific, which is the whole point: the
fix is not a constant but deriving the team size from the inner-loop length.

## Why x1 gains more than x2/x3

x1 launches 10368 teams and has a queue, so halving the team size puts more teams resident
AND stops two waves in four idling: 1.52x. x2 and x3 launch only **576 teams for 304 CUs**
-- 1.89 teams/CU, whole grid resident, nothing queued, occupancy 1.89*4 = 7.58 waves/CU
against 7.63 measured. They are starved of work, not badly shaped, so team size helps them
much less.

They launch over two indices and sweep the third serially so the reconstruction scratch can
be permuted (`wl`, `wl_jp1`, `wr` rotate; each face reconstructed once, reused by both
neighbours). Parallelising over the swept index means reconstructing every face twice.
Per direction they are only 17-21 % slower than x1, not the 2.4x an earlier framing of
"326 vs 137 ms" suggested -- that compared two kernels against one.

## If this is ever taken upstream

It is `src/athena.hpp`, shared by every problem generator, so it needs the full test suite.
Frame it as "AUTO is being asked a question it lacks the information to answer, and the
caller has it", not as a tuned constant.
