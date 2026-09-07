---
name: rt-chain-parallel-split
description: "2026-08-21: COMMITTED as 38311a8a on polar-average-perf (not pushed). Splitting the dhj RT into 3 kernels with the chain-block index as a parallel dimension gives 3.3x at 88 chains with table lookup (2087 -> 632 ms/100cyc) and 1.18x at production, bitwise identical at RT_NB=1. Behind problem/rt_split, default false."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T17:20:00.000Z
---

Fixes the starvation documented in [[rt-kernel-occupancy-limit]]. Follows
[[inflight-rt-band-scaling]].

## What it is

`problem/rt_split` (default false) in `src/pgen/deep_hot_jupiter_rt.cpp` selects a
three-kernel RT in place of the monolithic `par_for("2stream_rt", ...)`:

- **A `rt_pre`** `(m,k,j)` -- chain-independent: geometry, picket-fence coefficients, the
  3 V-band stellar down-sweep, `tau_down_r_f`, `B`, `Q_v`, written to global arrays.
  Body is the monolithic kernel's code SPLICED VERBATIM; the private arrays are replaced
  by `Kokkos::subview(tau_g, m, k, j, ALL)` etc. so the arithmetic cannot drift.
- **B `rt_chain`** `(m,blk,k,j)` -- one thread per column per block of RT_NB chains.
  `nblk = ceil(nchain/RT_NB)` times more threads than the monolithic kernel. Each block
  writes its own flux slot `rt_Fb(m,blk,k,j,i)`.
- **C `rt_apply`** `(m,k,j,i)` -- sums the block slots IN BLOCK ORDER starting from 0.0,
  exactly as the serial accumulator did, then applies the flux divergence.

New compile-time `RT_NNC` (default 72) sizes the private `I_ir_down_c` column in kernel B
instead of the monolithic `NN = 270`; there is a runtime fatal if `n1 > RT_NNC`. The
18 kB/thread private segment did not matter while occupancy was 0.5 waves/CU but does
once thousands of waves are resident.

## MEASURED (jobs 10978261, 10978313, 10978357, 10978399). RT ms per 100 cycles.

| config | 4 chains | 88 chains | 4 + table | 88 + table |
|---|---|---|---|---|
| monolithic RT_NB=4 (committed) | 124.8 | 1157.4 | 136.2 | 2087.1 |
| split RT_NB=1 | 108.3 | 879.9 | 119.9 | 915.4 |
| split RT_NB=2 | 123.7 | 630.3 | | |
| split RT_NB=4 | 158.1 | 562.7 | 171.5 | **632.5** |

**3.3x at 88 chains with the table lookup on** -- the case a real correlated-k scheme
pays. 2.06x without it. At production (4 chains) split RT_NB=1 is 1.15x FASTER than the
committed kernel and bitwise identical, so it is a free win independent of correlated-k.

Per-kernel at 88 chains, RT_NB=4: `rt_pre` 67.6, `rt_chain` 486.3, `rt_apply` 8.8.
**`rt_pre` is now the floor** -- it does not scale with chains but it is still the old
occupancy-starved column sweep.

**More parallelism past RT_NB=4 makes it WORSE** (RT_NB=1 gives 2816 workgroups and 880 ms
vs RT_NB=4's 704 workgroups and 563 ms). Per-thread ILP still beats extra waves, the same
conclusion as the block-size sweep in [[inflight-rt-band-scaling]]. RT_NB=4 stays optimal.

## Correctness

The whole configuration is reproducible: monolithic vs itself and split vs itself are both
bitwise identical at nlim=1 and nlim=200, so every difference below is real, not noise.

**Split at RT_NB=1 is bitwise identical to `fixed200.a` at 4 chains AND at 88 chains, with
the table lookup on.** That is the proof the refactor is faithful: RT_NB=1 sums the chains
left to right from 0.0, exactly the association the serial kernel used.

RT_NB=2 and RT_NB=4 split differ from the reference in the last bits (7.7e-7 relative on
bcc3 after 200 cycles, from a 1.8e-15 absolute difference on a near-zero field component;
`eint` is still bitwise identical at cycle 1). For RT_NB=2 this is PROVABLY just a
different summation association -- `(t0+t1)+(t2+t3)` instead of `((t0+t1)+t2)+t3`. The
RT_NB=4 case has the same association as the reference, so its residual is compiler
codegen (FMA contraction differing between the two kernel bodies), not a logic error.
Not worth chasing: a real correlated-k scheme has different physics anyway, and RT_NB=1
already establishes fidelity.

## Consequence for the correlated-k decision

88 chains with table lookup was 2087 ms/100cyc = wall ~2.3x, RT ~60 % of the run. With the
split it is 632 ms = **wall ~1.34x, RT ~31 %**. Add RT sub-cycling (RT runs once per RK
stage, 200 calls per 100 cycles -- halving that is standard practice) and it is ~1.15x.
Correlated-k is now a much cheaper proposition than the earlier measurement implied.

## Caching the layer coefficients: TRIED, IT LOSES (job 10978513)

`RT_CACHE` (compile-time, default 0) caches `e0`, `alp`, `bet` from the down-sweep so the
up-sweep needs no `expm1`, no divide, and with the table on no lookup and no two `log`s.
The identity is real -- down-sweep layer `i-1` IS up-sweep layer `i` -- and it is
**bitwise exact**: cached RT_NB=1 split matches `fixed200.a` at 4 and 88 chains with the
table on and off, and cached RT_NB=4 matches uncached RT_NB=4 exactly.

But it costs three more private `[NC][NNC]` arrays, and private memory here is
scratch, i.e. global memory. RT ms/100cyc at 88 chains:

| | ktab=false | ktab=true |
|---|---|---|
| RT_NB=4 uncached | **565.0** | **630.3** |
| RT_NB=4 cached | 1226.9 | 1243.1 |
| RT_NB=1 uncached | 879.9 | 915.4 |
| RT_NB=1 cached | 907.9 | 963.8 |

**4x the private footprint (2304 -> 9216 B/thread) costs 2x the time; the saved
transcendental is nowhere near enough to pay for the scratch traffic.** At RT_NB=1 the
footprint stays small and the trade is roughly neutral -- slightly better at 4 chains,
slightly worse at 88. Two independent points, both saying footprint drives this kernel,
not ALU. Consistent with RT_NB=8/16 losing in [[inflight-rt-band-scaling]].

**Do not retry this.** An `e0`-only variant (one extra array, recompute `x`) is untested
but on the same curve, so expect a loss too.

**Best configurations, both measured:**
- production, 4 chains: split + RT_NB=1 + RT_CACHE=1, **105.4 ms vs the committed
  124.8 (1.18x), bitwise identical.**
- 88 chains + table: split + RT_NB=4, uncached, **630.3 ms vs the monolithic 2087.1
  (3.3x).**

**Levers still unbuilt:** parallelizing `rt_pre` over radius (it is now the 67 ms floor);
FP32 for the `exp`/`expm1` path. NOT sub-cycling -- see below.
## Correlated-k access pattern: BOTH hypotheses refuted (job 10979033)

Two things were predicted to make the 630 ms/100cyc at 88 chains an overestimate. Neither
survived. The test added a `problem/rt_ktab_blk` mode modelling the realistic premixed
table access -- one (T,p) cell and one pair of interpolation fractions per CELL shared by
every chain in the block, with the chain index fastest so the block's k values are
contiguous. It returned bitwise identical numbers to the pessimistic per-chain lookup, so
it was a clean A/B of access pattern alone. **That code was DROPPED at the user's request
once it had answered the question -- it is not in the tree and not in git.** The numbers
below are the whole value of it; do not rebuild it.

RT ms/100cyc, split path, 88 chains (TOTAL, `rt_chain` in brackets):

| RT_NB | no table | table | table, blocked |
|---|---|---|---|
| 4 | 561.2 (484.9) | 625.2 (548.6) | 623.0 (546.9) |
| 8 | **514.8** (441.1) | 626.2 (552.7) | 641.1 (568.0) |
| 16 | 3671.0 (3596.6) | 3817.0 (3744.7) | 3768.5 (3696.1) |

1. **Band-blocking the lookup buys nothing** -- 0.4 % at RT_NB=4, and slightly NEGATIVE at
   RT_NB=8. The index arithmetic and the strided reads were never the cost.
2. **RT_NB=8 is only ~9 % better than 4 without the table and a wash with it.** RT_NB=16
   collapses: it loses parallelism (6 blocks -> 192 workgroups) AND quadruples the private
   footprint to 9216 B/thread. Stay at 4-8.

**But the real finding is better than either hypothesis:** the split kernel has already
absorbed the table cost. Adding the lookup costs **+11 %** here (561 -> 625) against
**+80 %** in the monolithic kernel (1157 -> 2087). The extra waves hide the table's memory
latency, which the occupancy-starved kernel could not. Marginal cost per chain with the
table is **5.4 ms/100cyc in the split kernel vs 23.4 in the monolithic**, a 4.3x
improvement.

**Design consequence:** there is no performance reason to align the chain blocking with
the band structure, so `(nband, ng)` can be chosen on physics grounds alone and the cost
just scales -- about 5.4 ms/100cyc per chain, table included. 11 x 8 = 88 costs ~623.

## RT sub-cycling: the "free 2x" claim was WRONG (checked 2026-08-21)

Measured t_rad = c_p p / (4 g sigma T^4/T) directly from `fixed200.a`
(c_p = 1.43e8, g = 942, dt = 10.14 s):

| | value |
|---|---|
| pressure range | 3.2e-6 to 238 bar |
| T range | 1613 to 12112 K |
| t_rad range | 0.088 to 8.99e4 s |
| **cells with t_rad < dt** | **21.3 %** |
| cells with t_rad < 1000*dt | 51.6 % |

**A fifth of the domain is already radiatively faster than the timestep**, so a blanket
sub-cycle factor is not safe -- it would only be defensible as a depth-dependent scheme,
and the deep half (t_rad > 1000 dt) is where the wasted work actually is.

Also, two separate things were being conflated. RT is called twice per cycle because it is
a SOURCE TERM evaluated in each SSP-RK stage (`user_srcs_func = SourceFunc` ->
`picket_fence_two_stream_RT(pm, bdt)`). Going to one call per cycle is operator
SPLITTING, not sub-cycling, and it trades integration accuracy for the 2x -- not free.

Caveat: t_rad above is the standard Newtonian-cooling estimate and assumes tau ~ 1; in the
optically thin top the true coupling is kappa-rho weighted and slower, which is presumably
why the fully explicit source is stable there at all. There is a commented-out
semi-implicit Newton-Raphson treatment of this source in the apply loop; the user has said
explicitly it is NOT wanted, so do not propose reviving it.

After the split kernel, RT is ~31 % of wall at 88 chains, so halving the calls buys ~15 %
of wall at the cost of accuracy in 21 % of the cells. Not recommended.


## Where things are

Committed as **38311a8a** on `polar-average-perf`, NOT pushed to the fork yet. The commit
also carries the blocked-band harness, which the split path's chain loop is written in
terms of and cannot be separated from. `bench/polar_ab/rt_split_kernel.patch` is the same
diff. Binaries `athena_split` (RT_NB=4), `athena_split1`, `athena_split2`.
Scripts `splitab.sh`, `splitdiag.sh`, `split2.sh`, `split3.sh`, `bandpmc.sh`.
`problem/rt_split` had to be added to `correctness.athinput` -- AthenaK overrides only
modify existing parameters. Build per [[viper-hip-build-recipe]].
Never write in [[never-write-in-run-dir]].
