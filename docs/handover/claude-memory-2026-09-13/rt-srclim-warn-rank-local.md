---
name: rt-srclim-warn-rank-local
description: BUG, unfixed - RTSourceLimiterWarn prints rank 0's local clip count as if it were the global one, and stays SILENT entirely when rank 0 itself never clips
metadata:
  type: project
---

`src/pgen/deep_hot_jupiter_rt.cpp:317` `RTSourceLimiterWarn(const int nclip)`:

```cpp
if (nclip <= 0 || rt_srclim_warned) return;
rt_srclim_warned = true;
if (global_variable::my_rank == 0) { ... << nclip << " cell(s) ... }
```

`nclip` comes straight out of `par_reduce_clip3/4`, so it is **rank-local and never
MPI-reduced**. The pgen has no MPI includes at all.

**Measured 2026-08-24.** The same bitwise-identical state reported **112 cells** at 16
ranks (`meshblock/nx2=nx3=8`, 128 blocks) and **10 cells** at 2 ranks (`nx2=nx3=64`, 2
blocks). Both are rank 0's share; neither is the domain total. `ea823e3a`'s "it clipped
10 cells -- the exact ten the dump analysis found" is therefore rank 0's ten.

**The worse half:** the early return fires on rank 0's OWN count, so when the offending
column does not live on rank 0 the warning **never prints**, however hard other ranks are
clipping. A run can clip every cycle and look clean. This is why "did the limiter fire on
viper?" cannot be answered from a log.

**Fix:** `MPI_Allreduce(MPI_IN_PLACE, &nclip, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD)` under
the tree's `MPI_PARALLEL_ENABLED` guard before the rank-0 print, so both the count and the
decision to warn are global. Needs a test. NOT YET DONE -- offered to the user, not
authorised. See [[dhj-general-eos-ck-blowup]].
