---
name: rt-ck-sweep-split-measured
description: 2026-09-12 rt_chain_ck SWEEP SPLIT (down- and up-sweep as separate threads, 2x waves, no private I_down column) MEASURED on the cs production restart: kernel 4.94 -> 4.39 ms/call (1.13x), rt_apply 0.65 -> 0.83, wall 25.0 -> 24.3 s (2.8%). NOT merged; patch kept in bench/prof/rt_new/. Kernel is latency-bound on the radial recurrence; 2x waves gave 1.13x not the 1.45x the wave^0.55 extrapolation predicted
metadata:
  type: project
---

Method: bench/prof/rt_new (apudev 11635026, 300-cycle restart of cs_mhd_prod3 dhj.00076.rst, 2 ranks,
plain timing + rocprofv3 --stats). Flag `problem/rt_ck_sweep_split` (default off) in patch
bench/prof/rt_new/0001-*.patch (commit db2b6892 of a deleted worktree; base 2e3ae96e). CPU: ck test
passes, 20-cycle hst identical to all printed digits, dumps 1 float32 ulp (reassociation).
Characterisation (agent, from the profile): one thread per (mb, chain block of 4, k, j); 816 serial
inner bodies per thread; VALUBusy 4.85%, MemUnitStalled 1.65% => latency on the recurrence, relieved
only by resident waves (load scan: 4.6x throughput from 176 to 2816 waves). Production = 1056 waves.
**Why not merged:** user's rule "leave 5 unless substantial"; 2.8% wall is not. Composable leftover:
RT_FP32 (already in tree, 1.42x on this kernel earlier) + this split might reach ~2x on the kernel,
still only ~10% of wall. The RT is 37% of GPU time; the remaining lever is fewer chains (ck_nquad,
bands) i.e. physics accuracy trade, not code. See [[rt-chain-parallel-split]], [[rt-kernel-occupancy-limit]].
