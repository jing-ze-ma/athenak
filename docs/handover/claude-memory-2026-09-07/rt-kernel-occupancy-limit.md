---
name: rt-kernel-occupancy-limit
description: "The dhj picket-fence RT kernel uses ~5% of the MI300A: VALUBusy 4-6%, occupancy ~0.5 waves/CU, because it is parallel only over (m,k,j) = 8192 columns. Chains/bands are a SERIAL loop inside each thread. Parallelizing over frequency band is the big unexploited lever."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T17:00:00.000Z
---

Measured on Viper MI300A with rocprofv3 counters, `bench/polar_ab/bandpmc.sh`
(job 10978180, `athena_band4`, correctness.athinput, meshblock 64x16x16, 32 blocks):

| | nchain=4 | nchain=88 |
|---|---|---|
| VALUBusy | 4.0 % | 6.0 % |
| MeanOccupancyPerCU | 0.44 | 0.53 |
| VALUUtilization | 98.8 % | 99.5 % |
| MemUnitStalled | 0.27 % | 0.10 % |

**The GPU is ~94 % idle even at 88 chains.** The kernel
(`picket_fence_two_stream_RT`, the `par_for("2stream_rt", ..., 0,nmb1, ks,ke, js,je)`
at ~line 2944) is parallel over (meshblock, k, j) ONLY = 32*16*16 = 8192 threads =
128 wavefronts on 1216 SIMDs -> 0.42 waves/CU, which is exactly the measured occupancy.
Radius is a serial recurrence and **frequency chains are a serial loop inside each
thread**, blocked RT_NB at a time.

**Consequence:** the linear 12.4 ms/chain in [[inflight-rt-band-scaling]] is NOT a
throughput limit. It is `exp()` latency on a dependency chain with no other wave to hide
it. VALUUtilization 99 % just means the lanes within the (few) active waves are not
diverged; it says nothing about the machine being busy.

**The lever:** make the chain index a parallel dimension. 8192 x 88 = 720k threads =
11.3k waves = ~9 waves/SIMD, i.e. full occupancy. 88 chains could then cost close to what
4 chains cost today. Requires (a) hoisting the chain-independent precompute (tau_down_r_f,
B, Q_v, the EOS/opacity lookups) into its own fully parallel kernel over (m,k,j,i) writing
tau and B to global arrays, (b) reducing F_ir over chains -- atomicAdd per cell, or a
per-chain buffer plus a reduce kernel.

**Prerequisite:** `constexpr int NN = 270` (three sites in the pgen, also in HEAD) sizes
the private per-column arrays. At 88 chains the private segment is 18212 B/thread. That
does not cost much today because occupancy is starved anyway, but once the chain axis is
parallel and 11k waves are resident, 18 KB/thread of scratch WILL bind. Shrink NN to the
real n1 (68 for this input) -- check the largest nx1 in production inputs first -- or move
the arrays to team scratch via `par_for_outer`.

Note the earlier `par_for_outer` + `ScrArray1D` version is still sitting commented out
directly above the live `par_for`. See [[inflight-rt-kernel-optimization]] for what was
tried and lost.
