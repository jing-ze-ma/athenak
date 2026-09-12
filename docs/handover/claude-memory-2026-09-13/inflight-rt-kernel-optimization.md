---
name: inflight-rt-kernel-optimization
description: "2026-08-21: RT kernel 160.1 -> 115.3 ms/100 cyc (1.39x, bitwise identical) by interleaving the four IR combinations for ILP; wall 39.09 -> 37.79 s, committed as 2be5d0d9. Team/LDS rewrite and all three launch-overhead fixes LOST."
metadata: 
  node_type: memory
  type: project
  originSessionId: 5e5b0978-be55-410d-9852-994d6627a6fd
  modified: 2026-08-21T14:25:00.000Z
---

**Task:** make `picket_fence_two_stream_RT` in `src/pgen/deep_hot_jupiter_rt.cpp` faster.
Branch `polar-average-perf`, head `c39b794c`. Working tree `src/athena.hpp` is back to HEAD.

## WHAT WORKED: interleave the four IR combinations (job 10976935)

`src/pgen/deep_hot_jupiter_rt.cpp`, `picket_fence_two_stream_RT`. The four (band,
quadrature) combinations are independent pairs of linear recurrences in radius; the
original ran them one after another, so each wavefront sat on a single dependency chain.
Stepping all four inside one radial loop (`NC = 4`, `I_ir_down_c[NC][NN]`, `I_ir_up_c[NC]`,
coefficient arrays `gamirc/fbc/muggc/wggc`) gives four independent strands to overlap.
Same work, same arithmetic, fluxes still summed in combination order 0,1,2,3 per interface.

**RT 160.1 -> 115.3 ms/100 cyc (1.39x). Wall 39.09 -> 37.79 s (3.3 %). BITWISE IDENTICAL.**
Binary `bench/polar_ab/athena_rtilp`, script `rtilp.sh`. **Committed as `2be5d0d9` and pushed to
`fork/polar-average-perf`.** See [[use-fork-not-origin]]. Follow-on work: [[inflight-rt-band-scaling]].

Why ILP is the only lever left: 8192 columns / 64 = **128 wavefronts for 912 SIMDs**. Six of
every seven SIMDs have no work at all and the mesh has no more columns to give, so adding
waves is impossible and hiding latency inside the existing waves is all that remains.

## WHAT LOST, with numbers (do not retry)

| attempt | RT ms/100 cyc | wall s | note |
|---|---|---|---|
| base `c39b794c` | 160.1 | 39.09 | reference |
| **interleaved combos** | **115.3** | **37.79** | bitwise identical, keep |
| team/LDS, `Kokkos::AUTO` | 440.7 | 48.20 | job 10976794 |
| team/LDS, 64-lane teams | 241.7 | 41.67 | job 10976857, `rt_team64.patch` |

The team-per-column rewrite (column state in LDS, transcendentals lane-parallel, one lane
for the recurrences) loses for a structural reason worth remembering: during the two
recurrences the ORIGINAL runs 64 columns per wavefront with every lane busy, while a team
runs one column per wavefront with 63 lanes parked at a barrier. The parallel phases cannot
pay that back. It also came out bitwise different -- a real bug in that rewrite, since the
interleaving change with the same discipline reproduces the reference exactly.

Also rejected on estimate, not measured: storing the per-interface coefficients to halve the
expm1 count. It removes 272 expm1 per column but adds ~13 kB of scratch traffic per column
(816 stores + 816 loads), which at the ~350 GB/s this kernel sustains costs more than the
transcendentals save.

## The constant-memory hypothesis is DEAD -- do not re-run it

Earlier reasoning (still in the traces) was: 3513 constant-memory launches per 100 cycles carry
408 ms of preceding gap because Kokkos' HIP ConstantMemory path event-synchronizes on the
previous constant-memory kernel. Two jobs tested it.

**Job 10976146** -- `Kokkos::Experimental::require(policy, HintLightWeight)` on all 5 `par_for`
RangePolicy sites *and* all 4 `par_for_outer` TeamPolicy sites (`athena_lw`,
`bench/polar_ab/lightweight.patch`): constant 3513 -> 907, but **53.39 s vs 39.07 s base**
(37 % SLOWER) and bitwise different. Hinting a TeamPolicy changes `Kokkos::AUTO` team sizing,
which is what wrecked the time.

**Job 10976233** -- hint on the 5 RangePolicy sites only (`athena_lwr`,
`bench/polar_ab/lightweight_range_only.patch`, script `lwrab.sh`): constant 3513 -> 1907,
time 39.40 s vs 39.02 s base -- **neutral**. The base binary was re-verified bitwise identical
to `fixed200.a` in the same job, so the reference is good.

**Job 10976411** -- the direct test: patch the kokkos submodule so the ConstantMemory launch
path stops serializing at all. 32 kB constant buffer split into 8 rotating 4 kB slots
(`HIPTraits::ConstantMemorySlots/ConstantMemorySlotSize`, `constantMemReusable` becomes
`std::array<hipEvent_t,8>`, new `constantMemNextSlot`, kernels take a `slot_offset` argument;
functors > 4 kB keep whole-buffer semantics and wait on all 8 events). Saved as
`bench/polar_ab/constmem_slots.patch`, binary `athena_cms`, script `cmslots.sh`. Transport and
codegen are otherwise untouched -- still 3513 constant launches. Result: **42.04 s vs 39.04 s
base, 7.7 % SLOWER.** Removing the sync does not help; letting the host run ahead and queue
deeper apparently costs more than the sync did. The kokkos submodule has been reverted to
clean and `build-gpu-bench` rebuilt from it. Do not re-try this either.

Gap-by-mechanism on the two traces (`mech_lwr.py`, reading `prof.lwr.base/st` and
`prof.lwr.new/st`) is the decisive number -- **total idle did not move**:

| | base | range-only hint |
|---|---|---|
| constant_memory | 3513 launches, 408.10 ms | 1907, 204.99 ms |
| global_memory | 1404, 40.65 ms | 1404, 70.59 ms |
| local_memory | 6418, 122.83 ms | 8024, 294.95 ms |
| other | 912, 20.33 ms | 912, 22.11 ms |
| **total gap** | **591.9 ms** | **592.6 ms** |

The gap followed the kernels to their new transport (local mean 19.1 -> 36.8 us). So the gap
belongs to *which kernel is being launched*, not to *how the functor travels*. The original
correlation was confounded: the big functors are simply the kernels that follow expensive
host-side work.

## Where the 592 ms actually is (`bench/polar_ab/gapsrc.py`, base trace)

Per 100 cycles, span 1532.6 ms, 12248 dispatches. Top contributors by preceding gap:

| kernel | n | busy ms | gap ms | mean us | median us |
|---|---|---|---|---|---|
| `MHD::CalculateFluxes<RSolver 3>` (outer) | 600 | 463.2 | 118.1 | 196.9 | 181.0 |
| `Kokkos::Impl::zero_with_...` | 995 | 2.6 | 67.7 | 68.1 | 71.0 |
| `ParallelReduce<...Combine...>` | 606 | 28.5 | 61.1 | 100.8 | 15.6 |
| `SourceFunc` | 200 | 42.2 | 36.4 | 181.8 | 181.4 |
| `picket_fence_two_stream_RT` | 200 | 160.5 | 36.1 | 180.3 | 179.6 |
| `Resistivity::SetResistivity` | 201 | 7.3 | 34.7 | 172.6 | 156.1 |
| `Coordinates::SrcTermsSphericalPolarMHD` | 200 | 18.0 | 33.2 | 165.9 | 165.6 |
| `HydrostaticEquilibrium` lambda #3 | 201 | 13.1 | 32.2 | 160.2 | 159.3 |

The signature to chase: a **flat ~160-200 us host gap, twice per cycle** (200 dispatches per
100 cycles = once per RK stage), in front of many unrelated kernels -- mean ~= median, so it is
a fixed cost, not a tail. Plus 995 near-empty `zero_with_...` kernels costing 67.7 ms of gap
for 2.6 ms of work. Next measurement should be an **API trace** (`rocprofv3 --hip-trace
--kernel-trace`) so the host calls inside those windows are visible; suspects are the host-side
fences behind `ParallelReduce` (median 15.6 us but mean 100.8 us) and whatever allocates/zeroes
per stage.

## Unexplained, possibly a real bug

`athena_lwr` differs from `fixed200.a` at the *same* byte (468596, line 1982) as `athena_lw`.
A pure `parallel_for` over a RangePolicy cannot depend on block size unless something races or
accumulates out of order, and there are **no float atomics** in the MHD path (`grep Kokkos::atomic`
finds only `outputs/track_prtcl.cpp`, `outputs/coarsened_binary.cpp`, `bvals/bvals_part.cpp`).
So this looks like a latent write-write race of the same family as the one fixed in `c39b794c`.
Job 10976611 ran the control and the answer is **the opposite of a build artifact**: a clean
`build-gpu-bench` binary (`athena_ctl2`, clean tree, clean submodule) is **bitwise identical to
`fixed200.a`**. So the build directory is fine, and `athena_lw`, `athena_lwr` and `athena_cms`
each genuinely change the answer -- all three self-reproducible, all three differing from the
reference at the *same* byte 468596 (line 1982 of `dhj.mhd_w_bcc.00001.bin`) despite sharing no
mechanism beyond "the launch transport changed". A pure transport change must not move a bit.
Prime suspect is a latent write-write race of the same family as `c39b794c`, exposed by a
different launch order. **Not yet done:** decode which variable/cell byte 468596 is, and check
whether it sits on the polar axis or the outer-x1 boundary.

## Concurrent session, 2026-08-21 14:20 -- resolved

Peer session `athenak-5a` began editing the **kokkos submodule** live: an 8-slot rotation of
the constant-memory buffer (`HIPTraits::ConstantMemorySlots = 8`, `constantMemReusable` becomes
`std::array<hipEvent_t,8>`, new `constantMemNextSlot`) in `Kokkos_HIP.cpp`,
`Kokkos_HIP_Instance.{cpp,hpp}`, `Kokkos_HIP_KernelLaunch.hpp`. That is the direct form of the
hypothesis refuted above. It was mid-edit and did not compile when I tried to rebuild.
**Both sessions share `build-gpu-bench/`** -- coordinate before building.

## API trace, job 10976611 (`prof.api/st`, 60 cycles, analysed by `bench/polar_ab/apigap.py`)

Steady state span 1039 ms, 8088 dispatches, **idle 471.5 ms in 6472 windows**; 82.7 % of that
idle is accounted for by host HIP calls overlapping the windows:

| host call | n | in-gap ms | % of idle |
|---|---|---|---|
| hipEventSynchronize | 2602 | 142.3 | 30.2 |
| hipMemcpyAsync | 4529 | 67.3 | 14.3 |
| hipLaunchKernel | 5286 | 42.0 | 8.9 |
| hipMemcpyToSymbolAsync | 1971 | 35.8 | 7.6 |
| hipStreamSynchronize | 3925 | 29.8 | 6.3 |
| hipFuncGetAttributes | 41 | 28.5 | 6.0 |
| hipDeviceSynchronize | 7369 | 19.6 | 4.2 |
| hipFreeAsync / hipMallocAsync | 1172 / 1429 | 7.4 / 7.2 | 3.1 |

Read it with care -- 30 % sitting in `hipEventSynchronize` is the constant-memory wait again
(~55 us of each call falls in a GPU-idle window), and job 10976411 already proved that deleting
that wait makes the run *slower*. Plausible reason: the wait is not pure waste, because the
`hipMemcpyToSymbolAsync` that follows it is stream-ordered work the GPU must execute between
kernels; running the host ahead just interleaves more copies.

New and not yet chased: **1429 `hipMallocAsync` + 1172 `hipFreeAsync` + 7369
`hipDeviceSynchronize` in 60 cycles** -- 24 allocations and 123 device fences *per cycle*.
Per-cycle View allocation inside the time loop is a plain bug-shaped cost; find the call sites
and hoist them. Worth ~34 ms of 471 ms directly, possibly more indirectly.
`hipFuncGetAttributes` (41 calls, 28.5 ms) is one-off per unique kernel -- ignore it.

## Still open, in order

1. DONE -- see the API trace section above. `zero_with_hip_kernel` is explained too: on this
   APU build (`KOKKOS_IMPL_HIP_UNIFIED_MEMORY`) `ZeroMemset<HIP>` replaces `hipMemsetAsync`
   with a kernel, so every `deep_copy(view, 0)` and every reduction scratch-flag clear is a
   dispatch. Next: hunt the per-cycle `hipMallocAsync`/`hipDeviceSynchronize` call sites.
2. Re-test the parked RT five-way split (`git stash` entry `rtsplit2 WIP`, duplicated as
   `bench/polar_ab/rtsplit2.patch`) only if a launch-overhead fix ever does land; on its own it
   was 152.1 -> 138.7 ms RT busy but 38.87 -> 40.61 s wall, because 8 extra launches per cycle
   cost ~88 us each.
3. The 472 us intensity sweep: split the four (quadrature, band) combos into their own threads
   -> 32768, 4x occupancy. Bitwise safe ONLY with per-combo partial `F` arrays summed in order
   c = 0,1,2,3; atomics would reorder and break identity.
4. `double_gray_two_stream_RT` and `double_gray_two_stream_RT_source` are dead copies, no call
   sites, still `NN = 270`. Left alone.

## Where things are

`bench/polar_ab/`: binaries `athena_lw`, `athena_lwr`, `athena_rtopt2`, `athena_rtsplit2`
(`athena_ctl` is a stale copy, ignore it). Scripts `lwab.sh`, `lwrab.sh`, `rtsplit2.sh`,
`rttrace.sh`; `rtcounters.sh` still FAILS ("Request exceeds the capabilities of the hardware").
Analysis: `gaps2.py`, `mech.py`, `mech_lwr.py`, `gapsrc.py` -- run from `bench/polar_ab/`.
Traces: `prof.lwr.base/st`, `prof.lwr.new/st`, `trace3/out`, `prof.split2/st`.
Bitwise reference `fixed200.a/bin/dhj.mhd_w_bcc.00001.bin`; base binary is
`bench/xe_long/b3_e13/athena`. Build per [[viper-hip-build-recipe]]. Never write in
[[never-write-in-run-dir]]. See [[dhj-run-to-run-nondeterminism]].
