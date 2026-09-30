---
name: meshblock-decomposition-gpu
description: "Meshblock decomposition on one MI300A. Its multi-GPU conclusion is SIZE-SPECIFIC and WRONG at production size -- see the 2026-08-26 correction at the end"
metadata: 
  node_type: memory
  type: project
  originSessionId: be8c40b1-14bc-4240-8f72-1fdd200e6113
  modified: 2026-08-19T15:12:48.119Z
---

**On one MI300A, MORE and SMALLER MeshBlocks are faster, up to 32.** Measured 2026-08-19, jobs
10955177-10955211, mesh 64x64x128, 3000 cycles each, identical work (all reached the same
simulated time), `<meshblock>/nx1` pinned to 64 by the RT column sweep:

| blocks | meshblock | 1 GPU | 2 GPU |
|---|---|---|---|
| 2 | 64x64x64 | 63.29 s | -- |
| 8 | 64x32x32 | 49.92 s | 45.52 s |
| 16 | 64x16x32 | 47.42 s | 45.39 s |
| **32** | **64x16x16** | **43.53 s** | 42.93 s |
| 64 | 64x8x16 | 45.75 s | 43.80 s |

The curve turns over after 32. **A second APU is not worth requesting** at THIS size: at 32 blocks it gains 1.4%.
**THAT CONCLUSION DOES NOT HOLD AT PRODUCTION SIZE -- see the correction at the end.** The apparent 2-GPU wins at 4 and 8 blocks were the finer split, not the hardware.

**This contradicts `docs/general_eos_gpu_porting.md`**, which advises "a few large blocks" for
accelerators -- that advice is about load balance across ranks, not single-GPU throughput.
Finer blocks cost MORE arithmetic (32 blocks carries 66% ghost overhead against 20% at 2) and
win anyway.

**Why, from rocprofv3:** the x2/x3 flux kernels and the polar average run 2.4-3.6x faster with
finer blocks, while the x1 flux kernel -- the one walking the contiguous `i` direction -- gets
20% SLOWER. The effect tracks the k-stride `ncells2*ncells1` (36.1 KB at 64x64x64 vs 10.6 KB at
64x16x16), i.e. it is transverse-stencil memory locality. Ghost overhead is real but smaller.

After the polar-average fix (commit 006bae07) the gap narrows to 1.35x (52.37 vs 38.88 s) but
32 blocks still wins. Cumulative: 63.29 -> 38.88 s, 1.63x.

**Profiling recipe:** `rocprofv3 --kernel-trace --stats --output-format csv -d out -o run -- athena ...`
works well at ~100 cycles. Counter collection (`--pmc`) SERIALISES kernels -- use nlim=3 or it
blows the 15-minute apudev limit. Kernel names in the CSV are Kokkos template soup; the enclosing
AthenaK function is recoverable with `re.search(r'ParallelFor<(?:par_for(?:_outer)?<)?(.*?)\{lambda', name)`.

Beware: verifying any of this against output files needs [[dhj-run-to-run-nondeterminism]] in mind.


---

## CORRECTION 2026-08-26: multi-GPU scaling is GOOD at production size

The 1.4% above was measured at **nx1 = 64** (64x64x128 = 524k cells). Re-measured on the
production grid **nx1 = 234** (234x64x128 = 1.92M cells, 3.7x larger), same 32 blocks,
600 cycles, jobs 11065811 / 11065815 / 11065830:

| GPUs | nodes | compute time | speedup | efficiency |
|---|---|---|---|---|
| 1 | 1 | 33.38 s | 1.00x | 100% |
| **2** | **1** | **18.39 s** | **1.81x** | **91%** |
| 4 | 2 | 12.08 s | 2.76x | 69% |

**The old problem was simply too small to saturate one MI300A**, so splitting it bought
nothing. Do not quote the 1.4% for production-sized runs.

**1, 2 and 4 ranks are all BITWISE IDENTICAL to the 1-rank run** (9,584,640 values, zero
differ, both dumps) -- worth knowing because this campaign has only ever run single-rank and
[[polar-mpi-host-mirror-bug]] killed every multi-rank polar run before 3882e37f.

**2 GPUs on one node is the sweet spot**: 91% efficient, and the second APU is already
allocated to the node (`Gres=gpu:2`) while every submit script in this campaign asks for
`--gres=gpu:1`. It is also the cheapest in node-hours. 4 GPUs buys wall-clock (59 h vs 90 h)
at 31% MORE node-hours, and needs the `apu` partition because **apu1 is MaxNodes=1**.

Untested: whether 4-rank efficiency improves with 64 blocks (meshblock nx2=8) instead of 8
blocks per rank.
