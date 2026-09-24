# ck-scratch: the ck chain kernel's column arrays in global memory (nx1 = 256 on MI300A)

Branch `ck-scratch` (from rt-integration e89954e2), worktree /viper/ptmp2/jinma/wt_ckscratch.
Work dir: /viper/ptmp2/jinma/ckscratch_0924 (binaries, scripts, logs, scratch tables).

## 1. Problem

Every ck run with n1 > 136 aborted in the first RT call with HSA_STATUS_ERROR_OUT_OF_RESOURCES
(/viper/ptmp2/jinma/cs_weak_0924/README.md). The cause is `rt_chain_ck` in
src/utils/two_stream_rt.hpp. It keeps up to ten whole-column arrays per thread, each
`[RT_NB][NN]` with NN the compile-time column tier (72/136/264/520):
I_down, Cmx (SPH), Krs (BSP), Kpc/Cc0/Cci/Cco (CCH >= 2), and Js0..2 (JAC).
All of them sit in the private (scratch) segment, which the runtime sizes for full occupancy.

## 2. Scratch per lane (bytes)

These are `.private_segment_fixed_size` values from the device code objects
(`tools/scratch.py <binary>`, `tools/chain.py` for the chain template tags).
Files: scratch_base.txt (e89954e2) and scratch_new3.txt (98544c67).

| rt_chain_ck instantiation (SPH=BSP=1, CCH=2, tm) | tier 72 | 136 | 264 | 520 | new (one instantiation, any n1) |
|---|---|---|---|---|---|
| production semi (FOP=0, JAC=0) | 17040 | 31376 | 60016 | 117376 | 976 |
| T4 / c2 (FOP=1, JAC=0) | 17072 | 31408 | 60096 | 117440 | 964 |
| Jacobian pass (FOP=1, JAC=1) | 24784 | 45264 | 86272 | (not built) | 1888 |
| all 34 chain variants, max | 24784 | 45264 | 86272 | 117616 | 1984 |

The other RT kernels were not changed. Their largest scratch values at 264 / 520 are:
grey monolithic reduce 17760 / 34144, grey chain 8832 / 17024, ck_lin 8464 / 16656.
These fit: the 136 chain at 45 kB runs in production.

## 3. The change (two_stream_rt.hpp only)

- The ten arrays are now slices of one device buffer, `rt_ckscr_ptr`. It is grown on demand
  (`CkScrEnsure`) and sized by the actual n1, thread count and number of groups.
  `CkScrCol[cc][i]` reads like the old array, so the kernel body is textually unchanged.
- The layout copies how the hardware lays out scratch. Threads form tiles of 64 (one
  wavefront, numbered by the par_for flattening of m, blk, k, j), and each tile owns one
  contiguous region. Element (cc, i) of lane l in group g is at
  `tile + ((g*RT_NB + cc)*n1 + i)*64 + l`.
- A switched-off array (its tag false) stays the one-element private array it was.
- Because the chain no longer uses NN, `launch_ck_tier` instantiates it once (tag 0) instead
  of four times. The JAC chain therefore now exists at every n1; before, it stopped at 264.
- Buffer size per GPU: for semi / T4 (7 groups) it is 4·n1·64·8 B·7 per tile.
  - Production (nx1 = 128, 12 blocks/GPU): 2.0 GB.
  - nx1 = 256 at 48 blocks/GPU: 15.6 GB.
  - A launch that would exceed `rt_ckscr_gb` (16 GB, a code constant) is split over chain
    blocks. This split never happened in the runs below.
- Registers are unchanged: 128 VGPRs, and the NC-deep register arrays stay.
- Still limited: ck_lin / ck_lin1 / grey keep their tiers, so the `n1 <= 520` guard stays.
  n1 up to 520 (nx1 up to 512) is fine. Beyond that, ck_lin's `Sc[NC][NN]` needs the same
  treatment.

A first layout, (group, cc, i, thread) with the thread index fastest (commit c72b3a2a /
1dadefa1, binary athena.gpu.new2, job 11956961), was also coalesced. But it put successive i
of one wave about 0.5 MB apart, and it was slower: semi 63.0 -> 82.6 ms/cycle (+31 %) and
c2 47.2 -> 52.1 (+10 %). The per-wave tiling (98544c67) fixed that.

## 4. Gates

CPU (serial, login node; `cpu/gates.sh new3`; athena.cpu.base vs athena.cpu.new3):

| test | semi | t4 | c2 |
|---|---|---|---|
| well-posed A2 (100 calls): rec.txt | bitwise | bitwise | bitwise |
| well-posed A2: hst / bin / rst data | bitwise | bitwise | bitwise |
| full hydro, ck RT per stage, 6 cycles from the wp IC: hst / bin / rst data | bitwise | bitwise | bitwise |

GPU (apudev, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, 2 ranks).
- Restart: bench/cs_hyd4_prod/rst/dhj.00160.rst, read in place: t = 2.43883e7 s = 79.96
  rotations, cycle 1256387.
- Input: /viper/ptmp2/jinma/ckfast_0923/prod_fast.athinput.

- c2, about 810 cycles (tlim = T0 + 3000 s): final rst, hst and bin are **bitwise**, base vs
  new (job 11957326).
- semi: the results agree to round-off but are not bitwise.
  - Both binaries are deterministic run to run: base vs base and new vs new are bitwise at
    1, 2, 3, 5 and 10 cycles (jobs 11957377, 11957453).
  - After 1 cycle, 61 370 of 11.4 M restart doubles differ: median relative difference
    9e-15, 280 above 1e-12. The maximum, 3.2e-11, is on a near-zero value (-4.36e-10,
    absolute difference 1.4e-20). The float32 bin output is still bitwise.
  - After that the differences grow chaotically in the thin top cells: after 10 cycles, bin
    velx differs in 1515 of 786 k cells, with max |d|/max|velx| = 5e-3 at i = 127.
  - Likely cause: code generation. With a private array, LLVM forwards a stored value to its
    later load and can then fuse the multiply into an FMA. With a global pointer that may
    alias Fb_g and the other Views, it has to reload, so a few operations round differently.
    The CPU build, which has no such fusion change, is bitwise.
- nx1 = 256: cs dhj hydro from the IC (prod_fast input with mesh/meshblock nx1 = 256 and
  mesh nx2 = nx3 = 64), 1 node, 48 blocks/GPU, 30 cycles (job 11957326).
  - The new binary completes all 30 cycles in every arm: semi, T4, c2.
  - T4 and c2 print ck_implicit lines with no NOT-CONVERGED. The first T4 call takes
    passes = 7 and reaches res 1e-8.
  - The base binary aborts in the first RT call with HSA_STATUS_ERROR_OUT_OF_RESOURCES
    (rc 134).

## 5. Timings (GPU, ms/cycle)

nx1 = 128 production state, job 11957326. Interleaved forward / reverse / forward, same
restart, about 810 cycles each, ms/cycle measured from the ndiag = 10 elapsed lines:

| arm | base r1 / r2 / r3 | new r1 / r2 / r3 | median base -> new |
|---|---|---|---|
| semi (production) | 61.34 / 63.03 / 63.20 | 50.45 / 50.68 / 50.87 | 63.03 -> 50.68 (-19.6 %) |
| c2 | 47.03 / 46.46 / 46.92 | 43.62 / 43.79 / 44.25 | 46.92 -> 43.79 (-6.7 %) |

nx1 = 256, 48 blocks/GPU, new binary, cycles 5-30, job 11957326:

| arm | ms/cycle |
|---|---|
| semi | 429.1 |
| c2 | 289.8 |
| T4 | 826.3 |

## 6. Files

- Binaries (build dirs deleted), each with a `.commit` file:
  - athena.{cpu,gpu}.base (e89954e2)
  - athena.{cpu,gpu}.new3 (98544c67)
  - athena.{cpu,gpu}.new2 (the first layout)
- Scripts:
  - build.sh
  - cpu/gates.sh
  - gpu/job.sub (128 A/B + 256)
  - gpu/p1c.sub (1-10 cycle agreement)
  - gpu/ana.py
  - gpu/cmpbin.py
  - tools/scratch.py
  - tools/summ.py
  - tools/chain.py

## HANDOVER

- ck-scratch is ready to merge into rt-integration.
  - It touches only src/utils/two_stream_rt.hpp: the kernel storage, the chain tier
    dispatch, and the JAC 264 limit.
  - It is faster on production (semi -20 %, c2 -7 %).
- After merging, any GPU production binary can run nx1 = 256. The weak-scaling test
  (cs_weak_0924) can then be rerun; hydro, or MHD after the split-hook fix.
- Open items:
  - `rt_ckscr_gb` is a code constant (16 GB), not an input key. It was kept that way to
    stay out of the pgen that other branches edit.
  - n1 > 520 still needs ck_lin's `Sc` moved the same way.
  - semi on GPU is round-off, not bitwise, vs base (see section 4). Any restart-continuity
    check that assumes bitwise GPU semi across this commit will fail.
