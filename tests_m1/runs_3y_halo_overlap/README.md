# Implicit-M1 Krylov halo overlapped with the interior operator (m1-int)

Date 2026-09-23, viper. Branch `m1-int` (worktree `/viper/ptmp2/jinma/wt_m1int`). Build and
run tree: `/viper/ptmp2/jinma/m1int_0924/`, which holds:
- `ref/` = rt-integration d80c84ca;
- `new/` = the merge 45e1f1b6;
- `ovl/` = this work;
- `bin/`, `cpu/`, `gpu/`, `runs/`.

Scripts are in `scripts/` here. `bash scripts/finish.sh` prints every table below.
This is TO-DO item 1 of `tests_m1/runs_3w_krylov/README.md` sect. 4.

## 0. The integration branch m1-int

m1-int is made of four branches:
- m1-time2b a55dea7b (m1-space2 + m1-vimp + H-ESDIRK2);
- m1-krylov c1ad9964;
- sc-scale bbf94802;
- rt-integration d80c84ca.

The merge commits are 0d61ff5c, e14a4c2b and 45e1f1b6 (tree eeab5735).
- **Only conflict:** `src/CMakeLists.txt`, where both new files are listed.
- **Bug found by reading the code and fixed in 0d61ff5c:** `implicit_halo_mpi` packed and
  unpacked component `c0` for every n of a multi-component exchange with a fixed first
  component (`c0 + n` is right). The only such exchange is implicit_vimp's
  M1_NVIMP_X = 12. The buffers were also sized for 14 > 12, which happened to be enough;
  they are now sized with M1_NVIMP_X explicitly.

**Merge gate** (`scripts/mklist.sh`, run on the login node with <= 12 workers; the CPU
partitions are closed to this account):

(a) CPU, defaults, bitwise vs rt-integration d80c84ca (hst plus bin data):
- 2-D He slab, t = 100: Eddington and vet_sc full, 1/2/4 ranks.
- 3-D box 84x32x32 in 4 blocks, 30 cycles: Eddington and vet_sc full, 1/2/4 ranks.
- radwave xy, Eddington (`rw2/c2`) and vet_sc full (`rw2/c3`): all 34 tab files
  identical.
  - The run_radwave.py inputs name keys that are new on m1-space2 (`radwave_eig*`,
    `vet_x1_periodic`, `vet_x1_npass`), and ref ignores them. Those keys are stripped
    here.
  - With them present, ref and new differ as intended.
- Deep hot Jupiter, 3 cycles on 4 ranks, run with the production input and the restart
  read in place. hst, the final bin and the final rst are identical:
  - hydro: `bench/cs_hyd4_prod` rst 00113, t = 1.723252e7 s, rotation 56.50,
    cycle 873580;
  - MHD: `bench/cs_mhd_prod4` rst 00081, t = 1.235251e7 s, rotation 40.50,
    cycle 821798.

(b) GPU, apudev job 11952430, defaults, bitwise vs rt-integration. hst, bin and rst
identical for:
- he3d_fast box (vet_sc, vet_mb_agg named), 1 and 2 GPUs, 30 cycles;
- the same box with Eddington, 2 GPUs;
- dhj hydro from rst 00113, 2 GPUs, 20 cycles.

(c) The switches on the merged code. NON-CONVERGED = 0 everywhere.
- plm + vimp + hesdirk2: 2-D on 1/2/4 ranks, 3-D on 4 ranks.
- plm + vimp: 2-D on 4 ranks.
- `implicit_halo_mpi` is **bitwise** vs off with vimp and hesdirk2 on 2-D 2 and 4
  ranks and 3-D 4 ranks. This checks the fix above.
- With `implicit_krylov_pipe` added, the result differs from off at round-off level:
  2-D totE 2.9e-13, KE1 5.1e-10.
- `vet_mb_agroup = 2` differs at round-off level (KE1 2.9e-10), as do halo_mpi + pipe
  with vet_sc on the 3-D box.

## 1. Which option: overlap or a 2-layer halo

Where the multi-rank overhead goes was measured from the rocprofv3 kernel trace of
runs_3w_krylov (`g2_ph`, 2 GPUs, `scripts/halogap.py`). GPU idle, classified by the
kernels around the gap, over the second 10 cycles:
- **Pack -> on-rank copy -> [host waits for the pack, Isend, Waitall] -> unpack:**
  1441 exchanges, 24 us median idle each, 41 ms in total.
- **After the unpack, the operator runs:** 25.4 us per launch.

The idle is message latency plus the host's round trip. It is as long as one operator
launch, so running the interior operator in it can hide nearly all of it.

The 2-layer halo would need `M^-1` applied on the overlap cells. The preconditioner is
the block-local red-black line GS, so a ghost cell's M^-1 value is the neighbour
block's, and computing it locally would change the preconditioner, i.e. the results.
The exchanges of P1 and P2 are also separated by the (omega, alpha) reductions. So the
2-layer halo either changes the solver or saves nothing. **Chosen: the overlap.**

## 2. What changed (commits a8689140, a98562fe, 55598fae)

**`<rad_m1>/implicit_halo_overlap`:**
- Default false, read only when named, so a run without it dumps the same parameters.
- Needs `implicit_halo_mpi`.
- Acts only when the exchange really goes through implicit_halo_mpi (neighbours on
  other ranks) and the operator is the stencil operator.

**New `ImplicitHaloOp(x, y, red, out)`** replaces every "Krylov halo then ImplicitOpX"
pair:
- in the fused loops, kf = 2 and 3;
- at the pipe's start-up and in P1/P2;
- in ImplicitApplyOp (the true residual).

With the switch off it calls exactly those two functions.

**With the switch on**, one exchange runs as:
- receives, pack, event, on-rank copy (`ImplicitHaloMPIPost`);
- **the operator on the interior cells**: every cell at least w from each face, w = 1
  (2 with implicit_vimp: M1VimpRow reaches +-2);
- then the host waits for the pack, sends, waits for the messages, and queues the
  unpack (`ImplicitHaloMPIFinish`);
- **the operator on the shell.**

Per-cell arithmetic is ImplicitStencilOp's, so y is bitwise. Where the operator reduces
(the fused loop's (rt,v) and (t,s), (t,t), and the pipe's start-up), the two partial
sums are reduced asynchronously into pinned memory and added after a fence. That is
round-off.

**Two corrections on the way:**
- **a98562fe:** splitting ImplicitHaloMPI into Post/Finish dropped the event wait
  before the sends.
  - Without it, the timestep collapsed at cycle 2 on the GPU (runs `s12` of job
    11952751).
  - The CPU missed it because the pack is synchronous there.
- **55598fae:** the first shell kernel ran over the whole block with an early return
  (13.8 us per launch against 22.6 us for the interior).
  - It ate the gain: w2 hmo 60.0 against hm 58.5 ms (`runs/s12_v1`, `runs/w2_v1`).
  - The shell cells are now indexed directly.

## 3. Gates of the switch (binary 55598fae; CPU `scripts/cpu_jobs2.txt`, GPU `scripts/ovl_gate.sh`)

- **Off, bitwise** (hst plus bin data):
  - CPU ovl vs merge: 2-D Eddington 2 ranks with halo_mpi, 4 ranks without it; 3-D
    Eddington 4 ranks with halo_mpi; plm + vimp + hesdirk2 4 ranks with halo_mpi.
    This was the a8689140 binary; the later commits do not touch the off path.
  - GPU (job 11953045): hm vs the merge binary on 2 GPUs, 30 cycles: hst identical.
- **On vs off.** Round-off only, NON-CONVERGED 0, fallbacks 0. Last row, relative:

| case | totE | KE1 | KE2 | F1top |
|---|---|---|---|---|
| 2-D Eddington, 2 ranks | 3.0e-13 | 3.5e-9 | 5.2e-9 | 6.4e-13 |
| 2-D Eddington, 4 ranks | 2.5e-13 | 1.2e-9 | 2.0e-8 | 3.9e-12 |
| 2-D vet_sc, 2 ranks | 1.6e-14 | 5.0e-9 | 2.7e-8 | 4.1e-12 |
| 3-D Eddington, 4 ranks | 5.6e-12 | 4.2e-10 | 1.9e-10 | 1.2e-11 |
| 3-D vet_sc, 4 ranks | 4.1e-12 | 1.8e-10 | 9.2e-11 | 4.4e-12 |
| 2-D plm+vimp+hesdirk2, 4 ranks | 2.2e-13 | 9.6e-10 | 3.1e-9 | 5.0e-13 |
| 2-D plm+vimp, 4 ranks | 8.8e-13 | 3.7e-9 | 4.0e-9 | 2.0e-12 |
| 3-D plm+vimp+hesdirk2, 4 ranks | 6.3e-13 | 9.2e-11 | 3.0e-11 | 3.0e-12 |
| 2-D Eddington + pipe, 4 ranks | 4.0e-13 | 4.5e-9 | 2.4e-8 | 7.0e-12 |
| 3-D plm+vimp+hesdirk2 + pipe, 4 ranks | 8.2e-13 | 9.8e-11 | 1.6e-11 | 1.3e-12 |
| GPU he3d vet_sc, 2 GPUs (hmo vs hm) | 2.1e-13 | 8.8e-10 | 7.9e-11 | 3.7e-12 |
| GPU he3d, pho vs ph | 6.9e-12 | 2.8e-10 | 5.5e-11 | 2.2e-12 |
| GPU 3-D plm+vimp+hesdirk2 (bhmo vs bhm) | 7.2e-13 | 2.8e-10 | 1.6e-11 | 2.1e-12 |

  - For scale, the runs_3w_krylov control (off with rad_flux_inner changed by 1 ulp) gave
    KE1 2.7e-9, totE 2.7e-12.
  - The max over all hst columns reaches 5e-3 on the 3-D box. That is on the net
    transverse momenta, which are round-off noise around zero (runs_3w_krylov
    sect. 2).
- **Restart bitwise, switch on** (2-D plm + vimp + hesdirk2, halo_mpi, 2 ranks):
  - Method: run to t = 40 with a restart at t = 20, then restart from it (`o_rsfull_2`,
    `o_rsrst_2`).
  - Final rst identical; all 20 hst rows after the restart identical.

## 4. GPU cost (apudev: s12 job 11953047, w2 job 11953048, prof job 11953050)

**Setup:**
- All jobs export `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.
- One binary, `bin/athena_ovl_boxgpu`, md5 08671642.
- Box and arms as in runs_3w_krylov (84x104x104 in 8 blocks, vet_sc full), 100 cycles.
- Repeat a runs in order, repeat b in reverse. Every run exited with rc 0.

**Arms:**
- hm = implicit_halo_mpi;
- ph = hm + pipe;
- `o` = + implicit_halo_overlap;
- hyd = the same box without `<rad_m1>`.

**Metric: ms/cycle over cycles 10-90.**
- A 9-17 s stall between cycles 90 and 100 hit random runs of every arm, including
  g1_off and hydro (for example g2_hm_b 90:13.0 s, 100:22.0 s).
- It is unrelated to this work (output or file system), so the window stops at 90.
- The runs_3w_krylov numbers used 10-100 and are not comparable. "Before" is the
  hm / ph arms of the same job.

rad = cycle minus hydro, from the a/b means.

| GPUs | hm (before) | **hmo** | ph (before) | **pho** | hydro | rad hm | **rad hmo** | rad ph | **rad pho** |
|---|---|---|---|---|---|---|---|---|---|
| 1 (off / pipe) | 50.1, 50.0 | - | 50.0, 50.3 | - | 21.1, 21.2 | 28.9 | - | 29.0 | - |
| 2 strong | 40.7, 40.8 | **39.3, 39.1** | 39.5, 39.7 | **39.2, 39.1** | 15.5, 15.6 | 25.2 | **23.7** | 24.1 | **23.6** |
| 2 weak (0.9 M/GPU) | 58.3, 58.6 | **58.3, 58.3** | 58.8, 58.7 | **58.8, 58.8** | 21.9, 22.0 | 36.5 | **36.4** | 36.8 | **36.9** |
| 4, 8 strong and weak | pending: apu jobs 11953052 (2 nodes), 11953054 (4 nodes) | | | | | | | | |

- Repeat spread: 0.3 ms or less in every arm.
- Picard passes per step: identical across arms.
- Inner iterations per solve: 15.26-15.36 (strong), 15.40-15.47 (weak).

**Host waits and GPU idle per cycle** (`runs/prof`, 2 GPUs, rank 0, (nlim20 - nlim10)/10;
`scripts/profsum.py`, `scripts/halogap.py`):

| arm | launches | hipStreamSync | hipEventSync | hipDeviceSync | GPU idle "copy -> unpack" (ms / 10 cycles) | total GPU idle (ms / 10 cycles) |
|---|---|---|---|---|---|---|
| g2_hm | 1073 | 265.6 | 191.0 | 83.0 | 41.3 | 257 |
| g2_hmo | 1187 | 167.0 | 192.2 | 83.0 | 4.9 | 192 |
| g2_ph | 1152 | 70.6 | 295.0 | 83.0 | 40.0 | 199 |
| g2_pho | 1253 | 68.4 | 293.0 | 83.0 | 5.4 | 179 |

- Each exchange still makes two host waits: the event after the pack, then MPI_Waitall.
  Neither is a queue drain, and hipEventSync is unchanged.
- The GPU now computes the interior during those waits.
- The "copy -> unpack" idle falls by 88 %, about 3.6 ms per profiled cycle.
- hipStreamSync drops by 99 per cycle in hmo. The fused loop's operator reductions used
  to block; now they are asynchronous and followed by one fence.
- There are ~110 more launches per cycle (2 operator kernels instead of 1).

**What this says:**
- **Within one node (2 GPUs):**
  - The overlap removes the halo idle on the GPU. The wall gain is 1.5 ms (hm -> hmo,
    -6 % of the radiation part) and 0.5 ms with the pipe (ph -> pho).
  - At the weak size, 0.9 M cells per GPU, it gains nothing.
  - The GPU sits idle between the host's launches anyway: 30 % idle with the halo idle
    gone, 179 of 624 ms in pho. The host's launch-and-wait chain sets the cycle, not
    the message latency.
- **Across nodes:** the message latency is larger there, so this is where the switch
  should matter. That is the pending f4/f8 data; `bash scripts/finish.sh` fills the
  table.
- **Radiation vs hydro:** best rad part is 23.6 ms (pho) against hydro 15.5 at 2 GPUs
  strong (1.5x), and 36.4 against 22.0 at 2 GPUs weak (1.65x). The radiation part is
  still not below hydro.

## 5. Next

1. **Launch count and host round trips (runs_3w_krylov TO-DO 3).**
   - With the halo idle hidden, what is left is host-bound: ~1200 launches per cycle
     and ~440 synchronising calls. 83 of those are hipDeviceSync outside the Krylov
     loop.
   - On one rank, a HIP graph of K1-P1-K2-P2 is the lever.
   - Across ranks, the pack event wait could be replaced by GPU-stream-triggered MPI or
     a persistent pack/unpack.
2. **Iteration count (runs_3w_krylov TO-DO 2).** Every inner iteration removed saves two
   exchanges plus two operators and two preconditioners.
3. **Production.** implicit_halo_overlap is exact to round-off and restarts bitwise. It is
   worth turning on together with implicit_halo_mpi once the multi-node numbers
   confirm a gain.
