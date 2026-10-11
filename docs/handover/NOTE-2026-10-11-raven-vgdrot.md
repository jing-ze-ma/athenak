# vgdrot-1011 (Raven): AG Car A vet_gd rotation spike and per-shell sweep loop

Work dir `/raven/ptmp/jinma/vgdrot_1011`. Branches (local only, not pushed):
- `vgdrot-1011` (worktree `/raven/u/jinma/ATHENAK/wt-vgdrot`), on rt-integration d28dacd8:
  - b1855fd8: batched exact-list builds, table prefetch, per-call allocations removed, opt-in vet_gd_halo_ipc.
  - a749711a: CUB scans in the batch, kernel-argument launches of the per-shell sweep kernels, IPC modes.
- `vgdrot2-1011` (worktree `/raven/u/jinma/ATHENAK/wt-vgdrot2`), on deepatm-1011 5cb2fec2:
  - 3e790a6c: merge of vgdrot-1011 into deepatm-1011 (this brings in vgdfuse-1010 / halo_exact as well, which deepatm-1011 did not have).

Every change is bitwise (hst every cycle + bin dumps, see gates). Defaults: all ON except vet_gd_halo_ipc (opt-in).

## What changed

### Phase 1: the rotation cycle
- Problem: the rotation cycle took 0.72 s against 0.27 s for a normal cycle (base). It had two causes.
  1. The exact-list build: 464 collective per-shell builds, each with fences, 2 handshakes and an Allreduce, run inside the sweep.
  2. VetGdTables: 0.11 s per rotation, mostly memset and pageable copies of about 200 MB of host mirrors.
- **VetGdHxBatch** (rad_m1_vetgd.cpp): at the start of VetGdSweep, every (pass, shell) whose exact-list record does not match is built in groups of G = 15 shells (G is set by the VGD_HXB_MB = 1024 MB pool per rank). Per group:
  - All marks, sources, flags, scans, bounds and offsets are enqueued on the device with no host sync, then ONE fence.
  - ONE count message per partner (it also carries the receive offsets), one wait for the offsets, ONE Allreduce vote.
  - The lists are made with the same expressions as VetGdHxBuild2.
  - Scratch lives in a pool allocated once.
  - List growth has 1/4 headroom: reallocs happen only at the first build.
  - Scans use cub::DeviceScan, which is exact integer arithmetic.
  - With this, `hx_mark` in the loop finds every record valid.
  - Env switches: `VGD_HXB=0` builds per shell as before; `VGD_HXB_CUB=0` uses the old VgdScan2.
- **VetGdTables**:
  - The host part goes into persistent pinned buffers.
  - The NEXT rotation's tables are computed on a host thread (std::async, same function, 8 threads) right after the current ones. They are uploaded into a second device set on their own stream once ready, then swapped in at the rotation, so the rotation spends 0 s on tables.
  - `VGD_TAB_PREFETCH=0` turns this off.
- **Quick wins**:
  - Persistent nan_report dummies in hydro_fluxes and conduction (were 4 + 2 device alloc/free per call).
  - fofc_report `ioff` is copied only when it changes.
  - Pinned persistent buffers in VetGdPost, VetGdTwin and VetGdHalfRange.
  - Under twin_det, the twin's two shell means and the 6 (deepatm: nhc) half-range shell means each use ONE gather instead of 2 or 6 rounds. These are the same rank-ordered sums.
- **Timing script**: `time.sh` + `ana.py`.
  - Wall time per cycle comes from the `elapsed` lines; cycles K0 = 2 .. NL-2 are used, so the LAST step (which carries the end-of-run output) is excluded.
  - It reports steady (c % 10 != 0), rotation (c % 10 == 0) and the average separately.
  - It writes hst every cycle and bin dumps at t = 0 and at the end, with rst off.
  - `bitcmp.sh` compares bin data after `<par_end>` and the hst rows.

### Phase 2: the host in the per-shell loop
- nsys (job 31064224, base path) shows the sweep window at 121 ms with only 64 ms of GPU work.
  - The host waits at the fence before each MPI post: cudaStreamSynchronize took 57 ms.
  - The shell, pack and scatter functors (0.5-4 kB) were launched through Kokkos' constant-memory path. That path does a cudaEventSynchronize on the previous such launch, so the host could never run ahead.
- **Kernel-argument launches** (`VgdLight` = `HintLightWeight`, as VetFor in rad_m1_vet.cpp) for the shell, pack, scatter and zero kernels. These are the same kernels and give the same bits. This is the main Phase 2 gain: -7 % steady.
- **vet_gd_halo_ipc** (opt-in: `rad_m1/vet_gd_halo_ipc = true` or env `VGD_IPC=1`).
  - What it does: a GPU-driven node-local CUDA IPC exchange of the batch-made exact lists, on the krylov implicit_halo_ipc protocol.
    - One ack/wait kernel, then the pack stores straight into the partners' receive buffers.
    - One flag/wait kernel, then the scatter reads the IPC buffer.
    - The host never waits in the sweep.
  - Results: it is bitwise. The profiled sweep is 88.5 ms against 106.9 ms on the MPI path, but the GPU time is now spin-waiting on rank imbalance (lockstep per shell).
  - Unprofiled, it gains only 0.7 % (0.2521 against 0.2538 s/cycle), so it stays opt-in.
  - Variant `VGD_IPC=2` (local pack + peer cudaMemcpyAsync) was slower: 0.294 against 0.276 before the launch fix.
  - The early version with a `__threadfence_system` after every remote store was slower too.
- Not done: splitting the shell into interior and boundary parts to overlap it with the halo. vgdfuse's `vet_gd_overlap` was already slower, and with IPC the GPU is no longer starved by the host.

## Gates (Raven gpudev, 1 node x 4 A100, 1 MeshBlock 480x64x64 per GPU)
- BASE = `bin/athena_he_a100_d28dacd8_base` (rt-integration d28dacd8, md5 1b8e4345), built with the same cmake as NEW.
- AG Car A input = `inputs/agc.athinput` (xthinfix agcar_A_xtf, production `all`, be, rotate_every 10). Runs are 32 cycles with rotations at cycles 10, 20 and 30; arms are interleaved.

| job | case | arm | steady s/cyc | rotation s/cyc | avg s/cyc | bitwise vs base |
|---|---|---|---|---|---|---|
| 31066702 | AG Car A | base (2 reps) | 0.2735 | 0.718 | 0.3195 | - |
| | | new5 (phase 1) | 0.2674 | 0.360 | 0.2770 | identical, 4 bins + 68 hst rows, both reps |
| | BSG (4 blocks) | base | 0.1714 | 0.516 | 0.2070 | - |
| | | new5 | 0.1656 | 0.224 | 0.1716 | identical, 8 bins + hst, both reps |
| 31067202 | AG Car A | base (2 reps) | 0.2734 | 0.717 | 0.3192 | - |
| | | **new6 (a749711a, default)** | **0.2538** | **0.346** | **0.2634** | identical, both reps; rep2 == rep1 |
| | | new6 + vet_gd_halo_ipc | 0.2521 | 0.344 | 0.2616 | identical |
| 31067404 | BSG (4 blocks) | base (2 reps) | 0.1714 | 0.517 | 0.2072 | - |
| | | **new6** | **0.1592** | **0.218** | **0.1653** (-20 %) | identical, 8 bins + 68 hst rows, both reps; rep2 == rep1 |

- Net on AG Car A: the rotation cycle drops from 0.717 to 0.346 s (normal-cycle cost + 0.09 s). The average drops from 0.319 to 0.263 s/cycle (**-17.5 %**). Steady cycles are -7.2 %.
- Fenced vet_gd timers over 32 cycles:
  - sweeps (including the batch): 4.39 s (base) / 3.28 s (new6) / 3.22 s (ipc).
  - tables: 0.51 s / 0 s.
  - moments + post: 0.41 s / 0.20 s.
- In the remaining rotation overhead (about 0.09 s), the batch is about 0.11 s per rotation, and it is GPU-bound: it is the marks, the scans and the lists of 464 shells.

## deepatm-1011 port (vgdrot2-1011, 3e790a6c)
- deepatm-1011 forks from 5b304cf9 and lacks the vgdfuse-1010 merge, so merging vgdrot-1011 also brings in vgdfuse-1010 (exact halo lists, the default halo_exact = 2).
- Conflicts were limited to rad_m1.hpp (members) and VetGdHalfRange stage 2. In stage 2 the batched gather is generalised to gdthin's nhc components (HREP skipped; C1 / J subtraction rules unchanged; same order).
- deepatm's sweep change (gdface typ-2 inner face) reads only ln chi / ln S, so the exact-list read marks are unaffected.
- DBASE = `bin/athena_he_a100_5cb2fec2_dbase` (deepatm-1011 5cb2fec2), NEW = `bin/athena_he_a100_m1` (3e790a6c). Input `inputs/agcva.athinput` (vetacc AG Car A input).

Jobs 31067401 / 31067402 / 31067403: 32 cycles, 2 interleaved reps, s/cycle.

| setup | arm | steady | rotation | avg | change in avg | bitwise (8 bins incl. m1_thsw / m1_fs + 68 hst rows) |
|---|---|---|---|---|---|---|
| production `all` (be) | dbase | 0.3103 | 0.924 | 0.3738 | - | - |
| | m1 | 0.2544 | 0.347 | 0.2639 | **-29 %** | identical in both reps; rep2 == rep1 |
| hesdirk2 + `time2_vet_col=predict implicit_vimp=true time2_vstage=true` (+ face default) | dbase | 0.4578 | 1.035 | 0.5175 | - | - |
| | m1 | 0.3962 | 0.482 | 0.4051 | **-22 %** | identical |
| the same + `implicit_blend_xthin=0 implicit_thin_switch=fse implicit_thin_corr=hc` | dbase | 0.4856 | 1.034 | 0.5423 | - | - |
| | m1 | 0.4220 | 0.477 | 0.4277 | **-21 %** | identical |

These gains include vgdfuse-1010 (exact halo lists), which deepatm-1011 lacked. On rt-integration, vgdfuse was already in the base and the gain was -17.5 %.

## Open
- BSG gate of new6: done on Raven, job 31067404 (bitwise, -20 %). The DeltaAI race can be cancelled.
- vet_gd_halo_ipc stays opt-in. The per-shell lockstep is now limited by rank imbalance in the shell kernel, not by the host.
- The remaining rotation overhead is the GPU work of the batch, about 0.09 s. It could move into the previous cycles together with the table prefetch, because the lists depend only on the angle, the cut and the depth.

