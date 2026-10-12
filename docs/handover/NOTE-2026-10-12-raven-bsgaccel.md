# NOTE bsgaccel-1011: speeding up the production rad_m1 implicit VET scheme (BSG true repro, AG Car A, He box)

Worker M (Raven). Work dir /raven/ptmp/jinma/bsgaccel_1011. Branch **bsgaccel-1011** (worktree
/raven/u/jinma/ATHENAK/wt-bsgaccel), based on prodfinal-1011 df483770. Local only: not pushed.

Commits:
- de543c12: flat 1-D reductions (rad_m1_parfor.hpp `par_reduce_flat`).
- 458afa91: vet_gd halo. Lean IPC protocol and merged main+twin MPI messages.
- 06c8a502: hydro x2/x3 flux sweeps chunked when there are too few pencils.

Every code change is **bitwise**: identical bin data and hst rows on GPU AG Car A, GPU reduced BSG, CPU reduced
AG Car A, the CPU rw_t10a hesdirk2 T/C ladders and the CPU sp rsw_T / rlat_T ladders. The only
result-changing recommendation is an existing opt-in key, `implicit_precond = mg_gf`, recommended for the BSG only.

Binaries (bin/):
- base = prodfinal_1011/bin/athena_he_a100_df483770 (md5 41db0c68).
- **athena_he_a100_c1** (md5 433867b7) = HEAD 06c8a502.
- athena_he_cpu_c1 (7d8aae5a) and athena_cpu_c1 (15936bcd).
- athena_box_a100_{base,c1} (PROBLEM=box_convection).

## 1. Profile (step 1). Production keys, 1 node x 4 A100, 1 MeshBlock per GPU, base binary

Cases:
- **AG Car A**: the production FRESH deck prodfinal_1011/prod/agcar_A_fresh.athinput. 4 blocks 480x64x64.
- **reduced BSG**: the bsgrec R deck (production scheme on bsg3d_truerepro2, original IC). Mesh 256x128x128, 4 blocks
  256x64x64, lowmem as in the deck.

Both decks are copied to inputs/ with the extra keys named at their defaults (bitwise vs the originals, checked).
Output: hst every cycle, bins at t = 0 and at the end, rst off.

Timing is the per-cycle wall from the `elapsed` lines over cycles 5..N-2 (jobs/ana.py):
- steady = cycles that are not multiples of 10;
- avg = all cycles, including the vet_gd rotation cycles.

| | AG Car A | reduced BSG |
|---|---|---|
| s/cycle steady / avg | 0.323 / 0.334 | 0.258 / 0.262 |
| Picard passes per solve (2 solves per step) | 2.63 (mean 2.71 incl. start) | 2.34 |
| Krylov iterations per pass | 2.7 (mg) | **21** (rbgs_fwd) |
| tensor (vet_gd formal solve, 1 build per step) | 0.138 s/cycle (40 %) | 0.076 (29 %) |
| Picard pass setup (assembly: asm, vimp_j, tcell, vlat_*) | 0.062 (18 %), about 12 ms per pass | 0.029 (11 %) |
| Krylov | 0.055 (16 %) | **0.121 (46 %)** |
| hydro + c2p + bvals + rest | about 0.06 | about 0.03 |

- From nsys (rank 0):
  - The GPU is busy 76-79 % of a cycle.
  - Most of the idle time is the vet_gd per-shell lateral halo: a host fence plus 12 MPI messages per shell, 460
    shells per cycle. The idle time before m1_vgd_hl_scatter is 48 ms/cycle (AG Car) and 31 (BSG).
  - The shell kernel itself is 58 / 29 ms per cycle, latency bound (117 registers, occupancy 25 %).
- From Kokkos kernel-timer and ncu:
  - The BSG Krylov is bandwidth-bound stencil + pcrx line sweeps at about 0.1 ms per call, 6 calls per iteration.
  - MDRange reductions run about 50x below bandwidth (m1_impl_vmmin 0.9 ms on 2M cells).
  - hflux_x2 / x3 launch only 66 teams with one block per GPU (occupancy 6 %, 2.2 ms each).
- Picard (implicit_picard_log):
  - BSG solves are at the 2-pass floor.
  - In AG Car about 35 % of solves need a third pass. The pass-1 residual is 5e-6..2e-5 at the top cells
    (i 455-482). Not attacked.

## 2. Candidates

Measured on the same node, interleaved, 2 reps (identical to 1e-3 s).
- Speed: steady s/cycle change vs base.
- Accuracy: vs base, last dump plus hst. tol = implicit_tol 1e-8.

### Code changes (bsgaccel-1011, all bitwise)

| # | change | AG Car A | red. BSG | accuracy | scope |
|---|---|---|---|---|---|
| C1 | flat reductions (emax, rhsmax, vmmin, gcnt, gms, eckq) | -3.0 % | -2.0 % (-3.1 % with mg_gf) | bitwise | every implicit solve |
| C2 | vet_gd: main+twin in one MPI message per partner (env VGD_MERGE=0 = old) | -2.8 % | 0 | bitwise | vet_gd MPI path |
| C3 | vet_gd_halo_ipc lean protocol: no ack kernel, 1 fence per flag kernel (env VGD_IPC_LEAN=0 = old) | IPC now -4.0 % vs MPI (old IPC: -0.5 %) | IPC -8.3 % vs MPI (old IPC: +3 %) | bitwise | node-local layouts only |
| C4 | hydro x2/x3 flux sweep chunks below 432 pencils | -1.6 % | -2.7 % | bitwise | all hydro on 1-2 blocks per GPU |

### Existing opt-in keys

| key | AG Car A | red. BSG | verdict |
|---|---|---|---|
| **implicit_precond = mg_gf** | +10.6 % (mg already 2.7 its) | **-39 %** (Krylov 21 -> 1.0 it/pass; Picard 2.68 -> 2.48) | **BSG only**. hst ≤ 6e-8 (transient), fields ≤ 8e-8 after 200 cycles, no growth. Needs periodic x2/x3 |
| vet_gd_halo_ipc = true | see C3 | see C3 | recommend wherever all lateral partners are node-local (Caltech 2 H200, DeltaAI 4 GH200, Raven 1 node). Bitwise |
| vet_gd_twin_lowmem = false | - | -5.6 % on top | where the memory fits (bitwise, DeltaAI mem2) |
| implicit_precond = mg | - | -23 % | **rejected**: implicit_vimp positivity fallback, stage fallback to BE, linear residual 2e-6. L_top changed by 1e-2 |
| implicit_precond_float | 0 | 0 | no |
| implicit_halo_ipc (Krylov) | +2.9 % | +3 % | no |
| implicit_op_team_red | 0 | 0 | no |
| implicit_eos_cache | +1.3 % | - | no |
| implicit_one_pass = 8 | +3 % | +1 % | no (BSG already at 2 passes) |
| implicit_accel = anderson | FATAL ("shared memory request is too large", A100) | same | broken on GPU here |
| implicit_gf_modes 1 / 3 | - | = 2 | keep 2 |
| LaunchBounds caps (shell kernel 128x5/6/8, par_for_lb 256x2/3) | 0 / +4 % | +1..+7 % | dropped (bitwise, no gain) |

### Calibration

implicit_tol 5e-9 vs 1e-8 on BSG gives hst differences ≤ 5e-13 and fields < 6e-8 (float32 dumps).

## 3. Order and bitwise checks

- CPU reduced AG Car A (480x32x32, 4 ranks x 3 threads, production FRESH deck, 30 cycles): c1 vs base **bitwise**.
  Also bitwise with VGD_MERGE=0.
- CPU rw_t10a ladders tJePvs_hesdirk2 T (1..16) and C (32..256): every file **bitwise** base vs c1, so the orders are
  unchanged (inputs from suiteC_1011/ladder).
- CPU sp ladders rsw_T and rlat_T (HFfM, 1..8): **bitwise** base vs r4 (= c1 code).
- mg_gf on the sp ladders: it needs periodic theta and nx2 > 2K.
  - Added rsw_Pt_T / rlat_Pt_T (periodic theta) with gf_modes 1.
  - The orders are identical to the default arm: rlat 1.54 / 1.33 (the known hesdirk2 lateral value); rsw eint
    1.99 / 1.85, E 1.76 / 1.40.
  - The solutions are identical to ≤ 1e-12 (L1 rel).
  - These ladders take ≤ 1 Krylov iteration either way, so this is a weak test of mg_gf. The BSG comparison above is
    the real one.
- GPU: AG Car A c1 and c1 + IPC vs base are bitwise (bin + 84 hst rows), and so is every intermediate build.
  BSG c1 vs base (200 cycles) is bitwise.

## 4. Combined results (A100, 2 reps; avg = including rotation cycles)

| case | base steady / avg | new steady / avg | speedup (avg) | new config |
|---|---|---|---|---|
| AG Car A, 4 x 1 block | 0.3229 / 0.3337 | 0.2995 / 0.3104 | 1.075x | c1, MPI path (multi-node layouts) |
| AG Car A | same | **0.2872 / 0.2978** | **1.12x** | c1 + vet_gd_halo_ipc |
| red. BSG, 4 x 1 block | 0.2579 / 0.2624 | 0.1351 / 0.1393 | 1.88x | c1 + mg_gf + ipc |
| red. BSG | same | **0.1293 / 0.1335** | **1.97x** | + vet_gd_twin_lowmem = false |
| red. BSG, 200 cycles | 0.2392 / 0.2457 | 0.1326 / 0.1376 | 1.79x | c1 + mg_gf + ipc (lowmem on) |
| BSG 2 blocks/GPU (256x256x128) | 0.4555 / 0.4644 | 0.2420 / 0.2483 | 1.87x | c1 + mg_gf + ipc (lowmem on) |
| He box (Cartesian vet_sc, hesdirk2, 4 x 84x52x52) | 0.0481 / 0.0481 | 0.0460 / 0.0460 | 1.045x (bitwise, hst) | c1 (mg) |
| He box | same | 0.0451 / 0.0450 | 1.07x | c1 + mg_gf (Krylov 6.2 -> 4.1 it/pass; Fres hst 6e-7, F1 4e-8) |

### He box (bench-2026-09-29-hebox deck)

Built with PROBLEM=box_convection, data vimphr_1009/tab. Timers over 55 cycles: Krylov 21 ms/cycle (44 %), vet_sc
tensor 4.7 ms (10 %), pass setup 5 ms (11 %), Picard 2.17 passes per solve. So the hot spot is again the Krylov
(mg, 6 it/pass). The vet_gd changes do not apply here; C1 and C4 do. mg_gf -2 %: optional, not recommended.

### Full BSG extrapolation

A linear fit per blocks-per-GPU n on A100 (avg s/cycle) gives:
- base: 0.060 + 0.202 n;
- new: 0.030 + 0.109 n.

The ratio is 1.86 at n = 4 and at n = 8, and 1.79 on the 200-cycle run. **Use 1.8x.** dt is about 97 s, so 57.4 d
is about 51k cycles.

| layout | est. base | est. new | 57 d wall |
|---|---|---|---|
| 2 H200 (Caltech), 8 blocks per GPU, lowmem, IPC | 1.9-2.1 s/cycle (main's estimate) | about 1.05-1.15 | **about 15-17 h** (from 27-30 h) |
| 4 GH200 (DeltaAI), 4 blocks per GPU, lowmem off (73.6 GB be/hr), IPC | about 1.25 | about 0.65-0.7 | **about 9-10 h** (from about 18 h) |
| Raven 16 A100 (4 nodes), 1 block per GPU, MPI path, lowmem off | about 0.26-0.3 | about 0.14-0.17 | about 2-2.5 h (from about 4 h) |

The Raven row is the reduced-BSG per-GPU timing. The inter-node MPI per-shell latency is not measured; the
Raven 4-node gpu queue wait dominates anyway.

## 5. Recommendation

1. Merge bsgaccel-1011 (3 commits, bitwise, no new keys; the env switches are for A/B only).
2. BSG decks: `implicit_precond = mg_gf`. Plus `vet_gd_halo_ipc = true` when single-node. Plus
   `vet_gd_twin_lowmem = false` when it fits (4 GH200: yes; 2 H200 8 blocks/GPU: no).
3. AG Car decks: `vet_gd_halo_ipc = true` (1 node x 4 A100: -4 %, bitwise). Keep mg. A restart keeps working:
   the key is read only when named.
4. Not done / open:
   - The vet_gd shell kernel (latency bound, about 20 % of a cycle after these changes).
   - The AG Car third Picard pass at the top cells (about 5 % potential).
   - A hybrid IPC/MPI halo for multi-node layouts.
   - GH200/H200 confirmation of the full-BSG numbers (TASK-2026-10-12-deltaai-bsgaccel-fullbsg.md, prepared, not
     posted).
