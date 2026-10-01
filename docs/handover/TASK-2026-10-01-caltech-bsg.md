# TASK for Caltech: BSG arm 2 (Ma, Bildsten & Jiang 2026 reproduction), port gate + production (duplicate of viper)

**User 10-01 18:00:** queue the BSG arm-2 run on Caltech (and DeltaAI) in addition to viper. Machine: Resnick HPC,
**NVIDIA H200, 4 per node, x86 host**, CUDA. The DeltaAI copy is TASK-2026-10-01-deltaai-bsg.md.

## 1. Context

- **What:** reproduction of the blue-supergiant (BSG) envelope of Ma, Bildsten & Jiang 2026 (arXiv 2609.28656):
  20 Msun, log L/Lsun 5.161, Teff 15835 K, 3-D spherical wedge 20-61 Rsun x (pi/3 x pi/3), implicit M1
  radiation (vet_col closure, TOPS opacities), ideal gas mu 0.62, no SMR, ~13.5 M cells.
- **Arm 2** (this task) = lower top (61 Rsun), hydrostatic top ghosts (`he_bc_outer = hse`, g(1-Gamma) from the M1
  force) + marshak, a grey radiative-equilibrium atmosphere relaxed in 1-D, `dfloor 1e-20`, radial sponge
  (outer 1.8 Rsun) + floor-gas sponge (`he_sponge_dmax 1e-17`). Arm 1 (80 Rsun top, outflow) runs on viper only.
- **viper work dir** `/viper/ptmp2/jinma/bsg_1001/arm2/` (not visible to you). Its RESULTS.md, summarised:
  the relaxed column has L_top/L_in 0.9936-1.0014 over 2e5 s, Picard 0 non-converged / 3.0 passes, 47x fewer
  Newton fallbacks than arm 1, FOFC only in the floor gas (56-61 Rsun); the 3-D IC is the k,j-mean of the
  relaxed column at 4e5 s, v = 0, re-balanced only below 25 Rsun. Arm 2 costs 0.60x arm 1 per cycle.
- **Duplicates are intended.** viper runs arm 1 (job 12050971) and arm 2 (job 12051492, link 1); both PENDING on
  viper at the time of writing. Caltech and DeltaAI each get arm 2 too. Whichever copy starts first is the
  reference; the others are a cross-machine check. **Do not cancel anything on another site**, and do not
  cancel your own copy because another one started: the user decides.

## 2. Bundle (`docs/handover/bsg_1001_bundle/`, 1.5 MB, on rt-integration with this note; the CODE to build is bsg-arm2, section 3)

| file | what |
|---|---|
| `bsg3d_arm2.athinput` | 3-D production input = viper `arm2/prod/bsg3d_arm2.athinput`, only the 3 file paths replaced by `@BUNDLE@/...` |
| `bsg_col_arm2.athinput` | 1-D port-gate column = viper `arm2/col/ver25.athinput` (production keys and grid B' on a 206 x 4 x 4 thin column at the equator, 1 MeshBlock, no seed, frozen MLT), tlim cut to 2e4 s (161 cycles) |
| `ic_bsg_arm2.txt` | 3-D / column IC (relaxed arm-2 column), md5 6f28bd8e16c7144d756700f38e366e83 |
| `rosseland_tops_x0.7_z0.008.txt`, `planck_tops_x0.7_z0.008.txt` | TOPS opacity tables (X 0.7, Z 0.008), md5 in `MD5SUMS` |
| `SETUP.sh <run dir>` | checks the md5s and writes both inputs into `<run dir>` with `@BUNDLE@` = the absolute bundle path |
| `gate_compare.py <A> <B> <athenak>/vis/python` | port-gate comparison (hst columns, last hydro_w / m1 dumps, FATAL/nan counts) |

The inputs keep absolute paths after SETUP.sh, and the restart file embeds them: **keep the bundle directory where it
is for the whole run** (the IC file is read again on restart for the a_ref / frozen-MLT targets). Copy the bundle
out of the git checkout (e.g. into your run area) if you want it independent of later checkouts.

## 3. Code: branch `bsg-arm2` (NOT rt-integration yet)

- Build **fork/bsg-arm2 = d0abe21e** (`git fetch fork bsg-arm2`). It is e415d47d + 463da97c (he_star_m1:
  `he_bc_outer = hse`, `he_sponge_mode = radial`, `he_sponge_dmax`) + d0abe21e (make_ic_mlt_star.py `--atm re`,
  not needed to run). This is exactly the source of the viper arm-2 smoke and production binary
  (`athena_gpu72_463da97c`, md5 f8368b30227270f175d29f298dcc4db6; 463da97c and d0abe21e have the same src/).
- It **lacks 57fcd2b8** (rad_signal_speed: bitwise restarts and a safe first dt). Consequence: the first step after
  each restart uses the gas signal speed instead of the radiation-modified one, so a restarted run is not a
  bitwise continuation. Harmless (one step, implicit M1). Move to rt-integration once it contains bsg-arm2; a
  NOTE will announce that. Do not switch binaries in the middle of a link.
- PROBLEM = **he_star_m1**; build one GPU binary and one CPU binary (for the port gate) from the same commit.
- Known CUDA traps (DeltaAI notes 09-29), check they do not bite here:
  - an `inline` function in a header that contains a `KOKKOS_LAMBDA` and is compiled into several .cpp files gives
    an uninitialised nvcc lambda-wrapper pointer (SIGSEGV at startup, `__nv_hdl_wrapper_t` in gdb;
    NOTE-2026-09-29-deltaai-c2p-track-cuda.md). **That fix (952abbfe, branch fix-c2p-track-cuda, merged in rt-integration) is NOT in
    bsg-arm2**: `utils/c2p_track.hpp` still has the inline Before/After lambdas. It hit only the MHD path (hydro got
    the working copy), and BSG is hydro, so it should not bite; if the GPU run segfaults before the first cycle with
    `c2p_track` / `__nv_hdl_wrapper_t` in the backtrace, build bsg-arm2 + cherry-pick 952abbfe on a local branch,
    rerun 4a-4b with it, and say so in the report.
  - MPI requests that are posted and never waited on leak request handles and abort under Cray MPICH after a few
    minutes (NOTE-2026-09-29-deltaai-seam-mpi-leak.md; that cubed-sphere fix is not in bsg-arm2 either, and is not
    needed: BSG is spherical-polar). A run that aborts with `req != NULL` after some thousand cycles is this class.
  - M1 itself has run clean on both NVIDIA machines (bench-2026-09-29-hebox: dt, Picard and Krylov counts identical
    to viper). `vet_col` (the column closure, team scratch per radial column of 206 + ghosts) is new on CUDA: if
    it fails with a scratch-size error, report the message; do not change the closure.

## 4. PORT GATE (before any production)

### 4a. 1-D column, CPU vs GPU (`bsg_col_arm2.athinput`, 161 cycles to t = 2e4 s, 1 rank, ~10 s on a GPU)

```
bash docs/handover/bsg_1001_bundle/SETUP.sh <run>/in
cd <run>/c1 && <cpu binary> -i ../in/bsg_col_arm2.athinput > log 2>&1        # CPU, 1 rank
cd <run>/g1 && srun -n 1 <gpu binary> -i ../in/bsg_col_arm2.athinput > log 2>&1  # 1 GPU; repeat as g2
python3 docs/handover/bsg_1001_bundle/gate_compare.py c1 g1 vis/python
python3 docs/handover/bsg_1001_bundle/gate_compare.py g1 g2 vis/python
```

Viper reference (`/viper/ptmp2/jinma/bsg_1001/arm2/port_gate`, CPU login node gcc 14 binary of the same source;
GPU job 12056122, 1 MI300A, binary md5 f8368b30):

| quantity | viper CPU (c1) | viper MI300A (g1 = g2) |
|---|---|---|
| hst rows / t_end | 21 / 2.0e4 s (161 cycles) | same |
| dt at the end | 1.243014e+02 s | 1.243014e+02 s |
| Picard passes, mean | 3.000 | 3.000 |
| L_top/L_in, end | 0.999466 | 0.999466 |
| L_top/L_in, range | [0.999373, 0.999467] | same |
| FATAL / nan | 0 / 0 | 0 / 0 |
| Newton fallbacks (end-of-run `<rad_m1>` line) | 0 | 1 |
| wall | 15 s (login node) | 8.0 s / 9.7 s |

CPU vs MI300A, `gate_compare.py c1 g1`: hst columns max rel diff dt 1.0e-8, L_top 6.8e-10, L_bot 2.6e-7,
E_rad 1.7e-8, e_gas 1.6e-8, KE_int 1.9e-5, v1sq_wall 1.8e-3, Min_top 2.3e-4; last dumps (scale = max |A|):
dens 6.8e-8, eint 9.2e-8, m1_e 1.2e-7, m1_f1 3.8e-5; velx 5.8e-2 over all cells but **4.0e-6 with rho > 1e-17**
(the difference sits in the floor gas at 54-61 Rsun). vely/velz/m1_f2/m1_f3 differ O(1): in a 4 x 4 column they are
round-off-seeded noise (|v_t| < 4e2 cm/s vs |v_r| 2.7e4 cm/s below the floor gas). The two GPU repeats were
**bitwise identical** (all hst columns and dumps).

**Pass criteria on your machine:** your CPU vs your GPU, and your GPU vs the viper numbers above: dt, Picard mean
and L_top/L_in to the printed digits; hst energies and dens/eint/m1_e at <= ~1e-6; velx (rho > 1e-17) <= ~1e-4;
0 FATAL/nan; GPU repeats bitwise (if not, report it: the viper GPU path is deterministic since the FOFC one-writer
fix). Anything at 1e-3 or worse in the energies is a port bug: stop and report.

### 4b. 3-D smoke (production input, 20 cycles)

Production is 4 MeshBlocks of 206 x 128 x 128 = **1 node x 4 GPUs** on your machine (1 rank per GPU, 1 block per
GPU): the smoke IS the production layout.

```
srun -n 4 <gpu binary> -i <run>/in/bsg3d_arm2.athinput -d <smoke dir> time/nlim=20 time/ndiag=1
```

Viper reference:
- smoke job 12051483 (same input and binary, 2 MI300A with 2 blocks each, 10 cycles, 2 repeats): rc 0,
  0 FATAL/NaN, 0 Newton fallbacks, Picard 3.0, **dt 88.18 s** (88.177 at cycle 9-10), hst of the repeats
  byte-identical;
- timing job 12051329 (apudev, 1 block of 206 x 128 x 128 per MI300A, the production per-GPU load, 60 cycles,
  2 repeats interleaved with arm 1): **1.147 / 1.140 s/cycle**, dt 88.16, Picard 6.3 mean (that run used the
  pre-relaxation IC; with the final IC the smoke gives 3.0, so expect 3-6).

Pass: rc 0, 0 FATAL/nan, dt 88.2 s (+-0.1 %), Picard ~3-6, no NOT-CONVERGED. Record s/cycle (median of cycles
5-20, from the `elapsed=` log lines) for the report; H200/GH200 ran the M1 He box at about MI300A speed, so expect
~1.1 s/cycle.

## 5. Production: arm 2 fresh start

- `-i <run>/in/bsg3d_arm2.athinput` unchanged: tlim 4.96e6 s (57.4 d), cfl 0.3, rst every 8.64e4 s (1 d), bin
  hydro_w / m1 / m1_vet every 4.32e4 s, hst 2000 s, log 1e4 s, fofc_report every 4.32e4 s (as viper).
  `implicit_resid_fatal 1e-2` stays on (FATAL on a diverged Picard solve).
- 1 node x 4 GPUs, wall limit = site maximum (24 h), run with `-t <limit - 20 min>` so the final rst is written.
- Expected length: 4.96e6 / 88.2 = 56,300 cycles; at ~1.1 s/cycle ~17-19 h, i.e. one link if dt holds. The
  convective velocities may lower dt later, then link 2+.
- Chain links with `afterany` and a guard at the top of every link (as viper `prod/link2.sh`):
  1. stop if the previous `run.log` tail contains nan / FATAL / error, or rc != 0 without the wall-limit message;
  2. stop if the last hst dt dropped below 0.5x the first link's dt (dt jump) and report;
  3. restart from the newest `rst/bsg3d.*.rst` with `-r <rst> -d <run dir>`; keep the previous `run.log`
     (rename it), single attempt per link (deterministic binary, no blind retries).
- Smoke rule: run 4b again (10 cycles) after ANY change to binary, input or command-line keys, and for the
  first restart link (`-r <rst> time/nlim=<ncycle+10>` in a scratch dir copy).

## 6. Caltech specifics

- Build: `docs/handover/caltech-2026-09-26/scripts/build_inc.sh bsg gpu <sha of fork/bsg-arm2> he_star_m1` on a
  compute node (header of the script: `sbatch -A carnegie_poc -p expansion -c 32 --mem=64G ...`), and the same with
  `cpu` for the gate binary. Stack as the hebox bench: gcc/13.2.0, cuda/12.9.0, hpcx/2.17.1; `Kokkos_ENABLE_CUDA=On`,
  `Kokkos_ARCH_HOPPER90=On`, nvcc_wrapper. Record both md5s. (bsg-arm2 is a branch, not rt-integration: fetch it
  explicitly, `git fetch fork bsg-arm2`, and pass the sha.)
- Run: `srun --mpi=pmix`, 1 rank per GPU, one GPU visible per rank, **H200 only** (`--gres=gpu:nvidia_h200:4`, 1
  node), exclude hpc-sm-01-09 and hpc-sm-02-16 as before. Run area under `/resnick/groups/carnegie_poc/jingze/` (not
  scratch: purged after 14 days). An H100 smoke is fine if labelled; production on H200.
- Report: `docs/handover/NOTE-2026-10-0x-caltech-bsg.md`.

## 7. Report back (NOTE pushed to the fork)

One NOTE as named in section 6, pushed to rt-integration on the fork (fetch, fast-forward on top of the fork tip, never
force; if that fails, push branch `bsg-results-caltech` and say so). Content:
1. binaries: commit sha, md5 (GPU and CPU), stack; any port fix you needed (e.g. the 952abbfe cherry-pick);
2. gate 4a: the summary lines of `gate_compare.py` (CPU vs GPU, GPU vs GPU) next to the viper numbers, pass/fail;
3. smoke 4b: job id, rc, FATAL/nan, dt, Picard, s/cycle (median cycles 5-20);
4. production: job ids per link, start time, s/cycle and dt per link, the guard results, final t, and the run dir;
5. anything unusual (NOT-CONVERGED, Newton fallbacks per cell-pass, FOFC below 56 Rsun in the fofc_report).
Update the NOTE when production starts and when it ends (or fails); no other interim reports.
