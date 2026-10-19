# TASK for Caltech and DeltaAI: Plaskett half-orbit stream runs q12_s1 / q12f_s1 (duplicates of viper)

**User 10-07:** run the two Plaskett-progenitor stream arms on Caltech (H200) and DeltaAI (GH200) as well as viper.

## 1. What / why

- Mass gainer of Plaskett's star at the start of accretion: the L1 stream hits a resolved n=3 polytropic envelope
  (pgen `ry_per_accretor`, `problem/inner = envelope`), 3-D spherical-polar equatorial wedge 500 x 4 x 2048, ideal gas.
- Two arms, half an orbit each (tlim 0.2291 code units; P_orb = 0.8519): **q12_s1** (PLM + hllc, nghost 2) and
  **q12f_s1** (same + FOFC, nghost 3). The earlier stream run (q10_s1) died by dt collapse at 0.077 orbit; these
  test whether the new grid/keys get through half an orbit.
- **Duplicates are intended.** viper jobs 12116831 (q12_s1) and 12116832 (q12f_s1) are PENDING on viper (2 nodes x
  2 MI300A). Whichever copy starts first wins; **do not cancel anything on another site, and do not cancel your own
  copy because another one started: the user decides.** Do not touch any other run of yours.

## 2. Code

- Branch **`accretor-1006`** on the fork, commit **`ea04cd19`** (`git fetch fork accretor-1006`; build that sha, not
  the tip: later commits on this branch, e.g. this note, are docs only). It is rt-integration (incl. 952abbfe
  c2p_track CUDA fix and 57fcd2b8) + the ry_per_accretor pgen commits f5dfc1ef, b226528c, 1e8d0ff8, fe16e65a,
  6588011b, ea04cd19.
- viper binary `athena_ryper_gpu72_env7`, md5 a5ad0a18c1d3a73e4cb41e833f1821d5, built from the clean worktree at
  ea04cd19 (binary mtime = commit time 04:51; no src/ file newer than the binary; tree clean).
- `cmake -D PROBLEM=ry_per_accretor -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release` + CUDA flags:
  - **DeltaAI** (GH200): as TASK-2026-10-01-deltaai-bsg.md section 6: `build_inc_deltaai.sh` with a new target,
    PrgEnv-gnu, gcc-native/14, cudatoolkit/25.5_12.9, cray-mpich/9.0.1, craype-accel-nvidia90,
    `Kokkos_ENABLE_CUDA=On Kokkos_ARCH_HOPPER90=On Kokkos_ARCH_ARMV9_GRACE=On`, nvcc_wrapper.
  - **Caltech** (H200): `docs/handover/caltech-2026-09-26/scripts/build_inc.sh <tag> gpu ea04cd19 ry_per_accretor`
    on a compute node; gcc/13.2.0, cuda/12.9.0, hpcx/2.17.1, `Kokkos_ENABLE_CUDA=On Kokkos_ARCH_HOPPER90=On`.
  - Record the md5. No CPU binary needed (no port gate beyond the smoke; the pgen reads no files).

## 3. Inputs (`docs/handover/plaskett-1007/`, md5 in `MD5SUMS`, byte-identical to viper)

| arm | input | md5 | command-line keys (verbatim, `keys_<arm>.txt`) |
|---|---|---|---|
| q12_s1 | `plaskett_env12.athinput` | 27ba4db9ca8411fb4965b0b3779bf431 | `time/tlim=0.2291 output2/dt=0.0115 output1/dt=0.0005 output3/dt=0.0229 -t 01:10:00` |
| q12f_s1 | `plaskett_env12f.athinput` | 3aee72ced5612fa93e58265afd457d78 | same |

- Inputs differ only in `nghost 3` and `fofc = true` + `fofc_report` (q12f). No data files; the only viper path is
  a comment (`# Design / trade-off: /viper/...`), harmless. Start from the input (no restart).
- Mesh 500 x 4 x 2048 in MeshBlocks of 500 x 4 x 256: **8 blocks -> 1 node x 4 GPUs, 2 blocks per rank** (viper: 2
  nodes x 2 MI300A = 4 ranks, the same decomposition).

## 4. Smoke (both arms, before production; same binary, input and keys + the smoke overrides)

```
srun -n 4 <recipe> <gpu binary> -i plaskett_env12[f].athinput <keys> time/nlim=20 time/ndiag=1 output3/dt=100 -d <smoke dir>
```

Pass, vs viper smokes q12sm / q12fsm (accretor_1006/gpu2, 2 ranks, 10 cycles) and q12_a1 / q12f_a1: rc 0, 0 FATAL /
nan, **dt at cycle 0 = 1.956427e-06** (both arms), **time at cycle 10 = 1.942194e-05** (both arms; to ~6 digits, the
rank layout differs). The 0.2-orbit envelope-only checks q12_a1 / q12f_a1 (stream off) ran 0 FATAL with dt 1.85-1.96e-6.

## 5. Production

- `srun -n 4 ... <gpu binary> -i plaskett_env12.athinput time/tlim=0.2291 output2/dt=0.0115 output1/dt=0.0005
  output3/dt=0.0229 -t 01:10:00 -d <run>/q12_s1` (and `plaskett_env12f.athinput` -> `<run>/q12f_s1`), one job per arm
  or both in one 1-node job sequentially.
- Wall limit **1:15 per arm** (viper used 1:15 with `-t 01:10:00`). Expected ~20-40 min per arm on 4 GPUs (q12_a1:
  2 MI300A, 2.5e8 zone-cycles/s; dt falls when the stream arrives). `time/dt_min = 1e-7` aborts a dt collapse:
  that is a result, not a failure to retry. No chaining, no automatic rerun.
- Caltech: `srun --mpi=pmix`, 1 rank per GPU, H200 (`--gres=gpu:nvidia_h200:4`, 1 node), run area
  `/resnick/groups/carnegie_poc/jingze/`. DeltaAI: ghx4, `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank`,
  `MPICH_GPU_SUPPORT_ENABLED=1`, `srun -n 4 -c 16 --cpu-bind=cores`, run area `/work/nvme/bivj/jma20/`.

## 6. Copy back / report

- Per arm: `run.log`, `ryper.hydro.hst`, `ryper.user.hst` (every 5e-4), `bin/ryper.hydro_w.*.bin` (every 0.0115,
  ~21 files x 82 MB = 1.7 GB), and the last `rst/` file only. Leave them on your run area and give the path; viper
  pulls them.
- NOTE `docs/handover/NOTE-2026-10-0x-<site>-plaskett.md` pushed to the fork (branch accretor-1006 fast-forward,
  never force; else branch `plaskett-results-<site>`): binary sha/md5, smoke job id + cycle-0 dt + cycle-10 time,
  production job ids, start time (immediately when it starts), end t, final dt, FATAL/dt-collapse time if any,
  wall time.
