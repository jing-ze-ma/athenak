# NOTE caltech -> raven: FIN5.1 build + smoke on H200 (TASK-2026-10-15-caltech-fin51-build, f4af5026)

Caltech, 2026-10-14 (US time). Work dir /resnick/groups/carnegie_poc/jingze/fin51_1015 (buildq.sh, smoke51.sh, decks/, smk_*.out).
No running or queued production was touched.

## 1. Builds (compute jobs 4537256 + 4537257, expansion debug QoS; NEW tree on /resnick/groups/carnegie_poc/jingze/builds_fin51)
Same script (bsg_1009/build_inc_nofma.sh, BUILD_ROOT moved off /resnick/home), toolchain and flags as FIN4 / 5b7e1bc5 / FIN5 nofma:
gcc 13.2.0 + cuda 12.9.0 + hpcx 2.17.1 ompi, Kokkos 4.6.2 (d8e9af03), `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON
-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off`,
default Athena_M1_VGD_NMAX, no -DM1_VGD_NCCL (0 ACC_NCCL strings). Fresh configure (driver/dt_retry.cpp picked up).

| binary | git rev-parse HEAD | md5 |
|---|---|---|
| athena_gpu_he_star_m1_21c4839d_nofma | 21c4839d7eeea9485b67da435524b89a1e49a2ef | f06db8e6c49a5014db07577289196d8b |
| athena_gpu_he_star_m1_65c249b3_nofma | 65c249b345ef50655df8e510ccc9dfe36967972a | edfbb474f8935b3eeb4e443b85260f99 |
| base: athena_gpu_he_star_m1_ff03587a_nofma (FIN5) | ff03587a | c05170d98488f1f8fc4978c00b763fbc |
| base: athena_gpu_he_star_m1_7299d0d7_nofma (FIN5 Linhao) | 7299d0d7 | b71887f49368b556cc816761780c184b |

Binaries in /resnick/home/jingze/ATHENAK/builds/bin (copies of builds_fin51/bin).

## 2. Smoke (gpu debug QoS, 1 node x 2 H200, 2 ranks, srun --mpi=pmix -n 2 -c 8, 1 GPU per rank, UCX_TLS=^xpmem)
Jobs: He giant LARGE 4537120 (M1 M2 M) + 4537121 (B F FP) + 4541305 (M1 M2 redo, see 2b); E1S 4537122; Linhao 4537123.
Decks: He giant LARGE in the 2-GPU layout (it fits: 2 MeshBlocks 1697x64x64 per H200) = hegiant_large_1014/hegiant_large_fin51.athinput
(= 09a3b446 + `time/dt_retry = true`; md5 8a490843 for B/M/F/FP/first M1, 8ee79d2c for the M1/M2 redo = same keys, inline comments
stripped, see 2b); E1S e1s_1014/agcA_E1_fin4dtr_s3_tlim1e9.athinput (680dadb0); Linhao bsgli_1014/files/bsgLI_bc_caltech.athinput
(d77f9a10). F decks = deck + `dt_retry_force_cycle = 5`, `dt_retry_force_fac = 1.0` under `<time>` (fin51_1015/decks/*_F.athinput).
Bases: B = ff03587a (HG, E1S), 7299d0d7 (Linhao); M/M1/M2/F/FP = 21c4839d (HG, E1S), 65c249b3 (Linhao).
Keys: time/nlim 12 (M1 6), time/ndiag=1, output1 (hst) dt=1e-6, every other output of the deck dt=1e30 (t=0 + end, rst at end).
rst compare = bytes after `<par_end>`; bins = every hydro_w / m1 / m1_vet (/ m1_face) dump, data after `<par_end>`; hst = data rows
(M2: every row verbatim in M).

| pair | He giant LARGE | E1S | Linhao bc |
|---|---|---|---|
| B = M (base vs FIN5.1) | BITWISE rst (20 302 097 308 B), all bins, both hst | BITWISE rst (5 101 725 932 B), all bins, both hst | BITWISE rst (9 043 275 188 B), all bins, both hst |
| M = M2 (6 -> `-r` -> 12) | BITWISE rst, final bins, 6/6 rows verbatim (redo 4541305) | BITWISE | BITWISE |
| M = F (forced retry) | BITWISE rst, bins, both hst | BITWISE rst, bins, both hst | rst + bins BITWISE, hydro.hst identical; user.hst: ONE value differs (below) |
| M = FP (+ VGD_RING_POISON=1) | BITWISE | BITWISE | as M = F (F = FP fully BITWISE incl. user.hst) |

Linhao non-bitwise line (diagnostic only, state bitwise): user.hst column 7 `L_in`, row 6 (t = 525, the step after the forced
retry at cycle 5): M 1.2882820390310911e+38, F 1.7963250591179169e+38 (x1.394); all other rows/columns identical. The L_in
accumulator (Linhao-bc member, 300b6ce8) keeps the rejected attempt's contribution across a retry; it is not reset by the
dt_retry rollback. Cosmetic for physics, but L_in is an energy-budget column: worth a reset in the retry path.

DT_RETRY: F and FP each have exactly one `DT_RETRY cycle=5 t=... dt X -> X reason=forced resid=0.000000e+00 cut 1`; M/M1/M2 none.
Log checks, all arms: FATAL 0 (except the first M2, 2b), NaN 0, NON-CONVERGED>0 0, `out-of-slab accesses main 0 twin 0`,
`dense band saves 0` (FP prints "(POISON test)"), `ic_rad: no fused twin intensities` 0 in EVERY arm (incl. the FIN5 base on the
large deck at 2 blocks/GPU).

### s/cycle (cycles 3..11, one 8-cycle sample, hpc-sm-02-10) and GPU peak (nvidia-smi max, both GPUs equal)
| deck | B s/cyc | M s/cyc | F (1 retry) | FP | B peak MiB | M peak MiB | M - B |
|---|---|---|---|---|---|---|---|
| He giant LARGE (2 blocks/GPU) | 1.742 | 1.753 (+0.6 %) | 1.975 | 2.121 | 91 633 | 106 591 | +14 958 MiB (+14.6 GiB) |
| E1S | 0.4535 | 0.4614 (+1.7 %) | 0.518 | 0.538 | 31 313 | 36 035 | +4 722 MiB (+4.6 GiB) |
| Linhao bc | 1.316 | 1.310 (-0.5 %) | 1.460 | 1.523 | 51 525 | 60 965 | +9 440 MiB (+9.2 GiB) |
The checkpoint memory is ~1.5-1.9x the per-rank rst data (rst/rank: HG 10.2 GB, E1S 2.55 GB, LI 4.5 GB), more than the
+2.5 GB/GPU expected for E1; all fit (max 106.6 GB of 141 GB per H200). M1 (6 cycles) peaks 108 119 MiB on HG LARGE.

## 2b. FOUND: restart of the He giant LARGE deck FAILS with FIN5.1 (rst header > 45 056 bytes)
First M2 (`-r` M1/rst/hegiant.00001.rst, job 4537120): `### FATAL ERROR in parameter_input.cpp at line 268  <par_end> is not found
in the first 40KBytes.` ParameterInput::LoadFromFile reads 4096-byte chunks and FATALs once > kBufSize*10 bytes are read without
`<par_end>`, so the effective limit is 11 x 4096 = 45 056 bytes. `<par_end>` line end offsets: FIN4 large gate A rst 44 845 (just
fit), FIN5 B 44 952, FIN5.1 M 45 282, F 45 363 -> FAIL. For the same deck FIN5.1 dumps ~330 B more parameters than FIN5, which pushed it over.
Workaround here (no key changed): hegiant_large_fin51.athinput now has the inline `# ...` comments stripped (408 keys parse identical,
`<comment>` block untouched), md5 8ee79d2c; header 39 186 -> redo M1/M2 4541305 BITWISE vs M. Other headers now: E1S 38.9 KB,
Linhao 40.0 KB (F), Caltech FIN4 He giant 200-Rsun production rst 41.6 KB, bsgLIw4 rst 39.5 KB: all under 45 056, but with
only 4-6 KB margin. Suggest raising the limit in parameter_input.cpp (e.g. kBufSize*64) in the next code drop.

## Not done / open
- No production switched to FIN5.1 (separate GO). He giant LARGE gate on FIN5.1 (NOTE aa27fbd6) runs next: job 4537065 (gate51).
- Smoke dirs hold ~200 GB of rst (fin51_1015/hgl, e1s, li); can be pruned once Raven has read this.
