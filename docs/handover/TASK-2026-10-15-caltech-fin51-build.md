# TASK raven -> Caltech: build FIN5.1 (fin51-1014, fin51-linhao-bc) on H200, nofma; short bitwise + forced-retry smoke — debug QoS

DRAFT from Raven worker F51 (fin51_1014) for the coordinator to push (edit hashes if they change before the push).
Background NOTE on Raven: /raven/ptmp/jinma/fin51_1014/NOTE.md.

FIN5.1 = FIN5 (fin5-1014 ff03587a: band ring + accel, bitwise to FIN4) + fin4-dtretry (5b7e1bc5 + 053bd996 log wording: <time>/dt_retry,
cycle checkpoint + dt cut-and-retry, default ON for fresh runs, never-triggered = bitwise) + an int-overflow fix:
rad_m1_icrad.cpp `vgd_itw.d.extent_int(0) > 1` -> `extent(0)` (with a fused twin and > 2^31 twin values per rank the fresh-start
radiation IC skipped the twin J correction: "ic_rad: no fused twin intensities"; e.g. FIN4 He giant LARGE 1697 smoke). Bitwise to FIN5
wherever no overflow occurs (restart-size products / rotavg rst copy hardened too, values unchanged).
The dt_retry checkpoint visitor was regenerated for the FIN5 members (ring bookkeeping, accel, NCCL state).
ONE code for ALL new production (He giant large, E1 remap/wedge, BSG). Do NOT switch, restart or touch any running/queued production.

## 1. Builds (report md5 + `git rev-parse HEAD`)
- `fin51-1014` @ `21c4839d` -> `athena_gpu_he_star_m1_21c4839d_nofma`
- `fin51-linhao-bc` @ `65c249b3` (= fin5-linhao-bc 7299d0d7 + merge fin51-1014) -> `athena_gpu_he_star_m1_65c249b3_nofma`
- EXACTLY the toolchain/flags of your FIN4 / 5b7e1bc5 / FIN5 nofma builds:
  `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`,
  default `Athena_M1_VGD_NMAX`, no -DM1_VGD_NCCL. Build in a compute job (new build dir; CMakeLists adds driver/dt_retry.cpp).

## 2. Smoke (one debug job, 1 node x 2 H200, your production layout; fresh, rst at the end)
Base = your FIN5 build if you have one (ff03587a / 7299d0d7 nofma), else your 5b7e1bc5 nofma (E1 deck) / FIN4 d505d21f nofma.
Arms per deck (He giant LARGE deck in the 2-GPU production layout if it fits, else your HG FIN4 deck; your E1S deck; Linhao-bc deck with 65c249b3):
- B base, `time/nlim=12`;  M FIN5.1 same;  M1 FIN5.1 nlim 6 then M2 `-r` -> 12.
- F = M with two keys ADDED UNDER `<time>` in a deck copy (cannot be added on the command line of a fresh run):
  `dt_retry_force_cycle = 5` and `dt_retry_force_fac = 1.0`  (one forced retry, same dt) -> run.log has ONE `DT_RETRY cycle=5 ... reason=forced`.
- FP = F with env `VGD_RING_POISON=1` (NaN in lent ring slabs).
- Checks: final rst data byte-equal after `<par_end>`: B = M = F = FP, M = M2; user.hst rows identical; `out-of-slab accesses main 0 twin 0`,
  `dense band saves 0`; FATAL / NaN 0; fresh-start log has NO `ic_rad: no fused twin intensities` line (fused-twin decks).
  If base is FIN4 / 5b7e1bc5 on a fused-twin deck that printed that line, B differs from t = 0 by design: then compare M = F = FP, M = M2 only.
- Record s/cycle (3..11) and GPU peak (nvidia-smi max) B vs M (expected: +~2.5 GB/GPU E1 checkpoint, +0.5..1 % s/cycle).

## 3. Report: docs/handover/NOTE-2026-10-15-caltech-fin51-build.md
md5s + HEADs, job id, the bitwise table, s/cycle + GB/GPU, any non-bitwise line in full.
Expected (Raven A100 / Viper MI300A results in the Raven NOTE): all BITWISE.

## Coordinator notes
- Priority: short (debug QoS). Do not delay running production. Production use of FIN5.1 = a separate GO.
- E1 decks: keep `time/dt_retry = true` explicit (a restart of an old rst without the key keeps it off). Link hooks: DT_RETRY lines carry
  no STOP words (053bd996); e1chain_hook.sh already ignores DT_RETRY lines.
