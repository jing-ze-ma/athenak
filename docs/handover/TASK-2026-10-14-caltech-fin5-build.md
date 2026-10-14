# TASK raven -> Caltech: build FIN5 (fin5-1014, fin5-linhao-bc) on H200, nofma; short bitwise smoke vs your FIN4 builds — debug QoS

DRAFT from Raven worker MERGE (fin5_1014) for the coordinator to push (edit hashes/branch names if they change before the push).
Background NOTE on Raven: /raven/ptmp/jinma/fin5_1014/NOTE.md.

USER adopted two BITWISE engineering branches on top of FIN4 (no method change, every result bitwise to FIN4):
- mem2-1014: vet_gd band ring (compact per-shell band intensities, lent wide slabs) + fin4-mem buffer items: -30..-60 % device memory.
- accel-1014: parallel gf_solve / gf_fac, twin/hr shell sums with the per-sum order kept, Anderson A100 fix,
  NCCL/RCCL list halo (COMPILED OUT by default; needs -DM1_VGD_NCCL + env ACC_NCCL=1; NOT for Caltech now).
FIN5 is for NEW runs only. Do NOT switch, restart or touch any running/queued production (BSG Linhao chain, He giant chains, E1 e1s_1014).

## 1. Builds (report md5 + `git rev-parse HEAD`)
- `fin5-1014` @ `ff03587a` -> `athena_gpu_he_star_m1_ff03587a_nofma`
- `fin5-linhao-bc` @ `7299d0d7` (= fin4-linhao-bc 300b6ce8 + the same merge) -> `athena_gpu_he_star_m1_7299d0d7_nofma`
- EXACTLY the toolchain and flags of your FIN4 / 8e018272 / 5b7e1bc5 nofma builds:
  `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`,
  default `Athena_M1_VGD_NMAX`. No -DM1_VGD_NCCL. Build in a compute job (incremental from your FIN4 build dir is fine).
- Note: E1 on Caltech runs fin4-dtretry (5b7e1bc5, dt_retry). FIN5 does NOT contain dt_retry (needs a fin5-dtretry branch with the
  regenerated checkpoint visitor; not made yet). Do not use FIN5 for E1 until such a branch exists.

## 2. Smoke (one debug job, 1 node x 2 H200, your production layout; keys NONE beyond nlim / outputs)
For each pair (FIN4 binary you already run, FIN5 binary): same deck, fresh, `time/nlim=12`, rst at the end; FIN5 also 6 -> `-r` -> 12.
- He giant: your `hegiant_fin4.athinput` (f266cd9b), base `athena_gpu_he_star_m1_d505d21f_nofma` (1440c053) vs `ff03587a_nofma`.
- BSG Linhao: your Linhao-bc deck, base `athena_gpu_he_star_m1_8e018272_nofma` vs `7299d0d7_nofma` -- NOTE 7299d0d7 contains 300b6ce8
  (restart-safe L_in history column, bb735538/c9ca4669): the rst/hst may differ from 8e018272 in the L_in history column / pgen
  trailer ONLY; compare the hydro + m1 rst data blocks and all hst columns except L_in. If you have a 300b6ce8 build, use it as base
  (then everything must be byte-equal).
- Extra arm FIN5 with env `VGD_RING_POISON=1` (fills lent band slabs with NaN): must be bitwise to FIN5.
- Checks: final rst data byte-equal after `<par_end>` (FIN4 vs FIN5, FIN5 vs FIN5-poison, FIN5 vs FIN5 restart); user.hst rows identical;
  run.log lines `out-of-slab accesses main 0 twin 0` and `vet_gd_band_ring ... dense band saves 0`; FATAL / NaN 0.
- Record s/cycle (cycles 3..11) and GPU peak memory (nvidia-smi max) for FIN4 and FIN5.

## 2b. Check your He giant FIN4 run.log (no rerun): FIN4 bug found by Raven (fin5_1014 NOTE)
- grep your HG FIN4 (and hegiant-large) fresh-start run.log for `ic_rad: no fused twin intensities`. With the fused twin and 2 blocks
  per rank the twin ragged array is > 2^31 doubles; rad_m1_icrad.cpp:208 `vgd_itw.d.extent_int(0) > 1` overflows and the fresh-start
  radiation IC skips the twin J correction (layout-dependent IC; Viper NT=2 vs NT=4 differ 5 % in hst at t=0). Restarts unaffected.
  FIN5 does not hit it at these sizes (compact band). Report whether the line is present; do NOT restart anything.
- In the He giant smoke of section 2 the FIN4 vs FIN5 pair is therefore EXPECTED to differ from t = 0 if that line is present in the
  FIN4 log (FIN5 = the 1-block-per-rank answer); then compare FIN5 vs FIN5-poison and FIN5 restart only, and report.

## 3. Report: docs/handover/NOTE-2026-10-14-caltech-fin5-build.md
md5s + HEADs, job id, the bitwise table, s/cycle + GB/GPU FIN4 vs FIN5, any non-bitwise line in full.
Expected (Raven A100 / Viper MI300A): all BITWISE; memory -30..-60 %; s/cycle within a few % (LI -7 %).

## Coordinator notes
- Priority: short (debug QoS). Do not delay running production.
