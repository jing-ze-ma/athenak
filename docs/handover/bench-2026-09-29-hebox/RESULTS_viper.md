# bench-2026-09-29-hebox: viper results (MPCDF viper, AMD MI300A APU)

Commit 401875f0 (fork/rt-integration), `builds/build_inc_viper.sh box_gpu72` (HIP GFX942_APU,
Release, MPI), binary `/viper/ptmp2/jinma/builds/bin/athena_box_gpu72_401875f0`,
md5 fe495394fa25b6ecc89f67c3e9384044. Stack: **gcc/16 rocm/7.2 (hipcc 7.2.4) openmpi_gpu/5.0**,
HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1. Fresh start, hesdirk2, cfl 0.3, 1200 cycles, window cycles
400-1200 (ana_bench_hebox.py), 1 MPI rank per GPU, 2 repeats per job (r1 / r2).
Jobs: 12018455 (1 GPU, apudev), 12018456 (2 GPUs, apudev), 12018457 (4 GPUs = 2 apu nodes x 2).
Logs: `/viper/ptmp2/jinma/bench_hebox_0929/timing/n<N>_r<r>_<job>/bench.log`. Every run rc 0,
FATAL 0, nan 0, NON-CONVERGED 0, fallbacks 0.

| machine | GPU | GPUs (ranks) | ms/cycle r1 / r2 | wall s / sim s r1 / r2 | dt | Picard/solve | Krylov it/linear solve | NC / fb | stack |
|---|---|---|---|---|---|---|---|---|---|
| viper | MI300A | 1 (1) | 111.78 / 111.86 | 0.694 / 0.694 | 0.16111 | 6.476 | 6.55 | 0 / 0 | gcc/16 rocm/7.2 openmpi_gpu/5.0 |
| viper | MI300A | 2 (2) | 93.17 / 93.80 | 0.578 / 0.582 | 0.16111 | 6.476 | 6.55 | 0 / 0 | gcc/16 rocm/7.2 openmpi_gpu/5.0 |
| viper | MI300A | 4 (4, 2 nodes) | 75.26 / 77.10 | 0.467 / 0.479 | 0.16111 | 6.476 | 6.55 | 0 / 0 | gcc/16 rocm/7.2 openmpi_gpu/5.0 |

- ms/cycle = median of the 10-cycle windows; the whole-window means are within 0.7 % of it.
- Picard/solve and Krylov iterations are whole-run (cycles 0-1200) from the `<rad_m1>` summary
  lines, identical at 1, 2 and 4 ranks (42.44 BiCGStab iterations per transport solve).
- Window check: cycles 100-1200 instead of 400-1200 give 112.96 / 94.35 / 75.73 ms (r1), within
  1.1 % of the 400-1200 values, so the start-up is short.
- Scaling per cycle 1 -> 2 -> 4 GPUs: 1.20x / 1.47x (r1). The box is small (4 blocks of
  84 x 52 x 52); strong scaling is limited.
- Wall per run incl. start-up (~8 s to the first cycle): 2.4 min at 1 GPU, 2.1 min at 2 GPUs,
  1.7 min at 4 GPUs.
- Smoke: job 12018403 (2 GPUs, 200 cycles, Picard log on), clean, see README.md.
- Reference (not comparable: evolved box restart, ROCm 6.3, whole-run average):
  hebox_cfl2_0927 R3 (cfl 0.3, 1 GPU) 92.7 ms/cycle, 5.06 Picard passes per solve.
