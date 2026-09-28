# bench-2026-09-28: viper results (MPCDF viper, AMD MI300A APU)

Commit 11c9a5be (git archive build, HIP GFX942_APU, gcc/14 rocm/6.3 openmpi_gpu/5.0, Release, MPI),
binary md5 b8e3fff04cc117252ec9e9a97a0e1e2c. Fresh start, 2000 cycles, timing window cycles 1000-2000,
median of 8-cycle windows (ana_bench.py). 1 MPI rank per GPU; 4 GPUs = 2 nodes x 2 GPUs.
HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1. Jobs 12016722 (1 GPU), 12016723 (2 GPUs), 12016724 (4 GPUs);
2 interleaved repeats each; FATAL 0 and NOT-CONVERGED 0 in every run.

| machine | GPU | GPUs (ranks) | 1x ms/cycle | 1x wall s / sim s | 10x ms/cycle | 10x wall s / sim s |
|---|---|---|---|---|---|---|
| viper | MI300A | 1 (1) | 25.91 / 26.00 | 1.84e-3 / 1.85e-3 | 34.17 / 34.10 | 2.83e-3 / 2.82e-3 |
| viper | MI300A | 2 (2) | 19.19 / 19.25 | 1.37e-3 / 1.38e-3 | 25.13 / 25.05 | 2.09e-3 / 2.09e-3 |
| viper | MI300A | 4 (4, 2 nodes) | 12.71 / 12.84 | 0.91e-3 / 0.92e-3 | 16.64 / 16.64 | 1.39e-3 / 1.38e-3 |

Scaling 1 -> 2 -> 4 GPUs: 1x 1.35x / 2.03x, 10x 1.36x / 2.05x (per cycle).
Window dt: 1x 13.892 s at every rank count; 10x 12.076 / 11.991 / 11.974 s at 1 / 2 / 4 ranks (the 10x dt
depends slightly on the rank count; compare wall per simulated second across machines).
One rotation (1.101535e5 s) at 4 GPUs: 1x ~100 s wall, 10x ~153 s wall (fresh-start dt; the production dt
after spin-up is smaller).
