# NOTE 2026-10-08 (Caltech): He giant port verified on H200; N445 fresh run to the remap point queued

Reply to NOTE-2026-10-08-viper-hegiant-files-caltech and NOTE-2026-10-08-viper-hegiant-remap-caltech. User 10-07: Caltech runs
the He giant itself (N445 fresh -> 10.42 d -> remap N897 with remap/he_remap_giant.py -> 30 d); the bundle could not be copied.
Leave viper hegiant897 / DeltaAI alone; the user decides which copies continue.

- **Binary:** f3a66907 (hegiant-opn-1007), PROBLEM=he_star_m1, CUDA 12.9, gcc 13.2, hpcx 2.17.1 OpenMPI, Kokkos HOPPER90, nvcc default
  FMA; `athena_gpu_he_star_m1_f3a66907` md5 f5d85053fbc8fa080e7ff0feeae949bb.
- **Fresh-start smoke** (job 4212175, 1 H200, 1 rank, 4 blocks, 60 cycles, input = viper's with the 3 paths changed): rc 0, FATAL 0,
  NaN 0, Picard 60 / mean 4.95 / max 10 / NON-CONV 0, he_ic_balance 0.00369707 at 2.15098e11, dt 20.04690 -> 20.04713,
  **0.955 s/cycle** (median cycles 10-60). Versus viper ref 12123795 (viper's comparison): last hst t 1.0e-14 rel, mass 5.3e-14 rel,
  Picard and he_ic_balance identical in the printed digits. **PORT PASS.** 1 H200 (0.955 s/cyc) ~ viper scout128 on 4 MI300A (0.862).
- **Production (queued, Priority, est. start 10-09):** N445 fresh to `time/tlim=9.004208e5` (XKEYS `time/cfl_number=0.3`), two
  copies racing, the first to start cancels the other: 2 H200 chain 4212559-62 (4 h links) and 1 H200 job 4212784 (24 h).
  Run dirs /resnick/groups/carnegie_poc/jingze/hegiant_1007/{n445,n445_1g}. The final rst at tlim is the remap source.
- Remap files staged (md5 equal to MD5SUMS): he_remap_giant.py, phi_code.py, grid/p2_897.npy. Remap on a CPU compute node
  (~11 GB RSS; login nodes cap 8 GB/proc), then a 60-cycle N897 smoke with the B''' keys before any N897 production.
