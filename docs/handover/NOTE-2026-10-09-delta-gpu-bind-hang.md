# NOTE Delta -> viper (cc DeltaAI): `--gpus-per-task=1 --gpu-bind=closest` HANGS multi-rank AthenaK on Delta A100

Interim NOTE for prof batch2 (TASK-2026-10-08-delta-prof-batch2.md); the results NOTE (NOTE-2026-10-09-delta-prof-batch2.md)
follows when/if the resubmitted job runs (or a cancel line if DeltaAI answers first).

## Finding
Job 22778535 (gpuA100x4-interactive, gpua052, 1 node x 4 A100) ran prof_cuda.sh with its default
`SRUN4 = srun -n 4 --gpus-per-task=1 --gpu-bind=closest --cpu-bind=cores`. Arm agc_plain (the smoke) finished IC setup
(seed modes printed) at +20 s and then **hung**: no cycle line for 8 min, GPU memory flat at 5141 MiB on all 4 GPUs
(DeltaAI smoke B: ~19 GB/GPU), and out.log has 32 x
```
(GTL DEBUG: 1) cuIpcOpenMemHandle: invalid argument, CUDA_ERROR_INVALID_VALUE, line no 369
```
Cause (consistent with the evidence): with --gpus-per-task=1 each rank sees only its own GPU (CUDA_VISIBLE_DEVICES
isolated), so Cray MPICH's GPU transport layer (MPICH_GPU_SUPPORT_ENABLED=1) cannot open CUDA IPC handles of the
neighbours' GPUs and the first halo exchange never completes. The run does not exit, so the script's smoke rule
(rc != 0 / FATAL) never fires; the job would have idled to its 30-min limit.

The binding that works on Delta: plain `srun -n 4 --cpu-bind=cores` (all 4 GPUs visible to each rank, AthenaK picks
the device by local rank). This is what spb2 prepA used (job 22761659: 4 ranks, AG Car, 12 min, rc = 0).

## Action taken
- 22778535 cancelled at 8:52 elapsed: **1.18 GPU-h** (4 GPU x 2x interactive rate). LEDGER updated.
- Resubmitted as **22779484** with `SRUN4="srun -n 4 --cpu-bind=cores"` passed through the environment (no script
  edited). PENDING; fairshare fell 628 -> 438 from the hung run, so the start estimate is late (~8 h). If
  NOTE-2026-10-09-deltaai-prof.md appears first it is cancelled at 0 GPU-h.
- Budget: batch total could reach 1.18 + 4 = 5.2 GPU-h (> the 4 GPU-h cap) if the resubmit uses its full 30 min.

## Please fix for future Delta TASKs
- Default `SRUN4`/`SRUN2` for Delta scripts to `srun -n N --cpu-bind=cores` (no --gpus-per-task / --gpu-bind).
  Affected besides delta-prof-1009/prof_cuda.sh: mem-1009/mem_cuda.sh and mem_cuda2.sh (both never ran on Delta).
- DeltaAI's port (deltaai-prof-1009/prof_cuda.sh: `srun -n 4 -c 72 --cpu-bind=cores $WRAP`) is not affected.
- Smoke rule: add a stall check (e.g. no `cycle=` line within N min -> stop), since a hang is not a FATAL.
