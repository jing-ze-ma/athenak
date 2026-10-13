# TASK raven -> DeltaAI: FIN3 final validation, 4 long Cartesian IDORT runs on GH200 (1 node, 4 GPUs, <= 1 h)

From Raven integrator IN (fin3-1012), via the coordinator. These runs race Caltech (CPU TASK) and Viper; the first complete result counts.

## 1. Build (report md5)
- Fork branch **`fin3-1012` @ `281dd686`**.
- PLAIN build (built-in pgens, NO `PROBLEM`) for GH200, as in your other builds:
  `-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off`
  (`build_inc_deltaai.sh`, with the plain / default-pgens target).
- `export MPICH_GPU_SUPPORT_ENABLED=1`.

## 1b. Partition
Use the regular ghx4 BATCH partition (1 node, <= 1 h), NOT ghx4-interactive, so that this job does not block your BSG link chain
(one interactive job per user). Do not touch the BSG run, its lock or its chain.

## 2. Files (docs/handover/fin3-1012/deltaai_cartlong/ on bsg-files-1009)
- `runs/{ddf,ddiff,ddn}/FIN3_d/` and `runs/pipe/FIN3_c220/`: `in.athinput` + `args.txt` (the command-line keys).
- `run_gpu.sh <binary>`: one run per GPU (CUDA_VISIBLE_DEVICES 0-3), each a single rank, in parallel.
- Wall limit 1 h. On Raven CPU (4 threads) these took 45-120 min; on a GPU they should take minutes.
  - If a run does not finish, report its last `cycle=` line.

## 3. Report: NOTE-2026-10-13-deltaai-fin3-cartlong.md
- STARTED: binary md5, job id.
- Per run: rc, wall, the last `cycle=` line, and `grep -c FATAL run.log`.
- Then `tar czf cartlong_gpu.tgz runs` (bin/, tab/, run.log). Post the path, or scp it to Raven /raven/ptmp/jinma/fin3_1012/deltaai/.
