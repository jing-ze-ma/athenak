# NOTE 2026-10-22 viper -> caltech (USER): GPU memory interim (re TASK 60f93337) + request: your RSG E H100 OOM details

Measured (viper MI300A, Kokkos-tools allocation log, RSG E det deck at N=64, 548x64x64 MB, 3 MB/GPU = 6.73 M cells/GPU, 2 ranks,
keys harmonic + extrap2 + accept_fallback + posfix=4, forced dt_retry):
                                  after restart read   after 1st dtr Save   peak (retry + rst write)
  old 4de8168d                    45.0 GB              48.6 GB              49.7 GB
  new 1976bac68 + 4 keys          46.2 GB              50.0 GB              52.7 GB
  -> the fix adds ~0.5 GB/MB steady, ~1.0 GB/MB peak; ~7.5-7.9 kB/cell overall.
Per-cell check of your numbers: HG 4xH200 17 M cells/GPU -> ~130 GB predicted (you: 131, consistent = cells/GPU, not a regression);
HG 6xH100 11.5 M -> ~86-91 GB (genuine OOM); HG 12xH100 5.75 M -> ~45 GB (fits); RSG E 8xH100 6.75 M -> ~53 GB: SHOULD FIT.
So the RSG E H100 OOM (> 80 GB at the first dtr Save after the restart) is NOT reproduced on MI300A: CUDA-specific or configuration.
Largest items per MB: implicit work array iw (92 comps) 2.0 GB, other ImplicitInit 2.0, dtr checkpoint 2.34, vet_gd ragged 2.23,
Time2Init 1.28, RotAvgSnap 1.15.
Cuts (branch fork/gpumem-1022 = 8969bc12, opt-in, default bitwise; gates in progress, not for production yet):
  <time>/dt_retry_host = true          ~2.3 GB/MB  (dtr checkpoint in pinned host memory)
  <rad_m1>/time2_posfix_host = true    ~0.1-0.2 GB/MB
  <rad_m1>/restart_pack_per_mb = true  ~0.8 GB/MB off the restart-write transient
Next on viper: CUDA measurement on Raven A100 (1 MB/GPU: restart read -> dtr Save -> retry -> rst write, nvidia-smi sampling).

REQUEST (to find the OOM): please push for the failed RSG E 8xH100 run:
  1. the job script + srun/mpirun line (ranks per node, --gpus-per-task / --gpu-bind, CUDA_VISIBLE_DEVICES handling);
  2. the full stdout/stderr around the OOM (Kokkos label, requested size, rank);
  3. nvidia-smi (or the Kokkos memory report) at/near the failure, per GPU: how many processes per GPU and their memory;
  4. binary md5 + cmake cache (Kokkos version, CUDA arch, Kokkos_ENABLE_CUDA_UVM / IMPL_CUDA_MALLOC_ASYNC settings), MPI flavour
     (CUDA-aware?), and the deck/keys + restart file used.
Two ranks on one GPU, or a CUDA malloc-async / UVM pool setting, would explain ~2x.
