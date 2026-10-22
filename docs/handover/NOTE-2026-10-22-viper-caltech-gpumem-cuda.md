# NOTE 2026-10-22 viper -> caltech (USER): GPU memory on CUDA (A100) -- RSG E 3 MB/GPU fits 80 GB even keys-off; your OOM has another cause

Raven gpudev 31415728, 4 x A100 40 GB, RSG E N=64, 548x64x64 MB (two GPUs with 2 MB, two with 1 MB), 1976bac68 + harmonic + extrap2
+ accept_fallback + posfix=4, Kokkos CudaSpace log per rank, nvidia-smi every 2 s:
                                              keys off (1 MB / 2 MB per GPU)        keys on, gpumem-1022 3617c794
  fresh, steady after 1st dtr Save            17.26 / 33.99 GB (16.7 GB/MB + 0.5)   --
  fresh, peak                                 18.07 / 35.60 GB (17.5 GB/MB)         --
  restart read -> dtr Save                    17.02 / 33.52 GB                      14.37 / 28.23 GB (13.9 GB/MB)
  restart, peak (forced retry + rst write)    17.91 / 35.31 GB (17.4 GB/MB)         15.26 / 29.21 GB (13.95 GB/MB)
  nvidia-smi max on the 2-MB GPU (incl. context): 33.3 GiB keys off, 27.96 GiB keys on. No CUDA-specific restart-read transient.
=> per MB the same as MI300A. RSG E at 3 MB/GPU (6.75 M cells) on H100: ~51 GB steady / ~53 GB peak keys off, ~42.4 GB peak keys on
   (+0.6 GB context). It should NOT OOM on 80 GB.
Likely causes of your RSG E 8xH100 OOM: two ranks on one GPU (2 x 53 = 106 GB) -- check nvidia-smi per GPU / CUDA_VISIBLE_DEVICES per
rank / --gpu-bind; or a different deck (dt_retry_views = all, vet_gd_band_ring off, MeshBlocks not 548x64x64). Please send the items
in the 8c5b12d request.
HG (~7.5 kB/cell keys off, ~6.2 keys on): 6xH100 11.5 M cells/GPU ~86-91 GB off (OOM) / ~72 GB on (fits, tight);
12xH100 5.75 M ~45 off / ~37 on.
Keys (fork/gpumem-1022 = 3617c794, opt-in, default bitwise; MI300A GPU gates + Raven CPU gates PASS: key-off bitwise vs 1976bac68,
keys-on == keys-off bitwise on a retry/posfix stretch, C2==C, R==C):
  <time>/dt_retry_host = true (dtr checkpoint in pinned host memory; unchanged vra_d skipped, <time>/dt_retry_host_skip_clean)
  <rad_m1>/time2_posfix_host = true      <rad_m1>/restart_pack_per_mb = true
Cost MI300A +3.9 % s/cycle (dtr save 45 ms/cycle); A100 timing still queued. Not in fin51-r-fix yet.
