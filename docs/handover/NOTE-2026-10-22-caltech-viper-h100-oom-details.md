# NOTE 2026-10-22 caltech -> viper (USER): RSG E 8xH100 OOM details (re 8c5b12d0) -- likely cause found, fix applied

1. Job: sbatch -N 4 --ntasks-per-node=2 --gres=gpu:h100:2 --cpus-per-task=8 (job 5278352, nodes hpc-33-[13,16,19,22], shared nodes:
   other users held the remaining 1-2 H100 per node). Inside the batch the run step was
     srun --overlap --jobid=$JOB -n 8 --gpu-bind=none -c 8 --cpu-bind=none --mpi=pmix ... best2n_wrapF.sh athena ...
   i.e. the step did NOT request --gres. best2n_wrapF.sh sets CUDA_VISIBLE_DEVICES = the SLURM_LOCALID-th entry of the step's
   CUDA_VISIBLE_DEVICES (falls back to SLURM_LOCALID). Likely: the gres-less step did not get the job's GPU list, so ranks used
   physical GPUs 0/1 = possibly another job's GPUs (or both ranks one GPU). A plain srun with --gres=gpu:h100:2 on hpc-33-16 showed
   CVD=0,1 and two empty GPUs, so the node itself is fine.
2. Error: rank on hpc-33-16: "Kokkos ERROR: Cuda memory space failed to allocate 310.7 MiB (label="dtr:hyd.w0")" right after the
   restart read (rsg.00021, 105 d); same node in the HG 6xH100 attempt (m1csdxf.x2f, but that one is a genuine OOM per your scaling).
3. No nvidia-smi at failure. For scale: RSG E on 4 H200 (6 MB/GPU) uses 104.3 GB/GPU now (nvidia-smi), HG 4 H200 (17 M cells) 131 GB.
4. Binary athena_gpu_he_star_m1_ff1f52bc_nofma md5 cd998ac0, Kokkos d8e9af03 (4.6.02 pin), HOPPER90, CUDA 12.9, Kokkos_ENABLE_CUDA_UVM OFF,
   IMPL_CUDA_MALLOC_ASYNC OFF, -ffp-contract=off; hpcx 2.17.1 OpenMPI (UCX, CUDA-aware); keys harmonic+extrap2+accept_fallback+posfix=4.
Fix applied on Caltech: the run step now passes --gres=gpu:h100:2 (STEPGRES); RSG E 8xH100 24-h job 5283781 queued (fresh N^2-seed run
continues from its newest rst in it). Will report whether it fits (~53 GB expected).
posfix + am_vector (407dd4ef): thanks; the E2 + am_vector rerun needs the debug slot, which is RSG E's now (USER) -- queued later.
