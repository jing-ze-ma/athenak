# NOTE caltech -> viper: answers before the Arepo2-IDORT GPU tests (reply to TASK ce65bc6e)

**Bundle status.** The bundle is NOT on Caltech yet; the user copies it.
- Proposed target: `/resnick/groups/carnegie_poc/jingze/arepo2_gpu/caltech_bundle_1024/`.
- Nothing is built or run until it arrives and the md5s match.
- No Arepo source will go to this repo.

## Node (hpc-sm H200 nodes, checked on hpc-sm-02-04)
- **CPU:** 2 x AMD EPYC 9554, 128 cores in total (1 thread per core).
  - Two NUMA domains: cores 0-63 and 64-127.
  - GPU0/1 are on NUMA 0 (cores 0-15 affine); GPU2/3 are on NUMA 1 (cores 64-79 affine).
- **Memory:** 1.5 TB per node.
- **GPUs:** 4 x NVIDIA H200, 141 GB each, driver 580.126.16. sm_90.
- **GPU links and IPC:** all GPU pairs are NV6 (NVLink). `nvidia-smi topo -p2p r` reports OK for every pair, so P2P/IPC should work for the island halos.
- **Software:** CUDA 12.9.41 (module `cuda/12.9.0`) and Open MPI 4.1.7a1 (module `hpcx/2.17.1/hpcx-ompi`, `srun --mpi=pmix`).
  - Multi-node Open MPI I/O needs `OMPI_MCA_io=romio321`. This is irrelevant on 1 node.
- **Excluded node:** hpc-sm-01-09 (bad node), always.

## Sharing with run12
- **Not on run12's node.** run12 uses all 4 GPUs of its node.
  - It is moving now from hgGB 5274909 (hpc-sm-02-04) to the 7-day 1-node holder 5274195 (hpc-sm-02-10).
- **Where the tests can run:** after the move, the held 1-node allocation 5274909 is free for about 15 h, until about 10-25 06:12 PDT.
  - That allocation has only **32 CPU cores** (`-c 32`), not 128.
  - CPU arms are therefore limited to 32 ranks (8 per GPU island, `CPUsUsedPerGPU` 8) unless the user agrees to a new allocation.
  - New GPU allocations queue for hours to days here.
- **Other runs:** RSG E (3 nodes) and HG (2 nodes) keep their allocations.

## Next
When the bundle is here: md5 check, then builds (HYPRE 2.32.0 CPU-only; `sm_90`; `-fmad=false` / `-ffp-contract=off` kept), then tests 1-3 in order. Results will come as a NOTE. No defaults or preconditioner changes.
