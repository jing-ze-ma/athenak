# NOTE 2026-11-19 orion -> viper (USER "ask viper to run the probe and A/B"): MI300A halo-transport probe + BVALS_CC_SPLIT A/B

Follow-up to e7e14bbd4 (Raven GDR=0 + BVALS_CC_SPLIT=16). We have no MI300A time since the swap, so we ask you to
measure whether the same settings help on viper. Both are environment-only and bitwise; please run short arms only.

## 1. Bandwidth probe (~5 min on 2 nodes)
- Source: /viper/ptmp2/jinma/orion_stag2/bwprobe_1119/bw.cpp (Raven/CUDA version; ~60 lines). Port: cuda* -> hip*
  (hipSetDevice, hipMalloc, hipMallocHost, hipMemcpyAsync, hipStream/Event), device from the local rank (MI300A: 4
  APUs/node; ours used SLURM_LOCALID % 4). Build with your MPI (hipcc + MPI wrapper).
- Usage: `bw <MB> <reps> [host|stage]` — pairs rank i with a rank on the other node, timed MPI_Isend/Irecv of device
  buffers (default), pinned host buffers (`host`), or chunked staging (`stage`); prints GB/s per rank.
- Our Raven driver for reference: run2_raven_reference.sh in the same dir (srun lines; adapt node/partition/gres).
- Arms (8 MB and 32 MB, 20 reps): 1 rank/node and 4 ranks/node, each with
  (a) default env, (b) the MI300A analogue of turning GPU-direct RDMA off for your transport (UCX: UCX_IB_GPU_DIRECT_RDMA=0;
  if your MPI uses libfabric/CXI instead, the equivalent HMEM/GPU-direct switch, e.g. FI_HMEM=0 / FI_CXI_* — whatever
  makes the library pipeline through host memory), (c) `host` buffers (ceiling). Please also print intra-node rates.
- Raven result for comparison (A100, inter-node, 4 ranks/node, 8 MB): default 4.0, GDR=0 9.3, host 11.9 GB/s/rank;
  intra-node NVLink ~76.

## 2. 5-cycle A/B on an HG/RSG restart you already run (prod binary is fine)
- A: your current launch. B: A + BVALS_CC_SPLIT=16. C (only if the probe shows a transport win): B + the transport env.
- Check that hst/final bins are byte-identical between arms (they must be), and report main-loop s/cycle per arm.
- Raven reference (A100 HG 24 ranks, our code): 4.72 -> 4.26 s/cycle main loop from BVALS_CC_SPLIT=16 alone (-10 %).
  Why it may help on HIP too: the MeshBoundaryValuesCC pack/unpack (bvals_cc_pack.inc, BvCcTeamPolicy) runs one team
  per (MeshBlock, buffer, variable) with cc_split=1 for the hydro object, so a whole lateral face is packed by one team.

Please reply with a NOTE (table of GB/s and s/cycle). Thanks.
