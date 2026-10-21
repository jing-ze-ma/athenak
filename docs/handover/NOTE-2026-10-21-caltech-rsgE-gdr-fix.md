# NOTE 2026-10-21 caltech -> raven: RSG E 2-node crash = GPUDirect RDMA in VetGdCsHalo; fixed by staging through host; chain requeued

Test 5223799 (2 x 4 H200, hpc-sm-01-10 + 02-03, prod deck a7e734c4 + 9a56b4a1, nlim 20):
- B UCX_IB_GPU_DIRECT_RDMA=no: rc 0, 20 cycles, clean (216 s incl. ic_rad).
- C UCX_RNDV_SCHEME=put_zcopy: abort after 38 s, "Local protection error on mlx5_3:1/IB".
- A default (GPUDirect): abort, "Transport retry count exceeded on mlx5_3:1/IB" (same as 5126654 on 02-[04,11]; earlier test arm A1: "failed to arm iface").
All failures in ic_rad -> VetGdCsHalo MPI_Waitall. Earlier 2-node decks without the cubed-sphere vet_gd exchange ran fine with GPUDirect on.
So on Caltech this exchange fails when sent straight from GPU memory; staging through the host works.

- Production: link.sh now exports UCX_IB_GPU_DIRECT_RDMA=no (all MPI goes through the host, some slowdown). New chain 5223807 -> 5223808 -> 5223809
  (2 x 4 H200, 7 d; 5126655 refused at start because of a log-name guard on our side, 5126656 cancelled). Binary and deck unchanged.
- Suggestion: stage only the VetGdCsHalo buffers through the host in code (or pinned host buffers), so the rest of the MPI traffic can use GPUDirect again.
