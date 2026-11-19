# NOTE 2026-11-19 orion -> viper (USER "send viper the note"): Raven multi-node A100 runs: set UCX_IB_GPU_DIRECT_RDMA=0 (+ BVALS_CC_SPLIT=16)

Measured on Raven pHG262 (A100 40 GB, 4 GPUs/node, CUDA-aware OpenMPI/UCX module), staggered clean code, same restarts:

**MPI bandwidth probe** (inter-node pairs, 4 ranks/node sending at once, 8 MB device buffers):
| setting | GB/s per rank |
|---|---|
| default (GPUDirect RDMA on; module sets it to 1) | 4.0 |
| `UCX_IB_GPU_DIRECT_RDMA=0` (UCX pipelines through host itself) | 9.3 |
| pinned host buffers (ceiling) | 11.9 |
Intra-node (NVLink cuda_ipc) is fine at ~76 GB/s. `UCX_PROTO_ENABLE=y`, rndv get/put_zcopy, MAX_RNDV_RAILS=1: no effect.

**Full runs** (s/cycle, 5 cycles; hst and final bins byte-identical, restart check bitwise):
| run | default | GDR=0 | GDR=0 + BVALS_CC_SPLIT=16 |
|---|---|---|---|
| HG 24 ranks / 6 nodes | 4.40 | 4.02 | 3.39 |
| RSG 12 ranks / 3 nodes | 6.28-7.31 | 5.61-6.37 | 5.44-6.09 |

- `UCX_IB_GPU_DIRECT_RDMA=0` is environment only: no code or key change, bitwise by construction.
- `BVALS_CC_SPLIT=16` is an existing env knob (more teams per buffer in the MeshBoundaryValuesCC pack/unpack); copies only, bitwise. The cs pack kernel (BvCcTeamPolicy) also lacks CUDA launch bounds (255 regs, 12 % occupancy); our commit adds them under KOKKOS_ENABLE_CUDA.
- Old prod's HG/RSG binaries pay the same transport cost on their (smaller) halo, so your Raven production runs should also gain (our estimate for prod HG: -0.15 to -0.2 s/cycle from GDR=0 alone; BVALS_CC_SPLIT=16 helps the hydro pack too). Worth a short A/B on your side before switching.
