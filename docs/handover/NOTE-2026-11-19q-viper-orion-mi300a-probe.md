# NOTE 2026-11-19 viper -> orion: MI300A halo-transport probe + BVALS_CC_SPLIT A/B (reply to 91fcab9eb)

Run on viper pHG 13255125, vipa1073/1081/1186 (MI300A, 2 APUs per node, not 4). Files:
/viper/ptmp2/jinma/ccfull_1117/bwprobe_1119 (bw.cpp HIP port: device = SLURM_LOCALID % 2; run_probe.sh, probe.out,
run_ab.sh, run_abc2.sh, ab.out); run dirs /viper/ptmp2/jinma/ccfull_1117/runs/ab_*.

## Transport

UCX 1.18.0 under Open MPI 5.0.10 (pml ucx), NICs mlx5_0 / mlx5_1; transports rc_mlx5, dc_mlx5, ud_mlx5, rc_verbs,
posix / sysv / cma, tcp. Our production env XENV = `UCX_RNDV_THRESH=inf UCX_TLS=^rocm_ipc` (UCX_PROTO_INFO, 8 MB):
inter-node "multi-frag eager zero-copy copy-out" on rc_mlx5/mlx5_0:1, intra-node "multi-frag eager copy-in copy-out":
it forces the eager copy path everywhere.

## 1. Bandwidth probe (device buffers unless noted; mean GB/s per rank, 20 reps)

| env | 1 rank/node 8 MB | 1 rank/node 32 MB | 2 ranks/node 8 MB | 2 ranks/node 32 MB |
|---|---|---|---|---|
| (a) XENV (production) | 8.8 | 4.8 | 8.9 | 5.0 |
| (b) XENV + UCX_IB_GPU_DIRECT_RDMA=0 | 1.0 | 1.0 | 1.0 | 1.0 |
| (c) XENV, pinned host buffers | 9.7 | 7.7 | 9.7 | 7.6 |
| default UCX env (rndv + rocm_ipc) | 33.7 | 44.1 | 22.3 | 24.5 |

Intra-node (2 ranks, one node): XENV device 8.9 / 5.0, default device 53.4 / 79.4, XENV host 9.7 / 7.7 GB/s (8 / 32 MB).
On MI300A, GDR=0 is the wrong direction (about 1 GB/s); the win is dropping XENV (4-9x inter-node, 6-16x intra-node).

## 2. A/B with the production RSG binary (athena_he_gpu72_q_864d120d_p1109) and keys, from rsg.01658, 3 nodes x 2 ranks

Main loop s/cycle excludes the first cycle.

| arm | env | cycles | s/cycle | hst + final bins |
|---|---|---|---|---|
| A | XENV | 5 | 6.44 | reference (5 cycles) |
| B | XENV + BVALS_CC_SPLIT=16 | 5 | 6.81 | identical to A |
| B20 | XENV + BVALS_CC_SPLIT=16 | 20 | 4.70 | reference (20 cycles) |
| C1 | default UCX + BVALS_CC_SPLIT=16 | 20 | 3.84 | identical to B20 |
| C2 | default UCX + BVALS_CC_SPLIT=16 | 20 | 3.83 | identical to B20 |

- BVALS_CC_SPLIT=16 gives no win here (B 6.81 vs A 6.44 over 5 cycles, within the noise of a short arm); bitwise.
- The default UCX transport cuts the main loop by 18 % (4.70 -> 3.84 s/cycle) and was bitwise identical in both
  repeats (hst and all final bins, md5).

## Caveat (why XENV is there)

XENV was added around 10-21 for run-to-run determinism: with the default rendezvous / rocm_ipc path, Caltech and
viper saw non-bitwise results and restart divergence (rsg-rst-fix-1021), possibly a completion/ordering race on device
buffers, i.e. a correctness issue and not only reproducibility. Two 20-cycle repeats that agree are not enough to drop
it: we do not recommend changing the production env until longer and repeated runs (and a restart continuation) are
bitwise as well.
