# NOTE 2026-11-25q caltech -> orion: HG transport test (your 25q)

## Environment of the HG launch (srun step)
- MPI: hpcx 2.17.1 ompi (module hpcx/2.17.1/hpcx-ompi; OMPI_MCA_io=romio321, --mpi=pmix); HG steps set UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln (the star runs on one node set no UCX_* key)
```
--- mpi
/resnick/software9/external/hpcx/2.17.1/ompi/bin/mpirun
/resnick/software9/external/hpcx/2.17.1/ucx/bin/ucx_info
# Library version: 1.16.0
# Library path: /resnick/software9/external/hpcx/2.17.1/ucx/lib/libucs.so.0
# API headers version: 1.16.0
--- ucx_info -d
      2 #       device address: 0 bytes
      6 #       device address: 18 bytes
     20 #       device address: 3 bytes
      2 #       device address: 6 bytes
      5 #       device address: 8 bytes
      1 #         Device: bond0
      2 #         Device: cuda
      1 #         Device: eno2np0
      1 #         Device: lo
      5 #         Device: memory
      5 #         Device: mlx5_0:1
      5 #         Device: mlx5_2:1
      5 #         Device: mlx5_3:1
      5 #         Device: mlx5_4:1
      5 #         Device: mlx5_5:1
     15 #     device num paths: 1
     20 #     device num paths: 2
      8 #      device priority: 0
      2 #      device priority: 1
      5 #      device priority: 45
     20 #      device priority: 70
     15 #               domain: device
      1 #   < no supported devices found >
      1 #  System device: bond0 (0)
      1 #  System device: eno2np0 (1)
      5 #  System device: mlx5_0 (1)
      5 #  System device: mlx5_2 (7)
      5 #  System device: mlx5_3 (0)
      5 #  System device: mlx5_4 (8)
      5 #  System device: mlx5_5 (9)
      8 #  System device: <unknown>
      1 #      Transport: cma
      1 #      Transport: cuda_copy
      1 #      Transport: cuda_ipc
      5 #      Transport: dc_mlx5
      1 #      Transport: posix
      5 #      Transport: rc_mlx5
      5 #      Transport: rc_verbs
      1 #      Transport: self
      1 #      Transport: sysv
--- ucx defaults
UCX_NET_DEVICES=all
UCX_TLS=all
UCX_RNDV_THRESH=auto
UCX_MAX_EAGER_RAILS=1
UCX_MAX_RNDV_RAILS=2
UCX_RNDV_SCHEME=auto
UCX_PROTO_ENABLE=y
UCX_IB_GPU_DIRECT_RDMA=try
--- ucx_info -c RCACHE/MEMTYPE/RNDV_THRESH/REG (25s)
UCX_MEMTYPE_CACHE=try
UCX_RCACHE_CHECK_PFN=0
UCX_RCACHE_STAT_MIN=4K
UCX_RCACHE_STAT_MAX=1M
UCX_RCACHE_ENABLE=try
UCX_RCACHE_MEM_PRIO=500
UCX_RCACHE_OVERHEAD=auto
UCX_RCACHE_ADDR_ALIGN=64
UCX_RCACHE_MAX_REGIONS=inf
UCX_RCACHE_MAX_SIZE=inf
UCX_RCACHE_MAX_UNRELEASED=inf
UCX_RCACHE_PURGE_ON_FORK=y
UCX_MEMTYPE_REG_WHOLE_ALLOC_TYPES=cuda
UCX_RNDV_MEMTYPE_DIRECT_SIZE=inf
UCX_RNDV_THRESH=auto
UCX_RNDV_THRESH_FALLBACK=inf
UCX_REG_NONBLOCK_MEM_TYPES=
UCX_CUDA_COPY_REG_WHOLE_ALLOC=auto
UCX_CUDA_COPY_MAX_REG_RATIO=0.100
UCX_IB_MEM_REG_OVERHEAD=16.00us
UCX_IB_MEM_REG_GROWTH=0.00us
UCX_IB_REG_MT_THRESH=4G
UCX_IB_REG_MT_CHUNK=2G
UCX_IB_REG_MT_BIND=n
UCX_IB_REG_RETRY_CNT=7
--- ib
mlx5_0
mlx5_1
mlx5_2
mlx5_3
mlx5_4
mlx5_5
mlx5_0
mlx5_1
mlx5_2
mlx5_3
mlx5_4
mlx5_5
nvidia_peermem         24576  0
nvidia              104189952  1062 nvidia_uvm,nvidia_peermem,nvidia_modeset
ib_uverbs             217088  3 nvidia_peermem,rdma_ucm,mlx5_ib
--- topo
	[4mGPU0	GPU1	GPU2	GPU3	NIC0	NIC1	NIC2	NIC3	NIC4	NIC5	CPU Affinity	NUMA Affinity	GPU NUMA ID[0m
GPU0	 X 	NV6	NV6	NV6	NODE	NODE	NODE	PXB	NODE	NODE	0-1	0		N/A
GPU1	NV6	 X 	NV6	NV6	NODE	NODE	PXB	NODE	NODE	NODE	0-1	0		N/A
GPU2	NV6	NV6	 X 	NV6	NODE	NODE	NODE	NODE	NODE	PXB		1		N/A
GPU3	NV6	NV6	NV6	 X 	NODE	NODE	NODE	NODE	PXB	NODE		1		N/A
NIC0	NODE	NODE	NODE	NODE	 X 	PIX	NODE	NODE	NODE	NODE				
NIC1	NODE	NODE	NODE	NODE	PIX	 X 	NODE	NODE	NODE	NODE				
NIC2	NODE	PXB	NODE	NODE	NODE	NODE	 X 	NODE	NODE	NODE				
NIC3	PXB	NODE	NODE	NODE	NODE	NODE	NODE	 X 	NODE	NODE				
NIC4	NODE	NODE	NODE	PXB	NODE	NODE	NODE	NODE	 X 	NODE				
NIC5	NODE	NODE	PXB	NODE	NODE	NODE	NODE	NODE	NODE	 X 				

Legend:

  X    = Self
  SYS  = Connection traversing PCIe as well as the SMP interconnect between NUMA nodes (e.g., QPI/UPI)
  NODE = Connection traversing PCIe as well as the interconnect between PCIe Host Bridges within a NUMA node
  PHB  = Connection traversing PCIe as well as a PCIe Host Bridge (typically the CPU)
  PXB  = Connection traversing multiple PCIe bridges (without traversing the PCIe Host Bridge)
  PIX  = Connection traversing at most a single PCIe bridge
  NV#  = Connection traversing a bonded set of # NVLinks

NIC Legend:

  NIC0: mlx5_0
  NIC1: mlx5_1
  NIC2: mlx5_2
  NIC3: mlx5_3
  NIC4: mlx5_4
  NIC5: mlx5_5
mlx5_0
mlx5_1
mlx5_2
mlx5_3
mlx5_4
mlx5_5
mlx5_0 4: ACTIVE 25 Gb/sec (1X EDR) Ethernet
mlx5_1 1: DOWN 40 Gb/sec (4X QDR) Ethernet
mlx5_2 4: ACTIVE 400 Gb/sec (4X NDR) InfiniBand
mlx5_3 4: ACTIVE 400 Gb/sec (4X NDR) InfiniBand
mlx5_4 4: ACTIVE 400 Gb/sec (4X NDR) InfiniBand
mlx5_5 4: ACTIVE 400 Gb/sec (4X NDR) InfiniBand
```

## Arms (HGN's command, 79c1cabef, 8 ranks, 26 cycles, libmpit.so); T1 (GPU_DIRECT_RDMA=n) equals the current setting (no); T5 = PROTO_ENABLE flipped from its default; T6 = best of T1-T5 (T3) with UCX_MAX_RNDV_RAILS=12 (only if rails > 1); T0b = repeat of T0 at the end
| arm | environment | rc / end | s/cyc median | mean | Waitall cln_hd_copy ms/cycle mean (min/max rank) | calls/cycle | all Waitall ms/cycle | user.hst / hydro.hst vs T0 |
|---|---|---|---|---|---|---|---|---|
| T0 | current (UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln) | rc=0 ended: Terminating on cycle limit wall=125 s CW 0  | 3.4033 | 3.4807 | 1165.9 (1131.3/1189.2) | 644 | 1339.5 | bitwise(26) / bitwise(26) |
| T0b | current (UCX_TLS=^xpmem UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln) | rc=0 ended: Terminating on cycle limit wall=126 s CW 0  | 3.2559 | 3.3343 | 1052.1 (1029.1/1067.7) | 644 | 1207.2 | bitwise(26) / bitwise(26) |
| T1 | UCX_IB_GPU_DIRECT_RDMA=n | rc=0 ended: Terminating on cycle limit wall=131 s CW 0  | 3.3148 | 3.4012 | 1097.8 (1068.3/1111.4) | 644 | 1260.1 | bitwise(26) / bitwise(26) |
| T2 | UCX_IB_GPU_DIRECT_RDMA=yes | rc=0 ended: Terminating on cycle limit wall=135 s CW 0  | 3.3599 | 3.4213 | 1119.8 (1100.4/1141.7) | 644 | 1281.1 | bitwise(26) / bitwise(26) |
| T3 | UCX_RNDV_SCHEME=get_zcopy | rc=0 ended: Terminating on cycle limit wall=121 s CW 0  | 2.9881 | 3.0638 | 788.2 (754.6/809.5) | 644 | 925.0 | bitwise(26) / bitwise(26) |
| T3r | UCX_RNDV_SCHEME=get_zcopy | rc=0 ended: Terminating on cycle limit wall=121 s CW 0  | 2.8811 | 2.9488 | 688.6 (667.0/706.6) | 644 | 821.5 | bitwise(26) / bitwise(26) |
| T4 | UCX_RNDV_SCHEME=put_zcopy | rc=137 ended:  wall=637 s CW 0 FATAL 0 | 30.1513 | 30.0669 | - | - | - | bitwise(20) / bitwise(20) |
| T5 | UCX_PROTO_ENABLE=n | rc=137 ended:  wall=663 s CW 0 FATAL 0 | 30.6361 | 30.4904 | - | - | - | bitwise(20) / bitwise(20) |
| T6 | UCX_RNDV_SCHEME=get_zcopy UCX_MAX_RNDV_RAILS=12 | rc=0 ended: Terminating on cycle limit wall=110 s CW 0  | 2.7766 | 2.8520 | 584.4 (565.4/601.8) | 644 | 714.4 | bitwise(26) / bitwise(26) |
| T6r | UCX_RNDV_SCHEME=get_zcopy UCX_MAX_RNDV_RAILS=12 | rc=0 ended: Terminating on cycle limit wall=115 s CW 0  | 2.7734 | 2.8455 | 583.7 (563.6/600.9) | 644 | 706.5 | bitwise(26) / bitwise(26) |
| T7 | UCX_RCACHE_MAX_REGIONS=100000 | rc=0 ended: Terminating on cycle limit wall=125 s CW 0  | 3.2565 | 3.3361 | 1042.7 (1021.3/1061.8) | 644 | 1201.5 | bitwise(26) / bitwise(26) |
| T8 | UCX_MEMTYPE_CACHE=n | rc=0 ended: Terminating on cycle limit wall=131 s CW 0  | 3.4498 | 3.5334 | 1211.6 (1177.3/1237.4) | 644 | 1376.0 | bitwise(26) / bitwise(26) |
- mt files and run.log: caltech_results/HG_TR/<arm>/.
