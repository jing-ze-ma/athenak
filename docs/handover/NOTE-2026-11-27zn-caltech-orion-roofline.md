
# NOTE 2026-11-27zn caltech -> orion: per-kernel roofline of a bsg cycle on H200 (your 27zn)

- Binary and state: rel6 01deec542 as is (no -lineinfo rebuild); bsg rst 00009 (cycle 29657), read-only hard link; outputs off, clean_verbose 1.
- Kokkos labels reach nsys and ncu through a 15-line NVTX Kokkos tool of mine (roofline_1127/kp_nvtx.cpp: one nvtxRangePush per parallel_for/reduce/scan, named by its label, no fencing).
- The .nsys-rep and .ncu-rep files are kept outside git in /resnick/groups/carnegie_poc/jingze/roofline_1127/ (ns2, ns7, ncu_<label>/).

## 1. nsys, 1 H200 (01-17, node otherwise idle)
- Method: two runs, 2 cycles and 7 cycles. ns7 minus ns2 = exactly 5 warm cycles.
- GPU kernel time is 2.055 s per cycle; the top-25 Kokkos labels and the host API / memcpy / MPI totals are below.
```
# GPU kernel time over 5 warm cycles (ns7 - ns2): 10.2744 s = 2.0549 s/cycle
rank | % | s/cycle | launches/cycle | mean us | Kokkos label
1 | 12.0 | 0.2468 | 250.0 | 987.3 | vet_vgd_shell
2 | 10.7 | 0.2205 | 26.4 | 8353.8 | cln_rhsaf
3 | 5.7 | 0.1166 | 12.4 | 9399.8 | cln_drag
4 | 5.6 | 0.1147 | 388.8 | 295.1 | cln_pcr2
5 | 4.5 | 0.0926 | 26.4 | 3506.3 | cln_fs0
6 | 3.9 | 0.0807 | 16.8 | 4802.9 | cln_facecoef
7 | 3.5 | 0.0709 | 78.4 | 904.9 | cln_apply_st
8 | 3.1 | 0.0630 | 1.0 | 63022.8 | cln_gd_fs
9 | 2.8 | 0.0569 | 31.6 | 1800.3 | cln_cslat
10 | 2.6 | 0.0529 | 129.6 | 408.5 | cls_restrict
11 | 2.5 | 0.0511 | 8.8 | 5806.8 | cln_rhs
12 | 2.4 | 0.0497 | 6.8 | 7315.1 | cln_velocity
13 | 1.8 | 0.0360 | 20.4 | 1763.1 | cln_fwc
14 | 1.6 | 0.0329 | 35.2 | 935.0 | cln_fsch
15 | 1.4 | 0.0283 | 26.4 | 1072.1 | cln_fwb
16 | 1.3 | 0.0266 | 6.8 | 3908.7 | cln_stencil
17 | 1.2 | 0.0255 | 20.4 | 1248.2 | cln_xdonor
18 | 1.2 | 0.0249 | 1.0 | 24915.7 | cln_gd
19 | 1.2 | 0.0239 | 2.0 | 11944.5 | cln_wbcell
20 | 1.0 | 0.0209 | 6.0 | 3482.1 | cln_wbxc
21 | 1.0 | 0.0205 | 52.8 | 387.9 | cln_fwbs
22 | 1.0 | 0.0199 | 150.0 | 132.6 | cln_scale
23 | 0.9 | 0.0182 | 26.4 | 689.6 | cln_fs1
24 | 0.9 | 0.0177 | 8.8 | 2009.6 | cln_newton
25 | 0.8 | 0.0172 | 26.4 | 651.8 | cln_fs3

# cuda_api_sum per cycle (5 warm cycles):
- cudaEventSynchronize: 0.8649 s/cycle, 835.4 calls/cycle
- cudaStreamSynchronize: 0.5562 s/cycle, 588.4 calls/cycle
- cudaDeviceSynchronize: 0.3989 s/cycle, 251.0 calls/cycle
- cudaMemcpyAsync: 0.2264 s/cycle, 356.4 calls/cycle
- cudaLaunchKernel: 0.0112 s/cycle, 3244.4 calls/cycle
- cudaMemcpyToSymbolAsync: 0.0033 s/cycle, 835.4 calls/cycle
- cudaEventRecord: 0.0023 s/cycle, 835.4 calls/cycle
- cudaFree: 0.0019 s/cycle, 28.8 calls/cycle

# cuda_gpu_mem_time_sum per cycle (5 warm cycles):
- [CUDA memcpy Device-to-Device]: 0.0037 s/cycle, 21.0 calls/cycle
- [CUDA memset]: 0.0019 s/cycle, 85.8 calls/cycle
- [CUDA memcpy Host-to-Device]: 0.0012 s/cycle, 1040.6 calls/cycle
- [CUDA memcpy Device-to-Host]: 0.0004 s/cycle, 175.2 calls/cycle

# mpi_event_sum per cycle (5 warm cycles):
- MPI_Allreduce: 0.0005 s/cycle, 0.0 calls/cycle
- MPI_File_close: 0.0003 s/cycle, 0.0 calls/cycle
- MPI_File_write: 0.0003 s/cycle, 0.0 calls/cycle
- MPI_Waitall: 0.0001 s/cycle, 0.0 calls/cycle
- MPI_Finalize: 0.0001 s/cycle, 0.0 calls/cycle
- MPI_Gather: 0.0000 s/cycle, 0.0 calls/cycle
- MPI_Bcast: 0.0000 s/cycle, 0.0 calls/cycle
- MPI_Alltoall: 0.0000 s/cycle, 0.0 calls/cycle
```
- Host round trips per warm cycle:
  - cudaEventSynchronize 0.86 s (835 calls);
  - cudaStreamSynchronize 0.56 s (588);
  - cudaDeviceSynchronize 0.40 s (251);
  - cudaMemcpyAsync 0.23 s (356);
  - H2D copies 1041 per cycle, but only 1.2 ms of copy time.
  - MPI is negligible on 1 rank.
  - So the host waits about as long as the GPU computes.

## 2. ncu --set full on the top 10 labels, one warm launch each
- The 1-GPU ncu replay ran out of device memory: bsg holds about 70 GB on one H200, and replay needs save/restore space.
- So ncu ran on 4 GPUs (the production layout), profiling rank 0 only (1 of the 4 MeshBlocks). The launch was skipped past 2 warm cycles; node 02-17 had only these runs.
- The % of cycle column is from the 1-GPU nsys.
| Kokkos label | % of cycle (nsys, 1 GPU) | dur us (1 launch, rank 0 of 4) | DRAM % peak | FP64 pipe % | FLOP/B (DRAM) | GFLOP/s | L1 / L2 hit % | sectors/req (global ld) | regs/thr | local ld+st B (spills) | occupancy achieved / theor. % | active threads/warp % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cln_apply_st | 3.5 | 247.7 | 52.8 | 4.7 | 0.220 | 570.9 | 28 / 30 | 9.0 | 44 | 0 | 59 / 62 | 100 |
| cln_cslat | 2.8 | 550.9 | 21.3 | 11.6 | 1.392 | 1458.9 | 30 / 72 | 7.4 | 109 | 0 | 24 / 25 | 100 |
| cln_drag | 5.7 | 4055.4 | 14.0 | 8.3 | 1.803 | 1241.1 | 67 / 73 | 2.9 | 80 | 6.38e+09 | 37 / 38 | 98 |
| cln_facecoef | 3.9 | 1193.0 | 65.3 | 7.8 | 0.404 | 1294.0 | 27 / 88 | 7.3 | 56 | 3.12e+09 | 52 / 56 | 100 |
| cln_fs0 | 4.5 | 928.9 | 18.1 | 6.6 | 1.009 | 896.5 | 60 / 56 | 2.7 | 80 | 3.38e+08 | 37 / 38 | 100 |
| cln_gd_fs | 3.1 | 20749.6 | 11.4 | 10.1 | 2.226 | 1247.6 | 57 / 79 | 5.4 | 80 | 5.53e+10 | 37 / 38 | 100 |
| cln_pcr2 | 5.6 | 68.4 | 11.7 | 13.9 | 4.216 | 2429.8 | 20 / 40 | 7.5 | 56 | 1.68e+07 | 47 / 50 | 99 |
| cln_rhsaf | 10.7 | 1462.7 | 23.9 | 5.5 | 0.575 | 675.5 | 58 / 71 | 2.5 | 80 | 2.72e+09 | 36 / 38 | 94 |
| cls_restrict | 2.6 | 60.8 | 26.2 | 2.5 | 0.212 | 273.5 | 27 / 16 | 8.1 | 72 | 0 | 34 / 44 | 100 |
| vet_vgd_shell | 12.0 | 278.3 | 4.1 | 16.1 | 10.524 | 2100.1 | 69 / 79 | 3.1 | 98 | 4.5e+06 | 22 / 25 | 82 |
- Columns:
  - FLOP counts double add + mul + 2 x fma (sass, predicated on).
  - DRAM bytes are read + write.
  - Spills are local-memory sectors x 32 B.
  - Active threads per warp = smsp__thread_inst_executed_per_inst_executed / 32.
  - Peaks: HBM3e 4.8 TB/s, FP64 33.5 TFLOP/s (non-tensor; data sheet).
- Reading:
  - Every kernel is below the bandwidth roof. The ridge is about 7 FLOP/B; only vet_vgd_shell (10.5) is right of it.
  - Most kernels sit at 1-2.4 TFLOP/s with DRAM at 4-26 % of peak. That points to latency / occupancy, not bandwidth, as the limit.
  - The 80-register kernels (cln_rhsaf, cln_drag, cln_fs0, cln_gd_fs) run at 37 % theoretical occupancy and spill to local memory. cln_gd_fs moves 55 GB of local traffic in one 20.7 ms launch; cln_drag 6.4 GB; cln_facecoef 3.1 GB; cln_rhsaf 2.7 GB.
  - cln_apply_st, cln_facecoef and cls_restrict are bandwidth-leaning: DRAM 53 / 65 / 26 %, at 7-9 sectors per request, so their loads are not coalesced (4 is ideal for 8-byte loads).
  - vet_vgd_shell: 98 registers, 22 % achieved occupancy, 82 % active threads per warp (divergence).
- Plot: docs/handover/caltech_results/roofline/roofline_bsg_1h200.png (mesa2 + idort_style). Table: docs/handover/caltech_results/roofline/ncu_table.md. nsys list: docs/handover/caltech_results/roofline/nsys_top25_and_api.txt.
