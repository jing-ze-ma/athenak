# NOTE 2026-11-27zn orion -> Caltech: per-kernel roofline profile of a bsg cycle (1 H200)

Goal: find which kernels are far from the memory/compute roofline (memory-bound kernels that move too much data, uncoalesced access, register spills, low occupancy, divergence), to steer the ray-kernel rewrite and the VET-solve work. About 0.5 node-hour.

Binary: rel6 01deec542 (your build, -lineinfo if a rebuild is cheap; else as is). State: the bsg restart from 27zc, 1 GPU if it fits (else 4 GPUs, profile rank 0), 5 cycles after 2 warm-up cycles, outputs off, rad_vet/clean_verbose=1.
1. nsys: kernel summary (time share, calls, mean) for the top 25 kernels, plus total time in cudaMemcpy/host syncs/MPI per cycle (to see host round trips).
2. ncu on the top 10 kernels by time (one launch each, --set full or at least SpeedOfLight + MemoryWorkloadAnalysis + Occupancy + LaunchStats + SourceCounters): achieved DRAM bandwidth (% of peak) and FP64 throughput (% of peak), arithmetic intensity (FLOP/byte), L1/L2 hit rates, global load efficiency / sectors per request (coalescing), registers per thread, local-memory spills, achieved vs theoretical occupancy, warp execution efficiency (divergence).
One table (kernel, % of cycle, DRAM %, FP64 %, FLOP/B, sectors/req, regs, spills, occupancy, divergence) + a roofline plot of those kernels (mesa2.mplstyle + idort_style.py). Keep the .nsys-rep / .ncu-rep files in pages/ (not in git if large); push table + plot to caltech_results/roofline/.
