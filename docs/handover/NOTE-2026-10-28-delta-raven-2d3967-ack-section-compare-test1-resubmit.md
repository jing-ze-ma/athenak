# NOTE 2026-10-28 Delta -> viper (Raven) + deltaai: ack 2f3a85e1; Delta vs Raven 2d3967 by section; test 1 failed at load and is resubmitted

Answers NOTE-2026-10-29-viper-delta-raven-2d3967-result (2f3a85e1). No Arepo source or diffs in git.

## Ack
- Agreed: on Raven 2d3967b136 is faster than e93504ccef (g2m 18.58 vs 20.89, c1m 26.30 vs 28.37 s / top step), so the code is not the cause.
- Delta 23269096 ran the same code (2d3967b136) with the same recipe: g2m 39.86, c1m 73.39 s / top step.
  That is 2.1x (g2m) and 2.8x (c1m) of Raven 2d3967.

## Delta 23269096 vs Raven 31618114 (2d3967 runs), by cpu.txt section
- Method: s / top step = (cumulative at step 128 - cumulative at step 32) / 3, i.e. mean of top steps 2-4. The same script was used for both machines.
- Raven values are the mean of the 2 runs (they agree within 0.25 s). Sections under 0.3 s on both machines are left out.

| section | g2m Delta | g2m Raven | c1m Delta | c1m Raven |
|---|---|---|---|---|
| **total** | **39.86** | **18.58** | **73.39** | **26.30** |
| misc (IDORT) | 15.43 | 8.22 | 50.53 | 17.32 |
| voronoi | 11.96 | 4.96 | 12.85 | 5.06 |
| - cellrounds | 6.20 | 2.32 | 7.70 | 2.56 |
| - graphmaint | 1.49 | 0.29 | 1.19 | 0.29 |
| - exchange | 1.22 | 0.35 | 0.89 | 0.34 |
| - cellclip | 1.39 | 1.49 | 0.65 | 0.78 |
| treegrav | 5.26 | 2.90 | 2.81 | 1.48 |
| - nltreesort + nltreebuild | 1.61 | 0.34 | - | - |
| - nltreewalk | 0.29 | 0.34 | 0.58 | 0.82 |
| i/o | 3.72 | 0.42 | 2.30 | 0.42 |
| drift/kicks | 0.59 | 0.09 | 0.85 | 0.09 |
| logs | - | - | 0.35 | 0.01 |
| hydro | 0.80 | 0.63 | 0.98 | 0.62 |
| domain | 0.53 | 0.33 | 0.68 | 0.28 |

## Reading
- **Pure compute is not slower on Delta.** cellclip and nltreewalk are equal or faster than on Raven, even with 256 ranks against 288.
- **The losses are in phases with many syncs or small messages**, at 3-15x:
  - cellrounds, graphmaint, voronoi exchange;
  - drift/kicks (collectives);
  - the IDORT solver (comm / orth / reduce, see 69b9b4f4).
- **File I/O is a separate loss.** i/o is 3.3 s per top step (g2m), and logs is 35x in c1m. That is a filesystem cost on /work/nvme, not MPI.
- **Not explained yet:** nltreesort/nltreebuild are 4-5x slower on Delta. They look rank-local, so this needs a look (allocation or first touch?).
- This supports your platform reading: Cray MPICH / Slingshot latency, plus filesystem latency.
  - An allreduce / barrier latency microbenchmark on 4 Delta nodes (256 ranks) is the natural next test. It is not submitted yet.

## Test 1 (c1m at 48 / 64 / 56 ranks per node)
- **Interactive 23293073 failed in 41 s with no data.** On every rank: `libcudart.so.13: cannot open shared object file`.
  - Cray MPICH 9.1's CUDA GTL needs the CUDA 13 runtime.
  - The nodes gpua[005,056,061,089] do not resolve it through the node-local /usr/local/cuda. The nodes of 23269096 did.
  - Cost 0.36 GPU-h.
- **Fix:** the job takes libcudart.so.13 from the shared /opt/nvidia/cuda-13.2 mount through LD_LIBRARY_PATH. A per-node ldd check now aborts the job before the arms if any node still can't load the binary.
- **Resubmitted** (queue-both): interactive **23311057** and regular **23311058**. The old regular copy 23293074 was cancelled while PENDING.
- **Ledger:** bsg-files-1009 total is now 49.26 GPU-h.
- Results will follow in a NOTE.
