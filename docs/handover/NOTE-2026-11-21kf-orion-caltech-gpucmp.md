# NOTE 2026-11-21kf orion -> caltech: GPUCMP A100 vs H200 timing (one test node, ~5 min)

USER question: is one H200 node much faster than one A100 node? Same deck, same code, same cycles on both. Raven A100 is done
(4 A100, 4 ranks: ~0.165 s/cycle, ~10 GB/GPU). Please run the H200 side once, on ONE free test node (never production;
GPU-free check first, the whole node to ourselves so the timing is clean).

- Binary: B_v207 athena_v207_5ba73ecab (md5 86d31021, CUDA 12.9 + g++ 13, gated). Same code as the A100 binary (5ba73ecab).
- Deck: docs/handover/gpucmp_1121/t9h_gpucmp.athinput (fresh t9h wedge 512x32x56, 4 MeshBlocks 512x16x28, no restart needed).
  Replace DATA_DIR in the 3 file keys with your copies. These are your FRESH_DIR production files; the md5s must match:
  - ic_t9l62_LIC.txt ea9164342f67e7343e4df9ab5e49c82a (FRESH_DIR/t9l62/)
  - rosseland_mesakap_gs98_x0.70_z0.02.txt d81b7a468b0d9470a09fd43ef03273d5 (FRESH_DIR/tables/)
  - planck_ext2_gs98_x0.70_z0.02.txt 99d668f69b16a502f50f1ec555fb0cd1 (FRESH_DIR/tables/)
- Run: 4 ranks, 1 per H200, KOKKOS_MAP_DEVICE_ID_BY=mpi_rank OMP_NUM_THREADS=1, your usual 1-node srun line. Command-line keys:
  `-i t9h_gpucmp.athinput time/nlim=100 job/final_restart=false output2/dcycle=1000000 output3/dcycle=1000000 output6/dcycle=1000000`
  (hst only, rst off, no clean_pass_damp key).
- While it runs, sample memory in a second step (overlap):
  `nvidia-smi --query-gpu=timestamp,index,memory.used,utilization.gpu --format=csv,noheader -lms 500 > mem.csv`
- Report to caltech_results/GPUCMP/h200_n4/: the log (gz), t9l62.user.hst, t9l62.hydro.hst, mem.csv (gz), and res.txt with:
  node, binary md5, the `elapsed=... cycle=` lines, `main loop time`, peak memory.used per GPU, md5 of both hst files.
- If the node had any other GPU process during the run, run it again.
- R136: skipped (no matching R136 restart on both sides).

-- orion (for the USER)
