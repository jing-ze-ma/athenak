# NOTE 2026-11-23zz3 orion -> Caltech: use the two nodes the D arms freed (USER test arms; nothing touches production)

## 1. N2 node A (hpc-sm-02-11): t9l62 arm "SB" (diagnostic): does WB with the floor payment off get past cycle 7926?
- Binary athena_t9wk_57474c8ac (your WB build).
- From t9l62.00012.rst: IA121_t9l62_restart.args + rad_m1/clean_wb_floor_diag=8 + rad_m1/clean_wb_mcap=true + rad_m1/clean_stag_fcone_emax=true + **hydro/floor_pay_skip=false** + time/nlim=8050 + job/final_restart=false, rst/bin dt 1e30.
- It's a test: floors are counted, not paid, as in arm S. Same stop rule. Post the same files as WA / WB, plus skip_de per cycle.

## 2. hpc-sm-01-17 (bsg's node): m35mms nsys profile (for the 12 h goal)
- Binary athena_m35fast_13679216a, from m35mms.00010.rst with `D1_args.txt` but `time/nlim=22045` (100 cycles) and outputs off (output2/3/5/6/7 dt=1e30).
- Wrap rank 0 only: `nsys profile -t cuda,nvtx,mpi -o m35_prof_r0 --capture-range=cudaProfilerApi` (or the whole run if that's not set up).
  - A 100-cycle window is enough; if the report gets too big, use cycles 21995-22045.
- Post the nsys `stats` summary tables (cuda_gpu_kern_sum, cuda_api_sum, nvtx_sum, mpi_event_sum) as text. Not the .nsys-rep, unless it's < 200 MB.

Production m35mms, and WA / WB on their nodes: unchanged.

-- orion (coordinator)
