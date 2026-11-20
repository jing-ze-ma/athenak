# NOTE 2026-11-20 caltech -> orion (USER): batches 1, 2, 3 + v20g gate, running results

Reply to 111be2f40 / b21808faa / 1997c7d5e / 4d9022943. All four received (35 arms). One running NOTE for all; one line per arm
below as each finishes; per-arm files in `docs/handover/caltech_results/<TAG>/`.

## Setup on our side
- 13 builds (batch 1 x6, batch 2/3 x6, v20g), one worktree per commit, each built on its own H200 node in a holder step.
  Your cmake line, except: **Unix Makefiles instead of Ninja** (no ninja on the Caltech nodes) and a shared ccache launcher.
  Toolchain: cuda/12.2.1 (spack, gcc-11.3.1) + system g++ 11.5 + hpcx 2.17.1 OpenMPI 4.1.7; Kokkos 4.6.02 (d8e9af031) from
  the submodule; `-ffp-contract=off` (the tree's CMakeLists adds `--fmad=false` for device code with it).
- `armc.sh` used as is, with two Caltech changes (copy `armc_cal.sh`): it loads the two modules above (runtime libs), and runs
  `res.py` as a step in the holder (our login-node rule). GPU binding: unchanged (KOKKOS_MAP_DEVICE_ID_BY=mpi_rank; every
  arm is checked for 4 distinct GPUs in mem.log).
- Queue order: P1, P2, P3 + v20g gate (restart lane, hpc-sm-02-16 first), P4, P5, then P6 > P7 (FMP/FMM long pair) > P8 > P9.
  Pairs start together (90 s apart) on comparable nodes, never on 02-16 (suspect for timing; only rt/poison arms go there).
- ETA (rough): builds ~20 min; P1-P3 + v20g within ~1.5 h of that; batch 1 done in ~3 h; batches 2/3 done in ~6-7 h.
  P10 (nsys): will check if nsys is on the nodes; last.

## Update after 21d / 21e / 20h
- v20g gate (21d, amended 21e: NH=2, V20GR_rsg added), CONV batch (20h, build conv 9e8234b9b) and conserv-ctl2 5ca370ffc
  (CS_*_ctl; the old-ctl arms had not started, so only ctl2 runs) are queued; builds conserv-ctl2 + conv queued.
- **Master order now: v20g gate > P1-P3 > CONV (CV_* RSG x4 together, then df2/k01, then HG x3) > P4/P5 (LAG) > batch 2/3**,
  unless you say otherwise. v20g and P3 rt/poison arms use the restart lane (02-16 first), so they run alongside P1/P2.
- LW*: real.py output instead of bins (bins deleted after); LV_*: full run.log.gz + npass.py output.

## Results (one line per arm)
- V20GP_rsg (hpc-sm-02-16): rc=0 wall=161 | cycles 589788 - 589792 | steady 2.906 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 |  
- V20GR_rsg (hpc-sm-02-11): rc=0 wall=158 | cycles 589788 - 589792 | steady 2.872 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 |  
- V20G_rsg (hpc-sm-01-17): rc=0 wall=153 rc=0 wall=151 rc=0 wall=148 | cycles 589788 - 589792 | steady 2.858 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF  xcheck vs V20G_rsg/a: V20GP_rsg/rsg.hydro.hst=same V20GP_rsg/rsg.user.hst=same V20GR_rsg/rsg.hydro.hst=same V20GR_rsg/rsg.user.hst=same
- TQK_rsg (hpc-sm-01-11): rc=0 wall=94 | cycles 589788 - 589794 | steady 3.021 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 |  
- V20G_hg (hpc-sm-02-16): rc=137 wall=126 rc=137 wall=116 rc=137 wall=25 | cycles 624961 - 624964 | steady 3.919 s/cycle | GPU MiB max 116703 min 115915 | CW a 0 FATAL 1 | rt: 0 same, 5 DIFF 
  FAILED: last 100 lines in caltech_results/V20G_hg/FAIL_tail100.txt
