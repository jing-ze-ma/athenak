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
  **V20G_hg needs a fix from you**: every leg dies at an HG restart WRITE: `### FATAL ERROR in ../../src/outputs/restart.cpp at line 1654:
  rad_m1 hesdirk2 slope not written correctly to rst file, restart file is broken.` (leg a at the final rst after reaching N0+3;
  leg b at its N0+2 dump, so leg b2 reads a broken rst). RSG rst writes of the same build are fine (V20G_rsg rt 5/5 same).
  V20F_hg (v20f, same form) is queued next on 02-16 and will show whether it is v20g-specific. Tail in caltech_results/V20G_hg/.
- Memory note: GPU peak is ~109 GB/GPU for RSG and ~117 GB/GPU for HG at 4 ranks (you expected 50 / 70); both fit in 141 GB.
- TQS_rsg_tie (hpc-sm-02-09): rc=0 wall=91 | cycles 589788 - 589794 | steady 3.011 s/cycle | GPU MiB max 109141 min 107955 | CW a 0 FATAL 0 |  
- CS_rsg_ctl (hpc-sm-01-11): rc=0 wall=101 | cycles 589788 - 589796 | steady 2.859 s/cycle | GPU MiB max 110487 min 109355 | CW a 0 FATAL 0 |  
- V20F_rsg (hpc-sm-02-11): rc=0 wall=154 rc=0 wall=147 rc=0 wall=152 | cycles 589788 - 589792 | steady 2.864 s/cycle | GPU MiB max 108763 min 107661 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF  xcheck vs V20F_rsg/a:
- CS_rsg_tie (hpc-sm-02-10): rc=0 wall=103 | cycles 589788 - 589796 | steady 3.082 s/cycle | GPU MiB max 110517 min 109385 | CW a 0 FATAL 0 |  
- CS_hg_ctl (hpc-sm-01-11): rc=0 wall=124 | cycles 624961 - 624969 | steady 5.322 s/cycle | GPU MiB max 121651 min 121095 | CW a 0 FATAL 0 |  
- (02:15) Got 03aa43e0b + bef7f30cc: v20g2 240d9eb32 building now; V20G2_hg / V20G2_rsg (rt NH=2, default romio321) are at the
  top of the queue, then V20G_hg_ompio (OMPI_MCA_io=ompio cross-check). V20F_hg, THF_hg_rt_v7 and batch-4 BR_hg_h200 run with
  ompio. Batch 4 is queued after everything (builds h200-p4/h200/h200mb2 under way; nsys 2023.2.3 is on the nodes; no Kokkos
  nvtx connector found, so PF_* run without KOKKOS_TOOLS_LIBS). Starts are now staggered 30 s (RSG restart read is ~24 s).
- CS_hg_tie (hpc-sm-02-09): rc=0 wall=135 | cycles 624961 - 624969 | steady 5.551 s/cycle | GPU MiB max 121663 min 121107 | CW a 0 FATAL 0 |  
- V20FP_rsg (hpc-sm-02-11): rc=0 wall=159 | cycles 589788 - 589792 | steady 2.884 s/cycle | GPU MiB max 108763 min 107575 | CW a 0 FATAL 0 |   xcheck vs V20F_rsg/a: V20FP_rsg/rsg.hydro.hst=same V20FP_rsg/rsg.user.hst=same
- V20G_hg_ompio (hpc-sm-01-11): rc=0 wall=239 rc=0 wall=232 rc=0 wall=231 | cycles 624961 - 624965 | steady 4.200 s/cycle | GPU MiB max 116703 min 115915 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF 
- V20F_hg (hpc-sm-02-09): rc=0 wall=252 rc=0 wall=229 rc=0 wall=251 | cycles 624961 - 624965 | steady 4.184 s/cycle | GPU MiB max 116703 min 115915 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF 
- CV_ref_rsg (hpc-sm-02-11): rc=0 wall=129 | cycles 589788 - 589808 | steady 2.622 s/cycle | GPU MiB max 109321 min 108187 | CW a 0 FATAL 0 |  
- CV_tight_rsg (hpc-sm-01-11): rc=0 wall=142 | cycles 589788 - 589808 | steady 3.328 s/cycle | GPU MiB max 109321 min 108187 | CW a 0 FATAL 0 |  
- CV_estdf_rsg (hpc-sm-01-17): rc=0 wall=133 | cycles 589788 - 589808 | steady 2.933 s/cycle | GPU MiB max 108889 min 107701 | CW a 0 FATAL 0 |  
- CV_est_rsg (hpc-sm-02-03): rc=0 wall=125 | cycles 589788 - 589808 | steady 2.763 s/cycle | GPU MiB max 109321 min 108187 | CW a 0 FATAL 0 |  
- CV_k01_rsg (hpc-sm-02-04): rc=0 wall=132 | cycles 589788 - 589808 | steady 2.920 s/cycle | GPU MiB max 109321 min 108187 | CW a 0 FATAL 0 |  
- CV_ref_hg (hpc-sm-02-11): rc=0 wall=191 | cycles 624961 - 624981 | steady 4.999 s/cycle | GPU MiB max 115921 min 115323 | CW a 0 FATAL 0 |  
- CV_k01_hg (hpc-sm-01-11): rc=0 wall=189 | cycles 624961 - 624981 | steady 5.212 s/cycle | GPU MiB max 116499 min 115951 | CW a 0 FATAL 0 |  
- CV_df2_rsg (hpc-sm-02-03): rc=0 wall=130 | cycles 589788 - 589808 | steady 3.058 s/cycle | GPU MiB max 109321 min 108187 | CW a 0 FATAL 0 |  
- CV_estdf_hg (hpc-sm-01-17): rc=0 wall=223 | cycles 624961 - 624981 | steady 6.834 s/cycle | GPU MiB max 116499 min 115951 | CW a 0 FATAL 0 |  
- CV_tight_hg (hpc-sm-02-09): rc=0 wall=276 | cycles 624961 - 624981 | steady 9.653 s/cycle | GPU MiB max 116499 min 115951 | CW a 4 FATAL 0 |  
- AB1_rsg (hpc-sm-01-17): rc=0 wall=90 | cycles 589788 - 589794 | steady 2.999 s/cycle | GPU MiB max 109015 min 107829 | CW a 0 FATAL 0 |  
- AB3_rsg (hpc-sm-02-09): rc=0 wall=93 | cycles 589788 - 589794 | steady 2.854 s/cycle | GPU MiB max 108391 min 107205 | CW a 0 FATAL 0 |  
- AB2_rsg (hpc-sm-01-17): rc=0 wall=114 | cycles 589788 - 589794 | steady 7.972 s/cycle | GPU MiB max 108943 min 107757 | CW a 1 FATAL 0 |  
