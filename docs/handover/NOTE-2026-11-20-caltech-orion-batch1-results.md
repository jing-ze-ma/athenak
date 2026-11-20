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
- (02:45) Got 7b574f6f9 / b2f7ba19c / dfafbf9f1: v20g3 35a5b37e6 and convdiag 0fbbe9f91 building; V20G3_hg / V20G3_rsg (rt NH=2,
  romio321) at the top, then CD_hg, CE_hg, CXS_hg, CEX_hg (convdiag v2, kappa 0.01; CD_hg log sent gzipped whatever its size).
  V20G2_* never started (v20g2 build was still running; stopped), so no v20g2 results. Old CD_hg v1 never ran, so CD2_hg is skipped.
  CV_estdf_hg had already run before the drop (result above). Ompio cross-check: V20G_hg_ompio and V20F_hg both rt 5/5 same,
  so the HG rst-write FATAL is the romio321 path.
- AB4_rsg (hpc-sm-02-09): rc=0 wall=94 | cycles 589788 - 589794 | steady 3.000 s/cycle | GPU MiB max 109447 min 108315 | CW a 9 FATAL 0 |  
- AB5_rsg (hpc-sm-01-17): rc=0 wall=86 | cycles 589788 - 589794 | steady 2.970 s/cycle | GPU MiB max 108889 min 107703 | CW a 0 FATAL 0 |  
- LW1170_hg (hpc-sm-02-11): rc=0 wall=594 | cycles 624961 - 625061 | steady 5.023 s/cycle | GPU MiB max 115761 min 115163 | CW a 0 FATAL 0 |  (real.py runs when both LW arms are done)
- LW260_hg (hpc-sm-02-04): rc=0 wall=668 | cycles 624961 - 625061 | steady 5.815 s/cycle | GPU MiB max 116343 min 115791 | CW a 0 FATAL 0 |  real.py -> caltech_results/LW_real.txt: LW1170_hg t50114000 5020740/4283389/3464277 max 26809180.00 | t50116000 5223392/4502883/3690951 max 64026220.00 | t50117900 5306397/4565043/3773074 max 77538000.00 | t50119800 5458910/4745484/3979947 max 149545088.00 | t50121800 5771095/5106550/4380722 max 132743576.00 | t50122400 5732220/5063369/43
- LV_loc300k10_rsg (hpc-sm-02-03): rc=0 wall=1037 | cycles 589788 - 590088 | steady 3.220 s/cycle | GPU MiB max 109385 min 108415 | CW a 2 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_loc300k10_rsg/a/run.log: stages 598 passes mean 5.04 max 200 sum 3016 CW 2 | t2s1 4.50/5 t2s2 6
- (02:50) 09cc73c0b real_diag: the LW bins were already deleted after real.py (per 21e item 3), before your NOTE arrived. LW1170_hg / LW260_hg
  took only ~10 min each here, so I queued exact reruns (LWr1170_hg / LWr260_hg, same build v20f + keys, bins kept) at the end of
  the main queue (before batch 4). real_diag.py runs on each in a holder step -> caltech_results/LW_real_diag.txt, then the bins go.
- LV_ref300k10_rsg (hpc-sm-01-11): rc=0 wall=1436 | cycles 589788 - 590088 | steady 4.533 s/cycle | GPU MiB max 109385 min 108415 | CW a 17 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_ref300k10_rsg/a/run.log: stages 601 passes mean 5.78 max 200 sum 3474 CW 17 | t2s1 102.50/200 t
- AB6_rsg (hpc-sm-01-11): rc=1 wall=38 |  | steady  s/cycle | GPU MiB max 77271 min 29151 | CW a 0 FATAL 4 |  
  FAILED: last 100 lines in caltech_results/AB6_rsg/FAIL_tail100.txt
  AB6_rsg: the build rejects the key at setup: `rad_m1/force_reference = wb_arad(_md) with scheme = clean needs the reference work given by the problem generator (force_reference_work = split)` (FATAL, 4 ranks). Not rerun; send a fix or a different key if you want it.
- LV_ref300k01_rsg (hpc-sm-02-09): rc=0 wall=1074 | cycles 589788 - 590088 | steady 3.317 s/cycle | GPU MiB max 109385 min 108415 | CW a 0 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_ref300k01_rsg/a/run.log: stages 600 passes mean 5.98 max 17 sum 3590 CW 0 | t2s2 7.80/17 t2s3 4
- LV_loc300k01_rsg (hpc-sm-01-17): rc=0 wall=1053 | cycles 589788 - 590088 | steady 3.257 s/cycle | GPU MiB max 109385 min 108415 | CW a 2 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_loc300k01_rsg/a/run.log: stages 601 passes mean 5.10 max 200 sum 3067 CW 2 | t2s1 200.00/200 t2
- (03:00) Got 8ed1a5170: v20j 6b5388315 building; V20J_rsg / V20JP_rsg / V20JR_rsg / V20J_hg queued right after V20G3 (off 02-16, romio321). Result lines carry your 3 xchecks (JP/JR vs J/a; user cols 1-26 + hydro vs V20G_rsg/a and V20G_hg_ompio/a; s/cycle side by side).
- HGV40_12 (hpc-sm-02-09): rc=0 wall=290 | cycles 624961 - 625001 | steady 4.920 s/cycle | GPU MiB max 115411 min 114865 | CW a 0 FATAL 0 |  
- FMP300_hg (hpc-sm-02-11): rc=0 wall=1509 | cycles 624961 - 625261 | steady 4.724 s/cycle | GPU MiB max 116469 min 115909 | CW a 0 FATAL 0 |  
- FMM300_hg (hpc-sm-02-04): rc=0 wall=1532 | cycles 624961 - 625261 | steady 4.820 s/cycle | GPU MiB max 116305 min 115765 | CW a 0 FATAL 0 |  
- AB6b_rsg (hpc-sm-02-11): rc=1 wall=28 |  | steady  s/cycle | GPU MiB max 29151 min 29151 | CW a 0 FATAL 4 |  
  FAILED: last 100 lines in caltech_results/AB6b_rsg/FAIL_tail100.txt
  AB6b_rsg also FATALs at setup: `he_star_m1.cpp:638: he_star_m1: needs <rad_m1>/force_reference = wb_arad (the well-balanced pair carries the reference force of the initial column)`. Skipped.
- (03:10) Got 6ea6452db + 2fc2377c5: AB8_rsg queued after AB7; xtloc 947c469d8 building; LV_loc8_300k01/k10 (pair) then LV_ref8_300k10
  queued after the v20j gate. LV_* result lines now also carry `grep -c -i retry run.log`.
- V20JP_rsg (hpc-sm-02-09): rc=0 wall=155 | cycles 589788 - 589792 | steady 2.957 s/cycle | GPU MiB max 108711 min 107577 | CW a 0 FATAL 0 |  xcheck:
- V20JR_rsg (hpc-sm-02-04): rc=0 wall=159 | cycles 589788 - 589792 | steady 2.938 s/cycle | GPU MiB max 108279 min 107091 | CW a 0 FATAL 0 |  xcheck:
- AB8_rsg (hpc-sm-02-04): rc=0 wall=94 | cycles 589788 - 589794 | steady 3.043 s/cycle | GPU MiB max 110447 min 109315 | CW a 0 FATAL 0 |  
- V20J_rsg (hpc-sm-02-11): rc=0 wall=160 rc=0 wall=151 rc=0 wall=153 | cycles 589788 - 589792 | steady 2.946 s/cycle | GPU MiB max 108711 min 107577 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF xcheck: V20JP_rsg/rsg.hydro.hst==V20J_rsg/a V20JP_rsg/rsg.user.hst==V20J_rsg/a V20JR_rsg/rsg.hydro.hst==V20J_rsg/a V20JR_rsg/rsg.user.hst==V20J_rsg/a | rsg user cols1-26==V20G_rsg/a, hydro== (s/cycle 2.946 vs V20G_rsg 2.858)
- (03:20) Got 7d87e5d57: lzcs 8ea56557c building; L1_rsg / L2_hg queued after AB8 (GPU peak in their result lines). v20j and xtloc built (8-9 min each); lzbook (11 files incl. a widely included header) is still compiling, so AB7 waits for it.
- CD_hg (hpc-sm-02-10): rc=0 wall=127 | cycles 624961 - 624963 | steady 33.763 s/cycle | GPU MiB max 120235 min 119677 | CW a 1 FATAL 0 |  
- CE_hg (hpc-sm-02-03): rc=0 wall=181 | cycles 624961 - 624967 | steady 21.521 s/cycle | GPU MiB max 120151 min 119595 | CW a 2 FATAL 0 |  
- CXS_hg (hpc-sm-02-04): rc=0 wall=176 | cycles 624961 - 624967 | steady 18.959 s/cycle | GPU MiB max 120311 min 119755 | CW a 2 FATAL 0 |  
- THF_hg_rt_v7 (hpc-sm-01-11): rc=0 wall=683 rc=0 wall=233 rc=0 wall=682 | cycles 624961 - 624965 | steady 4.688 s/cycle | GPU MiB max 113245 min 112647 | CW a 7 FATAL 0 | rt: 5 same, 0 DIFF 
- (03:30) Got 5df9ce5b5 + 533a80ee2: SPG runs as its own block after L1/L2 (and once lzbook is built): 4 comparable nodes claimed exclusively (nothing else, not even res.py of other arms, on them), round 1 g/p5/lz/j 90 s apart, round 2 rotated by one node; result lines carry node, steady s/cycle, GPU peak and cols1-26+hydro vs SPG_g_1. LEAN: lean 425f2548f building; LN_rsg / LNP_rsg / LN_hg (ompio) are the very last arms (after batch 4).
- CEX_hg (hpc-sm-02-10): rc=0 wall=186 | cycles 624961 - 624967 | steady 22.313 s/cycle | GPU MiB max 120311 min 119755 | CW a 2 FATAL 0 |  
- V20J_hg (hpc-sm-02-09): rc=0 wall=244 rc=0 wall=244 rc=0 wall=250 | cycles 624961 - 624965 | steady 4.377 s/cycle | GPU MiB max 116079 min 115291 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF xcheck: V20JP_rsg/rsg.hydro.hst==V20J_rsg/a V20JP_rsg/rsg.user.hst==V20J_rsg/a V20JR_rsg/rsg.hydro.hst==V20J_rsg/a V20JR_rsg/rsg.user.hst==V20J_rsg/a | rsg user cols1-26==V20G_rsg/a, hydro== (s/cycle 2.946 vs V20G_rsg 2.858) | hg user cols1-26==V20G_hg_ompio/a, hydro== (s/cycle 4.377 vs V20G_hg_ompio 4.200)
- V20G3_rsg (hpc-sm-02-11): rc=0 wall=157 rc=0 wall=145 rc=0 wall=142 | cycles 589788 - 589792 | steady 2.868 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF  xcheck V20G3_rsg/a vs V20G_rsg/a: V20G3_rsg/rsg.hydro.hst=same V20G3_rsg/rsg.user.hst=same
- V20G3_hg (hpc-sm-02-16): rc=0 wall=245 rc=0 wall=238 rc=0 wall=226 | cycles 624961 - 624965 | steady 4.193 s/cycle | GPU MiB max 116703 min 115915 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF 
- AB7_rsg (hpc-sm-02-09): rc=0 wall=278 | cycles 589788 - 589794 | steady 37.669 s/cycle | GPU MiB max 109605 min 108471 | CW a 0 FATAL 0 |  
- CK1_hg (hpc-sm-02-09): rc=0 wall=152 | cycles 624961 - 624967 | steady 11.850 s/cycle | GPU MiB max 120151 min 119595 | CW a 1 FATAL 0 |  
- CK1R_hg (hpc-sm-02-09): rc=0 wall=125 | cycles 624961 - 624967 | steady 5.019 s/cycle | GPU MiB max 120151 min 119595 | CW a 0 FATAL 0 |  
- LV_loc8_300k10_rsg (hpc-sm-01-17): rc=0 wall=1020 | cycles 589788 - 590088 | steady 3.147 s/cycle | GPU MiB max 109425 min 108455 | CW a 2 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_loc8_300k10_rsg/a/run.log: stages 598 passes mean 4.83 max 200 sum 2886 CW 2 | t2s1 4.50/5 t2s2 | grep -ci retry run.log: 7
- LV_loc8_300k01_rsg (hpc-sm-01-11): rc=0 wall=1145 | cycles 589788 - 590088 | steady 3.563 s/cycle | GPU MiB max 109425 min 108455 | CW a 5 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_loc8_300k01_rsg/a/run.log: stages 596 passes mean 6.72 max 200 sum 4007 CW 5 | t2s1 6.36/11 t2s | grep -ci retry run.log: 7
- LV_ref8_300k10_rsg (hpc-sm-02-10): rc=0 wall=1087 | cycles 589788 - 590088 | steady 3.368 s/cycle | GPU MiB max 109425 min 108455 | CW a 2 FATAL 0 |  npass: /resnick/groups/carnegie_poc/jingze/orion_runs_1121/LV_ref8_300k10_rsg/a/run.log: stages 601 passes mean 5.12 max 200 sum 3080 CW 2 | t2s1 200.00/200  | grep -ci retry run.log: 9
- THF300_ref_rsg (hpc-sm-02-04): rc=0 wall=1099 | cycles 589788 - 590088 | steady 3.410 s/cycle | GPU MiB max 94893 min 93577 | CW a 0 FATAL 0 |  
- THF100_ref_hg (hpc-sm-01-17): rc=0 wall=579 | cycles 624961 - 625061 | steady 4.944 s/cycle | GPU MiB max 112977 min 112453 | CW a 0 FATAL 0 |  
- THF300_rsg (hpc-sm-02-03): rc=0 wall=1640 | cycles 589788 - 590088 | steady 5.228 s/cycle | GPU MiB max 94893 min 93867 | CW a 14 FATAL 0 |  
- CXPR_hg (hpc-sm-01-17): rc=0 wall=132 | cycles 624961 - 624967 | steady 5.382 s/cycle | GPU MiB max 115761 min 114973 | CW a 0 FATAL 0 |  hesdirk2: k1_xpos(th=1.000000e-01) cell-stages(rank 0)=5.112400e+05 face_fct calls=1.200000e+01 skipped(no face flux)=0.000000e+00 K1 flux from U^n=1.000000e+00 limited cells=4.137338e+09 max rounds=1.612000e+0
- CXP_hg (hpc-sm-02-04): rc=1 wall=436 | cycles 624961 - 624965 | steady 17.627 s/cycle | GPU MiB max 115761 min 114973 | CW a 11 FATAL 1 |  hesdirk2: k1_xpos' in block 'rad_m1' is absent from the restart's input: added from the command line (= 0.1)
  FAILED: last 100 lines in caltech_results/CXP_hg/FAIL_tail100.txt
  CXP_hg (kappa 0.01 + time2_k1_xpos=0.1): `FATAL <time> dt_retry: cycle 624965 failed 9 times (max 8), dt 4.259657e-01` after 4 cycles (17.6 s/cycle, 11 CLEAN WARNINGs); no end-of-run hesdirk2 summary was printed (the quoted text above is the key-absent notice, ignore it). Full run.log.gz + tail sent.
- THF100_hg (hpc-sm-02-09): rc=0 wall=1132 | cycles 624961 - 625061 | steady 10.568 s/cycle | GPU MiB max 112969 min 112445 | CW a 11 FATAL 0 |  
- URF_rsg (hpc-sm-01-17): rc=0 wall=182 | cycles 589788 - 589818 | steady 3.359 s/cycle | GPU MiB max 108363 min 107257 | CW a 0 FATAL 0 |  
- (04:05) Got 1a2d534cd / 6984c31c3 / 8857717ed / ca6cf5814:
  * LZTEL withdrawn before any arm started; the lztel build was stopped. lzcs-tel 61a578136 building; T1_rsg (xcheck vs V20JR_rsg
    common rows/cols), T2_rsg, T3_hg, T4_rsg (rt 3+3) queued after SPG; L1_rsg / L2_hg moved behind them.
  * **lzcs 8ea56557c does NOT compile** (sm_90, CUDA 12.2.1): `rad_m1_clean.cpp(4965), (4970): error: identifier
    "radm1::RadiationM1::CleanWriteBack()::lzom" is undefined in device code` (3 errors, in rad_m1_clean_p2.cpp). So L1/L2 have not run;
    send a fixed commit if you still want them (if lzcs-tel shares that code, it will fail the same way; I will report).
  * xpos2 330120901 building; CXP_hg had already run (dt_retry FATAL above), CXPR_hg was running and is let finish. CONV2 order queued:
    CNX_hg, CXQ_hg, CXQR_hg, CXQN_hg (full run.log.gz, k1_xpos line quoted).
  * SPG A/B started now (v20j gate done): claiming 4 nodes.
- UR8_rsg (hpc-sm-02-03): rc=0 wall=206 | cycles 589788 - 589818 | steady 3.361 s/cycle | GPU MiB max 104911 min 103627 | CW a 0 FATAL 0 |  
- UR0_rsg (hpc-sm-02-04): rc=0 wall=178 | cycles 589788 - 589818 | steady 3.396 s/cycle | GPU MiB max 94751 min 93823 | CW a 0 FATAL 0 |  
- GCA300_rsg (hpc-sm-01-11): rc=0 wall=1338 | cycles 589788 - 590088 | steady 4.227 s/cycle | GPU MiB max 95661 min 94343 | CW a 15 FATAL 0 |  
- CNX_hg (hpc-sm-01-17): rc=0 wall=180 | cycles 624961 - 624967 | steady 17.710 s/cycle | GPU MiB max 120151 min 119595 | CW a 2 FATAL 0 |  k1_xpos: 
- (04:15) Got 1912abd32 / a07337754 / 6430ac473: KT HG trio, then KT RSG trio, now at the head of the queue (each trio starts together on 3 comparable nodes, not 02-16; full run.log.gz). SPG A/B moved behind the KT trios (its node claims released). lzcs2 b4ba38764 building; L1_rsg / L2_hg re-queued on it behind T1-T4.
- GCA300_ref_rsg (hpc-sm-02-10): rc=0 wall=1437 | cycles 589788 - 590088 | steady 4.554 s/cycle | GPU MiB max 95621 min 94595 | CW a 17 FATAL 0 |  
- KT1_hg (hpc-sm-01-11): rc=0 wall=193 | cycles 624961 - 624981 | steady 5.190 s/cycle | GPU MiB max 116353 min 115805 | CW a 0 FATAL 0 |  
- KT01_hg (hpc-sm-02-10): rc=0 wall=197 | cycles 624961 - 624981 | steady 5.577 s/cycle | GPU MiB max 116353 min 115805 | CW a 0 FATAL 0 |  
- RF100_hg_f1 (hpc-sm-02-03): rc=0 wall=609 | cycles 624961 - 625061 | steady 5.063 s/cycle | GPU MiB max 112977 min 112453 | CW a 0 FATAL 0 |  
- RF100_ref_hg (hpc-sm-02-04): rc=0 wall=595 | cycles 624961 - 625061 | steady 4.942 s/cycle | GPU MiB max 112977 min 112453 | CW a 0 FATAL 0 |  
- KT1_rsg (hpc-sm-01-11): rc=0 wall=128 | cycles 589788 - 589808 | steady 2.771 s/cycle | GPU MiB max 108777 min 107589 | CW a 0 FATAL 0 |  
- KT001_hg (hpc-sm-01-17): rc=1 wall=431 | cycles 624961 - 624965 | steady 17.544 s/cycle | GPU MiB max 115761 min 114973 | CW a 11 FATAL 1 |  
  FAILED: last 100 lines in caltech_results/KT001_hg/FAIL_tail100.txt
- KT01_rsg (hpc-sm-02-10): rc=0 wall=136 | cycles 589788 - 589808 | steady 3.218 s/cycle | GPU MiB max 109209 min 108075 | CW a 0 FATAL 0 |  
  KT001_hg: ### FATAL ERROR: <time> dt_retry: cycle 624965 failed 9 times (max 8), dt 4.259657e-01 after 4 cycles (11 CLEAN WARNINGs, 17.5 s/cycle); same failure as CXP_hg. Full run.log.gz sent.
- KT001_rsg (hpc-sm-02-03): rc=0 wall=161 | cycles 589788 - 589808 | steady 4.289 s/cycle | GPU MiB max 109209 min 108075 | CW a 0 FATAL 0 |  
- T1_rsg (hpc-sm-01-11): rc=0 wall=97 | cycles 589788 - 589794 | steady 2.832 s/cycle | GPU MiB max 108711 min 107577 | CW a 0 FATAL 0 |  xcheck: rsg.hydro.hst rows1-4 cols1-10 ==V20JR_rsg/a rsg.user.hst rows1-4 cols1-30 ==V20JR_rsg/a | CW 0
- RF300_rsg_f1 (hpc-sm-02-11): rc=0 wall=1112 | cycles 589788 - 590088 | steady 3.448 s/cycle | GPU MiB max 94893 min 93915 | CW a 0 FATAL 0 |  
- RF300_ref_rsg (hpc-sm-02-09): rc=0 wall=1093 | cycles 589788 - 590088 | steady 3.388 s/cycle | GPU MiB max 94893 min 93923 | CW a 0 FATAL 0 |  
- T4_rsg (hpc-sm-02-16): rc=0 wall=162 rc=0 wall=152 rc=0 wall=151 | cycles 589788 - 589794 | steady 2.828 s/cycle | GPU MiB max 108711 min 107577 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF 
- L1_rsg (hpc-sm-02-09): rc=0 wall=94 | cycles 589788 - 589794 | steady 2.824 s/cycle | GPU MiB max 109587 min 108453 | CW a 0 FATAL 0 |  
- L2_hg (hpc-sm-02-09): rc=0 wall=115 | cycles 624961 - 624967 | steady 5.003 s/cycle | GPU MiB max 121049 min 120491 | CW a 0 FATAL 0 |  
- T3_hg (hpc-sm-02-11): rc=0 wall=345 | cycles 624961 - 625011 | steady 5.176 s/cycle | GPU MiB max 115715 min 115169 | CW a 0 FATAL 0 |  
- (04:30) KC1c_hg / KC009_hg / KC1m_hg at the head of the queue (KT HG + RSG trios are done, above). ELEAK (20z + 20al): AB9_rsg / AB9o_rsg / AB9rt_rsg had not started; eleak built at a411520d7 and the three queued right after KC (result lines carry the floor_audit-cycle and CLEAN WARNING..unpaid counts).
- KC1c_hg (hpc-sm-02-09): rc=0 wall=198 | cycles 624961 - 624981 | steady 5.183 s/cycle | GPU MiB max 116353 min 115805 | CW a 0 FATAL 0 |  
- KC009_hg (hpc-sm-02-11): rc=0 wall=206 | cycles 624961 - 624981 | steady 5.616 s/cycle | GPU MiB max 116353 min 115805 | CW a 0 FATAL 0 |  
- SPG_g_1 (hpc-sm-02-04): rc=0 wall=106 | cycles 589788 - 589798 | steady 2.740 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_p5_1 (hpc-sm-02-10): rc=0 wall=98 | cycles 589788 - 589798 | steady 2.792 s/cycle | GPU MiB max 108265 min 107077 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_lz_1 (hpc-sm-01-17): rc=0 wall=101 | cycles 589788 - 589798 | steady 2.725 s/cycle | GPU MiB max 109227 min 108075 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_j_1 (hpc-sm-02-03): rc=0 wall=103 | cycles 589788 - 589798 | steady 2.790 s/cycle | GPU MiB max 108279 min 107091 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_g_2 (hpc-sm-02-10): rc=0 wall=100 | cycles 589788 - 589798 | steady 2.735 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_p5_2 (hpc-sm-01-17): rc=0 wall=108 | cycles 589788 - 589798 | steady 2.801 s/cycle | GPU MiB max 108697 min 107563 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_lz_2 (hpc-sm-02-03): rc=0 wall=99 | cycles 589788 - 589798 | steady 2.728 s/cycle | GPU MiB max 108777 min 107589 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- SPG_j_2 (hpc-sm-02-04): rc=0 wall=102 | cycles 589788 - 589798 | steady 2.790 s/cycle | GPU MiB max 108711 min 107577 | CW a 0 FATAL 0 |  xcheck: user cols1-26==SPG_g_1/a, hydro==SPG_g_1/a
- KC1m_hg (hpc-sm-02-09): rc=0 wall=196 | cycles 624961 - 624981 | steady 5.172 s/cycle | GPU MiB max 116353 min 115805 | CW a 0 FATAL 0 |  
- T2_rsg (hpc-sm-01-11): rc=1 wall=830 | cycles 589788 - 589864 | steady 3.492 s/cycle | GPU MiB max 108941 min 107971 | CW a 19 FATAL 1 |  
  FAILED: last 100 lines in caltech_results/T2_rsg/FAIL_tail100.txt
  T2_rsg (lzcs-tel, key on, 300 cyc): `FATAL <time> dt_retry: cycle 589864 failed 9 times (max 8), dt 3.681391e+00` at cycle N0+76 (19 CLEAN WARNINGs). Full log + tail sent.
- (04:40) Got ec649bc78: T3K_hg (key off, 50 cyc) queued right after KC / AB9, prefers 02-11; its result line carries the per-cycle dt list.
- CXQR_hg (hpc-sm-02-10): rc=0 wall=135 | cycles 624961 - 624967 | steady 5.051 s/cycle | GPU MiB max 116341 min 115785 | CW a 0 FATAL 0 |  k1_xpos: <rad_m1> time_scheme=hesdirk2: stage steps=6.000000e+00 backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00 vet clips=0.000000e+00 k1_xpos(th=1.000000e-01) cell-stages(rank 0)=3.243545e+06 face_fct calls=1.200000e+01 skipped(no face flux)=
- CXQ_hg (hpc-sm-02-11): rc=1 wall=436 | cycles 624961 - 624965 | steady 17.722 s/cycle | GPU MiB max 115761 min 114973 | CW a 11 FATAL 1 |  k1_xpos: ### NOTE: parameter 'time2_k1_xpos' in block 'rad_m1' is absent from the restart's input: added from the command line (= 0.1)
  FAILED: last 100 lines in caltech_results/CXQ_hg/FAIL_tail100.txt
  CXQ_hg (xpos2, kappa 0.01): same `dt_retry: cycle 624965 failed 9 times` FATAL as CXP_hg / KT001_hg (the quoted k1_xpos text is only the key-absent notice; no end-of-run summary).
- (04:50) Got 81f5d6818 + 8132e4882: T2K_rsg queued next to T3K_hg; v20l e067945c2 building (plain cmake line); V20L_rsg / V20L_hg (xcheck all columns vs V20J_*/a) next, then SPL_g/l rounds 1+2 on 2 exclusively claimed nodes (swapped in round 2).
- CXQN_hg (hpc-sm-01-17): rc=1 wall=466 | cycles 624961 - 624965 | steady 17.715 s/cycle | GPU MiB max 115761 min 114973 | CW a 11 FATAL 1 |  k1_xpos: ### NOTE: parameter 'time2_k1_xpos' in block 'rad_m1' is absent from the restart's input: added from the command line (= 0.1)
  FAILED: last 100 lines in caltech_results/CXQN_hg/FAIL_tail100.txt
- AB9_rsg (hpc-sm-02-11): rc=0 wall=93 | cycles 589788 - 589794 | steady 2.866 s/cycle | GPU MiB max 109531 min 108343 | CW a 0 FATAL 0 |  floor_audit cycle lines 18 | CLEAN WARNING..unpaid lines 0
  CXQN_hg (xpos2, no xt_local, kappa 0.01): same dt_retry FATAL at cycle 624965 (all kappa-0.01 HG arms from hegiant.00139 so far: CV_tight_hg CW, CXP, KT001, CXQ, CXQN).
- AB9o_rsg (hpc-sm-02-10): rc=0 wall=98 | cycles 589788 - 589794 | steady 2.862 s/cycle | GPU MiB max 109837 min 108703 | CW a 0 FATAL 0 |  floor_audit cycle lines 18 | CLEAN WARNING..unpaid lines 0
- T3K_hg (hpc-sm-01-11): rc=0 wall=371 | cycles 624961 - 625011 | steady 5.252 s/cycle | GPU MiB max 115717 min 115169 | CW a 0 FATAL 0 |  dt: 6.896552e+01 1.005172e+02 1.005172e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 1.063636e+02 5.000000e+01 
- LWr1170_hg (hpc-sm-02-09): rc=0 wall=607 | cycles 624961 - 625061 | steady 5.044 s/cycle | GPU MiB max 116347 min 115799 | CW a 0 FATAL 0 |  real_diag -> caltech_results/LW_real_diag.txt
- BA_rsg_v20g (hpc-sm-01-17): rc=0 wall=131 | cycles 589788 - 589808 | steady 2.695 s/cycle | GPU MiB max 109195 min 108061 | CW a 0 FATAL 0 |  
- (05:00) a7b20e2ef: CXQN_hg and CXQ_hg had both already run (results above), nothing to drop. 1508dd7a7: eleak2 = 8d8ed6e27 building; AB10_hg / AB10o_hg queued after the v20l gate (floor_audit count, unpaid CW count, EOS floors summary in the line).
- LWr260_hg (hpc-sm-02-03): rc=0 wall=697 | cycles 624961 - 625061 | steady 5.835 s/cycle | GPU MiB max 116343 min 115791 | CW a 0 FATAL 0 |  real_diag -> caltech_results/LW_real_diag.txt
- BA_rsg_p4 (hpc-sm-02-11): rc=0 wall=135 | cycles 589788 - 589808 | steady 2.685 s/cycle | GPU MiB max 108697 min 107563 | CW a 0 FATAL 0 |  
- BA_rsg_h200 (hpc-sm-01-11): rc=0 wall=129 | cycles 589788 - 589808 | steady 2.759 s/cycle | GPU MiB max 108265 min 107077 | CW a 0 FATAL 0 |  
- AB9rt_rsg (hpc-sm-02-16): rc=0 wall=175 rc=0 wall=158 rc=0 wall=131 | cycles 589788 - 589794 | steady 2.870 s/cycle | GPU MiB max 109963 min 108829 | CW a 0 FATAL 0 | rt: 5 same, 0 DIFF floor_audit cycle lines 18 | CLEAN WARNING..unpaid lines 0
- BA_rsg_mb2 (hpc-sm-02-11): rc=0 wall=137 | cycles 589788 - 589808 | steady 2.978 s/cycle | GPU MiB max 108545 min 107411 | CW a 0 FATAL 0 |  
- (05:05) Scheduler fix on our side: the main queue had drained at 04:43 and its runner exited, so batch 4 (BA_*) started and the later additions (T2K_rsg, V20L_*, AB10_*) waited ~20 min. Runner is persistent now and batch 4 yields to any waiting main-queue arm.
- BA_hg_v20g (hpc-sm-01-11): rc=0 wall=176 | cycles 624961 - 624981 | steady 4.904 s/cycle | GPU MiB max 120005 min 119457 | CW a 0 FATAL 0 |  
- (05:15) Got 0c0379c94: spd2b 170812c47 building (plain cmake); SD_rsg / SDR_rsg / SDP_rsg + the HG pair SDT_b_hg/SDT_s_hg queued after AB10; SDT_b/s RSG rounds 1+2 on 2 exclusively claimed nodes after the SPL rounds. Xchecks as asked in the result lines.
- AB10_hg (hpc-sm-02-09): rc=0 wall=122 | cycles 624961 - 624967 | steady 5.168 s/cycle | GPU MiB max 120991 min 120435 | CW a 0 FATAL 0 |  floor_audit cycle lines 18 | unpaid CW 0 | EOS floors:  hydro EOS floors energy created=2.721342e+33
- BA_hg_p4 (hpc-sm-02-11): rc=0 wall=179 | cycles 624961 - 624981 | steady 4.876 s/cycle | GPU MiB max 119391 min 118843 | CW a 0 FATAL 0 |  
- BA_hg_h200 (hpc-sm-01-11): rc=0 wall=177 | cycles 624961 - 624981 | steady 5.092 s/cycle | GPU MiB max 118789 min 118191 | CW a 0 FATAL 0 |  
- AB10o_hg (hpc-sm-02-10): rc=0 wall=109 | cycles 624961 - 624967 | steady 5.160 s/cycle | GPU MiB max 120831 min 120275 | CW a 0 FATAL 0 |  floor_audit cycle lines 18 | unpaid CW 0 | EOS floors:  hydro EOS floors energy created=2.721292e+33
- (05:25) Got 20af47b76: ksplit branch head = 1f9926b4982f (matches); building; KS1_hg then KS001_hg queued after the spd2b arms (KSPLIT/OLDVEC counts, full log). The LWr bins were already deleted after real_diag, so LWs1170_hg / LWs260_hg reruns are queued (low) for real_nb.py -> caltech_results/LW_real_nb.txt, then their bins go.
- (05:30) Got 810b64b31: ksplit2 = 6b18f5de5 building; KX1_hg / KX001_hg after KS1/KS001 (KS* use 1f9926b49). Result lines carry CW/FATAL, KSPLIT/OLDVEC counts, k1_xlim and face_fct left (Q<=0) text.
- (05:12 PST; the earlier (hh:mm) stamps in this NOTE ran ~20 min ahead, ignore them) Got 064b86841: v20m cf584d3a6 queued to build right after v20l (both on 02-16, ccache-warm); V20M_rsg / V20MP_rsg / V20M_hg queued right after the V20L gate, xchecks vs V20L_*/a (V20J_*/a if V20L is not done). Scheduler fix: from 04:48 to 05:09 the nodes sat idle (batch 4 yielded to main-queue arms that were waiting for builds); batch 4 now yields only to arms whose binary exists. T2K_rsg started 05:09 on 02-11.
- (05:13) Got fa5e9209d + c002f7759: KX1_hg / KX001_hg had not started, so they now run on ksplit3 = baf94d6e9 (the 6b18f5de5 build was stopped). LZTEL2: lztel2 a1874999f + lztel2-study 79406ea57 building; Z1..Z4, ZS_hg, ZS_rsg queued right after T2K_rsg (Z1 xcheck vs T1_rsg, Z3 dt list, CW counts).
- BA_hg_mb2 (hpc-sm-02-10): rc=0 wall=185 | cycles 624961 - 624981 | steady 5.252 s/cycle | GPU MiB max 119215 min 118667 | CW a 0 FATAL 0 |  
- BA_hg_v20g_m1 (hpc-sm-01-17): rc=0 wall=188 | cycles 624961 - 624981 | steady 4.977 s/cycle | GPU MiB max 116339 min 115791 | CW a 0 FATAL 0 |  
- Z1_rsg (hpc-sm-01-17): rc=0 wall=90 | cycles 589788 - 589794 | steady 2.839 s/cycle | GPU MiB max 108279 min 107091 | CW a 0 FATAL 0 |  xcheck: rsg.hydro.hst==T1_rsg/a rsg.user.hst==T1_rsg/a | CW 0
