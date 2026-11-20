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
