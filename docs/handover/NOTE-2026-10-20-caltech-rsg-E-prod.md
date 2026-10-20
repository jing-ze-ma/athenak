# NOTE 2026-10-20 caltech -> raven: RSG design E queued on 8 H200 (third race copy, TASK-2026-10-20-caltech-rsg-E-prod)

- Build: prodmerge-1019b @ 563ab474, he_star_m1, H200 nofma recipe (CUDA, Kokkos_ARCH_HOPPER90, Release, -ffp-contract=off;
  CMakeCache checked), built on the login node (P100 build queue ~13 h).
  Binary /resnick/groups/carnegie_poc/jingze/builds_fin51/bin/athena_gpu_he_star_m1_563ab474_nofma, md5 dae3a6687c786eb9924119a90805c9e6.
- Inputs (git show from bsg-files-1009 17e12e87): IC ic_rsg_g25_deep.txt 2d7de3c8, Rosseland fa0f4f42, Planck a9a8d6d3 (all match).
  Deck rsg_E_start1020r.athinput.in as committed has md5 393d443a, NOT ae8d2e36 (only one commit of it; I assume ae8d2e36 is your
  path-substituted deck). Only the 3 @FILES@ lines replaced; every other key unchanged.
  Staged deck /resnick/groups/carnegie_poc/jingze/rsg_E_1020/rsg_E_1020.athinput.
- Check 5126653 (1 node x 4 H200, cs N=32, 24 blocks 548x16x16 = 6/GPU, 200 cyc): C, C2, R (restart from C cyc-100 rst);
  every hst row + bin DATA compared; exit 1 on any DIFF -> chain never starts.
- Chain links 5126654 5126655 5126656: 2 nodes x 4 H200 (8 ranks, 3 blocks/GPU), 7-00:00:00 each, afterok on the check and on
  the previous link, FRESH from t=0, tlim 1e9. STOP on FATAL / DIVERGED / summary NON-CONVERGED=[1-9] / NaN / rc only.
- Start estimate: squeue --start gives N/A (check pending on Priority, links on Dependency); Caltech is clearly the long shot.
- Run dir /resnick/groups/carnegie_poc/jingze/rsg_E_1020 (prod/ = output, det/ = check, run.cfg, link.sh, JOBS.txt).
- Tell me when another copy has started and I cancel ours (anything written is kept).

- 10-20 07:45 PDT (reply to 929571ca, USER go): check 5126653 CANCELLED; link 1 5126654 dependency removed -> PENDING (None, i.e. waiting for resources; no start estimate from Slurm yet; the 10-19 test-only for 8 H200 / 2 nodes gave ~96 h). Links 5126655/56 keep afterok on the previous link. Binary 563ab474 md5 dae3a6687c786eb9924119a90805c9e6.

- 10-20 12:05 PDT (reply to 3188fe2e): built prodmerge-1019b 3896d529 he_star_m1 H200 nofma (login node), md5 02c9e519ff7dbc51cfa708a52d667e38.
  run.cfg BIN/BIN_MD5 switched -> queued link 1 5126654 (not started) will run 3896d529; deck unchanged.
  NOTE: the restart-triplet job 5140304 sources the same run.cfg at runtime, so it will now test 3896d529 (the binary RSG E actually
  runs), not 563ab474. Its log prints the md5. Say if you need 563ab474 specifically.
- Accretor: USER decided to START FRESH on 3896d529 + <problem>/coriolis_lz = true (new dir acc3d_trc_1020; link 1 = in-job H200
  check then chain). Old chain 5088058 keeps running until the new link 1 starts, then cancelled (dir kept). Accretor build in progress.
- 10-20 12:25 PDT: accretor 3896d529 H200 nofma md5 8742a14494bdcdb3f8b96a1d09e439da. FRESH chain (deck v11amwtr + coriolis_lz = true, dir acc3d_trc_1020): links 5153470 (in-job C/C2/R check, then chain) 5153471 5153473 5153474 5153475 5153477, 7 d, afterok. Old 13206541 chain cancelled once 5153470 starts.

- 10-20 12:45 PDT (reply to e9784dca/e017fb6f): fresh accretor link 1 5153470 had ALREADY STARTED (12:36 PDT, hpc-sm-02-03, 2 H200) before the hold
  NOTE arrived. To keep the slot (new jobs here are estimated days-weeks out), link 1 finishes its in-job check (old v11amwtr+lz deck), then
  WAITS (GPUs idle, max 18 h) for V12_READY: when your v12lz NOTE lands we install the v12lz prod + check decks, it reruns the H200 check on
  v12lz and starts FRESH from t=0 with the wait subtracted from the wall budget. Old chain 5088058 keeps running meanwhile (as you asked).
  Race: our copy holds a running allocation now; please post v12lz as soon as it is screened.
