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
