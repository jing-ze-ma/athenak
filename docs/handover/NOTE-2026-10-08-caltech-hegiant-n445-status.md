# NOTE 2026-10-08 Caltech: He giant N445 status + N897 remap/smoke queued

From: Caltech session. To: viper, DeltaAI. Docs only.

## N445 fresh run (to t 9.004208e5 s = 10.42 d)
- The 1-H200 copy (job 4212784, 1 rank, 4 MeshBlocks) started first and cancelled the 2-GPU chain (4212559-62), as planned.
- At 10.6 h: t = 7.40e5 s, cycle ~39150, dt ~17.8 s, no FATAL/NaN/NON-CONVERGED; restarts up to hegiant.00017.rst.
  About 0.97 s/cycle on 1 H200 (viper: 0.862 s/cycle on 4 MI300A). ETA about 2.5 h more, inside the 24 h limit.

## Queued after it (user OK 10-08)
- Remap 4234926 (CPU node, afterok:4212784): final N445 rst -> N897 with hegiant-files-1007/remap/he_remap_giant.py
  `--factor 1 --mbx2 64 --mbx3 64 --grid grid/p2_897.npy --nx1 897` (recipe 1307089d). Checks against viper's numbers:
  rad_cons_* ~1e-16, rad_Fr_change_rel ~0.7 %, rad_hyd_bad_after 0, rad_F_ratio_max 1.
- N897 smoke 4234927 (1 H200, afterok remap): 60 cycles with the first-link keys (cfl 0.3, restart_refill_ghosts=true,
  implicit_opac_newton_slope_max=3, output last_time from rst_info.py). Compare with viper smoke 12120235
  (Picard mean 18.5 max 42, dt 9.22 -> 18.59).
- No N897 production on Caltech until the user decides.

## Duplicates
Our 10.42 -> 30 d N897 leg would duplicate DeltaAI link 1 3337674 (8bcfab91) and the viper hegiant897 chain.
The user decides which copies continue. Do not cancel anything on account of this note.
