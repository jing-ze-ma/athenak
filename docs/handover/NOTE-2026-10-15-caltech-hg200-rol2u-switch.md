# NOTE Caltech -> Raven (10-14 21:15 PDT): He giant 200-Rsun PRODUCTION switched to ROL2U from day 30 (USER)

- USER decision: the queued production link hgF4w 4513229 (est. start 10-14 22:50, 4 H200) keeps its queue position but now runs
  ROL2U: hegiant_1014/wide4/run/run.cfg BIN = athena_gpu_he_star_m1_5ab8c87e_nofma (fin51-r-hdr = 1bd0c2ce + rst header search
  limit 1 MB; md5 0677cd3be292391bc64225823bbf6bb4), XKEYS = time/tlim=1e9 time/dt_retry=true + the 5 ROL2U rad_m1 keys.
  FIN4 config saved as run.cfg.fin4_d505d21f (revert = copy back). Restart from its own day-30 rst (hegiant.00060, t 2.592e6).
- Basis: hg gate 4552182 PASS (P = B bitwise, C = L2, C = F; ROL2U Picard 7.2 -> 3.1, ~20 % less s/cycle, 4 H200).
- The separate 2-GPU ROL2U He giant chain hgR2w (4552184-86) and its guard were removed (duplicate). hgGuard 4514966 (RLIM 150)
  still guards hgF4w. 1-day twin 4547436 and follow-on 4513230 read the same run.cfg.
- AG Car: agcE1r 4553459-61 queued (5e0f83d0); agcE1s production continues until agcE1r runs cleanly.
