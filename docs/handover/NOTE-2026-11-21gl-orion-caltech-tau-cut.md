# NOTE 2026-11-21gl orion -> caltech (USER via coordinator, TOP PRIORITY): BSG with a pure optical-depth ray cut

The USER's idea is to cut the sc rays by radial optical depth from the top instead of by the gradient test. v20.3 already supports this with existing keys (CleanRayCut: cut = min(gradient cut, tau_top cut)). The gradient test is switched off with `clean_ray_cut_r = 1e30`, and the depth is set with `clean_ray_cut_taumin`. The decks are committed here (keys in <rad_m1>, not on the CLI). Each is bsg_wngl_fresh with tlim 1.12e6 plus the two keys:
| arm | deck | md5 (8) | priority |
|---|---|---|---|
| **BI_bsgtau1000** | bsg_diag_tau1000.athinput | aaed058b | **TOP: next free node, ahead of BG / BH** |
| **BJ_bsgtau30** | bsg_diag_tau30.athinput | 0f595f02 | second; only if a node is free (tau 30 = the old scheme's vet_col_lat_taucut) |
Both: **B_v203**, fresh start, tlim 1.12e6, outputs as in the deck. Commit the hst, **raycut.txt** (we need the cut shell vs time), oldvec, log tail, plus the **m1prof (NOTE fy) at 0.5 d and at the end**. Refresh the hst about every 20 min, delete the arm's own rst at the end, and put the deck and binary md5 in info.txt. (The tau-10 arm from the coordinator's first message is dropped.)

-- orion (coordinator, for the USER)
