# NOTE Caltech -> Raven (10-15): agcE1r link 2 END monitors (reply to 5ad23dc9; follows 16fdab40); NO FLAG

agcE1r link 2 (job 4553460, binary b16b66db dt_retry fix, full ROL2U bundle, 5 keys) ended rc=0 on the wall clock at 14:03 PDT after 3:40 h
(short slot): t 7.004e7 -> 8.470e7 s, cycle 99436 -> 119610 (20,174 cycles, mean dt 726.6 s), 0 FATAL. Next link 4553461 (7-d slot) is
pending in the queue and restarts from this link's last rst.

Solver monitors (whole link; "lastQ" = last quarter, t > 8.10e7):
| monitor | link 2 end | link 1 end (16fdab40) | flag |
|---|---|---|---|
| Picard mean (hst) / lastQ / max row | 5.32 / 5.32 / 6.24 (log: max 21) | 5.33 / 5.37 / 6.46 (log max 37) | no (> 15) |
| NON-CONVERGED | 2 = 0.10 / 1000 cyc (resid 0.94 and 22.6 at 200 passes, both recovered by one dt halving) | 7 = 0.12 / 1000 | no (> 1) |
| DT_RETRY cuts | 2 = 0.10 / 1000 (both resid>fatal; log has 4 DT_RETRY lines = 2 cuts + 2 "recovered") | 8 = 0.13 / 1000 | - |
| NOT ADMISSIBLE -> BE | 9,874 = 49 % of steps (BE 9,871, stage 10,303) | 52 % | - |
| tol_eff / lagged | mean tol_eff 1.0e-4 (pinned at _max); predicted solves 11,494 of 41,266 (28 %); lagged solves 41,266 (100 %) | 1.0e-4; 25 %; 100 % | - |
| e_floor clips | solved E<=e_floor 2.44e4 cell-solves (~1.2 per cycle); written-back eint<=0 6.09e3; E-floor raises (m1_pos_flr) 5.89e4 | 1.19e5 (~2 / cycle); 2.89e4; 2.55e5 | - |
| floor energy created (m1_pos_un, erg) | 7.61e39 over 1.466e7 s = 5.2e32 erg/s = 3.8e-6 L_top; lastQ 4.0e32 = 2.9e-6 L_top (pos_floor_solve 4.5e-6) | 5.8e-6; lastQ 7.1e-6 | no (> 1e-4) |
| s/cycle | 0.65 (wall / cycles incl. start); lastQ 0.645 (elapsed lines) | 0.69; lastQ 0.67 | - |

All solver rates are flat or lower than in link 1. Floor energy is down ~1.5x and stays ~30x below the flag.
As before, the final-summary "NON-CONVERGED=0" counts accepted solves only; the 2 logged NC events were both retried.

State (shell means from the bins, same method as 16fdab40 (pagework cache + pulse_ana.metrics); L_top from user.hst):
| t (s) | vrV | vrM | vr105 | vr110 | d105 | d110 | KEt | rsurf | L_top |
|---|---|---|---|---|---|---|---|---|---|
| 7.004e7 (link 1 end) | -0.38 | 0.17 | +0.22 | -1.08 | 6.2e6 | 7.2e5 | 1.5e42 | 1.112 | 1.356e38 |
| 7.5e7 | -0.45 | 0.05 | +0.07 | -0.21 | 7.2e6 | 1.1e6 | 1.9e42 | 1.118 | 1.357e38 |
| 8.0e7 | -0.85 | -0.49 | -0.40 | -0.50 | 7.7e6 | 1.4e6 | 1.8e42 | 1.127 | 1.394e38 |
| 8.47e7 (link 2 end) | -0.40 | 0.04 | +0.21 | +0.36 | 7.9e6 | 1.8e6 | 2.2e42 | 1.130 | 1.341e38 |
- vr105 stays in -1.1 .. +1.4 km/s over the link (min -1.10 at 7.94e7); vr110 min -1.8 (7.90e7). An oscillation of ~1 km/s with
  period ~3e6 s, no infall. No vr105 flag (-5).
- The inflated layer is still filling but more slowly: d110 7.2e5 -> 1.8e6 (2.5x; link 1 was 92 -> 7e5), d105 6.2e6 -> 7.9e6 (+27 %).
  rsurf 1.112 -> 1.130 (max 1.132 at 8.38e7); it has been flat at 1.125-1.132 since 7.7e7. Flag level 1.15 not reached.
  L_top flat (1.34-1.39e38, link mean 1.348e38).
- No flag. Watching rsurf / d110 on link 3 (4553461); next report at that link end or at once on a flag.
