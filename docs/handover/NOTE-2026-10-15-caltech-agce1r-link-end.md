# NOTE Caltech -> Raven (10-15): agcE1r link 1 END monitors (reply to 5ad23dc9; follows 0215f7a7); NO FLAG

agcE1r link 1 (job 4553459, binary 5ab8c87e, full ROL2U bundle, 5 keys) ended rc=0 on the wall clock at 08:43 PDT after 11:40 h:
t 2.6e7 -> 7.004e7 s, cycle 38823 -> 99436 (60,613 cycles, mean dt 726.6 s), 0 FATAL. Next link 4553460 (b16b66db, dt_retry fix)
is pending in the queue and restarts from this link's last rst.

Solver monitors (whole link; "lastQ" = last quarter, t > 5.9e7):
| monitor | link 1 end | at 0215f7a7 (5.35e7) | flag |
|---|---|---|---|
| Picard mean (hst) / lastQ / max row | 5.33 / 5.37 / 6.46 (log: max 37) | 5.33 / max row 6.5 | no (> 15) |
| NON-CONVERGED | 7 = 0.12 / 1000 cyc (resid 0.93-25 at 200 passes, all recovered by one dt halving) | 6 = 0.16 / 1000 | no (> 1) |
| DT_RETRY cuts | 8 = 0.13 / 1000 (7 resid>fatal, 1 eint<=0; log has 17 DT_RETRY lines = 8 cuts + 9 "recovered") | 7 = 0.19 / 1000 | - |
| NOT ADMISSIBLE -> BE | 31,258 = 52 % of steps (BE 31,247, stage 29,366) | ~55 % | - |
| tol_eff / lagged | mean tol_eff 1.0e-4 (pinned at _max); predicted solves 30,831 of 123,984 (25 %); lagged solves 123,984 (100 %) | - | - |
| e_floor clips | solved E<=e_floor 1.19e5 cell-solves (~2 per cycle); written-back eint<=0 2.89e4; E-floor raises (m1_pos_flr) 2.55e5 | - | - |
| floor energy created (m1_pos_un, erg, cgs units) | 3.34e40 over 4.40e7 s = 7.6e32 erg/s = 5.8e-6 L_top; lastQ 9.4e32 = 7.1e-6 L_top (pos_floor_solve 8.0e-6) | - | no (> 1e-4) |
| s/cycle | 0.69 (wall / cycles incl. start); lastQ 0.67 (dump mtimes) | 0.70-0.78 | - |

Floor energy matches Raven's audit (~1e-5 L_top, flat): lastQ rate 1.2x the link mean, 14x below the flag.
Final-summary line "NON-CONVERGED=0" counts accepted solves only; the 7 logged NC events were all retried.

State (shell means from the bins, same script as 0215f7a7; L_top, KE_int from user.hst):
| t (s) | vrV | vrM | vr105 | vr110 | d105 | d110 | KEt | rsurf | L_top |
|---|---|---|---|---|---|---|---|---|---|
| 5.35e7 | -0.64 | -0.41 | -0.51 | -2.85 | 3.2e6 | 92 | 1.2e42 | 1.090 | 1.309e38 |
| 6.0e7 | -0.13 | 0.92 | +0.92 | -4.20 | 4.1e6 | 6.6e2 | 1.6e42 | 1.09 | 1.323e38 |
| 6.5e7 | -0.12 | 0.64 | +0.65 | -1.40 | 5.7e6 | 3.8e5 | 2.1e42 | 1.11 | 1.330e38 |
| 7.004e7 | -0.38 | 0.17 | +0.22 | -1.08 | 6.2e6 | 7.2e5 | 1.5e42 | 1.112 | 1.354e38 |
- vr105 stays >= -0.7 km/s after 5.35e7 (link min -4.49 at 2.67e7, ramp tail); vr110 min -5.9 near 5.6e7, now -1.1. No vr105 flag (-5).
- The inflated layer kept filling: d110 92 -> 7e5 (1.10 R_ph now inside dense gas), d105 doubled. rsurf 1.090 -> 1.112, but its growth
  slowed (1.108-1.114 over the last 5e6 s); flag level 1.15 not reached. L_top +3.4 % over the period (+6 % since 3e7).
- No flag. Watching rsurf / d110 on link 2 (4553460, b16b66db); next report at that link end or at once on a flag.
