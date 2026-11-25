# NOTE 2026-11-24zzzzg6 orion -> Caltech: adaptive-MLT t9: start (USER order, full run)

**USER go (11-25 ~05:00 PST, in the Orion chat): "let's just try a adaptive mlt t9 full run".** One idle node of holder 7333072 (production keeps its priority over the bsg jump pair zzzzg5: if no node is idle, the MLT run takes the next free one and NM2 waits). The preflight of step 2 is the gate; start step 3 only when it is clean. If your rule needs the USER's direct OK for a new production run, ask your coordinator, don't wait for me. Text by T9MLT:


This is a second, separate t9 run; it does not touch t9l62, bsg, m35mms or M35ZAMS. Package docs/handover/ia123rc4_decks/candidate_mlt_1125/: T9MLT_t9l62.athinput, ic_t9_mlt_a10.txt (2233 rows, text), watch_mlt.py, diff_vs_candidate_1124.txt, diff_vs_wedge_A10.txt, keydiff.py, MD5SUMS.

## What it is
- IC = a 1-D column with MLT convection, alpha 1.0: same GM, L, base state, opacity table and EOS as t9l62, M_env(> r_in) 0.3323 Msun, R_ph 56.5 Rsun (3.93e12 cm, T_eff 21.8 kK; the Jiang IC has 93.8 Rsun).
- The MLT flux (max 0.30 F) goes into the radiation as a conservative redistribution, weight w(t) = 1 -> 0 over 3.0e5 s (auto ramp from the column's N^2); after that the star is on its own.
- Deck = candidate_1124 IA123_t9l62 (md5 f2e3b57d) with 10 changed keys, all tagged "# T9MLT" (diff_vs_candidate_1124.txt): he_ic_file; hydro/wb_rmax 7.504e12 -> 4.5e12 (= 1.15 R_ph, as the production deck's ratio);
  mlt_flux_frozen false -> true; mlt_closure frozen -> adaptive; mlt_relax_time 1500; mlt_closure_taper 0.01; mlt_closure_dout 1.08e4; mlt_ramp_start / _time auto; mlt_src_time start (stage FATALs with sdirk2). No dt_fixed. Everything else the same, cfl 0.3, conv_max 1e-5, mg 2, GL 7x18, thinf on, vceil 3e8.
- Seed: the production rule (he_seed_band = n2, amplitude 1e-3, rng 1234), unchanged. It differs from t9l62's because the rule follows the column's own N^2: radial band 2.946e12-3.643e12 cm (Jiang: 4.324e12-5.767e12) and l = 16-64 (3969 modes; Jiang: l 7-26, 680 modes). Expect this; do not "fix" it.
- Binary athena_c4speed_d979d5268 (md5 3627c9a8, checked by the launcher before srun; a000d803c + A1 + A2, bitwise). Parameter names he_*, as in the deck.

## Checked on Orion (CPU, a000d803c)
- The IC file covers r = 1.252e12 to 8.5e13 cm, so x1min - 3 ghosts and x1max 5.632e13 + 5 ghost cells (5.83e13) are inside. R_ph 3.93e12; isothermal-ish atmosphere from there, reaching the density floor 3.6e-19 (= <hydro>/dfloor) at 4.09e12 (58.8 Rsun); T there 17.7 kK falling to 4.8 kK at x1max, above tfloor 1000 K. F_r(x1min) = 6.9466003350e13 (= implicit_flux_x1min to 4e-13), F_MLT = 0 at both x1 ends, so the closed wall and the energy flux BC are untouched.
- Key-by-key: the full deck differs from the wedge A10 deck only by geometry (nx2, nx3, x2, x3 limits, tile, meshblock), tlim, dt_fixed (not set here), seed l named, output cadence, and three named defaults (taper 0.01, src_time start, closure_dout already in A10): diff_vs_wedge_A10.txt. The full deck with only the wedge geometry and cadence put back gives a 20-cycle hst IDENTICAL (all columns) to the wedge A10 run. clean_init 3 passes, residual 4.3e-11, CW 0.
- 64 x 64 wedge arms to 8e5 s (viper): adaptive arm clean (CW 0, DT_RETRY 0, dt at its cap), L_top/L_in 0.95-1.76, no scaffold collapse, KE peak 2.6e43 (scaled) against 1.6e44 of the Jiang wedge. The wedge's control failed a pre-set rule (dt criterion), so this is not a validated prediction; 512 x 256 x 512 was never run. Restart of the closure and ramp bitwise (hst) on CPU slab and HIP wedge, adaptive and frozen.

## What to do (only on this NOTE's posting; new directory, no other run touched)
1. Copy the package to prod_mlt_1125/files (check MD5SUMS) and the deck to prod_mlt_1125/t9mlt; the deck's he_ic_file = /resnick/groups/carnegie_poc/jingze/prod_mlt_1125/files/ic_t9_mlt_a10.txt; tables unchanged (fresh_prod_1120/tables).
2. **Preflight (1 + 100 cycles), outputs and rst off, in prod_mlt_1125/pre:** report clean_init passes + residual + CW, the lines "mlt_flux_frozen" (max F_MLT/F, expect 0.298 at 3.28e12, 0 at the ends), "mlt ramp ... w(t = 0)", "seed ..." (l 16..64, 3969 terms), then dt (expect 96.85), s/cyc, L_top/L_in and L_bot/L_in at c1, 20, 50, 100 (wedge: c0 1.47, c20 1.20, ~1.01 by c40; L_bot 1.001), KE, floor_pay, clamped, BE, DT_RETRY, CW, FATAL, sum_esrc/L_in (hst P_esrc / L_in; round-off, < 1e-12). Post it. If anything is not as listed, or CW / FATAL / inf: STOP, do not continue.
3. If clean: start t9mlt fresh in prod_mlt_1125/t9mlt (4 H200, 4 ranks, one node, the same launcher and env as the other stars), outputs and rst cadence as the deck (hydro_w and vet 3 h, rst 2.592e5 s; first rst hard-linked in rst_keep), tlim 2.682e6 s, chain restarts as the other stars.
4. **STOP criteria (USER/COORD, watch_mlt.py RUNDIR every 10 min, its exit code 3 and RUNDIR/ALERT_MLT say why):** any CLEAN WARNING, FATAL, inf or nan; L_top/L_in < 0.8 or > 3 continuously for 2e4 s after t = 2e4 s (start-up excluded); dt < 20 s for 200 consecutive cycles. On a criterion: clean STOP by your usual mechanism (final rst), post, no relaunch. A slower dt without these is not a criterion: report it.
5. Status at c100, c500, c1000, then every 10 min as for the other stars, with w_mlt, L_MLT/L_in, KE and min dt since the last line. After the ramp (t > 3.0e5 s) watch L_top/L_in; a fall below 0.8 for 2e4 s is the AG Car signature.
6. **Restart rule.** `-r` the newest complete rst, the deck and the same binary, no key changes; the ramp (start 0, time 3.00635e5 s), F_sub and its clock travel in the rst. On every restart check the log for "mlt_ramp auto ... from restart file (computed at the fresh start)"; if it says "no stored numbers in the restart file: NO ramp", STOP and post. The first 5 hst rows must continue the last rows before.

## Cost
- Cycles: 2.682e6 / 96.85 = 27,700 if dt stays at 96.85 (the wedge stayed at its cap to 8e5 s; the Jiang run did not). s/cyc: t9l62 2.4-2.6 s; MLT arm was 1.12x J in the wedge, so 2.7-2.9. Total 21-23 h of stepping, about 22-24 h with outputs.
- Holder 7333072 ends 11-26 05:21 PST. Started by ~21:00 PST 11-24 it would end about 19:00-21:00 PST 11-25, 8-10 h before the wall. It fits only if dt stays near 96.85; if dt is halved for half of the run it does not (then ~35-45 h). At the wall: the launcher's final rst (or the newest periodic one, up to 2.6e5 s of model time lost) and a restart in another holder; a new holder is not requested by us.
- Disk as t9l62 (rst ~88 GB each; keep the newest two).

## Open
The USER chose this on a partial test. The MLT column is 1.7x more compact than Jiang's; it convects from t ~ 2e5 s and its thermal readjustment (3-4e6 s) is longer than tlim.

-- orion (coordinator)
