# tests_r12: relax_ic.py chases a seconds-long t=0 transient (see tests_r11/RESULTS-0919.md, RETRACTIONS).
# DO NOT USE the ic_passNN.txt profiles it writes.

# 2026-09-19 morning (see memory note he4-r11-floor-arms-verdict, last two paragraphs)
# fz0 = frozen-hydro probe; relax2.py + ic_rx1.txt + fr1 = T-pre-relaxed IC: NULL (dipole unchanged) -> do not use ic_rx1.
# fm0 (mlt off: dipole flips), ff0 (rad force off: same), fc3/fc07 (cfl 0.3/0.0375: dipole LINEAR in dt).
# athena_v11 + r12_arm11.sh: problem/mlt_split_deposit. fs0/fs3 = deposit AFTER the sweep (worse, superseded build);
# ft0/ft3 = deposit BEFORE the sweep: dipole 4.3 % -> 0.4 % of rho g.
# sd3 = 3-turnover recipe + mlt_split_deposit: momentum growth UNCHANGED (4.98e30 vs 4.76e30 at 0.25 turnover), KEr larger
# after 1 turnover, dt collapses 0.04 s at 2.1 turnovers (cancelled). The dipole was a real O(dt) bug but NOT the driver.
# athena_v12 + r12_arm12.sh: problem/e_ledger. el0 (0.5 turnover): energy gain 0.2-0.45 % L, entirely through the two-stream sweep;
# MLT deposit conservative; srcs and hydro cancel. Next: find the non-conservation inside the sweep.
