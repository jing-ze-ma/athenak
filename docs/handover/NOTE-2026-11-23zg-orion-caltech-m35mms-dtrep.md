# NOTE 2026-11-23zg orion -> Caltech: m35mms dt diagnostic (USER goal: m35mms to 10 d in ≤ 24 h on 1 node)

The USER wants m35mms to reach 10 days (8.64e5 s) in ≤ 24 h wall on one node (4 H200), without losing accuracy. Today it is on course for ~40 h, because dt dropped from 12 to ~7.4 after convection onset and late cycles cost ~1.7x the early ones. To find the cause we need two measurements on the REAL late state. Neither changes any state.

When a node frees (bsg's node when bsg finishes; ahead of the M35ZAMS smoke, since this takes only ~10 min):
1. Restart m35mms from its NEWEST rst into a separate scratch dir, same IA11 binary and deck:
   - `time/nlim=<rst cycle + 100>`, `job/final_restart=false`, no rst / bin outputs;
   - env `ATHENA_DT_REPORT=10`.
2. Post:
   - every DTREP / DTREG line: limiting cell, dir, r, rho, T, cgas, v, ceff, tauT_h, Erad, Ddd, cell_dt, cell_dt_gas, and the global cfl*dt_gasonly;
   - the CLEAN cycle / solves lines for those 100 cycles;
   - the same CLEAN cycle / solves lines for 20 cycles near cycle 2000 from the production log (early, for comparison).
3. Tell us the newest m35mms rst cycle / t and its size.

Do not stop or touch the production m35mms run.

-- orion (coordinator)
