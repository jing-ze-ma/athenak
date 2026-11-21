# NOTE 2026-11-21ff orion -> caltech: pc-sweep timing A/B on B_pcsym (LOWEST priority, free nodes only)

- Build caltech-runs/pcsym = 7c4f624f5 as **B_pcsym**. It is v20.3 d4b08dae4 + 15 lines adding the key rad_m1/clean_pc_sweep: rb (default, bitwise) or sym. Use the same recipe as B_v203 (he_star_m1); incremental is fine. Record the md5 and that the changed TUs compiled.
- Arm file: caltech_runs_1121/arms_pcsym.txt. It has 4 arms of 100 cycles each, MODE=speed, on restarts hegiant.00139 / rsg.01295 with the V203L recipe keys, rb vs sym per star. Run the rb/sym pair of each star **on the same node, back to back**.
- Report per arm: steady s/cycle, main loop time, total CLEAN iterations (the sum of the "(sum N)" values in the CLEAN cycle lines), CW, and the hst md5.
- Check: **rb must be hst-bitwise to V203L_*** on common rows. sym is expected to differ at solver tolerance.
- Order: after the H200SPD arms (NOTE ez / fe). Production and every launch item come first.

-- orion
