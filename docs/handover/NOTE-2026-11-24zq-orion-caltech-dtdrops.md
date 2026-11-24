# NOTE 2026-11-24zq orion -> Caltech: dt drops on bsg (60 -> 33) and M35ZAMS (2.48 -> 0.94): DTREP both; t9l62 status line is blank

From your 07:13 status:
- **bsg:** dt 33.4 at c37550 (it was 60 at c21700).
- **M35ZAMS:** dt 0.943 at c66760 (it was 2.48 for the whole run).
- **t9l62:** the status shows no cycle / dt for log.4. Is it running? Post its last log line.

**Please:**
1. For bsg and m35zams_z3: post dt every 1000 cycles over the last 15k cycles, and the latest bsg3d.log / user.hst floor columns (eos_efloor, any floor counters), so we can see whether the drop was sudden or gradual.
2. Turn on `hydro/dt_report=10` for each at its next rst (name every key, as before) for ~100 cycles, and post the DTREP / DTREG binder lines (r, rho, T, v, dir, ceff, tauT_h).
3. **bsg:** if the binder is floor-density gas (rho near dfloor) at large r with large |v|, that's the c777 mechanism. Then run DFLFIX's H3 arm (NOTE 24zm, caltech-dflfix-5b6159ad0, from the newest bsg rst, 1000 cycles) beside it straight away, and report.
4. **M35ZAMS:** if the binder is floor gas at the cap (rho 1e-17, |v| = floorcap_v), say so. That's the dt_floor_courant item (USER decision).
Don't revert or switch anything on a dt drop alone.

-- orion (coordinator)
