# NOTE 2026-11-24zzzzd7 orion -> Caltech: HG and RSG restarts with the ray-hole cure (merged code fixes)

The cure (a000d803c) is a blend close to the old `clean_stag_thinf_held`, which the USER dropped on 11-21. It has passed stars and the suite but **never an RSG or HG restart** (HTX = THINFCS's fix without the cure). One pair each on the CUDA build of `caltech-rel-a-c880e2703` (zzzzd6), thinf = true, time2_gtime_slope as in HTX:
| arm | start | args | cycles |
|---|---|---|---|
| HR | hegiant.00205 (as HTX) | HTX's args word for word, binary rel-a (he_star_vet) | +51 |
| RR0 | RSG newest production rst you hold for the RSG clean run (or rsg.01809 if present; say which, N0) | its own restart keys, thinf named true, vceil 3e8, vet_gd_allow_clamp = false, cfl 0.3, conv 1e-5, every stored key named | +51 |
| RR1 | same | RR0 with thinf = false (control) | +51 |
HR is scored against HTX (zzzzb1: e_floor sum, pos_gas, floor_pay, panel-edge |F|>cE at 1.1-1.5 and 1.5-2 R, ef<=0, L_top / L_mid / L_int / L_bot vs H0, dt, CW / FATAL / BE / fallback); the same table for RR0 / RR1. **Report CW first; a single CLEAN WARNING in HR or RR0 with none in HTX / RR1 is a FAIL and stops the series.** If the RSG rst or keys are not available, say so and skip RR.
Priority: after zzzzd6 (its build); ahead of zzzzd1's speed arms.

-- orion (coordinator)
