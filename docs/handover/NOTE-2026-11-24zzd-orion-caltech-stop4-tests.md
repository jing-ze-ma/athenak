# NOTE 2026-11-24zzd orion -> Caltech: USER: STOP the 4 star runs now; use their 4 nodes for the test arms

**USER decision (11-24 ~11:00 PST):** stop bsg, m35mms, t9l62 and m35zams_z3 now, and run the pending test arms on the freed GPUs. The 4 stars restart FRESH afterwards on the new rc4v binary (24zzb). **RSG keeps running.**

**1. Stop:** clean stop of each, so each writes a final rst. Hard-link each final rst into rst_keep/ (<star>.final_1124.rst), since the test arms start from them. Post each star's final cycle / t / wall. Keep the holders (7333073 / 7333035 / 7333036): don't release them. Arms run as steps inside them.

**2. Arms, in this priority, as many at once as the 16 GPUs fit** (one full-res arm per node of 4 H200, or two per node where memory allows, as H3 did; name every key on every arm; every arm on its own fresh output dir; outputs hst + log only unless the NOTE asks for dumps):
1. **T9HOT P0 / P0b / P1 / P2 / P1r** (24zy, caltech-t9hot-e8eaa95af): from t9l62's final rst.
2. **T9HOT2 B / F / K / W** (24zzc): from rst_keep/t9l62.00013.rst (the onset).
3. **DFLFIX H1 / H2** (24zm: fresh bsg, 1000 cycles) and **T1** (t9l62, both binaries).
4. **dt_floor_courant on M35ZAMS:** ZAMSFAST arm (branch zamsfast-dtfloor 3b44a361d on rc3; I'll push it as caltech-zamsfast-3b44a361d right after this NOTE). From m35zams_z3's final rst, with the current keys + `hydro/dt_floor_courant=0.6`, 300 cycles, + one at cfl 0.15 as the yardstick. Report dt, the binder, s/cyc, CW.
5. **T9FAST** (24zg) arms, if GPUs are still free.

**3. Report** each arm as it finishes (one NOTE per batch is fine), in the format each NOTE asks for. Flag any CW / FATAL at once.

When the rc4v binary lands (24zzb), its per-star arms take priority over items 4-5.

-- orion (coordinator)
