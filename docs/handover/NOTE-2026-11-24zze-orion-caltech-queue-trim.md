# NOTE 2026-11-24zze orion -> Caltech: trimmed test queue (USER: no unnecessary tests). This replaces 24zzd's list

The critical path is the fresh restart of the 4 stars on ONE binary (rc4v, being gated on Orion; push in a few hours). Only arms that change a decision for that start stay:

**KEEP, in this order (run concurrently on the 16 GPUs):**
1. **T9HOT2 bisect B / F / K / W** (24zzc), from t9l62.00013: it decides whether the fresh t9l62 needs a fix first.
2. **T9HOT P1 + P2 only** (24zy): key on, and key on at cfl 0.15 (accuracy yardstick), from t9l62's final rst. Drop P0b and P1r: the one binary's restart / key-off gates cover them. Use P0 = production's own last 300 cycles (no new run).
3. **ZAMSFAST dt_floor_courant arm + its cfl 0.15 yardstick** (24zzd item 4), from m35zams_z3's final rst.
4. **ROTALL rotation arms:** NOTE follows within ~1 h. They decide fixed vs rotating rays for the fresh start.

**DROP (answered already, or obsolete keys):**
- DFLFIX H1 / H2 / T1: H3 + the bsg production switch (dt 164, clean to 85 %) answered it. T1: t9l62's floor cells are 1e11x above dfloor, so the limiter can't act there.
- T9FAST arms (24zg): they're on mg 4 / the rc3 key set, now superseded (mg 2, rss_phase, dt_floor_courant). Redo later on the one binary only if t9l62 still misses 24 h.

**Idle GPUs:** if all of 1-4 run and GPUs are still free, leave them free for the one binary's per-star arms (priority 0 once it lands): fresh-start smoke + 300 cycles per star.

-- orion (coordinator)
