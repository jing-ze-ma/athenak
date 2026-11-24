# NOTE 2026-11-24zza orion -> Caltech: USER: clean_stag_mg_levels = 2 for ALL stars (reverses mg 4)

The USER decided (11-24 ~10:35 PST): `rad_vet/clean_stag_mg_levels=2` for every star (uniform settings). With 4: bsg +19 % and t9l62 +18 % per cycle; m35mms is ~20 % faster with 4, but stays inside 24 h with 2.

**Please, at each run's next rst** (no extra stop just for this; for t9l62 a clean stop now is worth it, since it's the slowest star):
- t9l62, m35mms, m35zams_z3, and bsg (if it hasn't reached tlim): resume with the same keys, but `rad_vet/clean_stag_mg_levels=2`. Name every key: the rst stores 4, so the key must be on the command line. Keep the pre-switch rst.
- No revert on a dt change. Revert only on CW / FATAL / inf / DT_RETRY.
- Post the 300-cycle watch for each: s/cyc before / after, dt, Krylov its per stage.
- RSG: keep as is (old-scheme hand-over; it doesn't name mg_levels; check its log echo and tell me what it runs).
- The orion deck dirs get mg 2 in the next sync; until then the command line wins.

-- orion (coordinator)
