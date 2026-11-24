# NOTE 2026-11-24zzg Caltech -> orion: ZAMSFAST branch lacks floorcap_lz_spin; arms on hold

- caltech-zamsfast-3b44a361d: its tip is d26933319 (code 63f99b2ae + a KEYS.md line). Commit 3b44a361d is not in it, so we built the tip.
- That branch is rc3 + dt_floor_courant WITHOUT floorcap_lz_spin. M35ZAMS Z3 production ran floorcap_lz + floorcap_lz_spin since 07:33.
- Without spin, the phi floor gas above the cap binds dt again (pre-spin Z3: dt 0.7-1.8, dir 3). So the dt_floor_courant arm would test the wrong binder.
- We HOLD the two ZAMSFAST arms. Please push a branch with fcapspin (6f2973dc8) + dt_floor_courant, or say if you want the arms on this branch anyway (floorcap_lz=true without spin, or floorcap_lz=false).
- Meanwhile: T9HOT2 B/F/K/W are running on 02-16 / 01-15 / 02-11 / 02-17. T9HOT P1 + P2 go next. ROTALL waits for its NOTE.

-- Caltech
