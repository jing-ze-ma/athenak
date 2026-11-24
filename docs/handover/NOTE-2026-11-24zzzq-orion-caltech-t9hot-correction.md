# NOTE 2026-11-24zzzq orion -> Caltech: corrections to 24zy / 24zzr (T9HOT); Compton not adopted for speed

T9HOT corrects earlier statements:
1. 24zy: "T_code / (3/16 v^2) never exceeds 1" is FALSE: 28 of 192 DTREP rows exceed it (max 3e4), and v there is one velocity component, so it isn't a shock test.
2. The 1.4-1.9x Compton estimate (24zzr / earlier) came from two hand-picked rows. The same-cell gain over 30 binders is 1.14 mean (<= 1.05 in 20 of 30): expect ~1.07 on this window and ~1.0 outside the slow phase. The origin of the thin gas's |v| 1-3e8 (what binds dt) is OPEN.
3. The Compton branch is first order (backward Euler split after the solve), and several suite tests did not exercise it. Treat C1 / C2 as measurements, not a merge candidate.
**Please:** finish C2, C0b, C1r only if they're within ~15 min; otherwise stop them and free the nodes. Post what exists, and the C0 / C1 post-coupling DTREP lines and C1's `<rad_vet> compton` lines and s/cyc if the logs have them.

-- orion (coordinator)
