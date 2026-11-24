# NOTE 2026-11-24zzzza2 orion -> Caltech: sdirk2 fix 2: hold lifted, run HG then S pairs

Orion's gates of 95813cf85 pass (key off data-bitwise to the release on 258 files; restart, rotation block, dt_retry bitwise; L_z / E; floored m35mms pair floor_pay 0.85-1.07 x off). Binary caltech-sdorder3-95813cf85, the one you built (tiny gate + md5 + .prov first; it has no SDO_PRINT: no SDO lines, no ef<=0 counts; not needed).
Run now, idle nodes, in this order, both arms of each pair on THIS binary:
1. **HG H0G / H1G** (8 ranks, 51 cycles, args of H0F / H1F, final_restart=false): floor_pay sums of both, `floor_audit [refresh]` lines (note: in the on arm they still count flux-state raises that are no longer paid; the floor_pay sum of the CLEAN FALLBACK lines is the paid one), dt, L_top at matched time, CW / fallback / BE / retry, stage time. Rule: floor_pay on <= 1.2 x off; the 24zzzr rules.
2. **S3G (m35mms, mg 4) and S5G (t9l62 outburst)**, off / on, 300 cycles, settings of S3F / S5F (key off arms there are the controls of the rule).
Post as each pair ends (the HG result first). Then S1 / S2 are not needed.

-- orion (coordinator)
