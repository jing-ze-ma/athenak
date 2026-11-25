# NOTE 2026-11-24zzzzb9 orion -> Caltech: what limits dt on the three stars (USER: fit 24 h on one node)

dt is constant at 12.007 (m35mms, then dips), 3.4805 (M35ZAMS), 96.85 (t9l62 quiet). We need to know which cell and which direction sets it. Short arms on idle nodes, copies of the decks, `hydro/dt_report` on (the existing key; a fresh deck needs it in the deck, a restart takes it on the command line), 25 cycles, hst + log only, binary 95813cf85 (md5 30d01765), T2 / c1 keys:
| arm | start | what |
|---|---|---|
| DZ | M35ZAMS c1 newest rst (or fresh) | DTREP lines: limiting cell (gid, kji, r, th, ph), which speed (radial / theta / phi, gas, rad-signal), h in each direction, share of cells within 5 % of the limit |
| DM | m35mms c1 newest rst (dt 12.0) AND the newest rst with dt < 11 if it exists (t > 5.4e4) | same, both |
| DT | t9l62 00015 (outburst) and the quiet stop rst | same |
Report: limiting direction and cell (is it a lateral r dtheta / r sin th dphi cell at the inner radius, a radial cell, floor gas, or the rad-signal speed), and what dt the SECOND limiting class would give (e.g. without the lateral limit: the radial-only dt). If dt_report prints only the worst cell, also give the min-dt per direction over the mesh if the log allows; say what it doesn't give.
Production has priority; after the HG repeat arms.

-- orion (coordinator)
