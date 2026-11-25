# NOTE 2026-11-24zzzzf1 orion -> Caltech: bsg without merged cells (USER)

Estimate without the merge: dt 169 -> ~84 s (the inner wall phi cell, 4.97e9 cm), 58.6k cycles, 12-13.5 h for a fresh bsg. Check it, one node, 4 ranks, binary c4 a000d803c (md5 3ae130e9), no outputs, 100 cycles per arm:
| arm | start | change |
|---|---|---|
| NM0 | fresh, candidate_1124 bsg deck | the merge lines removed from the deck (`merge_lat_rc/_f/_mom/_rad/_diag/_blockend`), `dt_report = 10` in the deck |
| MG0 | fresh, same deck | as the production deck + `dt_report = 10` (control; compare with c1's rows) |
| NM1 | newest c1 bsg rst (a COPY), developed state | `hydro/merge_lat_f=1,1,1` (the rst stores merge_lat_rc, so the module stays on; factor 1 merges nothing: say what the log prints: "MLAT merged lateral cells"), all other stored keys named |
| MG1 | same copy | as production |
Report: dt over 100 cycles and DTREP binders, s/cyc, passes / its, CW / FATAL / DT_RETRY, floor_pay, L_top, KE and L_z rows vs the control at equal time; s per 1000 s of model time; nfb; for NM1 whether switching the merge off at a restart gives a start-up signal (first 10 cycles apart). Don't put any of this into production. Priority: after HR L-arms and the rel-b suite; before G / H.

-- orion (coordinator)
