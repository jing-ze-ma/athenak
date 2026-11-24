# NOTE 2026-11-24zzzh orion -> Caltech: USER: stop RSG and the 4 stars; use the 7 nodes for a short test batch; stars restart FRESH afterwards

USER (11-24 ~13:10 PST): "you can stop the rsg and the 4 stars to make room for the tests". Nothing is restarted yet; the 4 stars will start fresh later on settings the USER picks after these tests.

**1. Stop** (clean; keep the final rsts, hard-linked <run>.stop1124b.rst): t9l62, m35mms, m35zams, RSG (bsg is already stopped). Post final cycle / t / wall for each. Keep all holders. Don't start the RSG 01809 run yet: hold it (the copied file stays); it needs the 3 nodes for tests first.

**2. Tests, all at once on the 7 nodes (28 H200). T2 binary 7d7c4a71 unless noted. Every key named. Hst + log only (plus what each asks). ≤ 300 cycles each. Stop rule CW / FATAL / inf / DT_RETRY. Report at each end.**
| node | test | from |
|---|---|---|
| 02-11 (bsg) | Compton arms C0 / C0b / C1 / C2 / C1r (24zzs, as running) | t9l62.00015 |
| 02-17 (m35mms) | M1: production keys, `rad_vet/clean_verbose=2`, 100 cycles; M2: M1 + `clean_stag_mg_levels=4`, 100 cycles; M3: M2 + `rad_vet/clean_stag_thinf=true`, 300 cycles | m35mms.00005 (c6278; or the final rst) |
| 02-16 (t9l62) | T1: production + verbose, 50 cycles; T2: + `clean_stag_thinf=true`, 300 cycles (L_top / L_in, dt, its, s/cyc against T1's 300) | t9l62's final rst |
| 01-15 (M35ZAMS) | Z1: thinf true vs false, 300 cycles each, same start | m35zams final rst |
| 01-17, 02-03, 02-09 (RSG's) | CSHERO A (K4), B (prod), C2 (fixed2 deck): 4 GPUs each, `time/nlim=0` first (24zzn, fixed in 24zzzg). If C2 FATALs, skip it. B first; A and C2 only if B keeps the X at cycle 600 | fresh |
**BSG thinf:** use bsg's own T2 fresh start (c0-300, in hand: 0.687 s/cyc, dt 165-169, L_top/L_in 1.064) as the baseline and run fresh bsg with `clean_stag_thinf=true` for 300 cycles on 02-11 after Compton; 8 cycles + restart bitwise is not needed.
**3. Not requested:** PE/PL/PO profiles, T9PASS, mg_gf, cfl ladders.
**4. Report per test:** s/cyc, dt, solves and its per stage, passes, CW / FATAL / inf / DT_RETRY, L_top / L_bot / L_in; for verbose arms the CLEAN PASS lines; for thinf: dt and s/cyc against the pair, L_top drift. Post each as it finishes; a partial table is fine.

-- orion (coordinator)
