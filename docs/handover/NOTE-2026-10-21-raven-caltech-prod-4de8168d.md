# NOTE 2026-10-21 raven -> caltech: GO — switch RSG E 5126654 to 4de8168d + the clip deck; accretor still waiting (v12lzc NO-GO)

fin51-r-fix = 4de8168d (pushed; all gates clean; IDORT cbg + cbh PASS) = 3896d529 + memcut (cs vet_gd band ring ON by default
where it applies + 2-bit read mask; BITWISE) + cs2clip (opt-in <rad_m1>/implicit_cs2_clip = rad: the F-limiter clip is not booked
as conservative2 gas momentum; fixes floor-gas kicks) + the cbg/cbh fallback fix (armed dt_retry_t2_be attempt: stage-1
NON-CONVERGED no longer skips the BE redo; default / no-fallback paths BITWISE).
Raven (R2, HG global R) and Viper (V2, HG global V) fresh chains now run 4de8168d + the *_clip decks. Nothing in production has started.

## RSG E (5126654, not started) — do now
1. run.cfg BIN -> your build builds_fin51/bin/athena_gpu_he_star_m1_4de8168d_nofma, BIN_MD5 9a56b4a109c575b7aa87a9a94d0124eb
   (no rebuild needed).
2. Deck -> docs/handover/rsg-E-1020/rsg_E_start1020r_clip.athinput.in (md5 4c442452, pushed in 23415519) with ONLY the 3 @FILES@
   lines substituted as before. It equals the committed rsg_E_start1020r.athinput.in + one <rad_m1> line `implicit_cs2_clip = rad`.
   That key is FATAL (unknown key) on a pre-4de8168d binary, so switch BIN and IN together.
3. Layout unchanged (2 nodes x 4 H200, FRESH t=0, tlim 1e9, same STOP rules). No new check (USER: no unneeded tests).
   The rdiff triplet 5140304 sources the same run.cfg; fine.
4. Post a NOTE with the run.cfg md5s; tell us when it starts (race unchanged).

## Accretor — still waiting
The v12lzc screening on Raven was NO-GO. Do not queue an accretor chain; no V12_READY yet. Raven run12 also stays held.
A NOTE with the fixed deck/binary follows once it passes a new screen. Keep the 4de8168d accretor build (0ff5fc54).
