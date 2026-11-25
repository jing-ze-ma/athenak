# NOTE 2026-11-24zzzzc6 orion -> Caltech: cubed-sphere interior pattern (CSHERO), and c4 at the next restart of the other stars

## 1. CSHERO arms (USER: fix the cs pattern; his candidate: radiation not using one flux at the seam)
Arm D (cs_vertex_3cell) was inert: that key is not read by the radiation or the main hydro update; my suggestion was wrong. Candidates with evidence: (a) the cs MeshBlock-face radiation flux: on Orion's small cs shell the tie moved E by 0.04-0.52 of E with thinf true and 0.03-0.04 with thinf FALSE (two flux values at lateral block faces in both; THINFCS's fix makes it one value but only under thinf = true), (b) merge_lat inside 82 Rsun (if the deck has it; the pattern starts in the innermost bins and spreads out).
CSHERO B's deck, 600 cycles, same node layout, imprint at c150 / 300 / 450 / 600 with cshero_imprint_phase.py (m4 and phase 200-400, sym_rms, 0-100, 100-200, 800-850), hst L_top, KE, CW / FATAL:
| arm | binary | change |
|---|---|---|
| Bp | 53290b5af (caltech-sdosdo) | B's deck unchanged (control on this binary; B itself was T2) |
| E | 87af1d65c (caltech-thinfcs) | thinf = true |
| G | 53290b5af | thinf = true (no fix) |
| H | 53290b5af | Bp without merge_lat (remove its keys; if the deck has none, say so and skip H) |
Run Bp first; if its m4 at c600 differs from B's (1.39e-4) by more than 10 %, say so and go on. Say what thinf is in B's deck. Both keys (thinf, fix) stay command-line / deck changes; nothing else moves. Idle nodes, after c4-related items, production first. Score at c600 (m4(200-400) ratio to Bp) and, if E or H drops it by more than a factor 3, keep that arm's dump c600 and stop the rest.

## 2. c4 for bsg, m35mms, M35ZAMS (USER yes)
Their c4 preflights are identical to c1 for 100 cycles, so **do not restart them fresh and do not stop them.** At the NEXT chain restart that happens anyway (wall limit or a crash), use athena_c4_a000d803c (md5 3ae130e9) instead of 95813cf85 (md5 30d01765), same args and deck keys, every key named; check that the first 3 restart cycles have CW 0, no fallback, and that a rst from 95813cf85 restarts on c4 with the hst continuing (report the first difference if any). If no restart comes before the run ends, nothing is done.

-- orion (coordinator)
