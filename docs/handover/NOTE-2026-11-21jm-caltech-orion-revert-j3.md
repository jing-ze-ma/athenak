# NOTE 2026-11-21jm Caltech -> orion: head at cfl 0.6 also fails in c33370-33460 -> reverted t9l62 to B_j3 from 00045 (ll item 4)

Head 462210d51 at cfl 0.6 (relaunched 23:27 from 00045, jk keys + cfl 0.6): c33433 t2s 2 stage fallback 1; c33434 BE not converged, dt_retry x2 (last change 8.3e-4, 2.5e-7); c33455 stage fallback 2; CW 8, dt_retry 3, no "last change inf"; hydro dt COLLAPSE x4, dt 41 -> 15.5 -> 0.38 at c33500; ~5 s/cycle. Stopped ~23:33 via STOP (wrote t9l62.00046, moved to t9l62/rst_head06/).
Reverted per ll item 4 (USER "go"): B_j3 athena_j3_a5257b7ee (md5 f71a22b6) from t9l62.00045 with the rst's own keys + output5/dt=21600, relaunched ~23:34 on testH1, END 00:20 (B_j3 has no graceful stop; last 6-h rst survives the 00:34 holder end). B_j3 got through this stretch earlier tonight with CW 5 / dt_retry 3 / FATAL 0.
So on this state the head fails at both cfl 0.9 and 0.6 in the same cycles where B_j3 only stumbles.

-- Caltech
