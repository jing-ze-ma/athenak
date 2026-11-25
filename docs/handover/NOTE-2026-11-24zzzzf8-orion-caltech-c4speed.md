# NOTE 2026-11-24zzzzf8 orion -> Caltech: c4 + speed (A1+A2): gates, then the switch of t9l62, m35mms, M35ZAMS (USER go)

Branch `caltech-c4speed-d979d5268` = c4 a000d803c + A1+A2 only (no key, no K1, no T4). Orion gates (258 files bitwise, restarts, suite) running. **Do now, on idle nodes, in this order; no star is touched before I post the go:**
1. Build (CUDA 12.9 + g++ 13, he_star_vet, own dir, gate 3.686e-11, md5 / .prov).
2. Identity: t9l62 (a COPY of a rst of the running c4 run, or the quiet stop rst), m35mms late, M35ZAMS: 100 cycles each on c4 and on c4speed, keys as production: hst printed-identical (user + hydro), CW / FATAL / BE / fallback 0, s/cyc. Any difference: first cycle and column, stop.
3. Continuation: for each of the three, from a production rst (a copy): c4 (or 95813cf85 for m35mms / M35ZAMS, whose runs are on 95813cf85) 5 cycles and c4speed 5 cycles from the same rst: printed-identical for c4; **for the two on 95813cf85 compare against 95813cf85** (their cure is inert: identical printed hst expected, as in your c4 preflights).
4. When 1-3 pass, I post the switch order (one star after another, right after a fresh rst, clean stop, restart on c4speed, first 5 rows and the log's key lines checked; bsg is NOT switched).
Priority: after HRL / HR150; ahead of K1.

-- orion (coordinator)
