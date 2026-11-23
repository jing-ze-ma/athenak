# NOTE 2026-11-23x Caltech -> orion + viper: IA11 m35mms bin cadence raised to 30 min (USER)

USER (11-23 ~14:15 PST, on Caltech): increase m35mms output cadence so the movie is smooth.

- The IA11 m35mms run was stopped cleanly and restarted from its own rst. Binary athena_ia11_bd9553537 and deck IA11_m35mms are unchanged.
- The restart passes on the command line only `output2/dt=1800 output3/dt=1800`: hydro_w and m1 bins every 30 min of model time (was 10800 s = 3 h). m1_vet / m1_face (43200), rst (21600), log and hst are unchanged.
- Restart: m35mms.00007.rst (md5 a5a8c142c8fe652c19cfbf7a979b0422), t ≈ 1.13e5 s, chunk 2 (log.2). It runs 0 FATAL / 0 CW, and dt ~7.8 continues smoothly. New bins are hydro_w/m1 00013, 00014, ... about 4 min of wall time apart.
- Cost: ~1.06 GB per frame, ~+500 GB over the run to tlim 8.64e5.
- An extra 14:18 restart (chunk 1, m35mms.00006.rst, ~7 min) ran WITHOUT the new dt: launch_ia11.sh did not forward extra keys. It is fixed now (shift 4 + "$@"; backup .pre1123x). The extra keys are recorded in the EX= field of chunks.txt.
- Viper: if you run m35mms anywhere, use the same output2/output3 dt 1800 for comparable movies.

-- Caltech
