# NOTE 2026-11-24zzzzl orion -> Caltech: uniform-settings pairs on the 4 sp stars (USER), on the nodes HG does not use

Binary: **T2** athena_ia123rc4v_262affd6f (md5 7d7c4a71), not the sdirk2 TEST binary. HG (02-09 + 02-03) has priority and is not disturbed. Use 02-11, 02-16, 02-17, 01-15, 01-17.

Two arms per star, same rst and args; the ONLY difference is three keys, both named in both arms:
- UA (as now): `rad_vet/clean_stag_thinf=false rad_vet/clean_stag_lat_sums=false rad_vet/clean_stag_tie_gid=false`
- UB (uniform): the same three `=true`
Both arms also name the USER's new uniform settings `rad_vet/vet_gd_allow_clamp=false hydro/vceil=3.0e8` (a clamped read is then FATAL: post it, do not turn the key on), `rad_vet/clean_stag_mg_levels` 4 for m35mms and 2 for the others, `time/cfl_number=0.3`, and every other production key as in the last restart args. (`job/final_restart=false`.) Restart only: the three keys are stored false in the rst, so they must be on the command line.

| pair | start | cycles |
|---|---|---|
| U1 bsg | bsg3d.pre-compton.rst (S1's) | 300 |
| U2 M35ZAMS | stop rst | 300 |
| U3 m35mms | stop1124b / m35mms.00006, mg 4 | 300 |
| U4 t9l62 outburst | t9l62.00015 | 300 |
(The quiet t9l62 pair is dropped: the outburst is the harder case.) First print the log lines that name the three keys in each arm (they must read true / false as intended). Order: U4, U3, U1, U2, one node each; ~3-16 min per arm. hst every cycle, log + CLEAN lines shipped; no rst, no dump.

Report per pair: s/cyc, dt (min / mean / max), its and passes, CW / FATAL / inf / DT_RETRY / BE / fallback counts (with the final summary lines of the log), L_top / L_bot / L_in and KE difference at matched time (end and max), M / E / Lz drift, log end = cycle limit + the run's own exit status.
HARMLESS, fixed before the run: no extra CW / FATAL / fallback / BE / retry / clamped read; dt mean within 1 % of UA's; seconds per simulated second <= 1.05 x UA; |L_top diff| <= 0.5 % (end and max); M / E drift not larger than 2 x UA's. (U4 diverges: thinf alone gave L_top up to 0.35 % there, so a 0.5 % bound is orientation only; report it, orion reads U4 with that in mind.)

-- orion (coordinator)
