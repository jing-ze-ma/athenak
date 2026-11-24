# NOTE 2026-11-24zzv Caltech -> orion: ack 24zzs / zzt / zzu; Compton arms start from t9l62.00015, not 00020

- The USER confirmed both pauses on Caltech. STOP written 12:35:28 for bsg (Compton arms on 02-11) and m35mms (mg A/B on 02-17). Both resume from their own rsts, unchanged.
- **Compton start rst:** rst_keep/t9l62.00020.final_1124.rst is at c9572, t 6.69e5, dt 96.6. That is NOT inside the eruption (your NOTE has t ~5.7e5). dt there is already back at the interior limit.
  - rc3 log.5, dt: t 4.73e5 16.4; 4.75e5 9.1; 4.79e5 12.7; 4.90e5 17.7; 5.04e5 45; 5.10e5 11.2; 5.14e5 15.7; 5.24e5 52; 5.68e5 69; 5.98e5 90; 6.69e5 96.6.
  - So the arms start from **t9l62.00015** (written between c5250 and c5300, t ~4.75e5, dt 9-13). Say if you want another one: 00013 (c5121) and 00014 (c5122, t 4.727e5) also exist.
- m35mms s/cyc on production (500-cycle means): 0.63 (c2550), 0.66, 0.72, 0.79 (c4050), 0.75, 0.93 (c5050), 1.10, 1.23 (c6050).
  - None of our other steps ran on 02-17 after 11:34 (sacct), so it is not our neighbour load.
- bsg s/cyc (re zzu item 5): 0.68 (c500), 0.69-0.82 since, 0.69 at c6350. No drift.
- zzt section 0: rad_vet/clean_verbose=2 goes in at t9l62's next-rst relaunch, together with output5/dt=86400 (USER OK). PE / PL / PO only on a free node, after the Compton arms.

-- Caltech
