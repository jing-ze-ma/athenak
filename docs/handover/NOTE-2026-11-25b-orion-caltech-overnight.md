# NOTE 2026-11-25b orion -> Caltech: overnight: keep the idle test GPUs busy (USER asleep)

The USER asks that idle GPUs do useful tests overnight. Rules as before: production stars untouched (stop / switch / restart only on the USER's word), no pushes by you except your result posts, no new production, nothing known-bad. Order of use for 02-03, 02-09, 02-16 (and any other idle holder node):
1. The HG profile of 25a (HGN nsys window, HGO is done): finish it.
2. As soon as I post a branch (next line, speed or memory branches: each NOTE names it), build star_vet + suite on free CPUs at once and run what the NOTE asks, first come first served; a node idle for more than 15 minutes without an order from me: say so in the node table.
3. If nothing is posted: the long-window check of the release on real states that is worth having (read-only copies, no production): for m35mms, bsg, t9l62, m35zams from their newest rst: 300 more cycles with rel-f (b09a7b3c1) and the running binary from the SAME rst (L_top, radius by rph.py at cycles 100, 200, 300, dt, s/cyc, its, CW / fallback / retry), one star per node; post one table per star. (This is the "longer window" for the continuity check; do not repeat it on a star already done.)
4. The bsg node (02-16, finished star): the 3rd arm of any of the above.
Post the node table as usual; stop all test arms at 08:00 PST unless told otherwise.

-- orion (coordinator)
