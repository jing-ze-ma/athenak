# NOTE 2026-11-25g orion -> Caltech: bsg KE arm, and star arms of a PROTOTYPE (caltech-fwbmem-a497b040f)

Thanks for the four 300-cycle pairs (25b item 3 is complete; do not repeat). Priority stays: CUDA suite of 25c, then HGV. Around them, on nodes free of those:

**1. bsg KE_int (02-03 or 02-16, ~4 min).** rel-f differs from the running binary in KE_int by -3.0e-4 at the first row and about -1.5e-3 later; L, R and M agree. First return (zero cost) the `MLAT`, `NOTE` and `added from the command line` / default lines at start-up of F300 and C300 (grep -E "^MLAT|NOTE|merge_lat|fofc|default"), and the value of KE_int at row 0 of both from the hst. Then one arm: rel-f from the same rst + `hydro/merge_lat_blockend=flat`, 20 cycles (or the value C300 actually ran, if it printed another). KE_int rows 0-20 against F300 and C300. This tells whether the bottom merge is the cause; do not chase further.

**2. Prototype arms (not a release, changes results; for the USER's decision).** caltech-fwbmem-a497b040f = 79c1cabef + one commit (thin-layer flux relaxes to the ray flux). Build star_vet only on idle CPUs (same recipe, tiny gate: post the value, do not judge it), .prov checked. Arms from the rst of your 300-cycle pairs, base = your next build 9bebf6a3, 300 cycles each, same keys, dumps as before, one star per node: t9l62, then m35mms, then bsg. Same tables plus: cone-cut counts (`CLEAN FCONE` lines), clamp / CW / fallback / retry, mass loss at the top if hst has it, L_mid against L_midE. If a run stops or CW > 0: keep the log, no retry, report.

Stop all by 08:00 PST as in 25b.

-- orion (coordinator)
