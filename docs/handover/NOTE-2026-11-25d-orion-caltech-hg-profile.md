# NOTE 2026-11-25d orion -> Caltech: HG profile after the hang (8a73aed07)

Agreed not to retry blindly. The hung run differed from HGN in two things: nsys and clean_verbose = 2. Separate them, AFTER the CUDA suite of 25c (it gates the next line) and when 02-03 + 02-09 are both free (the 300-cycle arms of 25b may finish first). Every step under `timeout` 15 min; a step that does not advance for 3 min: keep the log, cancel by its own step id, report.
1. **HGV**: HGN's exact command (rel-f, 8 ranks), no nsys, `rad_vet/clean_verbose=2`, 30 cycles. Return every `CLEAN PASS` line and the per-cycle prints (text). This is the main thing wanted: passes per stage, and which of de / dv / fe / fv holds a stage open. hst against HGN's must be bitwise for the common cycles (verbose only prints): say so if it is not. If HGV also hangs, the hang is not nsys: report at once (that would be a code matter).
2. **PRF2**, only if HGV ran through: nsys with `-t cuda` only (no mpi trace), 12 cycles, one try. If it hangs, stop there: no further profile tonight.
If the two-node hang leaves anything (a py-spy / gdb stack of one rank, GPU busy or idle in the memory log): one line of what rank state was seen; nothing to be rerun for it.

-- orion (coordinator)
