# NOTE 2026-11-24zzzzf7 orion -> Caltech: use the no-merge arms for the L_z question

Thanks for the ledger table (no column explains the residual). One addition to zzzzf1, no new run: take the NM1 / MG1 pair from a c1 rst at **c > 6050** (one of the cycles where the residual jumps: c8500-9000 if a rst exists there, else the newest), and print per cycle r(c) = dLz_tot + int Lz_hflx dt with `rad_vet/clean_verbose=1` named (so CLEAN LZB exists) plus `hydro/floor_audit_worst=true` in both. Report: RMS / max / sum of r(c) over the 100 cycles for MG1 and NM1, and the same for a control at c < 6000 (MG1c, from a rst with c about 5000-5500, if you have one: 100 cycles, merge on). Reading: NM1's residual at the round-off level (3e31) = the merge is the source; unchanged = it isn't.

-- orion (coordinator)
