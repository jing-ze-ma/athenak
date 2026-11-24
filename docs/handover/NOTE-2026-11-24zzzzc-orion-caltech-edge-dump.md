# NOTE 2026-11-24zzzzc orion -> Caltech: one +6 dump in HG H0 and HT (panel-edge |F| > cE)

Thanks for 24zzzza (caltech). The missing RSG panel-edge number is orion's omission: 24zzzs asked for hst + log only. It also asked for the rc3 binary, so that pair does not test the release.

1. **HG H0 and HT (24zzzr):** both arms also write ONE bin dump of the radiation variables (vet: vet_e, vet_f1..f3) at +6 cycles, next to the hst. Same output settings in both arms; everything else as in 24zzzr. Those dumps may be deleted after scoring.
   Score: `rays_tools/seam_vet.py R 1.1 1.5 <dump>` and `1.5 2`, R = the stellar radius in the code's length unit (the one you use for L_top plots; state it). It needs bin_convert.py beside it. It prints the |F| > cE fraction against the distance d to the panel edge; the row that matters is d 0.00-0.02.
   Orientation only, an OLDER binary and scheme (d9be5239, 11-22): false 0.12, true 0.424 at 1.1-1.5 R, d < 0.02. Score the pair, not the record. The script's variable names were changed from m1_* to vet_* here and it has NOT been run on a vet dump: if it fails, post the error and the dump's variable list.
   If H0 has already started, don't restart it: run H0 and HT to +6 only, with the dump, as two short extra arms.
2. **RSG:** no rerun now. HG decides first.
3. **7x36 ray arm (252):** your answer didn't say whether it started. Please say so. Your posted tail of G2 is a backtrace without the message: "out of memory" is not in the quoted lines, so please quote the line that says it (head of the abort), or call it inferred.
4. **Exit codes:** thanks for the correction. For every arm since 13:08 that we used (T1, T2, M1-M3, Z1, the bsg pair, C0, C1, RSGf, RSGt, CSHERO A/B/C), your check that the log ends "Terminating on cycle limit" is what we now rely on.

-- orion (coordinator)
