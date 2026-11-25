# NOTE 2026-11-24zzzzh6 orion -> Caltech: residuals received; MG3 (release code, radial merge) on 02-03

Thank you: the residual table is what we needed (E and M at round-off in all five production runs, key on, injection-checked). mpi_ranks 46/46 noted.

**MG3 on 02-03 (free):** the bsg jump window again, identical start, args and diagnostics as MG2 (rst bsg3d.00002, t 1.2099e6 to 1.50e6, 4 ranks, clean_verbose = 1, floor_audit_worst, same cap), but binary rel-e 08c48ac10 star_vet (md5 8737e5b0, your CUDA build) and `rad_vet/merge_lat_mom = radial` named explicitly. Differences from MG2 that come with the binary must be listed from the log (the NOTE lines of the restart rule; name on the command line every key whose stored value you keep). Report as for MG2 / NM2: r RMS / max / sum per window, the t 1.400-1.497e6 window, 10 largest cycles, L_top / KE / Lz_tot against NM2 at equal time, s/cyc, CW / FATAL / DT_RETRY / fallback.
Reading: NM2 4e33 and MG2 4e37; MG3 near NM2 (within a factor 3) = the radial merge cures the L_z error; near MG2 = it does not. No star is touched; production bsg stays as it is.
The LZ arms on 02-09 continue; the order of arms is yours.

-- orion (coordinator)
