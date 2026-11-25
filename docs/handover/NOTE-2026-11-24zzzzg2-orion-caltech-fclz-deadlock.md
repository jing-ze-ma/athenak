# NOTE 2026-11-24zzzzg2 orion -> Caltech: STOP: clean_stag_fcone_lz deadlocks on 2 or more ranks

KEYTRIM found: CleanStagFCone (rad_vet_clean.cpp, present since 286e27bb7, in 95813cf85 / c4 / rel-a / rel-b) calls an MPI_Allreduce inside a rank-0-only print when `clean_stag_fcone_lz = true`. Rank 0 waits in the Allreduce, the others in a halo exchange: the run hangs after the first cycles. So **do NOT run LZ1 or LZB1 of zzzzf4 (or any arm naming this key) on the present binaries**; if one is running and hung, kill your own step and post it. LZ0 / LZB0 stay valid as controls. The fix is being gated on Orion (vet-fclz-mpi, one commit); a new binary follows with a new NOTE. Nothing else changes.

-- orion (coordinator)
