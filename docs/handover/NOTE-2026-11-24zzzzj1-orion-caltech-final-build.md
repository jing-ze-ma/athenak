# NOTE 2026-11-24zzzzj1 orion -> Caltech: pre-build the release candidate (renamed) for its one CUDA pass

Branch **caltech-rel-f-b09a7b3c1** = all fixes of 03e874599 + the L_z wall fix (0b8ff1caa) + the Cartesian floor_from_ekin fix + docs, renamed M1 -> VET (code text only). Its CPU suite runs on Orion now; if that finds a fail I post it and the sha changes.
1. Build star_vet and the suite binary of it NOW on the free CPUs of 02-09 (after / beside the caltech-lzleak build, which is for LZB0z); gate deck as for rel-e.
2. The full CUDA suite (its README's count; 4 lists + mpi_ranks with Slurm unset) on the first of 02-03 / 02-09 that is free AFTER the Orion CPU suite is green (I post "GO"); LZB0z on 02-09 comes first (~2 min). The suite on 03e874599 (02-03) is an early warning: let it finish.
3. Report as usual: fails at once; star_smoke, floor_pay, conservation_lz / _cs, floorcap_spin, restart_fixes, mpi_ranks, eddington_cs, mhd_equil values against the tip's.

-- orion (coordinator)
