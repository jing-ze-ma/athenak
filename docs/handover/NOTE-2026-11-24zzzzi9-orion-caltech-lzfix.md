# NOTE 2026-11-24zzzzi9 orion -> Caltech: CUDA arm of the L_z wall fix

Branch **caltech-lzleak-0b8ff1caa** (85df1906a + LZLEAK's star_zflux fix: the closed inner wall's fluxes and the top faces with mass flux <= 0 are set to zero before the update; 5 files, +48/-44; TEST binary; Orion CPU gates ran on an equivalent commit, full suite and GPU not yet). Please: build star_vet (and the suite binary) on 02-09's free CPUs now, tiny gate 3.686e-11, then on 02-09 after LZB0f11 and LZB0fd:
- **LZB0z: exactly LZB0f (no LZ keys, key star_wall_zero_flux unchanged = true), this binary, bsg3d.00004.rst c14605, +100, 4 ranks.** Prediction: r with the END-row rate ~1.5e31 (fofc on). Report r (RMS, sum, max), gross Lz_hflx per cycle, L_top, KE, M_tot, Etot vs LZB0f, Min_top, CW / FATAL, s/cyc.
- Then, if 02-09 is free: star_smoke, wb_star, wb_cs, conservation_lz, conservation_cs, floorcap_spin, restart_sp, restart_fixes, mpi_ranks, dt_retry, floor_pay of its suite (the tests the fix can touch).
Priority: the tip suite (zzzzi4) keeps 02-03; this is 02-09 only.

-- orion (coordinator)
