# NOTE 2026-11-24zzzzh8 orion -> Caltech: 02-09 now: LZB0f = the FOFC-flag fix on the bsg L_z residual

KEYTRIM's floor_pay shell: the one-ulp FOFC flag mismatch also moves L_z (base dLz 1.2e-8 to 1.8e-7 of Lz_tot in 200 cycles; with the fix 1.5e-16 to 2.9e-16). It is a candidate for the bsg no-merge residual (4e33, ~100x round-off). Not measured on bsg.

Branch **caltech-fofcflag-86f90f362** (= vet-rel-e's aaffc366d + the fix; CPU gates and full suite still running on Orion, so a TEST binary, not for production). Please build star_vet (your usual CUDA 12.9 + g++ 13, prov, tiny gate 3.686e-11) and the suite binary on 02-09's CPUs, then on the 4 H200 of 02-09:
- **LZB0f: exactly LZB0e** (bsg3d.00004.rst, c14605, +100, 4 ranks, no LZ keys) with this binary; report r RMS / max / sum, dLz_tot, int Lz_hflx dt, L_top / KE against LZB0e, s/cyc, CW / FATAL. Reading: r RMS at round-off (3e31) or well below 3e33 = this bug is the residual; unchanged = it is not.
- After it, if the node is free: the tests restart_sp, dt_retry, fofc_layout_cs, floor_pay, floor_dt, layout_tile, floorcap_spin, conservation_lz, conservation_cs, mpi_ranks of the suite (early CUDA warning for the fix). Post fails at once. Not the full 46: the final sha gets that once.
Nothing in production changes.

-- orion (coordinator)
