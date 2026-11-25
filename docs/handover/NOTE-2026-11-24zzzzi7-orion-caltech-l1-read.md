# NOTE 2026-11-24zzzzi7 orion -> Caltech: L1 read; one more arm (fix binary, zero-flux off)

L1 (fofc off): r RMS 1.5e31 = round-off, and the L_z through the top falls from 2.9e35 to 1.4e32 per cycle (mass flux 8e12 g in 100 cycles in both): the L_z crossing the top is the FOFC fluxes of the floor gas with (nearly) no mass. Together with LZB0f (1.78e32, steady +): the FOFC flux at the top is the place. Hypothesis: the top take-back (star_zflux) acts on faces with `!(mass flux > 0)`: a face with exactly 0 mass flux but nonzero transverse momentum flux is taken back, and LZLEAK showed that take-back leaves a 2e-3 L_z error (Min_top is then 0, as you see).
**Arm LZB0f11, on 02-09, right after L2 ends (L4 stays cancelled; L2 only if it has started): the FIX binary 6a831ffe (86f90f362), LZB0f's args + problem/star_wall_zero_flux=false, +100.** Report r with the END-row rate (RMS, sum), gross Lz_hflx per cycle, M_tot, L_top, KE vs LZB0f. Expectation if right: r falls from 1.8e32 to ~2e31. LZLEAK's fix of the take-back (own branch) will follow for a CUDA check.

-- orion (coordinator)
