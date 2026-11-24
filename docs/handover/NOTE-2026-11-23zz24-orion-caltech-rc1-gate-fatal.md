# NOTE 2026-11-23zz24 orion -> Caltech: rc1 tiny-gate FATAL is expected (the gate deck is θ-periodic); rerun as below

The tiny gate deck (H200CHK/bsg_wngl_smoke_dtfac1e6.athinput) is θ-periodic (x2 π/3..2π/3, ix2_bc periodic). rc1 now FATALs on a θ-periodic ray set unless the seam keys are named (the USER's wedge rule). It also names keys removed in the release (implicit_pos_floor_solve, cs_cross_mono).

1. Please post the FATAL line itself (grep "FATAL" in the gate log) so we can confirm it's that.
2. Rerun the gate with: `rad_m1/vet_gd_theta_seam_sym=true rad_m1/vet_gd_wall_sym=true` on the command line.
   - If it FATALs on a removed key, drop that key from the gate deck copy.
   - The seam key changes the steady init slightly, so pass 1 may not be 3.686e-11. Report pass 0 / 1 / 2.
   - We'll compute rc1's reference on CPU and post it. Meanwhile, start WD and K0 if pass 1 is within 10 % of 3.686e-11 and pass 2 ≤ 1e-13; otherwise wait for our reference.
3. Use the same binary for WD / K0 (no rebuild).

-- orion (coordinator)
