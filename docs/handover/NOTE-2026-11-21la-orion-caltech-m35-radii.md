# NOTE 2026-11-21la orion -> caltech: m35mms lattice, shell extracts at more radii (CPU only, existing dumps, no new runs)

Thanks for CHK/m35mms_*.npz. Result (Orion CHECKER): the 0.75 R_ph pattern is NOT a 2-cell checkerboard (zero power at k >= 120).
It is a ~10-cell lattice (k ~ 24 per 256), aligned with the theta/phi grid axes. It is a v_r-dominated g-mode that grows from 1.5 d.
Isolated seeded modes on Orion/viper do not grow, so we need to see where it is driven.

Please extract, from the existing m35mms dumps 00012, 00016, 00020 and 00024 (CPU only):
- radial cells nearest r / R_ph,IC = 0.71 (first cell above the inner wall), 0.72, 0.74, 0.75, 0.78, 0.80, 0.82, 0.85;
- for each of those cells, all 256 x 256 (theta, phi) values of rho, v_r, v_theta, v_phi, e_int, E_rad, F_r, F_theta, F_phi (float32 is enough);
- the same format as CHK: one npz per dump -> caltech_results/CHK/m35mms_radii_<dump>.npz, arrays [radius, theta, phi] + r + time.

If you have a free H200 test GPU (not production), also say so. A full-width repro restarted from dump 00012 is the next step.

-- orion (for the USER)
