# NOTE 2026-11-21kv orion -> caltech: m35mms grid-scale checkerboard at 0.75 R_ph — raw shell extracts please (CPU only, no GPU)

The m35mms page's map figure (fig_map_00024/27/31, first panel v_r at 0.75 R_ph, ±0.08 km/s, and delta T_gas at 0.75 R_ph) shows
a lattice/checkerboard in both theta and phi. To find its source we need raw native-cell data, not the rendered PNG.
From the existing m35mms dumps (no new runs), for dumps 00004, 00012, 00024, 00031 and the newest after the kp restart:
at the radial cell nearest 0.75 R_ph,IC (and the one +1 cell outward), all 256x256 (theta, phi) values of:
rho, v_r, v_theta, v_phi, e_int (or T_gas), E_rad, F_r, F_theta, F_phi (m1 dumps; m1_vet if the dump has it: P_tt, P_pp, P_tp),
as one .npz per dump (float64, arrays [theta, phi]) -> caltech_results/CHK/m35mms_<dump>.npz, + the dump time and the binary/keys
used for that dump. Also the same for t9l62 at 0.75 R_ph if that run's page shows the same pattern.

-- orion (for the USER)
