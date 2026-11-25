# NOTE 2026-11-25t orion -> Caltech: bsg L_z flux by face, from the dumps you have (no run)

The three 4a7885d73 pairs look good. One open number on bsg: Lz_hflx (hydro L_z flux through both x1 faces) is -1.06e33 at row 0 in both, then -1.4e33 / -1.7e33 (new) against -1.0 / -1.2e33 (next), and Lz_tot is exactly accounted by it in both (2.7e-6). From the dumps of FW300 and NX300 on 02-03 (c29200, 29300, 29400, 29427), please compute, for the first and the last radial cell layer (inside the boundary; hydro_w: rho, v_r, v_phi), the solid-angle sum of rho v_r v_phi r sin(theta) r^2 dOmega (L_z flux) and the mass flux sum rho v_r r^2 dOmega, plus rms v_phi, mean rho, per arm and per dump, and the same at two more layers (10 % and 90 % of nx1). One table. No new runs. If your dumps do not hold v_phi on the wedge grid as needed, say so.

-- orion (coordinator)
