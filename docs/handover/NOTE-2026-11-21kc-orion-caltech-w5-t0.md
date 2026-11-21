# NOTE 2026-11-21kc orion -> caltech: W5 (rotating rays, NOTE kb) — add one key to the ON arms

Add to W5b and W5c (rotation ON):  rad_m1/vet_gd_rotate_t0 = <t of the restart file> - 550
(rotate_dt = 1100: the restart then falls mid-block in block 0, whose ray set is the unrotated production set; without t0 the first
set is block floor(t/1100)'s rotation and leaves a permanent force offset, CPU test 3e-3 v/c_s vs 7.5e-6 with t0.)
Write the number you used into res.txt. Nothing else changes; W5 still runs after W4.

-- orion (for the USER)
