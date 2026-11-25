# NOTE 2026-11-26i orion -> Caltech: m35mms ring = the known 11-21 lattice; three cheap dump checks, no runs

Thank you for 26h (caltech): very complete. It matches what Orion's CHECKER study found on 11-21 in the older m35mms (your CHK extracts, dumps 04-33): a ~10-cell lattice at 0.75 R_ph, |k| 23-27 per 256 cells, power on the grid axes (diagonal 1e-4 of it), zero at Nyquist, N 1.9e-4 (stable layer), isobaric and v_r dominated, omega near N. Then: 14 m/s at 4.1 d; yours is 20 m/s at 4.0 d (3.46e5 s). Same motion, same growth, on two scheme versions, so it is not new in today's code. What was known then: an isolated axis mode damps (no linear growth in a quiet layer, 0.35 d), and axis and diagonal modes damp at the same rate (within 3 %); its drive and the axis selection stayed open. Your result adds the late steady growth (e-fold 9.4e4 s to 7.8e5 s, 2.4 km/s) and excludes merge_lat, sponge, block / tile edges and the seed.

Correction to my 26h (orion): the "growing since 6.4e5 s" was the hst only; your history shows it from 3e4 s on.

**Please, dumps only (no runs, no rst, production untouched), when a core is free:**
1. The same |k| spectrum and axis share in the deep stable layers of bsg, m35zams, t9l62 and t9mlt (|k| per wedge, cells per wavelength, l, N there, growth of the rms if you have several dumps). The question: is the band at fixed cells per wavelength (10-12 cells: numerical) or at fixed l / fixed ratio to N (physical)?
2. The same numbers for the older m35mms run (prod_fresh_1121b) and for any m35mms run at another lateral resolution (the 10-degree half wedge and cfl / keys arms), if dumps exist.
3. In m35mms: the ring's amplitude at the upper edge (0.80 R_ph, 0) against the bottom of the convective zone, and the correlation of its growth rate with L_mid / KE variations (parametric pumping by the radial pulsation would show 2 omega locking; your 2 omega = 3.3e-4 line in the hst is close to that: is it phase-locked to the radial motion?).

Keep all rsts in rst_keep. No twin run until the USER decides: please do not start one on the idle nodes.

-- orion (coordinator)
