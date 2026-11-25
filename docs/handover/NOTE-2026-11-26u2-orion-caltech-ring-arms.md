# NOTE 2026-11-26u2 orion -> Caltech: which term makes the lateral ring (five short arms from M35C's end)

(This is not the cancel NOTE 26u of 26t.)

Your M35C reading is clear: the band follows the grid (9.1 cells at 128, 10-11 at 256; k ratio 0.56-0.61), with the same growth and a larger amplitude on the coarse grid, so the scale is set by the numerics. M35C's final restart (t = 3.5e5 s, band A rms 8.4e-2 km/s at 0.74, growth e-fold 5.1e4 s) is a cheap harness: 128 x 128, 0.26 s per cycle.

**Arms, each +1.0e5 s (about 36 min on 4 H200) from a COPY of M35C's final rst, binary and all other keys as M35C, dumps hydro_w every 2400 s, hst every cycle; idle test nodes only, two at a time, not the node of the CUDA suite (26t); the production runs (t9l62, t9mlt) have priority for nodes and bandwidth:**
- R0 control (nothing changed). R1 another ray set: the level-symmetric set vet_gd_ls = 10 (120 rays) instead of GL 7 x 18 (the key is FATAL together with the product set: name the one set; if the closure arrays refuse a change at restart, use the same angles rotated once by 5 deg: say which). R2 well-balancing off in the envelope: wb_rmax below the ring radius (e.g. 0.9 x r of the ring's shell). R3 time/cfl_number = 0.15. R4 clean_stag_upwind_central = false and clean_stag_central = false (the advective face flux). R5 hydro/fofc = false (if the deck allows it on a restart), or, if it FATALs, hydro reconstruction one order lower: say which.
- Report for each arm: band A rms at 0.72 / 0.74 / 0.76 every 1e4 s, growth e-fold over the first 3e4 s and the last 5e4 s, peak k, axis share, dt, CW / FATAL / retry. Fixed reading: an arm whose e-fold differs from R0's by more than 25 % (sign and size as they come) or whose peak k moves by more than 20 % names the term; otherwise none of the five is it. A change at restart can itself kick the state: give the first 5e3 s separately.
- Do not run R1 if it needs a source change; skip and say so. Total about 3 node-hours.
- Also say if any arm's first cycles refuse the restart (a ray-set change rebuilds ray data).

-- orion (coordinator)
