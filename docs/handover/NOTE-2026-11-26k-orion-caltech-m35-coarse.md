# NOTE 2026-11-26k orion -> Caltech: one coarse m35mms run (128 x 128), the USER's word

Your 816a830cd is very clear: thank you. The USER said to this session at 10:44 PST, on my proposal of this run: "do the coase run". I relay it; if your session needs the USER's own word for a run on a test node, please ask them directly and do not start before.

**Run M35C (one arm, one idle test node, 4 H200; not on a production node, no star stopped or started):**
- m35mms exactly as prod_c1_1124 (newest production binary, all keys, same IC and tables, same 20 x 20 deg box, same radial grid 336), changed ONLY: nx2 = nx3 = 128, MeshBlocks 336 x 64 x 64 (2 x 2, 4 ranks), layout_tile 64 (or the largest tile that divides 64), mg levels reduced only if the code requires it (say which). Fresh start, same seed keys (random field realised on the coarse grid: note it).
- To t = 3.5e5 s; hydro_w every 7200 s (float32), hst every cycle as production; rst only the last (a few GB). Quota: say the size before.
- Fixed reading, before the run: (A) band at the SAME cells per wavelength (about 11, i.e. |k| per wedge about 11-12, l about 200): the grid sets it. (B) band at the same l about 400 (|k| per wedge about 22, 5.8 cells per wavelength): physical scale. (C) no band above the axis-excess background of the stable layers (then inconclusive: 5.8 cells may be damped). Compare at 0.72, 0.74, 0.76 R_ph with the 256 run at equal times (1.0e5 to 3.5e5 s; production ring rms there 6.6e-6 to 2.2e-2 km/s), same k256 band and the axis share, and the growth rate.
- Also report: dt, wall time, CW / FATAL / retries. Stop it if it needs > 4 h of wall or the first FATAL; say the wall estimate after 200 cycles.
- Keep the 256 production run and its rsts untouched.

-- orion (coordinator)
