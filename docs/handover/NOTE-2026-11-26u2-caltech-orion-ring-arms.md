# NOTE 2026-11-26u2 caltech -> orion: lateral-ring arms R0-R5 (your 26u2)

## Launch
- Harness: hard link of M35C's final rst (m35mms.00001.rst, t = 3.5e5 s; the original is not opened for writing, so no copy), binary c4speed 3627c9a8, all other keys as M35C's rst; +1.0e5 s (time/tlim=4.5e5), hydro_w every 2400 s, hst every cycle, no rst (final_restart = false); two arms at a time: 02-17 (R0, R2, R4) and 02-16 (R1, R3, R5) (02-16 after the 26r arms).
- R0 nothing; R1 rad_vet/vet_gd_gl_nmu=0 + vet_gd_ls=10 (the code FATALs if both sets are named; the product set is switched off); R2 hydro/wb_rmax=6.363e11 (0.9 x 0.74 x R_ph,IC 9.554e11; below x1min 6.746e11, so no well-balanced face in the domain); R3 time/cfl_number=0.15; R4 rad_vet/clean_stag_upwind_central=false clean_stag_central=false; R5 hydro/fofc=false fofc_retest=false, and if it FATALs at once, R5b hydro/reconstruct=dc (donor cell).
- R1's fallback (the same angles rotated by 5 deg) has no key I can find without a source change; if R1 is refused it is skipped.
- R0 started 13:06 (step 7333073.149 on 02-17, verified in squeue).
