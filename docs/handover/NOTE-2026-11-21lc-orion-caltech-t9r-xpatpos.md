# NOTE 2026-11-21lc orion -> caltech: V214 D_x2 = t9r with the accuracy-preserving positivity fix (xpat = pos), TEST node

D_x showed stage-start positivity alone cures t9r at cfl 0.9 (0 CW, min dt 147). The old xpat cost 4-12x on radshock, so it is off.
The new form rad_m1/time2_k1_xpat=pos touches only cells whose exchange pair would take >= 99 % of what they hold (radshock/Marshak
bitwise to xpat off). Branch caltech-test-xpatpos = b86d79f71 (= B_p1 bf4f917d4 + 842fea905). Build B_xpos (CUDA 12.9 + g++13, smoke).
D_x2: exactly V214 D_x (copy of t9l62.00033, kp keys, cfl 0.9, clean_nl_print, run past 26520, rst off) but binary B_xpos and
rad_m1/time2_k1_xpat=pos (instead of true). Report as D_x (CW, DT_RETRY, min dt, OLDVEC eint<0/E<0 per stage) -> V214/D_x2.
Lower priority than NOTE lb (production restart on B_p1).

-- orion (for the USER)
