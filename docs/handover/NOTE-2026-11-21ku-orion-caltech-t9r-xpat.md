# NOTE 2026-11-21ku orion -> caltech: V214 D_x = t9r with the stage-start positivity fix (xpat), on the 4th test node

Orion analysis of t9r k9b_log: at cfl 0.9 the stage-2 start has a fast-growing count of eint<0 floor-gas cells left by the hydro
stage (1e3 -> 5.6e5 per stage by 26475) before the 26485 stall; the 26512 hydro dt collapse is one floor-density cell (rho 1e-20)
left at T 4e11 K. rad_m1/time2_k1_xpat (POS7, Patankar-form exchange pair) takes that count to ~0 on every deck at unchanged dt.
Branch caltech-test-fwbx = 75ce347cf (= f496d606b + xpat cherry-picked). Build B_fwbx (CUDA 12.9 + g++13, smoke gate), then on a free
TEST node: D_x = same as V214 D_a (copy of t9l62.00033, kp keys, cfl 0.9, clean_nl_print=true) + rad_m1/time2_k1_xpat=true,
run past cycle 26520, rst off. Report as V214 (CW, dt_retry, dt min, OLDVEC eint<0/E<0 per stage, floor_audit) -> V214/D_x.

-- orion (for the USER)
