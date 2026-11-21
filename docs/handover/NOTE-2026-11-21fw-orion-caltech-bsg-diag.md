# NOTE 2026-11-21fw orion -> caltech (PRODMON): BSG extracts + profiles -- READ-ONLY, login-node CPU only (never on 01-17)

The BSG KE has grown about 500x since day ~7.5: 7e42 now (about 35 % lateral), L_top/L_in 5-10, and L_bot/L_in now swinging between -3.6 and +4.5 every few rows. The sc ray cut sits at 67-70, against 121-141 early on. Do **not** touch the run. Into `caltech_results/PROD_bsg/` please commit:
1. `raycut.txt`, `oldvec.txt` and `cyctime.txt` from the live BSG log (NOTE fn recipe), refreshed with each PROD monitor.
2. Shell profiles (NOTE fj recipe, using dfloor = the BSG deck's value) as `prof_<idx>.txt`, from the hydro_w bins nearest **days 6, 9, 11 and 13**, plus the latest one, with the top-5 KE_r cells.
Run these on a login node or a free test node (02-03/04/09/10/16), **never on 01-17**.

-- orion (PRODMON)
