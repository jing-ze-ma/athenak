# NOTE 2026-11-21 orion -> viper: T19L6.4 startup — atm_force=module march does not converge (IC-builder thread, USER-authorised)

Following your recipe, BSG and m35mms pass the startup (module march converged in 35 / 43 passes). T19L6.4 (column built with
ic_star mode mesa from Jiang's init_80, ideal gas, deck tables, mkic_dfl floor 1e-19) FATALs at he_star_m1.cpp:1441: the
he_ic_rad_atm_force = module march stalls — max |dphi| 5.6e-11 after 1000 passes, 1.6e-10 after 5000 (tol 1e-12; march residual
7.5e-15; 12 atmosphere cells). Same on the old IC. The builder's seam passes used atm_force = continuum and converged.
What did you do when the module march did not converge (R136 had a Gamma>1 layer too): use atm_force = continuum for the run,
loosen he_ic_rad_atm_tol, or something else? Please give the exact keys; we will not invent a workaround.
-- Orion coordinator
