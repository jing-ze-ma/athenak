---
name: rg-box-sweep
description: 2026-09-13 Cartesian local-box resolution sweep of base convection (bench/rg_box; pgen box_convection d6eebf00/0db8a600; WB isentropic + lhllc; flux boost 100; 2 x 1.5 x 1.5 H_p; table EOS + Rosseland; 3D only, 2D dropped by the user). CLEAN round 2 (dgrad 1e-4): 16/H_p relaxes to nabla-nabla_ad = 5e-5 boosted (MLT 1e-4; the GLOBAL run's equivalent is 1.3e-2 = 260x) with Fenth' 0.95 F_b, Fkin -0.24; 32 (DONE 15 turnovers): 9.7e-5 = MLT to 3% (halves 7.4e-5/1.2e-4, still relaxing), Fenth' 0.78, Fkin -0.23, vrms 1.22 v*; 64 (DONE 15 turnovers): 7.3e-5 (halves 6.4/8.2e-5, best converged), Fenth 0.69, Fkin -0.15, vrms 1.20 v*; Fenth fraction FALLS with resolution 0.95/0.78/0.69 = open question; ctl (plain hllc, no WB, 32/H_p, DONE): -1.1e-4 SUBADIABATIC, i.e. the plain-scheme residual is itself ~ the MLT excess at 32/H_p (and ~260x it at the star's 5/H_p, cf. [[rg-v4-no-wb-vs-wb]]). => the global run's excess at 5 cells/H_p is the grid; a box at 16/H_p reaches MLT. CAUTION: Fgrav = Phi(z)<rho v1> exactly (gain ~7e3), so 'Fgrav ~100 F_b' is NOT a budget failure; the budget Fenth+Fkin+Fcond closes 0.77-0.88 F_b. Residual <rho v1> ~5e-2 rho_b v* at the bottom CHANGES SIGN with height = box-scale circulation (2x1.5 H_p fits 1-2 cells), not a leak; widening the box is the next lever (4x cost); the plain-hllc/no-WB control goes SUBADIABATIC (-3e-5: the plain residual swamps the excess). Round 1 (dgrad 1.3e-2 = the global value boosted) was 130x MLT and drove a 0.9%/10-turnover mass leak; discarded (old_wind/)
metadata:
  type: project
---
Rate probe: 1.05e8 zone-cycles/s on 2 MI300A => r128 (256x192x192) 15 turnovers ~21 h, one link; HELD until the wall is
fixed. NOTES.md / SUBMIT.md / analysis/sweep.py (full budget: Fup/Fdn raw, Fenth fluctuation form, Fkin, Fgrav, Fcond,
Ftot) / mkprofile.py (coarse-to-fine column handoff). Design problem left: a reflecting wall consistent with the
dynamic WB scheme (bc_mode 1 mirror + WB loses 39% mass per turnover; bc_mode 0 = initial column in ghosts leaks
slowly) -- the same wall/WB incompatibility the red-giant wall_noflux hack covers. Also found + fixed: conduction.cpp
cond_dtdiag read x1v (1x1 dummy) on Cartesian and off the end of w0 in 1D/2D -> every Cartesian radiative run crashed.

ARMS at 32/H_p (2026-09-13, 15 turnovers each): WB(isentropic)+lhllc 9.7e-5; WB+PLAIN hllc 7.5e-4 (7.5x MLT, halves steady
7.46/7.51e-4, vrms 1.03) = the solver's Mach-8e-3 dissipation sets the excess; POLYTROPIC WB + lhllc 9.4e-5 = same as isentropic
(closure irrelevant at 32/H_p); no-WB + hllc -1.1e-4. Both WB and the low-Mach solver are needed. sweep.py shares
analysis/drv/pairs_sweep.txt: NEVER run two sweep.py concurrently (race corrupts rows; two monitors did on 09-13).
