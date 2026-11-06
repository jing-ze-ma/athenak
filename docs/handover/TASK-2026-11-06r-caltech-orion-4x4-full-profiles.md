# TASK 2026-11-06 Caltech -> Orion (USER): full radial profiles from the 4 x 4 column (job 227389)

The USER asks why our AthenaK gas velocity differs from yours with luminosity matched. Your wall files stop at cell 23 (27 Rsun). Our `Lr` has a very different v_r further out at 0.5 d: at 26.41 Rsun it's +342 m/s in ours, -32 m/s in `his1d_g`. To see where the pulse starts we need all radii.

Please push `docs/handover/athenapp-bsg/repro3d/col4x4/orion/col4x4_prof_<t>.txt`, shell means of **all 256 radial cells**, at t = 0, 0.05, 0.1, 0.25, 0.375, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3 d, or at the nearest dumps (say which): r, rho, Pgas, Er(lab), Fr1(lab), Fr01(comoving), vr (rho-weighted), T, kappa_Ross (shell), and total Sigma_s + Sigma_a if it's in the dump. Double precision if possible. It should be small and cheap, from existing dumps. Short NOTE when pushed. Nothing else changes.
