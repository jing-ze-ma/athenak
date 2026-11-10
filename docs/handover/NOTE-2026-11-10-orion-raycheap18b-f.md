# NOTE 2026-11-10 ORION -> viper: raycheap18b addendum (f): the vet_gd FORMAL SOLUTION does not converge in space on pulse_u

Follows c1917d7c3. pulse_u SP-CMN X ladder (n 32..512, dt fixed) rerun with an extra `<output9> tab m1_lat` (same
slice as output1: x2 1.5208, x3 0.05) -> m1_jg, m1_hg = J, H_r of the current vet_gd formal solution at t = 1.
Restriction n -> 2n by r^2-weighted pair averages, d = L1_V(q_n - R q_2n)/L1_V(q_n) along the slice:
| quantity | d per pair (32/64 ... 256/512) | p |
|---|---|---|
| J of the vet_gd formal solution (m1_jg) | 1.06e-2 1.44e-2 2.21e-2 3.07e-2 | **-0.44 -0.62 -0.48 (diverging)** |
| H/J (m1_hg/m1_jg) | 2.01e-2 1.66e-2 1.73e-2 1.22e-2 | 0.27 -0.06 0.50 |
At n = 512: H/J in [-0.43, 0.66], 28 % of the slice cells in the fse blend band 0.3 < |H/J| < 0.6.
-> the closure input itself is not convergent in dx on this thin case (Kn 5): the differences between levels GROW
with resolution, as expected for ray effects of the fixed 42-direction geodesic set (vet_gd_geo 2) on a grid that
resolves the streaks better at each level. The M1 E / F1 space failure (p ~0.3) follows from it; it is upstream of
the half-range / thin-switch face path. Tests that would discriminate (not run): vet_gd_geo 3 / 4 (more directions:
d should drop at fixed n), and the rotation + EMA averaging (production) evaluated at a fixed rotation phase.
Decks: /orion/ptmp/jinma/raycheap18b_1110/hjdiag/X_<n>/ (in.athinput, tab/so.m1_lat.00001.tab).
