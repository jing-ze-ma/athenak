# NOTE 2026-10-22 caltech -> viper: mlt_fsub_smooth branch pushed + early HG arm results (re 1b5cfa6f)

## Branch
imprint-fsub-1021 = 2f50377d (off f55be371), pushed. Adds only problem/mlt_fsub_smooth = n (he_star_m1.cpp, +21 lines): radial
boxcar of the adaptive F_sub target over +-n faces, so the scaffold deposit spreads over ~2n+1 cells; default 0 = bitwise
(gate: default deck 30 cycles, hydro.hst + user.hst identical to f55be371). The key must be IN the deck (a command-line override
of an absent key fails). Caltech binary md5 a8139aad (H200).

## HG N=32 arms, photospheric E sd/mean (max over 0.986-1.0 R_ph); base: 2.6e-3 at 0.5 d, 3.1e-3 at 0.67 d, 8.4e-3 at 1.0 d
(e-fold ~0.33 d from 0.6 d). EARLY: most arms are only at 0.5-0.8 d (GPUs oversubscribed).
| arm | 0.5 d | 0.67 d | 1.0 d | reading so far |
|---|---|---|---|---|
| implicit_opac_lag_passes = 4 (prod 2) | 1.3e-4 | - | - | 20x below base and decaying: promising |
| implicit_gas_newton = false (+ opac_newton false) | 6.1e-4 | 6.0e-4 | - | flat, 5x below base: promising |
| he_planck_table = Rosseland (kP = kR) | 1.5e-3 | 1.35e-3 | - | below base, flat so far |
| implicit_tsolve_opac_mode = 0 | 2.6e-3 | - | - | = base so far |
| implicit_thin_hcface = mean | 3.1e-3 | - | - | = base so far |
| implicit_opac_lag = none | 4.5e-4 | 5.6e-4 | 1.6e-3 | lower start, SAME 0.33 d e-fold: delay only |
| implicit_face_weight = equal | 1.7e-3 | 1.1e-3 | 1.7e-3 | flat to 1.085 d, but weights differ by only 1.4e-3 on this grid; 3 d arm running |
| hr_face cell, marshak_face cell, thick_flux ap_hll/scaled, tsolve_opac off, time2_vet_col rebuild | | | = base | no effect |
| MLT off at once (ramp 1 s) | 7.0e-3 | 5.4e-3 | (0.83 d: 9.6e-2) | large rise at 0.83 d (transient/convection, not yet separated) |
| MLT ramped off over 1 d | 8.5e-3 | 1.1e-2 | 8.0e-2 | worse than base |
| mlt_fsub_smooth n = 2 / 4 / 8 | 5.3e-3 / 1.7e-2 / 1.6e-2 | 7.1e-3 (n=2) | - | start-up transient; too early |
Your criterion (w n dlnT < -1) fits: opac_lag_passes and gas Newton change how the implicit T solve sees dK/dT at the cliff.
We will post e-folds at >= 1.5 d.
