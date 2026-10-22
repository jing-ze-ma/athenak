# TASK 2026-10-22 caltech -> viper (USER: top priority): test the Kirchhoff face opacity on RSG E (restart from your growing 40 d state)

Branch kirch-1022 = 61580891 (on oprlx-1022 df41092b, on opaclag-all-1022 3028a788, on imprint-fsub-1021 2f50377d, on f55be371), pushed.
Key: <rad_m1>/implicit_face_opac = arith (default) | kirchhoff; optional implicit_face_opac_ng = 2..8 (default 8). Table opacity only.
On a fresh start it must be in the deck; on a restart it is accepted from the command line.
Default gate (Caltech): default deck 30 cycles, hydro.hst + user.hst IDENTICAL to f55be371.

## What it does
x1-face transport opacity (rho kappa_T)_f (was the arithmetic / distance-weighted mean of the cell values) -> Kirchhoff mean
1/(rho kappa)_f = int_{T_L}^{T_R} T^3/(rho_f kappa_T(rho_f,T)) dT / int T^3 dT, rho_f = sqrt(rho_L rho_R), 8-point Gauss-Legendre in T,
table lookups (incl. the vet_scatter EOS kappa_es); rebuilt every pass with the cell opacities. One helper replaces all 23 x1-face
uses (face-flux rows, inertia, vimp, write-back momentum, Krylov/off-diagonal operator, vet_lat, AP-HLL face tau). Thick-limit flux
~ (Phi(T_R) - Phi(T_L))/dr, Phi = int T^3/(rho kappa) dT -> d(outflow)/dT_own > 0 at any cliff slope (your criterion w n dlnT < -1
cannot be met). Lateral faces unchanged. Cost +15 %/cycle at N=32 (shared GPU), dt unchanged.

## HG N=32 restart from the growing state (seed 0.75 d)
| t (d) | panel L c-e kir | base | E sd kir | base |
|---|---|---|---|---|
| 0.88 | +1.3e-6 | +1.1e-5 | 4.7e-3 | 5.8e-3 |
| 1.00 | -6.2e-7 | +1.9e-5 | 5.5e-3 | 8.4e-3 |
| 1.13 | +3.0e-6 | +2.8e-5 | 6.0e-3 | 1.1e-2 |
Panel contrast ~10x lower; cliff disturbance e-fold ~1.1 d vs 0.3 d. First restart arm that differs from base at all.
N=128 HG real test (prod restart from 5 d, 12 H200) running at Caltech now.

## Asked of viper
1. RSG E N=32 restart from your 40 d growing state with implicit_face_opac = kirchhoff, to 60 d, vs your base continuation:
   panel contrast (primary, USER criterion), cliff-shell sd, dt/retries, s/cycle.
2. Reproducibility with the key ON (Caltech has no free GPUs while the N=128 test runs): C, C2 (rerun) and R (restart), bitwise.
