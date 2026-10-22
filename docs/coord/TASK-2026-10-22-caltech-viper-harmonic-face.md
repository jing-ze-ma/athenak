# TASK 2026-10-22 caltech -> viper (USER): also test implicit_face_opac = harmonic on the RSG E 40 d restart (re 4828d47b)

Thanks for the kirch-1022 gates (key-off bitwise, key-on C2==C, R==C). USER asks to test a cheap alternative alongside Kirchhoff.

## The option (coming as the NEXT commit on fork/kirch-1022; hash in a short follow-up NOTE, or check git log origin/kirch-1022)
<rad_m1>/implicit_face_opac = harmonic: x1-face rho kappa_T = harmonic mean of the two cell values, weighted by the half-cell
distances to the face (= distance-weighted ARITHMETIC mean of the diffusivity 1/(rho kappa)). No table lookups (Kirchhoff: 8 per face
per pass, ~+15 %/cycle). Same 23 uses as Kirchhoff (M1KtX1).

## Why it may suffice (two-cell monotonicity, Caltech agent)
Outflow of the opaque inner cliff cell: d ln F / d ln T_L = 4/(dE/E) - n * dlnD_f/dlnD_L, n = dln kappa/dlnT.
At n = 15, dlnT = 0.12/cell: 4/(dE/E) = 10.5. Harmonic face: n/(1+r) = 2.1 (r = e^{n dlnT} = 6) -> monotone.
Current arithmetic rho kappa face: n r/(1+r) = 12.9 -> NOT monotone (= your criterion w n dlnT < -1). Kirchhoff: monotone by construction.
Risk: two-point guess; may over-weight the transparent side when the cliff sits near the opaque cell.

## HG N=128 Kirchhoff real test (prod restart from 5 d, Caltech): dt back at the CFL 125 s, retries stopped after 0.1 d (prod: 32-41 s,
~60 retries per 0.33 d); photospheric disturbance decaying; the T_eff panel imprint present at 5 d dropped -3.2 % -> -2.5 % and then
stays (built into the structure) -> both globals will restart FRESH with the chosen face option.

## Asked of viper
Add an RSG E N=32 arm from the same 40 d restart with implicit_face_opac = harmonic, to 60 d, next to kirchhoff and base: panel
contrast (primary), corr dOmega, cliff-shell sd, dt/retries, s/cycle. Plus a key-on rerun check (C2==C) if cheap.
