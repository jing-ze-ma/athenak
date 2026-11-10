# NOTE 2026-11-10 ORION -> viper: raycheap18 (e) cs Part 2, interim (radial twins + all CS09 controls)

Part 2 of raycheap17 / (e) of raycheap18 on the cubed-sphere twins. d981e435c (CSP, CSP-FIXN), CS09 with the 09-25 keys.
Measures as st17.py (SHELL radial means Emean / F1mean + lateral spread; lateral: P_2 amplitude Amode; T: bin fields).
Lateral twins with CSP / CSP-FIXN and all CSP-CMN arms (bb58d80c2) still running (cs_rest, 4 nodes): follow-up NOTE.

## Part 0 (cs) drops
- **CS09 FATAL** on every cs twin: `on the cubed sphere ... supports only ... implicit_precond = rbgs_fwd | mg; this input
  has: implicit_precond != rbgs_fwd | mg` -> CS09 run with implicit_precond=mg (the only change).
- CSP on the twins whose sp case has vet_gd dropped (fs, pulse_edd, mvr, mvp, latmid): `implicit_flux_beam = halfrange on
  the cubed sphere needs vet_gd (the ray data of the cs sweep) or implicit_hr_model` -> implicit_hr_model=true there
  (sp silently runs those faces central; cs refuses).

## Controls: CS09 C (09-25 scheme on cs) -- 2nd order on every radial twin
| twin | Emean | F1mean | spread |
|---|---|---|---|
| pulse_cs | 1.93 2.06 2.08 | 1.87 2.05 2.07 | <= 1.6e-11 |
| fs_cs | 1.95 2.09 2.08 | 2.00 2.07 2.07 | <= 1.2e-11 |
| pulse_edd_cs | 2.01 1.99 2.00 | 1.91 1.98 2.00 | <= 4e-11 |
| mvr_cs | 1.98 1.99 2.05 | 2.39 2.28 1.98 | ~1e-14 |
| latatm_cs (z axis) | 1.97 1.99 | 1.98 1.97 | Amode **1.45 1.55** (thick: below 1.8) |
| latatmv_cs (vertex axis) | 1.97 1.99 | 1.98 1.97 | Amode 1.78 2.02 |
| latmid_cs | 2.16 2.09 | 1.96 1.99 | Amode 6.84 -3.20 (d 2.8e-3 2.4e-5 2.2e-4: erratic) |
(mvp_cs CS09 not finished yet.)

## CSP / CSP-FIXN, radial twins
| twin, mode | CSP E / F1 | CSP-FIXN E / F1 |
|---|---|---|
| pulse_cs T | -0.21 -0.05 -0.02 / -0.17 0.01 0.03 (flat, as sp) | **0.78 1.17 1.45 / 0.80 1.21 1.46** (converging) |
| pulse_cs C | 0.51 0.10 0.11 | 0.34 -0.09 0.11 |
| pulse_cs X | 0.58 0.29 0.24 | 0.29 -0.09 0.11 |
| pulse_cs spread C | **1.3e-1 1.5e-1 2.0e-1 3.4e-1 6.4e-1 (grows with n)** | 6.6e-2 5.6e-2 9.3e-2 1.5e-1 1.8e-1 (grows) |
| fs_cs T (hr_model) | 0.51 0.24 0.22 | **1.98 1.85 1.95 / 1.91 2.00 1.89** |
| fs_cs C | 0.36 0.26 0.28 | 1.08 0.65 0.37 / 1.46 0.95 0.77 |
| fs_cs X | 0.57 0.71 0.59 | 1.07 0.63 0.35 |
| pulse_edd_cs C / X | 2.01 1.99 2.00 / 1.99 2.00 2.00 | same |
| mvr_cs C / X / T | -0.52 4.75 2.06 / 2.03 2.06 2.10 / gas 2.0, E 0.89 7.76 1.99 | same (PASS) |
(T-mode m1_f2 / m1_f3 of the radial twins are ~0 fields: their "orders" are noise; spread of fs/mvr/edd ~1e-11..1e-14.)

## Reading
- The cs scheme itself (CS09) is 2nd order radially; production keys lose it as on sp: time fixed by dtc + x2 13
  (pulse_cs 1.45 rising, fs_cs ~1.9-2.0); space not (pulse_cs vet_gd as sp; fs_cs with implicit_hr_model -- on sp fs_u
  (central faces) passes space 2.07, so hr_model is the space suspect for fs_cs).
- **pulse_cs: the lateral spread of a spherically symmetric problem GROWS with resolution** under production keys
  (CSP 0.13 -> 0.64; FIXN 0.07 -> 0.18) while CS09 stays 1e-11: a symmetry break on the cs with vet_gd (+ rotation /
  EMA?) -- worth a look on its own; the radial-mean orders of pulse_cs are measured on top of it.
