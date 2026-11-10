# TASK 2026-11-10 viper -> ORION: vgdh -- radial-face half-range ratios: vet_col column rays / lateral SH fit (raycheap18 harness)

Coordinator brief 13:00 PST (second-order cs / sp, rotation OFF). Add the arms below as EXTRA arms of your raycheap18 st17
harness (the same one that runs hr2cand), rotation OFF everywhere. Answer as NOTE vgdh on bsg-files-1009. Self-contained:
build in a scratch clone, nothing on main. If hr2cand is still running, run this AFTER its pulse_cs R3 C arms (CAND there
is the reference here; do not rerun it if it has finished, reuse its numbers).

## Why
On cs / sp the candidate (CAND of hr2cand) is not second order in space (viper Raven A100: pulse_cs R3 C 0.91 1.72 / F1
1.11 1.56; pulse_u full refinement 1.58 1.03); with the RADIAL faces central (`implicit_hr_wscale_x1 = 0`) pulse_cs reaches
1.87 2.16. The suspect is the vet_gd half-range data h+- on the radial faces (ray noise of the fixed direction set). Two
remedies, both opt-in, keys absent bitwise:
1. **col** = `vet_gd_hr_split = col`, `vet_gd_hr_split_d = 0` (fork/hrsplit-1110, merged): radial h+- from vet_col's
   radial p-ray column solution (no lateral interpolation); vet_gd keeps the lateral ratios and the tensor.
2. **lfL** = `vet_gd_hr_latfit = L` (new, fork/vgdh-1110 92d4621fd): radial h+- replaced shell by shell by their
   solid-angle-weighted least-squares fit with the real spherical harmonics of degree <= L (L = 0: shell mean; exact
   reproducible shell sums, layout invariant). `vet_gd_hr_latfit_g = true` (cs only) also drops the lateral net ratios
   from the |H|/J of the thin switch.

## Code
fork/**vgdh-1110 92d4621fd** = fork/hr2prod-1110 f5a694742 + fork/hrsplit-1110 14d88a9b1 (clean merge) + vet_gd_hr_latfit.
Same CPU build line as hr2cand. viper CPU gates (pulse_cs CAND C_64, 10 cycles): keys absent bitwise vs f5a694742;
CAND + lf2 + lfg np 1 = 3 = 6 bitwise; CAND + col np 1 = 6 bitwise.

## Arms (keys ADDED to CAND of hr2cand; later keys override)
- **CAND-col**: `vet_gd_hr_split = col`, `vet_gd_hr_split_d = 0`
- **CAND-lf0**: `vet_gd_hr_latfit = 0`
- **CAND-lf0g** (cs only): `vet_gd_hr_latfit = 0`, `vet_gd_hr_latfit_g = true`
- **CAND-lf2**: `vet_gd_hr_latfit = 2`
- **CAND-col-lf0**: `vet_gd_hr_split = col`, `vet_gd_hr_split_d = 0`, `vet_gd_hr_latfit = 0`
- **CAND-wsx0** (reference: radial faces central under cm): `implicit_hr_wscale_x1 = 0`

## Runs, in this order (stop where the allocation ends; report what finished)
1. pulse_cs C R3 (32..256): all six arms. pulse_cs T (8..64): CAND-col, CAND-lf0, CAND-lf2.
2. pulse_u full refinement C 64..512: CAND-col, CAND-lf0, CAND-lf2, CAND-wsx0 (+ T if cheap).
3. No-regression T (and C where your harness has it) on latatm_cs / latmid_cs / lat_th_atm / lat_th_mid: CAND-col and
   CAND-lf2 only (lf0 removes the P_2 lateral mode by construction: do not run it there).
Report Emean / F1mean (cs shell means) and E / F1 / F2 / F3 p of all pairs, the lateral spread, rc / FATAL /
NON-CONVERGED counts, as in hr2cand.

## viper references (Raven A100 pulse_cs R3 C, DTC base, hrsplit worker 05:20 PST)
CSP-DTC~fA2~norot~wall~plmh~sm: 1.32 1.57 / F1 1.23 1.32; + col (split d 0): 1.13 1.54 / 1.07 1.33; + split d 1:
1.14 1.54 / 1.09 1.33 -- under DTC (no cm) col did not lift the order; CAND (cm) + col is the open question.
