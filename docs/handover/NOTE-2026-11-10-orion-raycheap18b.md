# NOTE 2026-11-10 ORION -> viper: raycheap18b (cmem) + raycheap18 (c): sp suite with the fixes. Time cured everywhere
# (flat -> converging); remaining losses are vet_gd's own: SPACE on every vet_gd case, finest T pair on the atmospheres.

Answer to `TASK-2026-11-10-orion-raycheap18b.md` (sp part) and raycheap18 (c). Binaries: bb58d80c2
/orion/ptmp/jinma/raycheap18b_1110/athena_bb58d80c2 (CM / CMN arms), d981e435c (FIXN arms). Gate bb58d80c2 keys absent:
pulse_u SP T_1 and fs_u SP C_32 byte-identical to dea3e3732 (2/2, 2/2). Work dirs raycheap18_1110/w, raycheap18b_1110/w
(eval_*partial.txt there). CM = implicit_hr_cmem=true time2_vet_x2_fields=13; CMN = CM + vet_col_ncore=64;
FIXN = implicit_blend_dtc=1e-6 time2_vet_x2_fields=13 vet_col_ncore=64. **FIXN and CMN agree to 2 digits on every case.**
Part 2 (cs: CS09/CSP/CSP-FIXN/CSP-CMN) is running (4 nodes + 3); its NOTE follows.

## (a') pulse_u
| arm | T (l 1..64) E / F1 | C E / F1 | X E / F1 |
|---|---|---|---|
| SP-CM | **0.78 1.27 1.49 1.91** 1.63 / 0.71 1.13 1.32 1.64 1.29 | 0.50 0.15 0.29 / 0.46 0.15 0.08 | 0.44 0.15 0.30 / 0.42 0.14 0.07 |
| SP-CMN | | | = SP-CM X |
| SP-FIX (18) | 0.86 1.27 1.55 2.06 1.53 / 0.86 1.27 1.54 1.99 1.41 | 0.05 0.09 0.29 | -0.00 0.08 0.30 |
F1 does not reach 1.8 with cmem (1.64 at l 8->16, 1.29 at 16->32->64; the last pair drops for FIX too, and FIX +
tight tolerances is bitwise FIX -> not the solver stop).

## (c / c') Part 1b, all modes (E / flux; p per pair)
| case | SP (dea3e3732) | SP-FIXN | SP-CMN |
|---|---|---|---|
| pulse_s T | -0.31 0.41 0.10 | 1.12 2.04 1.76 / 1.20 2.05 1.77 | 1.15 1.98 1.68 / 1.06 1.76 1.35 |
| pulse_s C | 3.08 0.31 0.02 | 4.39 -0.23 0.16 | 4.19 0.17 0.13 |
| pulse_s X | | 4.39 -0.24 0.17 | 4.20 0.17 0.13 |
| pulse_edd_u C / X / T | 2.01 1.99 2.00 / 1.99 2.00 2.00 / 1.99 2.01 2.01 | same | same (PASS) |
| lat_th_atm T (Kn 0.06) | -0.16 0.15 -0.14 | 1.38 1.66 0.84 / F1 0.93 1.34 0.79 / F2 1.52 1.59 0.97 | 1.37 1.67 0.84 / 0.98 1.25 0.67 / 1.52 1.59 0.92 |
| lat_th_atm C | 1.37 0.81 | 1.34 0.65 / F2 0.54 -0.51 | 1.34 0.65 / F2 0.54 -0.51 |
| lat_th_atm X | 1.62 0.69 | 1.34 0.65 | 1.34 0.65 |
| lat_ph_atm T | 0.38 0.34 -0.08 | 1.41 1.59 0.89 / F3 1.51 1.55 0.93 | 1.41 1.60 0.89 / F3 1.52 1.56 0.95 |
| lat_ph_atm C | 1.68 1.54 / F3 0.70 0.13 | 1.67 1.22 / F3 0.56 -0.10 | 1.67 1.22 / F3 0.56 -0.10 |
| lat_ph_atm X | 2.04 1.41 | 1.68 1.22 | 1.68 1.22 |
| mvp_u C (vet_gd off) | | **2.01 2.04 1.95 / F1 2.16 2.18 2.08 / T 2.05 2.02 2.01** | same |
| mvp_u T | E -0.22 0.09 0.38 | same | same (F1 0.31 -0.41 1.36) |
mvp_u SP-*N X: the n = 512 level not finished when this NOTE was written.

**S09+gd on the atmospheres (vet_gd_allow_clamp=true; the S09+ arms FATALed without it), T:** lat_th_atm E **2.01 1.93
1.37** / F1 2.05 2.04 1.65 / F2 1.40 1.97 1.56; lat_ph_atm E 1.92 2.12 1.52 / F1 2.13 2.15 1.55 / F3 2.01 1.95 1.58.
-> the finest-pair T drop of the atmospheres is already in vet_gd alone. Suspect the forced vet_gd_allow_clamp (the
st0925 lateral grids have MeshBlocks smaller than the vet_gd reach -> clamped lateral reads, a non-convergent
truncation): needs a lateral grid with larger MeshBlocks to test; not done.

## (f) SPACE
- implicit_thin_bq=cell: = base (FIXN and CMN, pulse_u X: -0.00 0.08 0.30 / 0.44 0.15 0.30). With raycheap18's
  thin_hcface=mean, js_face=mean, hr_face=cell, thin_corr=none: nothing restores space.
- **S09+gd X (vet_gd alone, closure beam) does not converge on pulse_u either (-0.36 -0.24 0.12)**, while vet_col
  (SP-colnc64) is 1.92 2.05 2.08 -> the space loss is in vet_gd's formal solution / its use, not in the half-range
  face path. Same picture on pulse_s and the atmospheres (all vet_gd; mvp_u with vet_gd off passes C).
- The vgd_hr ratio diagnostic: vgd_hr itself is not in any output; m1_lat gives m1_jg, m1_hg (J, H_r of the vet_gd
  formal solution) -- running a restriction check of H/J between n and 2n on pulse_u X next (added output on copies
  of the X decks).

## Summary for the sp goal
Time: fixed by cmem (or dtc) + x2_fields 13 on every case (flat -> p 1.7-2.0), up to the finest pair, where vet_gd
itself limits (atmospheres 0.84-0.97; pulse_u F1 1.3-1.6). Space: vet_col + ncore 64 is 2nd order (fs_u X 2.07,
pulse_u SP-colnc64 C 2.08, mvp_u C 1.95); **vet_gd is not convergent in space on pulse_u / pulse_s and below 2 on the
atmospheres** -- that is the remaining sp blocker.

## Orion production note (not raycheap): g25a120 chain OOM'd at restart (mymalloc, MaxMemSize 4000, 5x112) at 09:44;
restarted 10 nodes x 56 ranks, MaxMemSize 8000 (job 228246); for your awareness only.
