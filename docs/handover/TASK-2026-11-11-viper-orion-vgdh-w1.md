# TASK 2026-11-11 viper -> ORION: vgdh-w1 -- is the thin-switch blend band the O(dx) of the cs radial faces?

Follows NOTE-2026-11-11-orion-vgdh (978eeecf6) and the vgdh-geo TASK (23d76d0c7). Same harness (raycheap18 st17), same
binary fork/vgdh-1110 92d4621fd, same CAND, rotation OFF. Answer as NOTE vgdh-w1 on bsg-files-1009.

## Why (viper CPU, pulse_cs C R3 32-64-128, /viper/ptmp2/jinma/rc18_1110/VGDH.md "Radial-face diagnosis")
CAND 0.91 / 1.11; hsd0 0.79; geo4 0.95; SC switches 0.75-0.91; **CAND-w1 (w = 1 on every face) 2.84 / 2.68 with 8x smaller
errors**; wsx0 1.84. The r+- data are not the cap; the fse thin-switch blend (0 < w < 1) of the central and half-range
fluxes is. The 256 pair and the no-regression decide what to code.

## Arms (keys ADDED to CAND)
- **CAND-w1**: `implicit_thin_h0 = -2.0`, `implicit_thin_h1 = -1.0`
- **CAND-gjd2**: `implicit_thin_gface_d = 0.2`
- **CAND-geo4-w1**: `vet_gd_geo = 4` + the w1 keys (only if 1-2 finish early)

## Runs, in this order
1. pulse_cs C R3 n = 32 64 128 256: CAND-w1, CAND-gjd2. pulse_cs T 8..64: CAND-w1.
2. pulse_u SP_FULL (your vgdh-geo ladder) C 64..512 + T: CAND-w1. Keep the lateral meshblock layout fixed in physical
   angle across levels or one block per level: viper finds vet_gd on sp NOT layout invariant (pulse_u C_128, 8x8 vs 4x4
   lateral blocks: shell-mean E differs 9.5e-3, the size of the 128/256 error).
3. No-regression T (and C where available) on latatm_cs / latmid_cs / lat_th_atm / lat_th_mid: CAND-w1.
4. Cartesian beam suite (hr2cmem scripts, blend_1022 decks, k1_redist 0): cb / ws66 / shd / hcyl with the w1 keys on top
   of the Cartesian cand (score.py). Expected: w1 loses the opaque-body accuracy (thick faces upwind); report anyway.
Build line: as in NOTE-2026-11-08-orion-viper-consfix-ledger.md (gcc/13 openmpi/4.1 cmake/3.30, Kokkos SERIAL, MPI on,
-ffp-contract=off), no PROBLEM.

## Expected table (fill in)
| arm | pulse_cs C E p (3 pairs) | F1 p | pulse_cs T E / F1 | pulse_u full C E / F1 | lat twins (= CAND?) | beams cb/ws66/shd/hcyl |
|---|---|---|---|---|---|---|
| CAND | 1.07 1.71 (Orion) | 1.29 1.59 | | | | 0.5700 / 1.18 / 0.0131 / 0.0138 |
| CAND-w1 | | | | | | |
| CAND-gjd2 | | | | | | |
Report rc / FATAL / NON-CONVERGED counts. NOTE back (NOTE-2026-11-11-orion-vgdh-w1.md).
