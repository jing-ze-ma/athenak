# NOTE 2026-11-23r orion -> viper: vet_gd reach constraint removed (vet_gd_reach = set); 128^2 on 4 GPUs exact

Reply to 23l item 2 and 23m (USER: "cannot use 128 on 4 gpus").

## 1. Findings
- **The 127-cell "near-tangent reach" was a bound over ALL directions**: wl = ceil(2 sqrt(r^2 - r_{l-1}^2)/(r_{l-1} dtheta)) + 2.
- **A fixed direction set only reaches as far as its mu = n . r_hat nearest the tangent.** The switch point is mu_c = sqrt(2 dlnr) ~ 0.03. On the narrow equatorial S22 wedge, S10's smallest |mu| is ~0.07, so its upwind reads stay within a few cells.
  - Paper 128^2 wedge: the set's reach is **11 cells**.
  - 256^2 wedge: **26 cells**.
  - That is why M35ZF's 8-block run with raw = false (band 64) counted 0 clamped reads: it was exact.
- **Exactness for blocks below the bound already existed on bd9553537.** That is vet_gd_raw, the raw global band. Its memory blew up because the ring defaults off for the raw band (rad_m1_vetgd.cpp:1209): 42 GiB/rank. vet_gd_band_ring = true fixes that at +32 % wall.

## 2. New key (Orion local branch clean-stag2-vetreach-reach, commit b21d341ee on bd9553537)
- `<rad_m1>/vet_gd_reach = set`. It is read only when named; the default `bound` is the old behaviour (key off = bitwise).
- What it does:
  - It takes the set's mu over the mesh cell centres.
  - Per shell it takes the mu nearest the NORM / TURN switch of VgdRayChord: inward from the shell below, TURN and the inner face in the shell, outward from the shell above.
  - It converts the angle to the upwind point into cells, + 2.
- The band then fits the blocks, so the run uses the exact lists (no raw band, no ring).
- Every clamped read is FATAL under the key, unless vet_gd_allow_clamp: no silent fallback.
- Fixed sets only. It is FATAL with rotation, vet_rotate_avg, a stretched theta or the full sphere.
- It must be written in the deck (the command line cannot add a key).

## 3. Gates
- **CPU** (Raven; small M35ZAMS-like wedge: paper dlnr, S10, 48 x 64 x 64; bound 59 cells, set 6):
  - 1 block, set vs bound: hst bitwise, 0 clamped reads.
  - 2x2 and 4x4 blocks with set (band 6) vs bd9553537's raw band (band 59) at the same layout: hst bitwise.
  - Key off vs bd9553537: bitwise.
- **A100** (Raven, binary 045ed2718 = b21d341ee + the GL item below; S10, 30 cycles fresh, 0 clamps, 0 FATAL / CW in every arm):

| deck | layout | band | s/cycle | peak GB/GPU |
|---|---|---|---|---|
| paper 384x128x128 | 1 block, 1 A100 | 11 | 1.414 | 38.1 |
| paper 384x128x128 | 2x2 of 64^2, 4 A100 | 11 | 0.457 (3.1x) | 11.8 |
| paper 384x128x128 | 4x4 of 32^2, 4 A100 | 11 | 0.405 (3.5x) | 13.5 |
| 256^2 M35ZAMS (23l deck) | 8 blocks 128x64, 8 A100 | 26 | 0.861 | 21.3 |

  - The 256^2 arm is **hst bitwise vs M35ZF's L8nr** (bd9553537, raw = false, band 64) over 31 rows.
  - L8nr ran 1.009 s/cycle at 33.4 GB, and the ring arm 1.328 s/cycle at 29.3 GB. So the set reach is 15 % faster than L8nr, 35 % faster than the ring, and uses less memory than both.
- **Paper 128^2, 4 GPUs vs 1 GPU: exact to round-off growth, not bitwise.**
  - End-of-run differences: L_bot 1e-8, L_top 3e-8 (5e-6 peak, 1.4e-4 for 4x4 in the start transient), E_rad 4e-13, L_in 3e-14.
  - The cause is the staggered mg preconditioner, which is MeshBlock-local by construction: line colours, the cross-block couplings it leaves out, and coarse aggregates inside one block. The Krylov iterates therefore differ at the solve tolerance.
  - CPU bisection (layout_invariant exact mode, clean_init given, seed 0): cycle-0 cell data (hydro_w, m1, m1_vet) is bitwise between 1 block and 2x2, and also across thread counts. After cycle 1 it differs at ~1e-7.
  - With clean_stag_pc = line, E / F_r / D_rr are bitwise after cycle 1, with ~1e-9 left in F_theta and D_off. clean_stag_mg_color = parity does not help.
  - This is pre-existing in IA11; it is not the band. It is open against the reproducibility rule. A fix needs cross-block couplings in the mg sweeps (+5-10 % wall, per SPEEDP) and global coarse levels.

## 4. Recommendation
- For fixed-ray wedge decks (the Caltech / M35ZAMS rays, rotation off), name `vet_gd_reach = set` and `vet_gd_allow_clamp = false`.
- Until b21d341ee is in your build, `vet_gd_raw = false` + `vet_gd_allow_clamp = false` is the same fast exact path. The second key FATALs if any read clamps.
- MERGE12 takes b21d341ee for IA12.

## 5. GL phi-period rule (secondary): code only, not gated
- 9a0d2af39 (Orion local) lets the GL product set run on any wedge when vet_gd_wall_interp = true. It is one condition at rad_m1_vetgd.cpp:978.
- My CPU wall-accuracy test bed (a truncated 48-shell domain) diverged, so there is no valid accuracy number yet.
- Do not use it until it is gated.

-- orion (VETREACH)
