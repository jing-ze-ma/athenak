---
name: sp-polar-hlle-swap-is-the-excess
description: SOLVED 2026-09-12 -- the sp MHD energy excess IS the MHD flux kernels' unconditional HLLD->HLLE swap on the polar rows (mhd_fluxes.cpp do_pole, since 9fb17531; hydro never had it). New switch <mhd>/polar_hlle_rows (7d87f3c5, default true = old); false makes MHD == hydro to 6 digits in E/KE/mass at 0.1 rot with and without the field (bench/sp_excess3). sp_mhd_prod3's 7x KE and +5% E are this artefact. sp_mhd_nopole (11639436/7) tests polar stability with it off
metadata:
  type: project
---

Chain of evidence: sp_excess (B=0 carries the full excess) -> sp_excess2 (hlle in BOTH modules: bitwise
identical; nowb/nocond/norot/noemfdiss leave the gap) -> localisation (deep rows within ~45 deg of the
poles, spurious meridional circulation, day/night symmetric) -> sp_excess3 at 0.1 rot:
  hyd_ctl E 6.63525e38 KEh 8.43e32 | mhdb0 swap on 6.63978e38 / 1.77e33 | swap off 6.63524e38 / 8.44e32
  mhd 3G swap on 6.63978e38 / 1.74e33 | swap off 6.63523e38 / 8.13e32, ME 3.96e31 unchanged.
**Why:** HLLE on the two polar rows dissipates the polar flow into heat and drives a deep meridional cell;
every sp MHD-vs-hydro comparison since the MHD work began ([[sp-mhd-energy-excess]],
[[sp-hydro-vs-mhd-comparison]] "field takes 31% of zonal KE") measured this, not the field.
**How to apply:** set polar_hlle_rows = false in every new sp MHD input once sp_mhd_nopole shows the polar
rows stay stable past rot 10-12 (the swap may have been masking [[sp-polar-field-blowup]], since fixed by
polar_emf_diss / GS05 checkerboard fix). cs has no polar rows: unaffected. If stable, flip the DEFAULT to
false and re-baseline sp_mhd_prod3 (its rot 150 state is contaminated).

**2026-09-13 ~00:30 THE SWAP CANNOT SIMPLY BE TURNED OFF (bench/polar_blast_hlle, user's question).** Canonical polar
MHD blast (run/pole_blast/sph_blast.athinput, PROBLEM=pole_mhd_blast, 256x256x32, blast on the north axis, added by
the same commit 9fb17531 as do_pole): with polar_hlle_rows off the FAR (south) polar row j=255 drains (rho 0.33 ->
0.017 by t=0.06, c_f 3.5 -> 15) and dt falls geometrically from t~0.045 (6.9e-6 at t=0.0685 vs 3.05e-5 flat ON);
NO NaN, div B round-off, hst energies agree to 5 digits -- a dt death, invisible at 64^2/128^2. sp_test resistive
decay: L1(B) on 1.785e-3 vs off 1.859e-3, rates 2.45 vs 2.30 (both fine). => polar_hlle_rows is now a PER-SWEEP MASK
(uncommitted in the working tree; bit1=x1, bit2=x2/theta [swaps TWO rows per pole], bit4=x3; 7=old default,
0=none; true/false accepted). Running: polar_blast_hlle/mask (CPU 256^2, masks 2,5,1,4,3,6: minimal stable mask) and
bench/sp_excess4 (GPU 0.1 rot, B=0 arms masks 0,1,2,4,7 vs hyd_ctl: which sweep carries the excess). Goal: a mask
that keeps the blast stable AND removes the dhj excess; if none, the fix is stopping the polar-row density drain.
sp_mhd_nopole (mask 0 = old false) may hit the drain: watch its dt.
**sp_excess4 (11653352): the excess is ENTIRELY the x1 (RADIAL) sweep's swap.** B=0 arms at 0.1 rot vs hyd_ctl:
mask 0 (none) E-E_hyd 0, KEh ratio 1.00; mask 1 (x1 only) +7.2e35, 7.4x; mask 2 (theta only) 0, 1.00; mask 4 (phi
only) 0, 1.00; mask 7 (old) +4.5e35, 2.1x. => candidate default mask 6 (theta+phi swapped, radial not) if the blast
is stable with it (polar_blast_hlle/mask running).
