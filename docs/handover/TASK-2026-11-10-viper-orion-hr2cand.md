# TASK 2026-11-10 viper -> ORION: hr2cand -- the half-range production candidate on the cs / sp twins (raycheap18 harness)

Coordinator GO 12:30 PST. Add the arms below as EXTRA arms of your raycheap18 harness (st17-style composable arms), rotation
OFF everywhere. Answer as NOTE hr2cand on bsg-files-1009. viper runs the same pulse_cs / pulse_u arms on Raven A100 until
13:00 PST; whatever finishes first wins (report anyway).

## Code
fork/hr2prod-1110 **02798a1b2** (all keys opt-in, read only when named; keys absent bitwise vs 382e8cef3 on CPU; the cs
candidate key set np 1 = 3 = 4 bitwise, 12 cycles). New since 2aeefd3df:
- `implicit_hr_cmem_cap = k` (f99c21f48): cmem's F_c clamped to |F_c| <= k c max(E'_L, E'_R) per face (ws66 at
  implicit_cfl 256: cm 2.81 -> 1.1782 = base 1.1771; inactive in every order test: digits unchanged).
- `implicit_hr_recon_sfgate` (02798a1b2): plm on h only off the beams (cells in the vet_sc sf state Jh > J and their
  neighbours keep donor cell). Inert on cs / sp (no sf there); needed on Cartesian vet_sc.
- merged: hrhc-1110 (`implicit_thin_corr_plm`), aphll-time-1110 (cons2 cmem fix a1b8b62c7, implicit_aphll_fcimp).
Same CPU build line as hr2cmem.

## Cartesian result that motivates it (plane-parallel beam-free vet_sc pulse, your hr2cmem decks, viper CPU)
Candidate = cb keys + cm + cap 1 + gjump + plm q = h smooth + sfgate + thin_corr_plm (hc and jse kept, no DTC):
pk12 C 2.06 2.00 / F1 2.03 2.01, T 2.04 2.04 2.05; pk1 C 2.08 2.03, T 1.90 1.97 1.99; pk128 C 1.98 1.96, T 1.99 1.82 2.14.

## Arms (keys ADDED to your production arm; later keys override earlier ones)
- **CAND** (cs and sp):
  `implicit_hr_cmem = true`, `implicit_hr_cmem_cap = 1.0`, `implicit_thin_gface = jump`,
  `time2_vet_col = extrap2`, `time2_vet_x2_fields = 15`, `vet_col_ncore = 64`,
  `implicit_hr_recon = plm`, `implicit_hr_recon_q = h`, `implicit_hr_recon_lim = smooth`,
  `implicit_thin_corr = hc`, `implicit_thin_corr_plm = true`, `implicit_thin_recon = dc`,
  `vet_gd_upwind = shellmean`, `vet_gd_upwind_nbin = 8`, `vet_gd_upwind_wall = true`,
  rotation off: `vet_gd_rotate_every = 0`, `vet_gd_rotate_dt = 0`, `vet_rotate_avg = 1`, `vet_rotate_avg_ema = 0`;
  NO `implicit_blend_dtc`. (viper's st17 name: CSP-CM~fA2~norot~wall~plmh~sm~x2f15~gjump~cap~hcp.)
- **CAND-nogj**: CAND without `implicit_thin_gface`.
- **PRODN** (reference): production keys + rotation off (+ Fix A + wall as above), nothing else.
- **CENN** (reference): PRODN + `implicit_flux = central`, `implicit_flux_beam = closure`, `implicit_thin_switch = none`,
  `implicit_thin_corr = none`, `implicit_flux_faces = x1`, `implicit_blend = tau_f`, `vet_col_ncore = 64`.

## Runs
pulse_cs C with R3 (3-D refinement, 32..256) and T (8..64 or 128); pulse_u full refinement C 64..512 (+ T if cheap);
latmid_cs / latatm_cs and the sp lat twins T (no-regression). Report Emean / F1mean (cs shell means) and E / F1 / F2 / F3
p of all pairs, lateral spread, rc / FATAL / NON-CONVERGED counts.

viper reference numbers (pulse_cs C R3, Raven A100): PRODN-like CSP-DTC~fA2~norot~wall~plmh~sm~pred 1.35 1.57 /
F1 1.20 1.30; + gjump 1.84 1.72 / 1.77 1.47; central 1.87 2.18 / 1.67 1.99; CSP-CM~...~x2f15~gjump (old binary,
no cons2 fix, no cap) 0.91 1.73 / 1.11 1.57. pulse_u full refinement (A8norot_DTC + wall + predict): base 0.88 1.05,
+ plmh 1.13 1.65, + plmh + gjump 1.10 1.48.
