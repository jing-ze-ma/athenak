# TASK addendum 2026-11-11 viper -> ORION: the exact CAND key list + a CORRECTION of the vgdh-w1 table

**Correction.** The vgdh-w1 TASK (64170b860) expected-table row "CAND | ... | 0.5700 / 1.18 / 0.0131 / 0.0138" was WRONG:
those are the BASE beam scores (HR2ND_PROD 5), not CAND's. CAND contains plm-on-h, i.e. the known opaque-body failure;
your 693c4fae7 numbers (cand shd 0.0543, hcyl 0.1656 rF 0.77) are the expected ones.

**CAND (cs arm CSP-CM~fA2~norot~wall~plmh~sm~x2f15~gjump~cap~hcp) = the production sp/cs deck (st17 SP base) plus,
later keys overriding:**
- implicit_hr_cmem = true; vet_col_ncore = 64 (CM)
- vet_gd_upwind = shellmean; vet_gd_upwind_nbin = 8 (fA2)
- vet_gd_rotate_every = 0; vet_gd_rotate_dt = 0; vet_rotate_avg = 1; vet_rotate_avg_ema = 0 (norot)
- vet_gd_upwind_wall = true (wall)
- implicit_hr_recon = plm; implicit_hr_recon_q = h; implicit_thin_recon = dc (plmh; its implicit_thin_corr = none is
  overridden by hcp below)
- implicit_hr_recon_lim = smooth (sm)
- time2_vet_x2_fields = 15 (x2f15, overrides CM's 13)
- implicit_thin_gface = jump (gjump)
- implicit_hr_cmem_cap = 1.0 (cap)
- implicit_thin_corr = hc; implicit_thin_corr_plm = true (hcp)
The Cartesian "cand" of HR2ND_PROD 2b (cm + cap + gjump + plm q=h smooth + sfgate + thin_corr_plm) is the Cartesian
analogue (vet_sc decks; no vet_gd keys) and is what your beam rows should use. Beam-safe recipe = cm + cap + gjump only.
