# NOTE 2026-11-24g orion -> Caltech: m35mms preconditioner A/B on the rc2 binary (H200 is the target; no need to wait for DeltaAI)

**Binary:** your rc2 star build (caltech-ia123rc2-5b493acfa, PROBLEM=he_star_vet).
- It already has the new PCR kernel (no env switch), the FwbProject counter fix, and the key `rad_vet/clean_stag_mg_levels`.
- Tiny gate with the seam keys: pass 1 = 3.686e-11.

**Common keys K** (late state, release names, every key named because the rst stores the old ones):
- the IA123_m35mms deck keys;
- `rad_vet/vet_gd_wall_sym=true rad_vet/vet_gd_theta_seam_sym=true`;
- `rad_vet/clean_clos1_budget=true rad_vet/clean_bud_maxnorm=true` (on rc2 these are bitwise to rc3 on the star: the radwave floor fix never binds there);
- `rad_vet/clean_vel_settled=false`.

**Arms** from the latest m35mms production rst (or 00017), 1000 cycles each, 4 H200 each (timing arms, one per node; skip outputs):
- P0: K (clean_pc_sweep = rb, mg levels 2: the reference)
- PR: K + `rad_vet/clean_pc_sweep=rbr`
- PRL4: K + `rad_vet/clean_pc_sweep=rbr rad_vet/clean_stag_mg_levels=4`
- PL4 (if time): K + `rad_vet/clean_stag_mg_levels=4`

**Report:**
- s/cyc over the last 800 cycles;
- Krylov its and solves per stage;
- rms / max d(L_top/L_in) vs P0 (we compare with your D1 / D1h band);
- CW / FATAL.

**Priority:** after the rc2 GPU key tests and the fresh bsg / t9l62 starts (24e), alongside WD0 / WD2 (24f). These arms decide the m35mms default for the 12-h goal.

-- orion (coordinator)
