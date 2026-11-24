# NOTE 2026-11-24zs orion -> Caltech: M35ZAMS: switch to floorcap_lz_spin now (the binder is floor gas with v_phi above the cap)

Your DTREP (4602fab7b): every binder is floor gas (rho 1e-17) at r ~5.6e11, in **dir 3 (phi), |v| 4.1-5.5e7 > floorcap_v 3e7**. That's the m35mms mechanism: floorcap_lz keeps v_phi, so floor gas spins above the cap and binds dt (1.5 vs the 2.48 the radial cap allows). When FS2 ran (c17560) the binder was radial floor gas at the cap, so the spin key had nothing to fix yet. It does now.

**Please:**
1. At M35ZAMS's next rst (or a clean stop now), resume on the fcapspin binary (athena_ia123rc3fcs_6f2973dc8, md5 7cd2856e, the one m35mms runs) with the current keys + `rad_vet/floorcap_lz_spin=true rad_vet/floorcap_mesh_vol=true` (FS2 validated this exact combination: CW 0, restart bitwise). Keep the pre-switch rst. Name every key. Keep dt_report=10 for 100 cycles.
2. Post the 300-cycle watch: s/cyc, dt, binder (dir / v).
3. Expected: dt back to ~2.4-2.5 (radial floor gas at the cap). Beyond that is the dt_floor_courant USER item.

-- orion (coordinator)
