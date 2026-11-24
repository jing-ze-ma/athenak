# NOTE 2026-11-24zh orion -> Caltech: m35mms dt halved after floorcap_lz (10.55 -> 5.6): post the dt binder

Your floorcap_lz switch shows dt 10.55 (pre) -> 5.61 (after), i.e. ~2x more cycles, while s/cyc fell 1.108 -> 1.005.
- Please run a 50-cycle dt report on m35mms: at the next convenient restart, or as a side arm from its newest rst (`hydro/dt_report=10`, nlim +50, no outputs), on a free node or a GPU slice if possible. Post the DTREP / DTREG lines.
  - We expect the binder to be floor gas (rho < floorcap_rho) that now keeps its v_phi above floorcap_v, as on M35ZAMS.
- Keep m35mms running with floorcap_lz = true meanwhile: exact L_z is an absolute rule. The dt fix (a floor-gas dt key) is the USER's decision when they're back.
- M35ZAMS (Z3, floorcap_lz on): post dt every ~1000 cycles so we see if the same happens there.
- t9l62 layout: T1 2.867 vs T0 2.846 s/cyc (no gain, md5 identical). Wait for T2, then resume t9l62 production with the faster of T0 / T2 per the 24v rule (most likely T0, the current layout).

-- orion (coordinator)
