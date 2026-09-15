# Magnetized Gresho vortex: HLLD vs LHLLD Mach scan

Setup: `tst/inputs/gresho_mhd.athinput` (`pgen_name = gresho_mhd`), 64^2, PLM + RK2,
CFL 0.4, one turnover (`tlim = 2 pi/5 = 1.2566370614`). The uniform x1-field exerts no
force, so the vortex is a stationary solution and every bit of kinetic energy it loses
is numerical dissipation or work wound into the field.

Metric: `KE(t_end)/KE(0)` from `1-KE + 2-KE` in `gresho.mhd.hst`; `dME` is
`(1-ME + 2-ME)` end minus start. All 12 runs reached `tlim` with a finite history.

## beta = 100 (the usable column)

| mach | hlld ret. | lhlld ret. | gain | hlld dME | lhlld dME |
| --- | --- | --- | --- | --- | --- |
| 1e-1 | 0.4651 | 0.4792 | 1.030 | +3.375e-02 | +3.713e-02 |
| 1e-2 | 0.2745 | 0.3020 | 1.100 | +2.811e-02 | +4.053e-02 |
| 1e-3 | 0.2783 | 0.3187 | 1.145 | +6.934e-03 | +1.189e-02 |

## beta = 1e4 (NOT usable below mach 1e-1)

| mach | hlld ret. | lhlld ret. | gain | max KE/KE(0) hlld / lhlld |
| --- | --- | --- | --- | --- |
| 1e-1 | 0.4597 | 0.4890 | 1.064 | 1.00 / 1.00 |
| 1e-2 | 1.2339 | 1.0245 | 0.830 | 6.52 / 3.83 |
| 1e-3 | 0.5621 | 1.8685 | 3.324 | 6.57 / 2.74 |

At beta = 1e4 and mach <= 1e-2 the vortex is visibly unstable: kinetic energy *grows*
by factors of 2.7-6.6 during the turnover, drawn out of the enormous thermal reservoir
(`p0 = d0/(gamma M^2)`), and `dME` grows with it. Retention there measures a break-up,
not dissipation, so **beta = 1e4 is dropped** and only the beta = 100 column is read as
a dissipation measurement. (mach = 1e-1, beta = 1e4 is still quiescent, but a single
stable point does not make the column usable.)

## Commentary

- Leidi et al. 2022 (Balsara advected vortex) report LHLLD retention roughly
  *independent* of Mach while HLLD degrades roughly *proportionally* to Mach.
  **That is not what is measured here.** Both solvers degrade sharply from mach 1e-1 to
  1e-2 (hlld 0.465 -> 0.275, lhlld 0.479 -> 0.302) and then both saturate: from 1e-2 to
  1e-3 hlld is flat (0.275 -> 0.278) and lhlld barely moves (0.302 -> 0.319). Neither
  curve is Mach-independent, and neither is proportional to Mach.
- The LHLLD gain does grow monotonically as the Mach number falls (1.030, 1.100, 1.145),
  so the low-dissipation fix points the right way and helps more where it should. The
  effect is a few tens of percent, not the order-of-magnitude separation of the
  hydrodynamic Balsara vortex.
- The likely reason the scan does not reproduce the literature trend is that this setup
  couples the Mach number to the field strength. At fixed beta the Alfven speed relative
  to the peak rotation speed is `sqrt(2/(gamma beta))/M`, i.e. 1.1, 11, 110 at
  mach 1e-1, 1e-2, 1e-3. Below mach 1e-2 the vortex is magnetically dominated, the fast
  speed (not the sound speed) sets the numerical dissipation, and the chi magnetic floor
  `chi = min(1, max(c_uL,c_uR)/max(c_fL,c_fR))` correctly stops LHLLD from reducing the
  magnetosonic dissipation any further. A clean Mach-independence test would need beta
  scaled as 1/M^2 to hold the Alfven Mach number fixed.
- `dME > 0` in every beta = 100 run: the field is wound up in all cases, and LHLLD winds
  up more of it (it dissipates less of the flow), which is consistent with the higher
  retention rather than an artifact.
