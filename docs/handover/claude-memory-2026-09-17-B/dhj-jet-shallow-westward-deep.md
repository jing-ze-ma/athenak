---
name: dhj-jet-shallow-westward-deep
description: The dhj equatorial super-rotating jet is SHALLOW (~1 km/s above 30 mbar) on every grid and solver; below 0.1 bar the equator is WESTWARD and strengthens to rot 164 on sp. The deep jet never spins up. Measured 2026-09-09
metadata:
  type: project
---

Measured 2026-09-09 (viper) with docs/handover/scripts/dhjcs.py (pressure from the run's own
EOS table; 6-deg bins; equatorial = mean of +-3 deg bins). Runs: cs_hyd_rs/{hllc,lhllc,
ausmpup,ppmx} (cs, 09-07 defaults, pure hydro) and sp_dhj_hyd (sp, d25d8792, no rot_potential/WB,
dumps to rot 165).

**Every solver on cs gives the same jet**: eastward equatorial core 1.0-1.1 km/s at rot 4 centred
at p = 1-6 mbar, westward poleward of +-20 deg; solver-to-solver rms 11-15 % of peak at rot 2,
~20 % at rot 4 (chaotic spread). lhllc and ppmx tracked hllc/ausm+up in horizontal KE to 2-6 %
until their vertex dt collapse. cs matches sp in structure, position and magnitude (rms 12-18 %).
Time-mean over rot >= 6 of u_eq at 1e-2 bar: hllc 0.38, ausm+up 0.69, sp 0.73 km/s -- cs-hllc
runs at HALF strength in the upper jet (more snapshot variability, same total KE); ausm+up
matches sp best.

**The deep jet never spins up.** On sp out to rot 164 the equatorial zonal mean at 0.1-10 bar is
westward at every sampled time and strengthens monotonically: 1 bar -0.19 (rot 14) -> -0.37 km/s
(rot 164); 10 bar +0.12 -> -0.23. The eastward jet stays confined to p < 3e-2 bar, saturating
by rot ~40 at 1.0-1.55 km/s. Deep eastward flow appears only OFF-equator (0.2-0.7 km/s at
+-40-57 deg, 0.1-0.6 bar). cs at rot 14 matches all of this.

**How to apply:** "no deep jet at rot 15 on cs" is NOT a cs shortfall; it is what sp does at 10x
the integration. But it is also NOT the canonical hot-Jupiter super-rotation (GCMs: eastward
equatorial jet of several km/s down to ~1 bar). Whether the shallow, deep-westward state is the
physics of this model (H2 dissociation heat pipe, general EOS, drag/sponge, radiative timescales)
or a model defect is OPEN and is a science question, not a grid question. Figure:
session scratchpad jet/u_lat_p.png (regenerate with dhjcs.py). See [[cs-dhj-diagnostics-rot20]],
[[cs-hyd-rs-run]], [[sp-hydro-vs-mhd-comparison]].
