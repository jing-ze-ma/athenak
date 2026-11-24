# NOTE 2026-11-24zzzzs orion -> Caltech: t9l62 late phase: why 8-10 passes per cycle (one node per arm)

Goal (USER): t9l62 in 24 h on one node. Most of its run is the late phase (dt ~96, t 7e5 -> 2.68e6), at 8-10 clean passes per cycle and ~3.3 s; the early phase has ~5. If late passes fall to the early count, late s/cyc falls by ~1/3. Please run, T2 binary (7d7c4a71), settings of the U arms + opac_drag true, thinf=false:
| arm | rst | cycles | what |
|---|---|---|---|
| PL | prod_fresh_1121b/t9l62/rst/t9l62.00031.rst (t ~1.34e6) | +40 (leave the scheme change behind), then +25 | `rad_vet/clean_verbose=2`, every key named as for t9l62 restarts (24zzp), `mesh/layout_tile=128 time/cfl_number=0.3`, hst + log only |
| PO | t9l62.00015 (outburst) | +25 | same |
| PQ | quiet T2 t9l62 rst (stop rst, t ~2.5e5) | +25 | same |
Report per arm: s/cyc, its and passes per stage line, every `CLEAN PASS` line of 10 cycles (it, de, dv, tol, fe, fv, clos / opac), CW / retries. No nsys, no dumps. These decide where the late passes sit (opacity refresh, closure refresh, xdonor freeze at pass 4, gas Newton switch at 6, or convergence of dv). On idle nodes now, after item 4 (projections) and HG / S arms.

-- orion (coordinator)
