---
name: red-giant-dfloor-energy-diagnosis-wrong
description: "09-11 11:10: my 'density floor keeps the energy -> T 1e16' diagnosis was WRONG (raising rho at fixed e LOWERS T). The 1e16 K dt-dip cells are VACUUM cells: rho_old <= 0 after the update, floored to 1e-18 with v = -0 (the fv<=0 branch), and their energy GROWS cycle to cycle at v=0 (3.3e16 -> 6.1e16 K in 32 cycles) = energy re-supplied without mass. dfloor_keep_temperature (e *= fv) was implemented per my spec, is vacuous/wrong-direction, UNCOMMITTED in the tree; decision pending: vacuum reset (rho_old<=0 -> e at tfloor) + instrument the energy source"
metadata:
  type: project
---
Agent gate data: flag off bit-identical (regress/run_dflT_off); flag on 0 differing cells in
150 cycles (floor fires only in ghosts there; min active rho 2.5e-18); reproduction from the V9f
rst with the flag (regress/run_dflT_repro, job 196303): 8 collapses vs I8's 7, seven with
rho=1e-18, T 1.4e16-8.4e16, v=(-0,0,0); one new at 5.808e5 (rho 2e-13, at the ceiling).
The agent caught the spec error: keep T => e *= 1/fv (creates energy), not e *= fv.
Proposed instead: (1) VACUUM RESET in the density floor: if rho_old <= 0, set the internal
energy to e(dfloor, tfloor) (or the pressure floor) and m = 0 -- a cell created from nothing has
no meaningful energy; (2) instrument which operator feeds energy into a v=0 rho=dfloor cell
(runaway_scan ledger on that cell class, or the hydro energy flux at its faces): pressure work
from neighbours vs RT vs conduction. Files touched by the vacuous switch (uncommitted): eos.hpp,
eos.cpp, ideal/general_c2p_hyd.hpp, ideal/general_hyd.cpp, hydro.cpp, coordinates.cpp.
See [[red-giant-vceil-r4-passes-then-hangs]].
