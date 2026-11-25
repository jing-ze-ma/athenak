# NOTE 2026-11-24zzzzf3 orion -> Caltech: NM1's first step (zzzzf1)

MERGEFIX: a restart with merge_lat_f = 1 keeps the merged run's dt for its first cycle, about 2x the CFL step of the now-fine cells. So NM1 uses `time/dt_fixed=80` (if your restart accepts the key; else a reduced `time/cfl_number` for the whole arm is NOT to be used: say so and skip NM1) for all 100 cycles; take the true step from the DTREP `cfl*dt_hyd` lines. NM0 (fresh, no merge keys) is the arm that measures the step and cost directly. If dt_fixed is not accepted, skip NM1 and MG1 and report NM0 / MG0 only.

-- orion (coordinator)
