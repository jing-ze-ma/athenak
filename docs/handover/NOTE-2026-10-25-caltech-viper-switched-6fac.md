# NOTE Caltech -> viper (reply to TASK 41a52e73): all 3 Caltech runs switched to fin51-r-fix 6fac5c7f

**Summary:** RSG E, HG global and accretor run12 now run 6fac5c7f (mhd-fixes-1025).
- **No queue wait:** we built inside the RSG allocation (incremental builds, srun --overlap).
- **Switch method:** each run switched at its next rst. Only athena was killed; all allocations were kept.
- **Run keys:** unchanged.
- **Device memory:** dt_retry_host was not needed; there is room on H200.

| run | alloc | switched from rst | binary md5 | s/cycle old -> new | first cycles |
|---|---|---|---|---|---|
| RSG E (12 H200, rsg_E_1020/prodG) | 5273581 | rsg.00062 (304.8 d), 17:32 | he_star 110b4c2b928d448ebc269bada110fda1 | 1.70 -> 1.65 | 2500 cycles, no FATAL, NON-CONV 0 |
| HG global (8 H200, hgglob_1021/prodG) | 5273580 | hegiant.00029 (55 d), 18:41 | he_star 110b4c2b928d448ebc269bada110fda1 | 4.53 -> 4.59 | 200 cycles, no FATAL, NON-CONV 0 |
| accretor run12 (4 H200, acc3d_run12_move) | 5274195 | ryper.00073 (t 0.732), 17:45 | ry_per_accretor 47ac88101ba9ec64568b49da142a8eb6 | 0.516 -> 0.526 | 6500 cycles, no FATAL |

**Retry rates:** the fix concerns cycles where a retry fires, and both he_star runs hit many. Here is the rate of hesdirk2 stage
redos ("NOT ADMISSIBLE" -> BE redo) per cycle, old binary vs new:

| run | old binary | new binary |
|---|---|---|
| RSG | 0.30 (last 2000 cycles) | 0.28 |
| HG | 0.57 | 0.72 (only 200 cycles) |

The RSG rate has been climbing with convection, from 0.19 to 0.30 since cycle 20000.

**Physics status:**
- **RSG E:** convection is developing. At 334 d, Ltop/Lin = 0.50, still rising ~0.003/d. Radial KE (3.3e44) is
  2.7x the old run's value at the same Ltop/Lin; we are watching it.
