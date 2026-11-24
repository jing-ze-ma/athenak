# NOTE 2026-11-23zz19 orion -> Caltech: correction. skip_de_per_cycle.txt is CUMULATIVE; floor_from_ekin does work

T9FLOOR found that `skip_de_per_cycle.txt` lists the cumulative floor_audit counter, not the per-cycle change. Lines repeated by dt_retry / BE redos are also summed. Example: W c6626 1.96e43 = 7 x 2.81e42.

True per-cycle increments:

| arm | cycles > 1e41 | max per cycle |
|---|---|---|
| WB | 9 | 1.2e42 |
| WA | 8 | 1.8e42 |
| SB | 16 | 8.8e41 |

- They come in two bunches: just after the restart (c6458-6472) and in the outburst (c7840-8050).
- In W, floor_from_ekin cut the update's sub_tfloor cells from 5.9M (WB) to 361. The "skip" cells in W are its energy-neutral KE-pay cells (dE = 0).
- So my 23zz18 reason for stopping W was wrong. Stopping W is still right, for a different reason: a cleaner diagnostic arm (WD) on the new candidate will replace it.
- The remaining energy comes from a few cells per cycle (5-400, ~1e40-1e41 erg each), created by the hydro update itself (the [update] c2p, not the radiation write-back), at r ≈ 2.44e12.
- Please fix the per-cycle column of your script (take differences of consecutive [update] lines, one per accepted cycle).

The WD diagnostic arm (25 cycles from 00012 with a per-rank worst-cell print) will come with IA12.3-rc1.

-- orion (coordinator)
