# UPDATE 11-13 viper->Caltech (USER): T9L6.2 at Caltech = a SHORT TEST first, in a short 2-GPU backfill holder (< 4 h)

This replaces "run to tlim" in TASK e72c57b05: do not queue the 800 h production. The files are the same: docs/handover/t9l62-1113/.

## Test
A full-resolution wedge, one arm per H200, both arms at once in one 2-GPU slot. Expected wall 3-3.5 h.

| arm | deck | GPU |
|---|---|---|
| OW | t9l62_ourwall.athinput | 0 |
| LW | t9l62_final.athinput | 1 |

Command-line overrides for both arms:
  mesh/nx2=128 mesh/x2min=1.1780972450961724 mesh/x2max=1.9634954084936207 mesh/nx3=128 mesh/x3max=0.7853981633974483
  meshblock/nx1=512 meshblock/nx2=128 meshblock/nx3=128
That is 512 x 128 x 128 (8.4 M cells, ~42 GB) over theta 3pi/8..5pi/8, phi 0..pi/4, with the paper's cell size and periodic angular BCs.

## Steps
1. Smoke: 10 cycles, then restart and run 10 cycles. Check rc 0, FATAL 0, NaN 0, and that the restart is bitwise.
2. Run each arm with -t (slot end - 15 min). Write a restart at the end and hst every 1000 s.

## Report by NOTE (Caltech->viper)
- s/cycle, mean dt, cycles and simulated time reached
- GPU memory, stage retries, BE steps, NC
- history columns: L_in at the wall vs target 6.9466e13 x 4 pi r_in^2 x (wedge solid angle / 4 pi); L_out; v1sq_wall; Mdot_bot; total energy
- the rst paths, in case a later slot continues them

Viper runs the same 1-day wall comparison on a coarser grid (dayOW512r vs dayF512rB). Your full-cell-size H200 numbers decide the production layout.

## UPDATE2 (USER: "with what GPU we have now")
Use the short2g holder 6854764 that is running now (empty, ends ~16:55 PST) at once.
- Do the smoke plus restart check for both arms, then run both arms until 10 min before the end, with a restart at the end.
- Report what you got (s/cycle, dt, memory, wall L_in) even if it is only ~20-30 min.
- Continue from those restarts in the next short slot (test2g / any backfill).

## UPDATE3 (USER: "use 2 H200"): ONE run on both H200, replacing the two-arm split above
- Deck: t9l62_ourwall.athinput only. No LW arm.
- 2 MPI ranks, 1 H200 each, in the short2g holder now (or the next 2-GPU slot).
- Overrides:
    mesh/nx2=128 mesh/x2min=1.1780972450961724 mesh/x2max=1.9634954084936207 mesh/nx3=256 mesh/x3max=1.5707963267948966
    meshblock/nx1=512 meshblock/nx2=128 meshblock/nx3=128
  That is 512 x 128 x 256 = 16.8 M cells over theta 3pi/8..5pi/8, phi 0..pi/2, at the paper cell size: 2 meshblocks, ~42 GB per GPU.
  The full grid (67 M cells, ~336 GB) does not fit in 2 x 141 GB.
- Steps and report: as above (smoke + bitwise restart, run to the slot end - 10 min, restart at the end).
  Scale the wall L_in target by the wedge solid angle.
