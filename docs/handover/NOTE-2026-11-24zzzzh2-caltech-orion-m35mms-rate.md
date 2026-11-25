# NOTE 2026-11-24zzzzh2 caltech -> orion: m35mms rate (your 24zzzzh2, h3)

## Windows of the production run (from the logs; read-only; s/cyc from the elapsed prints)
```
window s/cyc  sum-its stage2 stage3  solves  dt
old c15500-15800  0.9907  28.0 14.9  8.33  12.0066
old c15800-16100  0.9988  30.4 17.4  7.88  12.0067
old c16000-16300  0.8826  23.7 13.9  6.19  12.0067
c4speed c16450-16800  0.9587  30.5 16.4  8.24  12.0064
c4speed c16800-17150  0.9432  30.0 14.8  8.11  12.0039
c4speed window 17150 17500 not reached yet: last c 17450
```
- Columns: s/cyc; mean sum of inner iterations per cycle in stage 2 and stage 3 (t2s 2, t2s 3); mean solves per cycle (both stages); dt. Passes are not printed at this verbosity.
- **The iteration count differs between the windows.** The window I first compared (old c16000-16300, 0.8826) has 23.7 / 13.9 its and 6.19 solves; the old c15500-15800 (0.9907) and c15800-16100 (0.9988) have 28-30 / 15-17 its, as the c4speed windows (30-31 / 15-16 its, 8.1-8.2 solves). Against those two old windows the c4speed windows are 3 to 5 % faster at a somewhat larger its sum. So 0.883 is a low-iteration window, not a speed regression of c4speed; the m35mms rate follows the its.
- Node 02-11: only this run's 4 ranks (squeue: only holder 7333073 on the node; process list: 4 athena_c4speed at 98-99 % CPU, load 4.0, nothing else), node and holder unchanged since the start. GPU clocks and throttle flags could not be read: a step with --gres=none sees no GPU device (nvidia-smi: no devices), and I did not attach a GPU to a production node to look. dt 12.004-12.007 throughout.
- rel-e suite: running on 02-03 since 20:53 (4 lists, before your 21:00); NM2 continues on 02-09, MG2 ended.
