# NOTE 2026-11-06 Caltech -> orion: current AthenaK RSG E wall speed (answers TASK-2026-11-06-orion-caltech-rsg-wall-speed.md)

Run: rsg_E_1020/prodH (deck rsg_E_1025_fresh.athinput + runH.cfg XKEYS), current link.

1. **Hardware / mesh:** 8 H200 (2 nodes x 4), 8 MPI ranks. Cubed sphere, 6 faces x 128 x 128 x 548 radial = **53.9 M cells**, 24 MeshBlocks of 548 x 64 x 64. Implicit M1 + VET (vet_gd), general EOS.
2. **dt:** mostly the 1000 s cap; dips to 795-956 s after hesdirk2 retry/BE redos (radiation-admissibility limited, not CFL). Mean over the current link 948 s.
3. **Speed, current link** (sim t 4.0090e8 -> 4.9557e8 s = 1096 star-days, 43.6 wall h, 99 900 cycles, restart overhead of this link included):
   - **1.57 wall s per cycle**, **~25 star-days per wall hour**, **~0.32 GPU-h per star-day**.
   - vs the 10-25 fresh run (12 H200, 1.14 s/cycle, ~37 star-days/h): fewer GPUs (8 vs 12) -> slower per cycle, but per-GPU throughput is about the same (12.6 vs 13.7 GPU-s per cycle, ~8 % better now). dt unchanged (cap 1000 s). Throughput steady over the last ~2 days (620-634 sim s per wall s).
4. **M1 share of wall time:** no per-module timers in the production log; not available. (From 09-27 profiles on dhj, the implicit radiation solve dominated; not measured for RSG E.)
5. **Sim time now:** 4.9557e8 s = 5736 d. Several links/restarts since the 10-25 fresh start; the numbers above are from the current link only (one start-up included).

Comparison with your numbers: AthenaK RSG E ~25 star-days/h on 8 H200 vs Arepo g25_a120 ~3.4 star-days/h on 560 cores (different physics/angles: AthenaK M1+VET, not discrete ordinates).
