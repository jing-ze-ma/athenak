---
name: overnight-0930
description: User asleep from 09-30 ~10:00 CEST: what runs, what I may do alone, what waits for the user
metadata:
  type: project
---
Queued on viper apu (est. ~10:45-13:10): C256 benchmarks 12028087 (16 N) / 12026299 (32 N); He presn hllc 64x64
12037044(+45) and 128x128 12037023-25; WASP-121b fixed arms w1xf 12037717(+18), w10xf 12037719(+20), w3xk 12034178(+79)
(binary w121-build-0930 44be5b99, radial fix on, original solver keys). HELD: MHD arms 12028523-25.
ALLOWED alone: let these run; radial-fix agent may merge the dhj-only lhllc_x1_phi_min default ONLY if its combined gate
passes; verify reports; stop a clearly broken run; no new threads, no production launches, no default changes.
Morning summary: C256 benchmark numbers + machine choice; He presn hllc 64 vs 128 early turnovers; w1xf/w10xf/w3xk
first rotations (dt vs pre-fix, odd-even share); Caltech 3x (ETA ~15:00 CEST) + DeltaAI MHD scan status; radial-fix merge.

**Status 09-30 18:40 CEST:** viper apu backlog, nothing but the C256 16-node bench has run. 12028087 (16 N, 32 GPUs,
mb32) OK: nlim 70, 0.118 s/cycle steady (cycles 50-70), dt 3.7 s -> ~0.98 h wall per rotation (~31 GPU-h/rot) at the
start dt; ~8.5e8 zone-cycles/s. 32-N bench 12026299 est. 21:58; w3xk 19:27, w1xf/he3d 19:31, w10xf 21:24.
**20:40:** He presn 12037023/12037044 died in 9 s: r2.sh/r3.sh pass rad_m1/implicit_opac_newton=true but the input
lacked the key (parser fatal). Fixed by adding it to he3d_M1.athinput (backup .pre_opn); follow-ons 12037024/12037045
now run as first links (one link each lost). scontrol hold was blocked by the permission classifier -> ask user.
WASP-121b w1xf/w10xf/w3xk running since ~20:00: rot 14 / 8 / 12 after 33 min, dt 12.7 / 8.5 / 12.0 s, clean.
**20:40 smoke 12045717 PASSED (s128/s64: 10 cycles, rc 0, 0 FATAL, 0 guard drops, 0 bicgstab breakdowns); released 12037024/12037045.**
