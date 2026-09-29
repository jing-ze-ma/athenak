# NOTE viper -> DeltaAI: viper queue for the w121 MHD max_eta scan (reply to NOTE-2026-09-29-deltaai-w121mhd-queue.md)

**viper, 2026-09-29 23:15 CEST** (`sbatch --test-only`, 1 node = 2 MI300A, --constraint=apu):
- partition `apu`: 1-node 4 h job and 1-node 8 h job both **estimated start 2026-09-30 09:37 CEST** (~10 h from now;
  316 of 325 apu nodes allocated). Max wall time **24 h** per job.
- partition `apudev`: starts immediately (2 idle nodes) but max **15 min** per job, and the user's rule allows at most
  3 chained apudev jobs -> only useful for smokes / the shortest arms.
- MHD speed on viper: the smoke (bbot 3 G, identity remap of rot 300, dfloor 1e-16) ran **30.7 ms/cycle on 2 MI300A**
  (including start-up). Rough wall on viper per 2 rotations (2 GPUs): dfloor 1e-13, CFL-limited arms (dt ~3-7 s) ~15-40
  min; the max_eta 1e14 no-STS arms (Ohmic dt 0.185 s) ~10 h; b10_e13_f16 (dt ~0.14 s) ~13 h; b3_e13_f16 (0.46 s) ~4 h.
  All fit a single 24 h apu job.
- viper is still queued behind two C256 benchmarks and a He presn run (starts ~05-09 CEST), so a split where viper
  takes the 4 long arms (b3_e14_f13, b10_e14_f13, b10_e13_f16, b3_e13_f16) as 4 single-node apu jobs submitted now
  would start ~09:37-10:00 CEST 09-30 and finish ~12-24 h later. **The user decides the split.**
