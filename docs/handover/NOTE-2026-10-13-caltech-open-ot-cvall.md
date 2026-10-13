# NOTE 2026-10-13 from Caltech: open-ot CVALL (re TASK-2026-10-13-caltech-open-ot-cvall.md)

## QUEUED (10-12 ~22:45 PDT)
- OLD = athena_cpu_omp_281dd686, md5 9ae9f09aeb81841a275e8d046a0eebe4; NEW = 281dd686 + patch/0001 (CV binary), md5 dc58375a9206e2c034fc5d57cfd8ac9d.
- Debug QoS chain (1 node x 32 cores, 28-min cutoff each, finished steps skipped): CV rounds 4429421 -> 4429434 -> 4429435, then
  CVALL rounds 4429436 -> 4429437 -> 4429438 -> 4429439 -> 4429440 -> 4429441.
- run_cvall.sh unchanged except /usr/bin/time -> stand-in; P = 6 (each dir runs OLD + NEW side by side at 2-4 threads; 16 would oversubscribe 32 cores ~4x).
  m1_shock_ref in the 5 rs decks -> our copies of caltech_od/shockref/ref_m*.txt.
- FULL steps longer than ~28 min cannot finish in debug (restart from scratch); they will be listed. Interim bit_all counts here every few hours.

## INTERIM after round 1 (job 4429436, 10-13 00:29 PDT)
- BIT done for 355 / 359: **131 BITWISE, 210 DIFF, 14 FAIL**; 59 FULL runs finished (rc=0) in the same round.
- **All 14 FAIL = a runner bug, not the patch:** OLD and NEW both FATAL with `More than one history, event log, or restart output block found
  in input file` -- these decks already have an rst <output> block and run_cvall.sh appends <output99> rst. Dirs: cart/hcyl/FIN3_n512,
  cart/pipe/FIN3_c2200, lad/pulse_k0.128/FINs_hesdirk2_C_512, sp/patm_L16_s/HFfZ_1841da_{64,128,256}, sp/patmT_L16_s/HFfZ_1841da_{64,128,256},
  sp/patmT_L16_e2_s/HFfZ_1841da_{64,256}, sp/patmT_L16_e4_s/HFfZ_1841da_256, sp/nlte2_r1000_tb/HFfSgZ_74a24b_320, sp/nlte2_tb_cart/HCFfZ_d842e3_320.
- Caltech fix (same for OLD and NEW): rstonce.py -- if the deck has an rst block, set its dt = 1e30 there; else append <output99> as before.
  The 14 bit.txt were cleared and are redone in the next rounds (the chain continues: NP2 4433856, then CVALL 4429437..41).

## FINAL BIT (10-13 01:01 PDT): 359 / 359 -- 145 BITWISE, 214 DIFF, 0 FAIL
- The 14 duplicate-rst FAILs re-done with the rstonce.py fix: all now BITWISE or DIFF.
- Note: 3 dirs carry time/nlim in args.txt, which overrides the BIT step's nlim=10 (e.g. cart/hcyl/FIN3_n512 ran 128 cycles for OLD and NEW; BITWISE).
- FULL step: 59 FULL runs finished (all rc=0) in round 1 before NPALL superseded FULL; FULL stopped after that (runner returns after BIT).
- docs/handover/caltech_cvall/bit_all.txt (359 lines) + cvall_light.tgz (29 MB: bit.txt, run.log, rc/time, tab/ incl. the 59 full/, no bin/, no rst).
- Full incl. bin/ is 350 MB gz (no rst): kept on Caltech /resnick/groups/carnegie_poc/jingze/cvall_1013/cvall_runs_full.tgz (md5 846eec3b87cb88f154e7e7299742123d);
  ask if you need specific bin/ files (we push them in < 45 MB pieces).
