# implicit_predictor = step in restart files, and made the default (0923)

Patch: /viper/ptmp2/jinma/m1_predrst_0923/predrst.patch (relative to 66d5c59a, not committed).

## What changed
- `radm1::kM1PredRstMagic = "M1PRED01"` (rad_m1_implicit.hpp): a marked header block
  (int32 have = 1, int32 ncomp = 3, Real pred_dt) written behind the RTWARM01 header,
  only when `impl_pred && pred_ok`; otherwise not one byte is written.
- The 3 channels of `ipred` are appended as per-MeshBlock slabs behind the warm-start
  slabs (gid order, same loop as wtemp/wder/warm), so restarts work across rank counts.
- Reader (pgen.cpp): the same 8-byte peek. The header gives the tail length. The slabs go
  straight into `ipred`, which the RadiationM1 constructor has already allocated, and set
  `pred_ok` and `pred_dt`. A file without the block loads as before: pred_ok = false, the
  first step starts cold, and a warning is printed if the run wants a predictor.
- Default: `implicit_predictor = fixcl ? step : none` (eddington, vet_sc, tau).
- Caveat: a restart file written before this change has `implicit_predictor = none` in its
  embedded parameters, so it keeps running with none. Add
  `rad_m1/implicit_predictor=step` on the command line to switch it on.

## Gates (CPU, gcc/14 openmpi/5.0; runs in /viper/ptmp2/jinma/m1_predrst_0923/runs)
2-D slab he_slab_m1_2d_V3edd with meshblock/nx2=16 (2 MBs). N = 20 steps, restart at 10.
The md5 is taken over the payload after `<par_end>` (first 10 hex digits: hydro_w, m1).
- Gate 1 (HEAD vs new). closure=m1: straight 5927458cbe/1bb5803a89 in both builds,
  restarted 6b9d5f1cbe in both. implicit_predictor=none: 3f54e97292/2e0cb14447 and
  b49612c704 in both. The restart files at the half-way point and at the end are
  byte-identical (cmp).
- Gate 2 (predictor on, the new default). Straight vs restarted:
  - edd, 1 rank: b25422f7f2 = b25422f7f2
  - edd, 2 ranks: 02e31291bd = 02e31291bd
  - vet_sc, 1 rank: cfad9d3de9/33cef1a236 = the same
  - vet_sc, 2 ranks: fcb0afb1f3/7e153919d6 = the same

  1 rank vs 2 ranks is not bitwise even in a straight run, with or without the predictor
  (s1 != s2; the closure=m1 control also differs), so the 1->2 check is a round trip:
  the 1-rank file was read on 2 ranks and written again after 0 steps. The payload md5 is
  identical: edd 71f7dd46..., vet 9a343269... The block is present (M1PRED01 count 1).
- Gate 3: a HEAD-written Edd restart read by the new binary gives the same final dumps
  (3f54e97292) as HEAD, and the restart file is byte-identical. With the override =step it
  runs, prints the warning and starts cold.
- Picard mean over 20 steps on this slab: predictor 2.95, none 2.90 (log.txt). The gain
  is measured on the 3-D box, README_PICARD.md.

## Gate 4 (GPU, apudev)
/viper/ptmp2/jinma/m1_predrst_0923/gpu/gpu_gate.sh runs he_slab_m1_3d on 2 ranks: 40 steps
straight vs 20 + restart + 20. The job is submitted automatically once the HIP build
finishes (job id in gpu/submit.out). The result is in gpu/gpu_gate.<jobid>: the s2 and r2
md5 lines must match. PENDING at the time of writing.
