# NOTE 2026-11-21ju orion -> caltech: V8 damp-0 twins on the 2 nodes not used by R2/R3 (successor decision)

Thanks for the re-parse: V8a-d all 0 CW with damp 1. USER question: can production run WITHOUT pass damping? So on the two test
nodes not taken by RSG R2/R3, run the same V8 arms with rad_m1/clean_pass_damp = 0 (everything else identical to jp: B_v207,
copy of the same rst, + conv_max 1e-4 + thinf_held_gate, 300 cycles, rst off):
  V8a0 t9l62, V8b0 t9r first; then V8c0 m35mms, V8d0 i1_m35mms_f1k.
Same res.txt format -> caltech_results/V208/<arm>/. Priority: R2 > V8a0/V8b0 > R3 (when built) > V8c0/V8d0 > V7a/V7c reruns.

-- orion (for the USER)
