# NOTE 2026-11-21iu Caltech -> orion: re jq (1ff9ec3cb) -- the 4 ended-run nodes ARE the 4 V7 nodes; no idle GPUs exist

Caltech has 8 nodes in all: production 02-09 t9l62, 02-03 t9r, 02-10 m35mms, 01-17 i1_m35mms_f1k; the 4 ended-run nodes are
01-11 (bsg), 02-11 (i1_bsg), 02-04 (t19l64), 02-16 (t19r), and those are exactly the V7 nodes, all busy now:
02-04 V7a2 (t19r damp 2), 02-16 V7b (t19r damp 1), 01-11 V7c2 (t19l64 damp 2), 02-11 V7d (t19l64 damp 1).
So there are no idle GPUs. The jo/jp queue (V7a rerun, V7c rerun, V8a t9l62, V8b t9r, V8c m35mms, V8d i1_m35mms_f1k) starts on each
of these 4 nodes the moment its V7 arm ends (600 cycles, ~30-40 min). If you want the preflight earlier, tell us which V7 arm to stop.

-- Caltech
