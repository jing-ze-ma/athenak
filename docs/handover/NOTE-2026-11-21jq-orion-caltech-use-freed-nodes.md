# NOTE 2026-11-21jq orion -> caltech: also use the 4 idle nodes of the finished/stopped runs for the jo/jp queue

The nodes whose production runs have ENDED (bsg and i1_bsg finished at tlim; t19l64 stopped by the USER; t19r dead) are idle
inside their holders until 17:46 PST. USER: no idle GPUs. Start the queued jo/jp items there now (one GPU job per node,
GPU-free check first), in the same order: V7a rerun, V7c rerun, V8a t9l62, V8b t9r, then V8c/V8d on the next free node.
Never a node whose production run is still running (t9l62, t9r, m35mms, i1_m35mms_f1k). No other change.

-- orion (for the USER)
