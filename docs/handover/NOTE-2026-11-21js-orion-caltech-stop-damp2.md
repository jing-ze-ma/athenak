# NOTE 2026-11-21js orion -> caltech: stop V7a2 and V7c2 (damp 2) now; start V8 on those two nodes

USER: clean_pass_damp = 2 is decided against (fails on the RSG stall case), so no damp-2 tests.
1. Stop V7a2 (t19r damp 2, 02-04) and V7c2 (t19l64 damp 2, 01-11) now — our own test jobs only. Post what they have (res.txt) as is.
2. Start V8a t9l62 and V8b t9r on 02-04 and 01-11 right away (as in jp/jr).
3. V7b / V7d (damp 1) run to their end; then V8c m35mms, V8d i1_m35mms_f1k on those nodes; then the V7a/V7c damp-0 reruns.

-- orion (for the USER)
