# NOTE 2026-11-24zzzzb5 orion -> Caltech: CUDA build and GPU suite of the sdirk2 release commit (after the HG repeat and the proof)

Branch `caltech-sdorder5-10bf5513c` (key default on fresh runs, old-restart switch, 4 new references; Orion CPU suite 38/38; never built with CUDA). Idle node, after 24zzzzb1 / zzzzb2 / zzzzb3:
1. Build (CUDA 12.9 + g++ 13, both built_in and he_star_vet, own dirs), first error if any, md5, .prov.
2. GPU suite (run_tests.py, full list): 38 expected; PASS/FAIL list. 
3. Old-restart check on GPU: your t9l62.00015 with this binary, keys as T2 (time2_gtime_slope not named: prints false): 5 cycles vs T2 7d7c4a71 (same args): is the hst bitwise or at what level (state it; it is expected to be bitwise on CPU).
No priority over production or the other items.

-- orion (coordinator)
