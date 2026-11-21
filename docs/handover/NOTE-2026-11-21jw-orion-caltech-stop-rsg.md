# NOTE 2026-11-21jw orion -> caltech: R2 shows rsg.01295 does not stall -> stop the RSG arms (R3, and E1 from jv)

R2 (v20.7, damp 0, rsg.01295): 0 CW, no 200-solve stage -> this restart does not reproduce the stall, so wb_local / ew_max arms on it
cannot decide anything. Unnecessary tests (USER):
1. Stop R3 now (our test job). Do NOT run E1. (B_wbl / the c2eb74fac rebuild are not needed.)
2. Use the freed nodes, in order: V8c0 m35mms, V8d0 i1_m35mms_f1k (damp 0 twins, as ju); then the V7a / V7c T19 damp-0 reruns.
The RSG stall tests stay on viper (rsg.01510).

-- orion (for the USER)
