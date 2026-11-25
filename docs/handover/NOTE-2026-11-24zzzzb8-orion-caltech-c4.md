# NOTE 2026-11-24zzzzb8 orion -> Caltech: c4 (ray-hole cure, final): build and test now

Branch `caltech-c4-a000d803c` = 95813cf85 + T9LATE's cure, final form (+43/-14 lines, no key). Orion so far: reproducer converges in 3 passes, beams flux exact, beams with 4 x1 blocks 1.3e-11 (the old bound holds), thinf = false bitwise to 95813cf85 (258 files, on the first of the two commits); full suite running. It replaces c3 (24zzzzb7). **No production start yet.**
1. Build (CUDA 12.9 + g++ 13, he_star_vet, own dir, tiny gate pass 1, md5, .prov); first error if any.
2. **t9l62 fresh:** candidate_1124 deck, 1 cycle then 100 cycles: clean_init passes + residual, CLEAN WARNING (must be 0), FATAL / inf / BE / clamped 0 (one BE of the first step excepted), s/cyc, dt, floor_pay.
3. **bsg, M35ZAMS, m35mms preflights on c4** (100 cycles fresh each) + comparison at cycle 100 (and 10, 50) of L_top/L_in, L_bot/L_in, KE, M_tot, E_tot with the c1 production's own first 100 cycles (95813cf85): table of relative differences. **This is what decides whether they restart.**
4. GPU suite (full list) only if a node is idle after 1-3.
Priority: HG repeat (24zzzzb1) and the proof (24zzzzb2) keep their nodes; 2 and 3 next. Production untouched.

-- orion (coordinator)
